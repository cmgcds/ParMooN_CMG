// =======================================================================
//
// Purpose:     TNSE2D Code with the implementaiton  of DO equations - COEFFICIENT PART ALONE
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 03.09.2014

// =======================================================================

#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <LinAlg.h>
#include <FESpace2D.h>
#include <SystemTNSE2D.h>
#include <SquareStructure2D.h>
#include <Structure2D.h>
#include <Output2D.h>
#include <CD2DErrorEstimator.h>
#include <MainUtilities.h>
// #include <TimeUtilities.h>
#include <TNSE2D_ParamRout.h>
#include <TimeDiscRout.h>

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>


#include<QuadAffin.h>
#include<QuadBilinear.h>

#ifdef INTELMKLBLAS

#endif
#include<mkl.h>
#include<omp.h>
#include "../Examples/TNSE_2D/DrivenCavity.h" // smooth sol in unit square

int main(int argc, char *argv[])
{
    //------- Initialise all the Essential Variables  ----- //
    TDomain *Domain;
    TDatabase *Database = new TDatabase();
    TFEDatabase2D *FEDatabase = new TFEDatabase2D();

    // Collection of Cells
    TCollection *coll, *mortarcoll = NULL;

    // -- finite Element Space fo Velocity and pressure
    TFESpace2D *Velocity_FeSpace, *Pressure_FeSpace, *fesp[2];
        // -- Fe Functions for Velocity and pressure --- //
    TFEVectFunct2D *Velocity_mean;
    TFEFunction2D *u1_mean, *u2_mean, *Pressure_mean, *fefct[2];


    // Class to Write the System output to VTK File
    TOutput2D *Output;

    // -- Parameter to pass older variables into the Assembly function --
    TAuxParam2D *aux, *NSEaux_error;

    // -- Multi Index to be used while defining the Assembly function
    MultiIndex2D AllDerivatives[3] = {D00, D10, D01};


    // Variables initialisation for output and Definitions
    const char vtkdir[] = "VTK";
    char *PsBaseName, *VtkBaseName, *GEO;
    char UString[] = "u";
    char PString[] = "p";
    char NameString[] = "VMS";

    // Set the output stream to write VTK files
    std::ostringstream os;
    os << " ";


     mkdir(vtkdir, 0777);

    // Write the accessories output file and check parameter consistencies for NSE
    OpenFiles();
    OutFile.setf(std::ios::scientific);

    Database->CheckParameterConsistencyNSE();
    Database->WriteParamDB(argv[0]);
    Database->WriteTimeDB();
    ExampleFile();

    // ------ DOMAIN CREATION ------------------ //
    // Create a Domain based on the input mesh File provided.
    Domain = new TDomain(argv[1]);
    GEO = TDatabase::ParamDB->GEOFILE;
    Domain->Init(NULL, GEO);

    // refine grid up to the coarsest level
    for (int i = 0; i < TDatabase::ParamDB->UNIFORM_STEPS; i++)
        Domain->RegRefineAll();

    int ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
    int NSEType = TDatabase::ParamDB->NSTYPE;
    int Disctype = TDatabase::ParamDB->DISCTYPE;
    int pressure_space_code, velocity_space_code;
    coll = Domain->GetCollection(It_Finest, 0);
    int N_Cells = coll->GetN_Cells();
    OutPut("N_Cells : " << N_Cells << endl);

    // fespaces for velocity and pressure
    GetVelocityAndPressureSpace(coll, BoundCondition, mortarcoll, Velocity_FeSpace,
                                Pressure_FeSpace, &pressure_space_code,
                                TDatabase::ParamDB->VELOCITY_SPACE,
                                TDatabase::ParamDB->PRESSURE_SPACE);

    // defaulty inf-sup pressure space will be selected based on the velocity space, so update it in database
    TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE = pressure_space_code;
    velocity_space_code = TDatabase::ParamDB->VELOCITY_SPACE;

    int N_U = Velocity_FeSpace->GetN_DegreesOfFreedom();
    int N_P = Pressure_FeSpace->GetN_DegreesOfFreedom();
    int N_TotalDOF = 2 * N_U + N_P;

    OutPut("Dof Velocity : " << setw(10) <<  N_U << endl);
    OutPut("Dof Pressure : " << setw(10) << N_P << endl);
    OutPut("Total Dof all: " << setw(10) << N_TotalDOF << endl);



    //--- Setup arrays and Functions for System Matrices and Arrays ( AX=b )
    double* sol     = new double[N_TotalDOF]; 
    double* rhs     = new double[N_TotalDOF];
    double* oldrhs  = new double[N_TotalDOF];
     
    memset(sol, 0, N_TotalDOF*SizeOfDouble);
    memset(rhs, 0, N_TotalDOF*SizeOfDouble);

    Velocity_mean = new TFEVectFunct2D(Velocity_FeSpace, UString,  UString,  sol, N_U, 2);
    u1_mean = Velocity_mean->GetComponent(0);
    u2_mean = Velocity_mean->GetComponent(1);
    Pressure_mean = new TFEFunction2D(Pressure_FeSpace, PString,  PString,  sol+2*N_U, N_P);
    
    //  interpolate the initial solution based on the initial condition of the problem
    u1_mean->Interpolate(InitialU1);
    u2_mean->Interpolate(InitialU2);
    Pressure_mean->Interpolate(InitialP);

   

    // parameters for time stepping scheme
	int m = 0;
	int N_SubSteps = GetN_SubSteps();
	double oldtau = 1.;
	double end_time = TDatabase::TimeDB->ENDTIME;
	int limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE;
	int Max_It = TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE;
	double tau;
	double* defect;
	double residual;
	double impuls_residual;

    int N_S =  50;
    int N_R  = 100;


    // Declare the variables here 
    double* p_tilde = new double[N_P*N_S];
    for( int i = 0 ; i < N_P*N_S ; i++)   p_tilde[i] = 1.0;

    // u_tilde      -- is a FEVectFunction with N_S Components 
    double* u_tilde  = new double[N_U * N_S];
    for( int i = 0 ; i < N_U*N_S ; i++)   u_tilde[i] = 1.0;   
    // v_tilde      -- is a vectfunction with N_S Components
    double* v_tilde  = new double[N_U * N_S];
    for( int i = 0 ; i < N_U*N_S ; i++)   v_tilde[i] = 1.0;
    // u bar
    double* u_bar = new double[N_U];
    for( int i = 0 ; i < N_U ; i++)   u_bar[i] = 1.0;
    //v bar
    double* v_bar = new double[N_U];
    for( int i = 0 ; i < N_U ; i++)   v_bar[i] = 1.0;
    //Co-efficient
    double* CoeffVector = new double[N_R*N_S];
    for( int i = 0 ; i < N_R*N_S ; i++)   CoeffVector[i] = 1.0;
    //ans
    double* ans = new double[N_R*N_S];
    for( int i = 0 ; i < N_R*N_S ; i++)   ans[i] = 1.0;
    //cov 
    double* covariance = new double[N_S*N_S];
    for( int i = 0 ; i < N_S*N_S ; i++)   covariance[i] = 1.0;

    TFEVectFunct2D* p_tilde_vectFunction = new TFEVectFunct2D(Pressure_FeSpace,PString,PString,p_tilde,N_P,N_S);
    TFEVectFunct2D* u_tildeVectfunction = new TFEVectFunct2D(Velocity_FeSpace,UString,UString,u_tilde,N_U,N_S);
    TFEVectFunct2D* v_tildeVectfunction = new TFEVectFunct2D(Velocity_FeSpace,UString,UString,v_tilde,N_U,N_S);

    // Setup Data Structures for Cells 
    int* GlobalNumbers = Velocity_FeSpace->GetGlobalNumbers();
	int* BeginIndex = Velocity_FeSpace->GetBeginIndex();

    // Setup Data Structures for Cells 
    int* GlobalNumbers_pressure = Pressure_FeSpace->GetGlobalNumbers();
	int* BeginIndex_pressure = Pressure_FeSpace->GetBeginIndex();	

    int LocN_BF[N_BaseFuncts2D];
    BaseFunct2D LocBF[N_BaseFuncts2D];
    BaseFunct2D *BaseFuncts;
    int *N_BaseFunct;

    // Get the Basis Funciton Details
    BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
    N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();


	// --- Quadrature Formula Arrays  ------------------//
	int N_Points2;
	double *Weights2, *t1, *t2;      // Weights - Quadrature Weights array ; t1  -- Quadrature point ( xi ) in ref coordinate ; t2 -  Quadrature Point ( eta ) in Ref Coordinate 
	bool Needs2ndDer[1]; 
	Needs2ndDer[0] = TRUE;
    bool Needs2ndDerPressure[1]; 
	Needs2ndDerPressure[0] = FALSE;
	double AbsDetjk[MaxN_PointsForNodal2D];
	double X[MaxN_PointsForNodal2D];
	double Y[MaxN_PointsForNodal2D];


	// FE Values Arrays 
	double** origvaluesD00;           // Shape function values at quadrature Points
	double** origvaluesD10;           // Shape Function Derivatives ( x ) at Quadrature Points
	double** origvaluesD01;           // Shape Function Derivatives ( y ) at Quadrature Points
    double** origvaluesD20;           // Shape Function 2nd Derivatives ( x ) at Quadrature Points
	double** origvaluesD02;           // Shape Function 2nd Derivatives ( y ) at Quadrature Points

    double** origValuesD00_press;
    double** origValuesD10_press;     // Shape Function Derivative Pressure ( x )
    double** origValuesD01_press;     // Shape Function Derivative Pressure ( y )

    cout << "  Start of Cell Looping " <<endl;

    double start = omp_get_wtime();

    for ( int cellId = 0 ; cellId < 1 ; cellId++)
    {
        // cout << " Cell - id : " << cellId <<endl;
        // Get the Cell object from the Cells Collection based on cell id 
        TBaseCell* currentCell = coll->GetCell(cellId);
        // Get the "ID" of Finite Element for the given 2D Element ( Conforming/NonConforming-Order Finite Element : eg : it could be Conforming-2nd order Finite Element )
		FE2D elementId = Velocity_FeSpace->GetFE2D(cellId, currentCell);
        // Get the Class object for that 2d FEM Element , which has all the details like Shape functions , Nodal point locations for that location, Reference Transformation ( Affine , Bilinear )
        TFE2D *element = TFEDatabase2D::GetFE2D(elementId);
        TFEDesc2D *fedesc = element->GetFEDesc2D();
        // Class for basis functions in 2D ( Number of basis functions ), basis function values and Derivatives
        TBaseFunct2D* bf = element->GetBaseFunct2D();
        // Get the Reference Elemet
		BF2DRefElements RefElement = TFEDatabase2D::GetRefElementFromFE2D(elementId);
		// Get the reference Transformation -- Affine Mapping / Bilnea Mapping of Triangle or Quadrilateral
		RefTrans2D referenceTransformation =  TFEDatabase2D::GetRefTrans2D_IDFromFE2D(elementId);
        // Get the number of basis functions in the Current Cell ( Number of Local DOF)
		int N_BaseFunct = element->GetN_DOF();
		// Type of Basis Function in 2D
		BaseFunct2D BaseFunct_ID = element->GetBaseFunct2D_ID();

        // Get the pressure values and DOF's 
        FE2D elementId_pressure = Pressure_FeSpace->GetFE2D(cellId, currentCell);
        TFE2D *elementPressure = TFEDatabase2D::GetFE2D(elementId_pressure);
        TBaseFunct2D* bf_pressure = elementPressure->GetBaseFunct2D();
        BF2DRefElements RefElementPressure = TFEDatabase2D::GetRefElementFromFE2D(elementId_pressure);
        RefTrans2D referenceTransformationPressure =  TFEDatabase2D::GetRefTrans2D_IDFromFE2D(elementId_pressure);
        int N_BaseFunct_pressure = elementPressure->GetN_DOF();
        BaseFunct2D BaseFunct_ID_pressure = elementPressure->GetBaseFunct2D_ID();

        // get cell measure
        double hK = currentCell->GetDiameter();

        switch (referenceTransformation)
		{
			case QuadBilinear:
			{
				int l = bf->GetPolynomialDegree();  												// Get the Polynomial Degreee  of the basis functions
				QuadFormula2D QF2 = TFEDatabase2D::GetQFQuadFromDegree(3*l);						// Get te ID of Quadrature Formula
				TQuadFormula2D* QuadratureRule = TFEDatabase2D::GetQuadFormula2D(QF2);              // Get the Quadrature Rule Objetc based on Quadrature ID
				QuadratureRule->GetFormulaData(N_Points2, Weights2, t1, t2);                        // get the Quadrature points , Weights 

				// Set the values on the Reference Cell
				TRefTrans2D* F_K = TFEDatabase2D::GetRefTrans2D(QuadBilinear);
				TFEDatabase2D::SetCellForRefTrans(currentCell, QuadBilinear);                              	// Set the Cell for Current reference Transformation

				// Get Original Coordinates from reference Coordinates and the Determinant of jacobian
				TFEDatabase2D::GetOrigFromRef(QuadBilinear,N_Points2, t1,t2,  X,Y, AbsDetjk);               // Get the Original Co-orinates for the cell from xi values

				// Get all the original Values from the Referece cell values.
				TFEDatabase2D::GetOrigValues(QuadBilinear, 1, &BaseFunct_ID, N_Points2, t1, t2, QF2, Needs2ndDer);
				

				// The below are 2D arrays in the form 
				// Values[QuadraturePointLocation][ShapeFunction]  i.e, the Value of Shapefunction at all quadrature points for each shape functions
				origvaluesD00       = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D00);            	// Shape Function Values at Quadrature Points 
				origvaluesD10       = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D10);               // Shape Function Derivative Values at Quadrature Points
				origvaluesD01       = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D01);				// Shape Function Derivative Values at Quadrature Point
                origvaluesD20       = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D20);				// Shape Function 2nd Derivative Values at Quadrature Point
                origvaluesD02       = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D02);				// Shape Function 2nd Derivative Values at Quadrature Point
				break;
			}


			case QuadAffin:
			{
				int l = bf->GetPolynomialDegree();  												// Get the Polynomial Degreee  of the basis functions
				QuadFormula2D QF2 = TFEDatabase2D::GetQFQuadFromDegree(3*l);						// Get te ID of Quadrature Formula
				TQuadFormula2D* QuadratureRule = TFEDatabase2D::GetQuadFormula2D(QF2);              // Get the Quadrature Rule Objetc based on Quadrature ID
				QuadratureRule->GetFormulaData(N_Points2, Weights2, t1, t2);                        // get the Quadrature points , Weights 

				// Set the values on the Reference Cell
				TRefTrans2D* F_K = TFEDatabase2D::GetRefTrans2D(QuadAffin);
				TFEDatabase2D::SetCellForRefTrans(currentCell, QuadAffin);                              	// Set the Cell for Current reference Transformation

				// Get Original Coordinates from reference Coordinates and the Determinant of jacobian
				TFEDatabase2D::GetOrigFromRef(QuadAffin,N_Points2, t1,t2,  X,Y, AbsDetjk);   // Get the Original Co-orinates for the cell from xi values

				// Get all the original Values from the Referece cell values.
				TFEDatabase2D::GetOrigValues(QuadAffin, 1, &BaseFunct_ID, N_Points2, t1, t2, QF2, Needs2ndDer);
				

				// The below are 2D arrays in the form 
				// Values[QuadraturePointLocation][ShapeFunction]  i.e, the Value of Shapefunction at all quadrature points for each shape functions
				origvaluesD00 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D00);            	// Shape Function Values at Quadrature Points 
				origvaluesD10 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D10);             // Shape Function Derivative Values at Quadrature Points
				origvaluesD01 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D01);				// Shape Function Derivative Values at Quadrature Point
                origvaluesD20 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D20);             // Shape Function 2nd Derivative Values at Quadrature Points
				origvaluesD02 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D02);				// Shape Function 2nd Derivative Values at Quadrature Point

				break;
			}
					
			default:
			{
                cout << " [ERROR] - Error in File : CoeffEqn_DO.C " <<endl;
				cout << " Unknown Reftype " <<endl;
				cout << " REF TYPE : " << referenceTransformation <<endl;
				exit(0);
				break;

			}
		}

        int* DOF = GlobalNumbers + BeginIndex[cellId];

        ////   ========================================== VELOCITY FE SPACE LOOP ========================================== /////

        // `i` Loop is for the columns of the Co-efficient Matrix  
        // This picks the ith column of coefficient matrix, which is of Size (N_R * 1)
        // The total size of coeff(phi) is [N_R * N_S] 
        // Phi stored in Column major order (i.e) Each colums of phi matrix are laid flat and appended at end of prev col
        for ( int i=0 ; i < N_S ; i++)
        {
            // cout << "i : " << i <<endl;
            // Get the ith components of all the arrays
            double* localCoeff = ans + N_R*i;     // Answer vector
            double* u_tilde_i  = u_tilde + N_R*i;
            double* v_tilde_i  = v_tilde + N_R*i;

            //`g` loop is for each component in the ith column of the Coefficient vector
            // This has size of (N_R)
            for (int g = 0 ; g < N_R ; g++)
            {
                // cout << "i : " << i << "   g : " << g  <<endl;
                // This `a` loops runs for summations in the ODE terms
                for ( int a=0 ; a < N_S ; a++)
                {
                    // cout << "i : " << i  << "   g : " << g << "      a : " << a <<endl;
                    double* phi_a           = CoeffVector + N_R*a;
                    double* u_tilde_a       = u_tilde + N_U*a;
                    double* v_tilde_a       = v_tilde + N_U*a;
                    
                    double val = 0;
                    // Begin Quadrature Integration 
                    for ( int quadPt = 0 ; quadPt < N_Points2; quadPt++)
                    {
                        double Mult = Weights2[quadPt] * AbsDetjk[quadPt];
                        double* orgD00 = origvaluesD00[quadPt];
                        double* orgD10 = origvaluesD10[quadPt];
                        double* orgD01 = origvaluesD01[quadPt];
                        double* orgD20 = origvaluesD20[quadPt];
                        double* orgD02 = origvaluesD02[quadPt];

                        for (int j = 0 ; j  < N_BaseFunct ;  j++)
                        {
                            int GlobalDOF   = DOF[j];
                            // < DDx_u_tilde_a + DDy_u_tilde_a , u_tilde_i > 
                            val += ((u_tilde_a[GlobalDOF] * orgD20[j]) + (u_tilde_a[GlobalDOF] * orgD02[j])) * u_tilde_i[GlobalDOF]*orgD00[j] ;

                            // < DDx_v_tilde_a + DDy_v_tilde_a , v_tilde_i > 
                            val += ((v_tilde_a[GlobalDOF] * orgD20[j]) + (v_tilde_a[GlobalDOF] * orgD02[j])) * v_tilde_i[GlobalDOF]*orgD00[j] ;

                            // < ( u_bar*Dx_u_tilde_a + v_bar*Dy_u_tilde_a + u_tilde_a*Dx_u_bar + v_tilde_a*Dy_u_bar ) , u_tilde_i
                            double ubar  = u_bar[GlobalDOF]*orgD00[j];
                            double vbar  = v_bar[GlobalDOF]*orgD00[j];

                            double dxubar = u_bar[GlobalDOF]*orgD10[j];
                            double dyubar = u_bar[GlobalDOF]*orgD01[j];

                            double utilde_a = u_tilde_a[GlobalDOF]*orgD00[j];
                            double vtilde_a = v_tilde_a[GlobalDOF]*orgD00[j];

                            double dx_utilde_a  = u_tilde_a[GlobalDOF]*orgD10[j];
                            double dy_utilde_a  = u_tilde_a[GlobalDOF]*orgD01[j];

                            val += -(ubar*dx_utilde_a + vbar*dy_utilde_a + utilde_a*dxubar + vtilde_a*dyubar)*(u_tilde_i[GlobalDOF]*orgD00[j]);

                            // < ( v_bar*Dy_v_tilde_a + u_bar*Dx_v_tilde_a + v_tilde_a*Dy_v_bar + u_tilde_a*Dx_v_bar ) , u_tilde_i
                            double dxvbar = v_bar[GlobalDOF]*orgD10[j];
                            double dyvbar = v_bar[GlobalDOF]*orgD01[j];

                            double dx_vtilde_a  = v_tilde_a[GlobalDOF]*orgD10[j];
                            double dy_vtilde_a  = v_tilde_a[GlobalDOF]*orgD01[j];

                            val += -(vbar*dy_vtilde_a + ubar*dx_vtilde_a + vtilde_a*dyvbar + utilde_a*dxvbar)*(v_tilde_i[GlobalDOF]*orgD00[j]);

                            val *= Mult;

                        } // Local DOF loop
                    } // Quadraure loop


                    // Add the current contributiuon to `g`th component of the localCoeff vector
                    localCoeff[g]  += phi_a[g]*val;

                    // Runs the second loop of summations
                    for ( int b=0 ; b < N_S ; b++)
                    {
                        double* phi_b           = CoeffVector + N_R*b;
                        double* u_tilde_b       = u_tilde + N_U*b;
                        double* v_tilde_b       = v_tilde + N_U*b;

                        double val = 0;
                        //Quadrature Loop
                        for ( int quadPt = 0 ; quadPt < N_Points2; quadPt++)
                        {
                            double Mult = Weights2[quadPt] * AbsDetjk[quadPt];
                            double* orgD00 = origvaluesD00[quadPt];
                            double* orgD10 = origvaluesD10[quadPt];
                            double* orgD01 = origvaluesD01[quadPt];
                            double* orgD20 = origvaluesD20[quadPt];
                            double* orgD02 = origvaluesD02[quadPt];

                            for (int j = 0 ; j  < N_BaseFunct ;  j++)
                            {
                                int GlobalDOF   = DOF[j];
                                // < u_tilde_a*Dx_u_tilde_b + v_tilde_a*Dy_u_tilde_b > u_tilde_i
                                double utilde_a = u_tilde_a[GlobalDOF]*orgD00[j];
                                double vtilde_a = v_tilde_a[GlobalDOF]*orgD00[j];

                                double dx_utilde_b  = u_tilde_b[GlobalDOF]*orgD10[j];
                                double dy_utilde_b  = u_tilde_b[GlobalDOF]*orgD01[j];

                                val += -(utilde_a*dx_utilde_b + vtilde_a*dy_utilde_b)*u_tilde[GlobalDOF]*orgD00[j];

                                // < v_tilde_a*Dy_v_tilde_b + u_tilde_a*Dx_v_tilde_b > v_tilde_i
                                double dx_vtilde_b  = v_tilde_b[GlobalDOF]*orgD10[j];
                                double dy_vtilde_b  = v_tilde_b[GlobalDOF]*orgD01[j];

                                val += -(vtilde_a*dy_vtilde_b + utilde_a*dx_vtilde_b)*v_tilde[GlobalDOF]*orgD00[j];


                                val *= Mult;

                            }

                            localCoeff[g] +=  phi_b[g]*phi_a[g] - covariance[a*N_S + b]*val;
                            

                        }

                    }


                }

            }


        }




        ////   ========================================== PRESSURE FE SPACE LOOP ========================================== /////

        DOF = GlobalNumbers_pressure + BeginIndex_pressure[cellId];

        switch (referenceTransformation)
		{
			case QuadBilinear:
			{
				int l = bf->GetPolynomialDegree();  												// Get the Polynomial Degreee  of the basis functions
				QuadFormula2D QF2 = TFEDatabase2D::GetQFQuadFromDegree(3*l);						// Get te ID of Quadrature Formula
				TQuadFormula2D* QuadratureRule = TFEDatabase2D::GetQuadFormula2D(QF2);              // Get the Quadrature Rule Objetc based on Quadrature ID
				QuadratureRule->GetFormulaData(N_Points2, Weights2, t1, t2);                        // get the Quadrature points , Weights 

				// Set the values on the Reference Cell
				TRefTrans2D* F_K = TFEDatabase2D::GetRefTrans2D(QuadBilinear);
				TFEDatabase2D::SetCellForRefTrans(currentCell, QuadBilinear);                              	// Set the Cell for Current reference Transformation

				// Get Original Coordinates from reference Coordinates and the Determinant of jacobian
				TFEDatabase2D::GetOrigFromRef(QuadBilinear,N_Points2, t1,t2,  X,Y, AbsDetjk);               // Get the Original Co-orinates for the cell from xi values

				// Get all the original Values from the Referece cell values.
				TFEDatabase2D::GetOrigValues(QuadBilinear, 1, &BaseFunct_ID_pressure, N_Points2, t1, t2, QF2, Needs2ndDerPressure);

				// The below are 2D arrays in the form 
				// Values[QuadraturePointLocation][ShapeFunction]  i.e, the Value of Shapefunction at all quadrature points for each shape functions
				origValuesD00_press       = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID_pressure, D00);
                origValuesD10_press       = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID_pressure, D10);             // Shape Function Derivative Values at Quadrature Points
				origValuesD01_press       = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID_pressure, D01);				// Shape Function Derivative Values at Quadrature Point
				break;
			}


			case QuadAffin:
			{
				int l = bf->GetPolynomialDegree();  												// Get the Polynomial Degreee  of the basis functions
				QuadFormula2D QF2 = TFEDatabase2D::GetQFQuadFromDegree(3*l);						// Get te ID of Quadrature Formula
				TQuadFormula2D* QuadratureRule = TFEDatabase2D::GetQuadFormula2D(QF2);              // Get the Quadrature Rule Objetc based on Quadrature ID
				QuadratureRule->GetFormulaData(N_Points2, Weights2, t1, t2);                        // get the Quadrature points , Weights 
				// Set the values on the Reference Cell
				TRefTrans2D* F_K = TFEDatabase2D::GetRefTrans2D(QuadAffin);
				TFEDatabase2D::SetCellForRefTrans(currentCell, QuadAffin);                              	// Set the Cell for Current reference Transformation

				// Get Original Coordinates from reference Coordinates and the Determinant of jacobian
				TFEDatabase2D::GetOrigFromRef(QuadAffin,N_Points2, t1,t2,  X,Y, AbsDetjk);   // Get the Original Co-orinates for the cell from xi values

				TFEDatabase2D::GetOrigValues(QuadBilinear, 1, &BaseFunct_ID_pressure, N_Points2, t1, t2, QF2, Needs2ndDerPressure);
                // cout << " AFFINE MAPPING1 " <<endl;

				// The below are 2D arrays in the form 
				// Values[QuadraturePointLocation][ShapeFunction]  i.e, the Value of Shapefunction at all quadrature points for each shape functions
				origValuesD00_press       = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID_pressure, D00);
                origValuesD10_press       = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID_pressure, D10);                // Shape Function Derivative Values at Quadrature Points
				origValuesD01_press       = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID_pressure, D01);				// Shape Function Derivative Values at Quadrature Point
				break;
			}
					
			default:
			{
                cout << " [ERROR] - Error in File : CoeffEqn_DO.C " <<endl;
				cout << " Unknown Reftype " <<endl;
				cout << " REF TYPE : " << referenceTransformation <<endl;
				exit(0);
				break;

			}
		}


        // cout << " REFTRANS FOR PRESSURE Done " <<endl;

        for ( int i=0 ; i < N_S ; i++)
        {
            // Get the ith components of all the arrays
            double* localCoeff = ans + N_R*i;     // Answer vector
            double* u_tilde_i  = u_tilde + N_R*i;
            double* v_tilde_i  = v_tilde + N_R*i;

            //`g` loop is for each component in the ith column of the Coefficient vector
            // This has size of (N_R)
            for (int g = 0 ; g < N_R ; g++)
            {
                // This `a` loops runs for summations in the ODE terms
                for ( int a=0 ; a < N_S ; a++)
                {
                    double* phi_a           = CoeffVector + N_R*a;
                    double* u_tilde_a       = u_tilde + N_U*a;
                    double* v_tilde_a       = v_tilde + N_U*a;
                    
                    double val = 0;
                    // Begin Quadrature Integration 
                    for ( int quadPt = 0 ; quadPt < N_Points2; quadPt++)
                    {
                        double Mult = Weights2[quadPt] * AbsDetjk[quadPt];
                        double* orgD00 = origValuesD00_press[quadPt];
                        double* orgD10 = origValuesD10_press[quadPt];
                        double* orgD01 = origValuesD01_press[quadPt];

                        for (int j = 0 ; j  < N_BaseFunct ;  j++)
                        {
                            int GlobalDOF   = DOF[j];
                            // < Dx_p_tilde_a, u_tilde_i > 
                            val += -(p_tilde[GlobalDOF] * orgD10[j]) * u_tilde_i[GlobalDOF]*orgD00[j] ;

                            // < Dy_p_tilde_a  , v_tilde_i > 
                            val += -(p_tilde[GlobalDOF] * orgD01[j]) * v_tilde_i[GlobalDOF]*orgD00[j] ;
                            
                            val *= Mult;
                        } // Local DOF loop
                    } // Quadraure loop

                    // Add the current contributiuon to `g`th component of the localCoeff vector
                    localCoeff[g]  += phi_a[g]*val;

                }   // a loop ( N_S ) 

            }  // g loop ( N_R )


        } // i-loop ( N_S )




    }


    double end = omp_get_wtime();

    cout << " TIME TAKEN : " << end - start <<endl;
    
    return 0;




}