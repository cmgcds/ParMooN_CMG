// ======================================================================
// instationary problem
// ======================================================================

#include <MacroCell.h>
#include <IsoBoundEdge.h>
#include <IsoInterfaceJoint.h>

extern "C"
{
#include <gridgen.h>

	void triangulate(char *, struct triangulateio *,
					 struct triangulateio *, struct triangulateio *);
}

/// Include files related to the cell looping for the RHS filligng part.

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

#include <QuadAffin.h>
#include <QuadBilinear.h>

/// ========================================================================
// example file
// ========================================================================

#define __SIN3__

void ExampleFile()
{
	OutPut("Example: advection.h" << endl);
}

// exact solution
void Exact(double x, double y, double *values)
{
	double t;

	t = TDatabase::TimeDB->CURRENTTIME;

	values[0] = exp(t) * (sin(2 * Pi * x) * sin(2 * Pi * y));
	values[1] = exp(t) * 2 * Pi * cos(2 * Pi * x) * sin(2 * Pi * y);
	values[2] = exp(t) * 2 * Pi * sin(2 * Pi * x) * cos(2 * Pi * y);
	values[3] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int BdComp, double t, BoundCond &cond)
{
	//  cond = NEUMANN;

	//  if(BdComp == 1 || BdComp == 3 || BdComp == 2 )
	//  {
	cond = DIRICHLET;
	//  }
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
	value = 0;
}

// initial conditon
void InitialCondition(double x, double y, double *values)
{
	double t,val;

	t = TDatabase::TimeDB->CURRENTTIME;
	// values[0] = exp(t)*(sin(2*Pi*x)*sin(2*Pi*y));
	double val_x = exp(-1.0* ( 1.0 / ( 1 - pow((2*x-1),4) ) ));
	double val_y = exp(-1.0* ( 1.0 / ( 1 - pow((2*y-1),4) ) ));
	values[0] = val_x*val_y;
	
	double r = 0.1;
	double centerx = 0.5;	
	double centery = 0.5;
	x = x - centerx;
	y = y - centery;	
	double len = pow((x*x + y*y),0.5);
	if( abs(x) < r && abs(y) < r)
	{
		val = 1;
	}

	else
		val = 0;

	values[0] = val;

}


void DO_Mean_Equation_Coefficients(int n_points, double *X, double *Y,
								   double **parameters, double **coeffs)
{
	double eps = 1 / TDatabase::ParamDB->PE_NR;
	double b1 = 1, b2 = -1, c = 1;
	int i;
	double *coeff;
	double x, y;
	double t = TDatabase::TimeDB->CURRENTTIME;

	for (i = 0; i < n_points; i++)
	{
		coeff = coeffs[i];

		x = X[i];
		y = Y[i];

		coeff[0] = 1;
		coeff[1] = 0;
		coeff[2] = 0;
		coeff[3] = 0;

		coeff[4] = 0.0;
	}
}

void DO_Mode_Equation_Coefficients(int n_points, double *X, double *Y,
								   double **parameters, double **coeffs)
{
	double eps = 1 / TDatabase::ParamDB->PE_NR;
	double b1 = 1, b2 = -1, c = 1;
	int i;
	double *coeff;
	double x, y;
	double t = TDatabase::TimeDB->CURRENTTIME;

	for (i = 0; i < n_points; i++)
	{
		coeff = coeffs[i];

		x = X[i];
		y = Y[i];

		coeff[0] = 1;
		coeff[1] = 0;
		coeff[2] = 0;
		coeff[3] = 0;

		coeff[4] = 0.0;
	}
}

// ======================================================================
// ASSEMBLY FUNCTION
// This fucntion will be called for Every Quadrature Point inside a Cell for Local Assembly
// ======================================================================
void DO_Mean_Equation_Assembly(double quad_wt, double *coeff, double *param,
							   double hK, double **derivatives, int *N_BaseFuncts, double ***LocMatrices, double **LocRhs)
{

	// The below values N_x, N_y , etc are ARRAY of Values
	// Which provides Value of a Particular Shape function or its derivative ( Based on a DOF )
	// at the given Quadrature POint.
	// Here the Size of the Array is equal to the NUmber of DOF in the cell

	// For Eg : Orig0[1]  gives the derivative of Shape Function number 1 ( Shape function of DOF 1 ) w.r.t x
	// at the given Quadrature Point.
	// Take these values based on the MULTIINDEX DERIVATIVE
	double *N = derivatives[0];
	double *Nx = derivatives[1];
	double *Ny = derivatives[2];

	double **A11, *F1;
	double val = 0.;

	A11 = LocMatrices[0]; // Local Stiffenss matrix

	F1 = LocRhs[0];

	double c0 = coeff[0]; // nu
	double b1 = coeff[1]; // b1
	double b2 = coeff[2]; // b2
	double c = coeff[3];  // c
	double f = coeff[4];  // f

	int N_DOF_perCell = N_BaseFuncts[0];

	for (int i = 0; i < N_DOF_perCell; i++) // Test
	{
		// Assemble RHS
		F1[i] = 0.0; // TO DO
		double val = 0;
		for (int j = 0; j < N_DOF_perCell; j++) // Ansatz
		{

			val += (b1 * Nx[j] + b2 * Ny[j]) * N[i]; //   TO DO
			val  +=  c0*((Nx[j] * Nx[i])   + (Ny[j] * Ny[i]) );
			// val  +=  c*N[i]*N[j];

			A11[i][j] += val * quad_wt;
		} // endfor j
	}	  // endfor i
}

// ======================================================================
// ASSEMBLY FUNCTION
// This fucntion will be called for Every Quadrature Point inside a Cell for Local Assembly
// ======================================================================
void DO_Mode_Equation_Assembly(double quad_wt, double *coeff, double *param,
							   double hK, double **derivatives, int *N_BaseFuncts, double ***LocMatrices, double **LocRhs)
{

	// The below values N_x, N_y , etc are ARRAY of Values
	// Which provides Value of a Particular Shape function or its derivative ( Based on a DOF )
	// at the given Quadrature POint.
	// Here the Size of the Array is equal to the NUmber of DOF in the cell

	// For Eg : Orig0[1]  gives the derivative of Shape Function number 1 ( Shape function of DOF 1 ) w.r.t x
	// at the given Quadrature Point.
	// Take these values based on the MULTIINDEX DERIVATIVE
	double *N = derivatives[0];
	double *Nx = derivatives[1];
	double *Ny = derivatives[2];

	double **A11, *F1;
	double val = 0.;

	A11 = LocMatrices[0]; // Local Stiffenss matrix

	F1 = LocRhs[0];

	double c0 = coeff[0]; // nu
	double b1 = coeff[1]; // b1
	double b2 = coeff[2]; // b2
	double c = coeff[3];  // c
	double f = coeff[4];  // f

	int N_DOF_perCell = N_BaseFuncts[0];

	for (int i = 0; i < N_DOF_perCell; i++) // Test
	{
		// Assemble RHS
		F1[i] = 0.0; // TO DO
		double val = 0;
		for (int j = 0; j < N_DOF_perCell; j++) // Ansatz
		{

			val += (b1 * Nx[j] + b2 * Ny[j]) * N[i]; //   TO DO
			// val += c0 * ((Nx[j] * Nx[i]) + (Ny[j] * Ny[i]));
			// val += c * N[i] * N[j];

			A11[i][j] += val * quad_wt;
		} // endfor j
	}	  // endfor i
}

void DO_Mode_RHS_Aux_Param(double *in, double *out)
{
	out = in;
}

/// ---- RHS CELL ASSEMBLY LOOP FOR THE MODE EQUATION --- ///

// This function will take care of the Additional terms that go into both the co-efficient equation and the mode equation
// The terms in the Coefficient equation and the Mode equation are repeating , So the local RHS matrix will be reused in both the
//

void DO_Mode_RHS(TFESpace2D *Fespace, TFEVectFunct2D *FeVector_C, int N_S,double* GlobalRhs_mode , int i_index)
{

	int N_Cells = Fespace->GetN_Cells();
	TCollection* coll = Fespace->GetCollection();

	// Get the Global DOF arrays INdex from the FE Space. 
    int* GlobalNumbers = Fespace->GetGlobalNumbers();
	int* BeginIndex = Fespace->GetBeginIndex();

	// --- Quadrature Formula Arrays  ------------------//
	int N_Points2;
	double *Weights2, *t1, *t2;      // Weights - Quadrature Weights array ; t1  -- Quadrature point ( xi ) in ref coordinate ; t2 -  Quadrature Point ( eta ) in Ref Coordinate 
	bool Needs2ndDer[1]; 
	Needs2ndDer[0] = TRUE;
	double AbsDetjk[MaxN_PointsForNodal2D];
	double X[MaxN_PointsForNodal2D];
	double Y[MaxN_PointsForNodal2D];


	// FE Values Arrays 
	double** origvaluesD00;           // Shape function values at quadrature Points
	double** origvaluesD10;           // Shape Function Derivatives ( x ) at Quadrature Points
	double** origvaluesD01;           // Shape Function Derivatives ( y ) at Quadrature Points
    double** origvaluesD20;           // Shape Function 2nd Derivatives ( x ) at Quadrature Points
	double** origvaluesD02;           // Shape Function 2nd Derivatives ( y ) at Quadrature Points


	for (int cellId = 0; cellId < 1; cellId++)
	{
		TBaseCell* currentCell = coll->GetCell(cellId);
        // Get the "ID" of Finite Element for the given 2D Element ( Conforming/NonConforming-Order Finite Element : eg : it could be Conforming-2nd order Finite Element )
		FE2D elementId = Fespace->GetFE2D(cellId, currentCell);
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
		double val = 0;
		
		double* C_Array = FeVector_C->GetValues();
		int len         = FeVector_C->GetLength();

		

		double* C_Array_i = C_Array + i_index*len;

		
		// Save Values of C at all quadrature points for I component
		double C_i[N_Points2];
		double C_x_i[N_Points2];
		double C_y_i[N_Points2];

		for ( int quadPt = 0 ; quadPt < N_Points2; quadPt++) C_i[quadPt] = 0;
		for ( int quadPt = 0 ; quadPt < N_Points2; quadPt++) C_x_i[quadPt] = 0;
		for ( int quadPt = 0 ; quadPt < N_Points2; quadPt++) C_y_i[quadPt] = 0;
		
		
		for ( int quadPt = 0 ; quadPt < N_Points2; quadPt++)
		{
			for (int j = 0 ; j  < N_BaseFunct ;  j++)
			{
				int globDOF = DOF[j];
				C_i[quadPt] += origvaluesD00[quadPt][j] * C_Array_i[globDOF];
				C_x_i[quadPt] += origvaluesD10[quadPt][j] * C_Array_i[globDOF];
				C_y_i[quadPt] += origvaluesD01[quadPt][j] * C_Array_i[globDOF];
			}
		}
		
		double rhs[N_BaseFunct];
		for ( int j = 0 ; j < N_BaseFunct; j++) rhs[j] = 0;

		for ( int a = 0 ; a < N_S ; a++)
		{
			double* C_Array_a = C_Array + i_index*len;

			double C_a[N_Points2];
			double C_x_a[N_Points2];
			double C_y_a[N_Points2];

			for ( int quadPt = 0 ; quadPt < N_Points2; quadPt++) C_a[quadPt] = 0;
			for ( int quadPt = 0 ; quadPt < N_Points2; quadPt++) C_x_a[quadPt] = 0;
			for ( int quadPt = 0 ; quadPt < N_Points2; quadPt++) C_y_a[quadPt] = 0;

			// Obtain all values for C_a
			for ( int quadPt = 0 ; quadPt < N_Points2; quadPt++)
			{
				for (int j = 0 ; j  < N_BaseFunct ;  j++)
				{
					int globDOF = DOF[j];
					C_a[quadPt] += origvaluesD00[quadPt][j] * C_Array_a[globDOF];
					// C_x_a[quadPt] += origvaluesD10[quadPt][j] * C_Array_a[globDOF];
					// C_y_a[quadPt] += origvaluesD01[quadPt][j] * C_Array_a[globDOF];
				}
			}
			


			double val = 0;
			//Get Coefficients b1 and b2
			double *Param[MaxN_QuadPoints_2D];
			double **Coeffs =  new double*[MaxN_QuadPoints_2D];
			for ( int i = 0 ; i < MaxN_QuadPoints_2D; i++)
			{
				Coeffs[i] = new double[10]();
			}
			
			

			DO_Mode_Equation_Coefficients(N_Points2,X, Y,Param, Coeffs);
			
			// INner Quadrature Loop 
			for ( int quadPt = 0 ; quadPt < N_Points2; quadPt++)
			{
				double Mult = Weights2[quadPt] * AbsDetjk[quadPt];
				double* orgD00 = origvaluesD00[quadPt];
				double* orgD10 = origvaluesD10[quadPt];
				double* orgD01 = origvaluesD01[quadPt];
				double* orgD20 = origvaluesD20[quadPt];
				double* orgD02 = origvaluesD02[quadPt];

				double b1 = Coeffs[quadPt][1]; // b1
				double b2 = Coeffs[quadPt][2]; // b1

				
				for ( int quadPt_1 = 0 ; quadPt_1 < N_Points2; quadPt_1++)
				{
					val += (b1*C_x_i[quadPt_1] + b2*C_y_i[quadPt_1])*C_a[quadPt_1];
					val *= Mult;
				}

				val *= C_a[quadPt] ;  // This is Final "f"

				for ( int j = 0 ; j < N_BaseFunct ; j++)
				{
					rhs[j] += val * orgD00[j] * Mult;
					
				}

			}

		}

		for ( int j = 0 ; j < N_BaseFunct ; j++)
		{
			int GlobalDOF   				= DOF[j];
			GlobalRhs_mode[GlobalDOF]   	+= rhs[j];
		}
		// -- 
	}
}



void DO_CoEfficient(TFESpace2D *Fespace, TFEVectFunct2D *FeVector_C_Mode,TFEVectFunct2D* FEVector_Phi, int N_S, int i_index)
{

	int N_Cells = Fespace->GetN_Cells();
	TCollection* coll = Fespace->GetCollection();

	// Get the Global DOF arrays INdex from the FE Space. 
    int* GlobalNumbers = Fespace->GetGlobalNumbers();
	int* BeginIndex = Fespace->GetBeginIndex();

	// --- Quadrature Formula Arrays  ------------------//
	int N_Points2;
	double *Weights2, *t1, *t2;      // Weights - Quadrature Weights array ; t1  -- Quadrature point ( xi ) in ref coordinate ; t2 -  Quadrature Point ( eta ) in Ref Coordinate 
	bool Needs2ndDer[1]; 
	Needs2ndDer[0] = TRUE;
	double AbsDetjk[MaxN_PointsForNodal2D];
	double X[MaxN_PointsForNodal2D];
	double Y[MaxN_PointsForNodal2D];


	// FE Values Arrays 
	double** origvaluesD00;           // Shape function values at quadrature Points
	double** origvaluesD10;           // Shape Function Derivatives ( x ) at Quadrature Points
	double** origvaluesD01;           // Shape Function Derivatives ( y ) at Quadrature Points
    double** origvaluesD20;           // Shape Function 2nd Derivatives ( x ) at Quadrature Points
	double** origvaluesD02;           // Shape Function 2nd Derivatives ( y ) at Quadrature Points

	double val = 0;
	double* C_Array = FeVector_C_Mode->GetValues();
	int len         = FeVector_C_Mode->GetLength();


	double* Phi_Array = FEVector_Phi->GetValues();
	int lenPhi         = FEVector_Phi->GetLength();
	double* phi_New = new double[lenPhi]();
	double* C_Array_i = C_Array + i_index*len;
	double* phi_Array_i = Phi_Array + i_index*len;



	for (int cellId = 0; cellId < 1; cellId++)
	{
		TBaseCell* currentCell = coll->GetCell(cellId);
        // Get the "ID" of Finite Element for the given 2D Element ( Conforming/NonConforming-Order Finite Element : eg : it could be Conforming-2nd order Finite Element )
		FE2D elementId = Fespace->GetFE2D(cellId, currentCell);
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
		double val = 0;

		// Save Values of C at all quadrature points for I component
		double C_i[N_Points2];
		double C_x_i[N_Points2];
		double C_y_i[N_Points2];


		for ( int quadPt = 0 ; quadPt < N_Points2; quadPt++) C_i[quadPt] = 0;
		for ( int quadPt = 0 ; quadPt < N_Points2; quadPt++) C_x_i[quadPt] = 0;
		for ( int quadPt = 0 ; quadPt < N_Points2; quadPt++) C_y_i[quadPt] = 0;

		for ( int quadPt = 0 ; quadPt < N_Points2; quadPt++)
		{
			for (int j = 0 ; j  < N_BaseFunct ;  j++)
			{
				int globDOF = DOF[j];
				C_i[quadPt] += origvaluesD00[quadPt][j] * C_Array_i[globDOF];
				// C_x_i[quadPt] += origvaluesD10[quadPt][j] * C_Array_i[globDOF];
				// C_y_i[quadPt] += origvaluesD01[quadPt][j] * C_Array_i[globDOF];
			}
		}


		// for ( int j = 0 ; j < N_BaseFunct; j++) rhs[j] = 0;

		for ( int a = 0 ; a < N_S ; a++)
		{
			double* C_Array_a = C_Array + a*len;
			double* phi_Array_a = Phi_Array + a*len;

			double C_a[N_Points2];
			double C_x_a[N_Points2];
			double C_y_a[N_Points2];

			for ( int quadPt = 0 ; quadPt < N_Points2; quadPt++) C_a[quadPt] = 0;
			for ( int quadPt = 0 ; quadPt < N_Points2; quadPt++) C_x_a[quadPt] = 0;
			for ( int quadPt = 0 ; quadPt < N_Points2; quadPt++) C_y_a[quadPt] = 0;


			// Obtain all values for C_a
			for ( int quadPt = 0 ; quadPt < N_Points2; quadPt++)
			{
				for (int j = 0 ; j  < N_BaseFunct ;  j++)
				{
					int globDOF = DOF[j];
					C_a[quadPt] += origvaluesD00[quadPt][j] * C_Array_a[globDOF];
					C_x_a[quadPt] += origvaluesD10[quadPt][j] * C_Array_a[globDOF];
					C_y_a[quadPt] += origvaluesD01[quadPt][j] * C_Array_a[globDOF];
				}
			}
			


			//Get Coefficients b1 and b2
			double *Param[MaxN_QuadPoints_2D];
			double **Coeffs =  new double*[MaxN_QuadPoints_2D];
			for ( int i = 0 ; i < MaxN_QuadPoints_2D; i++)
			{
				Coeffs[i] = new double[10]();
			}
			

			DO_Mode_Equation_Coefficients(N_Points2,X, Y,Param, Coeffs);
			
			// INner Quadrature Loop 
			for ( int quadPt = 0 ; quadPt < N_Points2; quadPt++)
			{
				double Mult = Weights2[quadPt] * AbsDetjk[quadPt];
				double* orgD00 = origvaluesD00[quadPt];
				double* orgD10 = origvaluesD10[quadPt];
				double* orgD01 = origvaluesD01[quadPt];
				double* orgD20 = origvaluesD20[quadPt];
				double* orgD02 = origvaluesD02[quadPt];

				// double b1 = Coeffs[quadPt][1]; // b1
				// double b2 = Coeffs[quadPt][2]; // b1
				double b1 = 1;
				double b2 = 1;

				val += (b1*C_x_a[quadPt] + b2*C_y_a[quadPt])*C_i[quadPt];
				val *= Mult;

			}
	

			for ( int i = 0 ; i < lenPhi ; i++)
			{
				phi_New[i] += val*phi_Array_a[i];
			}

		}

			
		
	}

	double timeStep = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
	for ( int i = 0 ; i < lenPhi; i++)
	{
		phi_Array_i[i] -= timeStep*phi_New[i];
	}

	delete[] phi_New;

}

