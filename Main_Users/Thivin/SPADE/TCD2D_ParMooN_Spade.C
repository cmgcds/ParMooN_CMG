// =======================================================================
//
// Purpose:     main program for scalar equations with new kernels of ParMooN
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 08.08.2014

// =======================================================================

#include <Domain.h>
#include <Database.h>
#include <SystemTCD2D.h>
#include <SystemTCD2D_ALE.h>
#include <FEDatabase2D.h>
#include <FESpace2D.h>
#include <SquareStructure2D.h>
#include <Structure2D.h>
#include <QuadAffin.h>
#include <DirectSolver.h>
#include <Assemble2D.h>
#include <Output2D.h>
#include <LinAlg.h>
#include <CD2DErrorEstimator.h>
#include <MainUtilities.h>
#include <TimeDiscRout.h>
#include<fstream>
#include<string> 
#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <mkl.h>
#include<cmath>
#include<random>


// =======================================================================
// include current example
// =======================================================================
#include "../Examples/TCD_2D/exp.h"
// #include "../Examples/TCD_2D/SinCos1.h"
// #include "../Examples_All/TCD_2D/Time3.h"
// #include "../Examples/TCD_2D/exp_0.h"
//    #include "../Examples/TCD_2D/exp_2.h"
// #include "../Examples_All/TCD_2D/exp_1.h"
// #include "../Main_Users/Sashi/TCD_2D/Hemker.h"

int main(int argc, char *argv[])
{
    int i, j, l, m, N_SubSteps, ORDER, N_Cells, N_DOF, img = 1, N_G;
    int N_Active;

    double *sol, *rhs, *oldrhs, t1, t2, errors[5], Linfty;
    double tau, end_time, *defect, olderror, olderror1, hmin, hmax;

    bool UpdateStiffnessMat, UpdateRhs, ConvectionFirstTime;
    char *VtkBaseName;
    const char vtkdir[] = "VTK";

    TDomain *Domain;
    TDatabase *Database = new TDatabase();
    TFEDatabase2D *FEDatabase = new TFEDatabase2D();
    TCollection *coll;
    TFESpace2D *Scalar_FeSpace, *fesp[1];
    TFEFunction2D *Scalar_FeFunction;
    TOutput2D *Output;
    TSystemTCD2D *SystemMatrix;
    TAuxParam2D *aux;
    MultiIndex2D AllDerivatives[3] = {D00, D10, D01};

    std::ostringstream os;
    os << " ";

    // ======================================================================
    // set the database values and generate mesh
    // ======================================================================
    // set variables' value in TDatabase using argv[1] (*.dat file), and generate the MESH based
    Domain = new TDomain(argv[1]);

    if (TDatabase::ParamDB->PROBLEM_TYPE == 0)
        TDatabase::ParamDB->PROBLEM_TYPE = 2;
    OpenFiles();

    Database->WriteParamDB(argv[0]);
    Database->WriteTimeDB();
    ExampleFile();

    /* include the mesh from a mesh generator, for a standard mesh use the
   * build-in function. The GEOFILE describes the boundary of the domain. */
    if (TDatabase::ParamDB->MESH_TYPE == 0)
    {
        Domain->ReadGeo(TDatabase::ParamDB->GEOFILE);
        OutPut("PRM-GEO used for meshing !!!" << endl);
    } // ParMooN  build-in Geo mesh
    else if (TDatabase::ParamDB->MESH_TYPE == 1)
    {
        Domain->GmshGen(TDatabase::ParamDB->GEOFILE);
        OutPut("GMSH used for meshing !!!" << endl);
    }                                            //gmsh mesh
    else if (TDatabase::ParamDB->MESH_TYPE == 2) //triangle mesh
    {
        OutPut("Triangle.h used for meshing !!!" << endl);
        TriaReMeshGen(Domain);
    }
    else
    {
        OutPut("Mesh Type not known, set MESH_TYPE correctly!!!" << endl);
        exit(0);
    }

#if defined(__HEMKER__) || defined(__BEAM__)
    TriaReMeshGen(Domain);
    TDatabase::ParamDB->UNIFORM_STEPS = 0;
#endif

    // refine grid up to the coarsest level
    for (i = 0; i < TDatabase::ParamDB->UNIFORM_STEPS; i++)
        Domain->RegRefineAll();

    // write grid into an Postscript file
    if (TDatabase::ParamDB->WRITE_PS)
        Domain->PS("Domain.ps", It_Finest, 0);

    // create output directory, if not already existing
    if (TDatabase::ParamDB->WRITE_VTK)
        mkdir(vtkdir, 0777);

    //=========================================================================
    // construct all finite element spaces
    //=========================================================================
    ORDER = TDatabase::ParamDB->ANSATZ_ORDER;

    coll = Domain->GetCollection(It_Finest, 0);
    N_Cells = coll->GetN_Cells();
    OutPut("N_Cells (space) : " << N_Cells << endl);

    // fespaces for scalar equation
    Scalar_FeSpace = new TFESpace2D(coll, (char *)"fe space", (char *)"solution space",
                                    BoundCondition, 1, NULL);


    N_DOF = Scalar_FeSpace->GetN_DegreesOfFreedom();
    N_Active = Scalar_FeSpace->GetActiveBound();
    OutPut("dof all      : " << setw(10) << N_DOF << endl);
    OutPut("dof active   : " << setw(10) << N_Active << endl);

    //======================================================================
    // construct all finite element functions
    //======================================================================
    sol = new double[N_DOF];
    rhs = new double[N_DOF];
    oldrhs = new double[N_DOF];

    memset(sol, 0, N_DOF * SizeOfDouble);
    memset(rhs, 0, N_DOF * SizeOfDouble);

    Scalar_FeFunction = new TFEFunction2D(Scalar_FeSpace, (char *)"sol", (char *)"sol", sol, N_DOF);

    //     // --- Get the list from the co-ordinates
    double *org_x_coord     = new double[N_DOF];
    double *org_y_coord     = new double[N_DOF];
    double *x_coord         = new double[N_DOF];
    double *y_coord         = new double[N_DOF];
    int *mappingArray       = new int[N_DOF];


    i=0;
    int N = pow(2,TDatabase::ParamDB->UNIFORM_STEPS  ) + 1;
    for ( int i = 0 ; i < N_DOF; i++)
    {
        int local_i = i/N;
        int local_j = i%N;
       
        //  i=0, 0, i=1, 0.3, i=2, 0.6 i=3 = 0.99
        x_coord[i] =  double(1.0/(N-1)) * local_i;
        y_coord[i] =  double(1.0/(N-1)) * local_j;
    }


    cout << " End File Read" <<endl;

    Scalar_FeSpace->GetDOFPosition(org_x_coord,org_y_coord);

    for ( int i=0 ; i < N_DOF; i++)
    {  
        // get the generated Value
        double xx = x_coord[i];
        double yy = y_coord[i];
        bool foundFlag = false;

        for ( int j=0 ; j<N_DOF;j++)
        {  
            if(abs(xx - org_x_coord[j]) < 1e-10 &&  abs(yy - org_y_coord[j]) < 1e-10 )
            {
                mappingArray[i] = j;
                foundFlag = true;
            }
        }

        if(!foundFlag) cerr<< " DOF NOT FOUND FOR " << i << " position : " << setw(8) << org_x_coord[i]<<setw(8) <<org_y_coord[i] <<endl;
    }
     

    // for (int i = 0; i < N_DOF; i++)
    // {
    //     cout << i << setw(8)<< mappingArray[i]<<endl;
    // }
   
   
    // ---------------------- GENERATE DATA FOR ALL REALISATIONS ------------------------------------ //
    // int N = pow(2,TDatabase::ParamDB->UNIFORM_STEPS  ) + 1;
    int N_Realisations = 4000;
    // int N_DOF =  N * N;
    double* x  =  new double[N_DOF];
    double* y  =  new double[N_DOF];

    for ( int i = 0 ; i < N_DOF; i++ )
    {
        int local_i = i/N + 1;
        int local_j = i%N + 1;

        x[((local_i-1) * N) + (local_j-1)] = local_j;
        y[((local_i-1) * N) + (local_j-1)] = local_i;
    }

   
    double *C = new double[N_DOF*N_DOF];  //MATRIX

    double norm = 0;
    for( int i =0  ; i < N_DOF ; i++ )
    {
        int local_i = i/N;
        int local_j = i%N;

        double actual_x = x[i];
        double actual_y = y[i];

        for ( int j=0 ; j < N_DOF ; j++)
        {
            double local_x = x[j];
            double local_y = y[j];

            double r = sqrt( pow((actual_x - local_x),2 ) + pow((actual_y - local_y),2 ));
            C[j*N_DOF + i] = exp(-1 * r );
            norm += C[j*N + i]*C[j*N + i];

        }

    }

    // Declare SVD parameters
    MKL_INT m1 = N_DOF, n = N_DOF, lda = N_DOF, ldu = N_DOF, ldvt = N_DOF, info;
    double superb[std::min(N_DOF,N_DOF)-1];

    double* S = new double[N_DOF];
    double* U = new double[N_DOF*N_DOF];
    double* Vt = new double[N_DOF*N_DOF];

    info = LAPACKE_dgesvd( LAPACK_ROW_MAJOR, 'A', 'A', m1, n, C, lda,
                        S, U, ldu, Vt, ldvt, superb );

    cout << endl <<endl;

    if( info > 0 ) {
                printf( "The algorithm computing SVD failed to converge.\n" );
                exit( 1 );
    }
   
    int energyVal = 0;

    double sumSingularVal = 0;
    for( int i=0;i<N_DOF;i++) sumSingularVal += S[i];

    double val = 0;
    for( energyVal =0 ; energyVal< N_DOF; energyVal++)
    {
        val += S[energyVal];
        if(val/sumSingularVal > 0.8) break;
    }

    cout << energyVal <<endl;

    int modDim = energyVal+1;
       
    double* Ut = new double[N_DOF*modDim]();
    double* Z  = new double[N_Realisations*modDim]();

    double* SolutionVector = new double[N_DOF * N_Realisations]();

    // -------------- Generate Random Number Based on Normal Distribution -------------------------//
    int k=0;
    int skip = N_DOF - modDim;
    int count =0;
    for ( int i = 0 ; i < N_DOF*N_DOF ; i++ )
    {  
        // cout << "i val " << i <<endl;
        if(count < modDim )
        {
            Ut[k] =  U[i];

            count++;
            k++;
        }
        else
        {
            i += skip;
            count = 0;
            i--;
        }
       
    }


    for( int k = 0 ; k < modDim ; k++)
    {  
        std::random_device rd{};
        std::mt19937 gen{rd()};
        std::normal_distribution<> d{0,1};

        double* norm1 = new double[N_Realisations];

        for(int n=0; n<N_Realisations; ++n) {
            Z[k*N_Realisations + n] =  S[k] * d(gen);
        }
    }



    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,N_DOF,N_Realisations, modDim , 1.0, Ut,modDim,Z,N_Realisations,0.0,SolutionVector,N_Realisations);

    // printMatrix(SolutionVector, N_DOF,N_Realisations);

    mkl_dimatcopy('R','T', N_DOF,N_Realisations,1.0,SolutionVector,N_DOF,N_Realisations);


    cout << " REALISATIONS COMPUTED " <<endl;

    /////////////////////////////////////// -------- END OF REALISATION DATA SETS ------------ ////////////////////////////////////////////////////////////////


    // -------- Output parameters------------//
    VtkBaseName = TDatabase::ParamDB->VTKBASENAME;
    Output = new TOutput2D(2, 2, 1, 1, Domain);
    Output->AddFEFunction(Scalar_FeFunction);



    std::ofstream fileout;
    std::string name = "Sol1" ;
    fileout.open(name);

    //======================================================================
    // SystemMatrix construction and solution
    //======================================================================
    // Disc type: GALERKIN (or) SDFEM  (or) UPWIND (or) SUPG (or) LOCAL_PROJECTION
    // Solver: AMG_SOLVE (or) GMG  (or) DIRECT
    SystemMatrix = new TSystemTCD2D(Scalar_FeSpace, GALERKIN, DIRECT);

    // initilize the system matrix with the functions defined in Example file
    SystemMatrix->Init(BilinearCoeffs, BoundCondition, BoundValue);

    for ( int RealNo=0 ; RealNo < 1; RealNo++)
    {
        // cout << " ============================================================================================================= " <<endl;
        cout << RealNo  <<endl;
        // cout << " ============================================================================================================= " <<endl;

        // Scalar_FeFunction->Interpolate(InitialCondition);
        std::string str= std::to_string(RealNo);
        for ( int i=0 ; i < N_DOF; i++)
            sol[mappingArray[i]] = SolutionVector[RealNo*N_DOF + i];
        

        if (TDatabase::ParamDB->WRITE_VTK)
        {
            os.seekp(std::ios::beg);
            if (img < 10)
                os << "VTK/" << VtkBaseName << ".0000" << img << ".vtk" << ends;
            else if (img < 100)
                os << "VTK/" << VtkBaseName << ".000" << img << ".vtk" << ends;
            else if (img < 1000)
                os << "VTK/" << VtkBaseName << ".00" << img << ".vtk" << ends;
            else if (img < 10000)
                os << "VTK/" << VtkBaseName << ".0" << img << ".vtk" << ends;
            else
                os << "VTK/" << VtkBaseName << "." << img << ".vtk" << ends;
            Output->WriteVtk(os.str().c_str());
            img++;
        }
        
        double origSolVal = sol[500];
        double origSolVal1 = sol[200];
        double origSolVal2 = sol[900];
        // cout << " Solution Norm : " << Ddot(N_DOF,sol,sol) <<endl;



        // assemble the system matrix with given aux, sol and rhs
        // aux is used to pass  addition fe functions (eg. mesh velocity) that is nedded for assembling,
        // otherwise, just pass with NULL
        SystemMatrix->AssembleMRhs(NULL, sol, rhs);

        //======================================================================
        // time disc loop
        //======================================================================
        // parameters for time stepping scheme
        m = 0;
        N_SubSteps = GetN_SubSteps();
        end_time = TDatabase::TimeDB->ENDTIME;

        UpdateStiffnessMat = TRUE; //check BilinearCoeffs in example file
        UpdateRhs = TRUE;           //check BilinearCoeffs in example file
        ConvectionFirstTime = TRUE;

        // time loop starts
        while (TDatabase::TimeDB->CURRENTTIME < end_time)
        {
            m++;
            TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;

            for (l = 0; l < N_SubSteps; l++) // sub steps of fractional step theta
            {
                SetTimeDiscParameters(1);

                if (m == 1)
                {
                    // OutPut("Theta1: " << TDatabase::TimeDB->THETA1 << endl);
                    // OutPut("Theta2: " << TDatabase::TimeDB->THETA2 << endl);
                    // OutPut("Theta3: " << TDatabase::TimeDB->THETA3 << endl);
                    // OutPut("Theta4: " << TDatabase::TimeDB->THETA4 << endl);
                }

                tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
                TDatabase::TimeDB->CURRENTTIME += tau;

                // OutPut(endl<< "CURRENT TIME: ");
                // OutPut(TDatabase::TimeDB->CURRENTTIME << endl);

                //copy rhs to oldrhs
                memcpy(oldrhs, rhs, N_DOF * SizeOfDouble);

                // unless the stiffness matrix or rhs change in time, it is enough to
                // assemble only once at the begning
                if (UpdateStiffnessMat || UpdateRhs || ConvectionFirstTime)
                {
                    SystemMatrix->AssembleARhs(NULL, sol, rhs);

                    // M:= M + (tau*THETA1)*A
                    // rhs: =(tau*THETA4)*rhs +(tau*THETA3)*oldrhs +[M-(tau*THETA2)A]*oldsol
                    // note! sol contains only the previous time step value, so just pass
                    // sol for oldsol
                    SystemMatrix->AssembleSystMat(oldrhs, sol, rhs, sol);
                    ConvectionFirstTime = FALSE;
                }

                // solve the system matrix
                SystemMatrix->Solve(sol, rhs);

                // restore the mass matrix for the next time step
                // unless the stiffness matrix or rhs change in time, it is not necessary to assemble the system matrix in every time step
                if (UpdateStiffnessMat || UpdateRhs)
                {
                    SystemMatrix->RestoreMassMat();
                }

                      if (TDatabase::ParamDB->WRITE_VTK)
        {
            os.seekp(std::ios::beg);
            if (img < 10)
                os << "VTK/" << VtkBaseName << ".0000" << img << ".vtk" << ends;
            else if (img < 100)
                os << "VTK/" << VtkBaseName << ".000" << img << ".vtk" << ends;
            else if (img < 1000)
                os << "VTK/" << VtkBaseName << ".00" << img << ".vtk" << ends;
            else if (img < 10000)
                os << "VTK/" << VtkBaseName << ".0" << img << ".vtk" << ends;
            else
                os << "VTK/" << VtkBaseName << "." << img << ".vtk" << ends;
            Output->WriteVtk(os.str().c_str());
            img++;
        }

            } // for(l=0;l<N_SubSteps;l++)

        } // while(TDatabase::TimeDB->CURRENTTIME< end_time)

        //======================================================================
        // produce final outout
        //======================================================================

        // if (TDatabase::ParamDB->WRITE_VTK)
        // {
        //     os.seekp(std::ios::beg);
        //     if (img < 10)
        //         os << "VTK/" << VtkBaseName << ".0000" << img << ".vtk" << ends;
        //     else if (img < 100)
        //         os << "VTK/" << VtkBaseName << ".000" << img << ".vtk" << ends;
        //     else if (img < 1000)
        //         os << "VTK/" << VtkBaseName << ".00" << img << ".vtk" << ends;
        //     else if (img < 10000)
        //         os << "VTK/" << VtkBaseName << ".0" << img << ".vtk" << ends;
        //     else
        //         os << "VTK/" << VtkBaseName << "." << img << ".vtk" << ends;
        //     Output->WriteVtk(os.str().c_str());
        //     img++;
        // }

        // cout << " Solution Norm After: " << Ddot(N_DOF,sol,sol) <<endl;

        fileout << sol[500]  << "," << origSolVal<< "," << sol[200] << "," <<origSolVal1 <<"," << sol[900] << "," <<origSolVal2 << endl;
        

        // set Current Time as Zero
        TDatabase::TimeDB->CURRENTTIME = 0;
        // delete SystemMatrix;

    }
    fileout.close();

    CloseFiles();

    return 0;
} // end main
