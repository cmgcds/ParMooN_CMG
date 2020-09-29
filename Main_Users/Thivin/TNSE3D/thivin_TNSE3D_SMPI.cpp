#include <Domain.h>
#include <Database.h>
#include <FEDatabase3D.h>
#include <LinAlg.h>
#include <FESpace3D.h>
#include <SquareStructure3D.h>
#include <Structure3D.h>
#include <Output3D.h>
#include <LinAlg.h>
#include <MainUtilities.h>
#include <TimeDiscRout.h>
#include <tetgen.h>
#include <SystemTNSE3D_ALE.h>
#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "DeformMesh3D.h"
#include <FE3D_ALE.h>   
#include <omp.h>
#include "mkl.h"
#include "mkl_spblas.h"
#include "mkl_types.h"

#ifdef _SMPI
#include "mpi.h"
 #include "MumpsSolver.h"
#endif

#include "../Examples/TNSE_3D/climate_2.h"
#include "thivin.h"

using namespace std;

int main(int argc, char *argv[])
{
    omp_set_num_threads(4);

    TDatabase *Database = new TDatabase();
    TDomain *Domain = new TDomain(argv[1]);
    TFEDatabase3D *FEDatabase = new TFEDatabase3D();
    TCollection *coll;

    TFEFunction3D *p, *u1, *u2, *u3, **Pressure, *fefct[2];
    TOutput3D *Output;

    TSystemTNSE3D_ALE *SystemMatrix_ALE;
    TAuxParam3D *aux;

    MultiIndex3D AllDerivatives[4] = {D000, D100, D010, D001};

    const char vtkdir[] = "VTK";
    char *PsBaseName, *VtkBaseName, *GEO, *PRM, *SMESH;

    char Name[] = "name";
    char UString[] = "u";
    char PString[] = "p";
    char NameString[] = "VMS";
    char SubID[] = "";
    std::ostringstream os;
    double stime;

    os << " ";
    mkdir(vtkdir, 0777);

    Domain = new TDomain(argv[1]);



    //THIVIN - added SMPI implementation
    #if defined(_MPI) || defined(_SMPI)
    int rank, size, out_rank;
    int MaxCpV, MaxSubDomainPerDof;
    
    MPI_Comm Comm = MPI_COMM_WORLD;
    MPI_Init(&argc, &argv);

    double time1, time2;

    MPI_Comm_rank(Comm, &rank);
    MPI_Comm_size(Comm, &size);
    
    TDatabase::ParamDB->Comm = Comm;
    
    int Refine;
    int metisType[2] = {0,0};  
    #endif 


    #ifdef _SMPI
    TMumpsSolver * solmumps;
    #endif



    /** needed in the new implementation */
    if (TDatabase::ParamDB->FLOW_PROBLEM_TYPE == STOKES)
        TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE = 1;

    /* ------ meshgenerator */
    if (TDatabase::ParamDB->MESH_TYPE == 0)
    {
        Domain->Init(TDatabase::ParamDB->BNDFILE, TDatabase::ParamDB->GEOFILE); // ParMooN  build-in Geo mesh
    }

    else if (TDatabase::ParamDB->MESH_TYPE == 1){
        cout << " GMSH " <<endl;
        Domain->GmshGen(TDatabase::ParamDB->GEOFILE);
    }
    else
    {
        OutPut("Mesh Type not known, set MESH_TYPE correctly!!!" << endl);
    }

    // Refine the Mesh
    for (int i = 0; i < TDatabase::ParamDB->UNIFORM_STEPS; i++)
        Domain->RegRefineAll();
    
    coll=Domain->GetCollection(It_Finest, 0);
    
    // --------------------------- FE Space  for TNSE3D ------------------------------------ //
    // Declare FE3D Objects to be passed inside the FESpace  Constructor

    // Get the Shape of the Collection
    Shapes CellShape;
    CellShape = coll->GetCell(0)->GetType();

    FE3D *fe3d_velo = new FE3D[coll->GetN_Cells()];
    FE3D *fe3d_presure = new FE3D[coll->GetN_Cells()];
    FE3D *fe3d_mesh = new FE3D[coll->GetN_Cells()];

    if ( CellShape == Hexahedron || CellShape == Brick)
    {
        for (int i = 0; i < coll->GetN_Cells(); i++)
            fe3d_velo[i] = C_Q2_3D_H_M;

        
        for (int i = 0; i < coll->GetN_Cells(); i++)
            fe3d_presure[i] = C_Q1_3D_H_M;
        
        
        for (int i = 0; i < coll->GetN_Cells(); i++)
            fe3d_mesh[i] = C_Q1_3D_H_M;
    }

    else if (CellShape == Tetrahedron)
    {
        cout << " Tetra Headral Cells " <<endl;
        
        for (int i = 0; i < coll->GetN_Cells(); i++)
            fe3d_velo[i] = C_P2_3D_T_A;

        for (int i = 0; i < coll->GetN_Cells(); i++)
            fe3d_presure[i] = C_P1_3D_T_A;
        
        for (int i = 0; i < coll->GetN_Cells(); i++)
            fe3d_mesh[i] = C_P1_3D_T_A;

    }

    else
    {
        cout << " ERROR : only shapetype Hexaheadral and Tetraheadral cells are supported " <<endl;
        cout << " INFO  : Error in Main program " <<endl;
        cout << " Cell Shape : " << CellShape <<endl;
        exit(0);   // THIVIN - EXIT Statement.
    }
    
    TFESpace3D *Velocity_FeSpace = new TFESpace3D(coll, Name, Name, BoundCondition, fe3d_velo);
    TFESpace3D *Pressure_FeSpace = new TFESpace3D(coll, Name, Name, BoundCondition, fe3d_presure);
    cout << "FESPACE Velocity - Declared " << endl;
    cout << "FESPACE Pressure - Declared " << endl;



    int N_DOF = Velocity_FeSpace->GetN_DegreesOfFreedom();
    int N_Active = Velocity_FeSpace->GetActiveBound();



    // ------- FE Space  for Mesh --------------------- //

    // for (int i = 0; i < coll->GetN_Cells(); i++)
    //     fe3d_mesh[i] = C_Q1_3D_H_M;

    TFESpace3D *fespace_mesh = new TFESpace3D(coll, NameString, UString, Grid_BoundCondition, fe3d_mesh);

    int N_DOF_mesh = fespace_mesh->GetN_DegreesOfFreedom();
    int N_Active_mesh = fespace_mesh->GetActiveBound();


    // --------------------- END OF Declaration of FESpace for Mesh and Velocity ---------------------- //

    // ------- Declaration of Solution and RHS Arrays -------------------- //
    int N_U = Velocity_FeSpace->GetN_DegreesOfFreedom();
    int N_P = Pressure_FeSpace->GetN_DegreesOfFreedom();
    int N_TotalDOF = 3 * N_U + N_P;


    #ifdef _SMPI
    if(rank==0){
    #endif
    cout << " N_DOF :  " << N_DOF << endl;
    cout << " N_Active :  " << N_Active << endl;

    cout << " ------------------------ MESH PARAMETERS ----------------- "<< endl;
    cout << " N_DOF Mesh :  " << N_DOF_mesh << endl;
    cout << " N_Active Mesh:  " << N_Active_mesh << endl;
    cout << " N_U :  " << N_U << endl;
    cout << " N_P:  " << N_P << endl;
    cout << " N_TotalDOF:  " << N_TotalDOF << endl;
    cout << " ------------------------ VELOCITY  PARAMETERS ----------------- "<< endl;

    #ifdef _SMPI
    }
    #endif

    double* solution, *rhs, *oldrhs;

    solution = new double[N_TotalDOF]();
    rhs = new double[N_TotalDOF]();
    oldrhs = new double[N_TotalDOF]();

    double *meshVelocity = new double[3 * N_DOF_mesh]();
    double *meshVelocityOld = new double[3 * N_DOF_mesh]();
    // -------  - END - Declaration of Solution and RHS Arrays -------------------- //


    // --- FEVectFunction Declaration ---------------------- //
    TFEVectFunct3D *MeshVelo_FEvect = new TFEVectFunct3D(fespace_mesh, (char *)"MeshVelo", (char *)"MeshVelo", meshVelocity, N_DOF_mesh, 3);
    TFEVectFunct3D *Velocity_FEvect = new TFEVectFunct3D(Velocity_FeSpace, (char *)"Velocity", (char *)"Velocity", solution, N_DOF, 3);
    TFEFunction3D *Pressure_FEvect  = new TFEFunction3D(Pressure_FeSpace, (char *)"Pressure", (char *)"Pressure", solution + 3*N_U,N_P);
    
    // THIVIN - Creating a new velocity vect function for external boundary Condition
    double *externalVelocity    = new double [3* N_U]();
    TFEVectFunct3D *External_Velocity_FEvect = new TFEVectFunct3D(Velocity_FeSpace, (char *)"Velocity", (char *)"Velocity", externalVelocity, N_DOF, 3);    
    // --- - END -  FEVectFunction Declaration ---------------------- //
    // ------------------- VTK File - Declaration  ------------------------ //
   VtkBaseName = TDatabase::ParamDB->VTKBASENAME;    
   Output = new TOutput3D(2, 2, 1, 1, Domain);
   Output->AddFEVectFunct(Velocity_FEvect);
   Output->AddFEVectFunct(MeshVelo_FEvect);   
   Output->AddFEFunction(Pressure_FEvect);
    // os.seekp(std::ios::beg);
    int output_write_counter = 0;

    // ---------------- - END - VTK File - Declaration  ------------------------ //

  
    // ----------- ALE - Class Declaration --------------------------- //
    BoundCondFunct3D *ALEBoundaryConditions[3];
    BoundValueFunct3D *ALEBoundValues[3];

    ALEBoundaryConditions[0] = ALE_BoundCondition;
    ALEBoundaryConditions[1] = ALE_BoundCondition;
    ALEBoundaryConditions[2] = ALE_BoundCondition;
    ALEBoundValues[0] = ALE_BoundValue_X;
    ALEBoundValues[1] = ALE_BoundValue_Y;
    ALEBoundValues[2] = ALE_BoundValue_Z;

    FE3D_ALE *mesh_deform = new FE3D_ALE(coll, ALEBoundaryConditions, ALEBoundValues,fespace_mesh,MeshVelo_FEvect);
    TFEFunction3D *feFunction_u1,*feFunction_u2,*feFunction_u3 ;
    std::vector<int> BoundIds;
    std::vector<int> freeSlipBoundIds;

    #ifdef _SMPI
    if(rank == 0) {
    #endif // DEBUG
    // Pick the Free surface Bdid from example file
    getFreeSurfaceBoundaryIds(BoundIds);

    getFreeSlipBoundaryIds(freeSlipBoundIds);

    //Function to Pick Free Surface DOF
	mesh_deform->Pick_free_surface_DOFs(fespace_mesh, BoundIds,coll);

	// Function to pick up FreeSli Boundary DOF's
	mesh_deform->pickDOFsOfFreeSlipBoundaries(fespace_mesh,freeSlipBoundIds,BoundIds);

    // -------------- END OF ALE MESH MOVEMENT INIT  --------------------- //

 
    
    feFunction_u1 = Velocity_FEvect->GetComponent(0);
    feFunction_u2 = Velocity_FEvect->GetComponent(1);
    feFunction_u3 = Velocity_FEvect->GetComponent(2);  

    int N_Cells = coll->GetN_Cells();
    OutPut("N_Cells      : "<< setw(10) << N_Cells <<endl);
    OutPut("Dof velocity : "<< setw(10) << 3*N_U << endl);
    OutPut("Dof pressure : "<< setw(10) << N_P << endl);
    OutPut("Dof all      : "<< setw(10) << N_TotalDOF  << endl);  
    
    #ifdef _SMPI
      }
    #endif

    int NSEType = 4 ;// TDatabase::ParamDB->NSTYPE; 

    #ifdef _SMPI
    if(rank == 0) {
    #endif
    ///////////// -------------------- MOVE THE INITIAL BLOCK OF MESH ------------------------- //////////////////////////////
    BoundCondFunct3D *GridBoundaryConditions[3];
	BoundValueFunct3D *GridBoundValues[3];

	GridBoundaryConditions[0] = Grid_BoundCondition;
	GridBoundaryConditions[1] = Grid_BoundCondition;
	GridBoundaryConditions[2] = Grid_BoundCondition;
	GridBoundValues[0] = Grid_BoundValue_X;
	GridBoundValues[1] = Grid_BoundValue_Y;
	GridBoundValues[2] = Grid_BoundValue_Z;
	
	deformMesh3D* defMesh = new deformMesh3D(coll,GridBoundaryConditions,GridBoundValues);
	defMesh->moveMesh();


	os.seekp(std::ios::beg);
    os <<  "VTK/"<<VtkBaseName<<".0000" <<output_write_counter <<".vtk" << ends;
	Output->WriteVtk(os.str().c_str());
	output_write_counter =  output_write_counter + 1;

    cout << " -------------------------- Initial Mesh Movement Done --------------------------------     " <<endl;
    
    #ifdef _SMPI
    }
    #endif
    ///////////// ------- MOVE THE INITIAL BLOCK OF MESH ------------------------- //////////////////////////////

    // --------------- SYSTEM MATRIX DECLARATION ------------------------------------------------ //

    // double* defect = new double[3*N_U+N_P]();
     
    // Setting up Arrays for providing values to existing ParMooN dataStructure //
    TFESpace3D **Velocity_FeSpaceArray, **Pressure_FeSpaceArray,**gridFeSpaceArray;
    Velocity_FeSpaceArray       = new TFESpace3D*[1];
    Pressure_FeSpaceArray       = new TFESpace3D*[1];
    gridFeSpaceArray            = new TFESpace3D*[1];

    
    // FE- Spaces Array
    Velocity_FeSpaceArray[0]    = Velocity_FeSpace;
    Pressure_FeSpaceArray[0]    = Pressure_FeSpace;
    gridFeSpaceArray[0]         = fespace_mesh;

    // FE- Vect Function Array
    TFEVectFunct3D **Velocity_FEvectArray;
    TFEVectFunct3D **MeshTemp_FEvectArray;
    TFEFunction3D **Pressure_FEvectArray ;
    Velocity_FEvectArray        = new TFEVectFunct3D*[1];
    Pressure_FEvectArray        = new TFEFunction3D*[1];
    MeshTemp_FEvectArray        = new TFEVectFunct3D*[1];
    Velocity_FEvectArray[0]     = Velocity_FEvect;
    Pressure_FEvectArray[0]     = Pressure_FEvect;
    MeshTemp_FEvectArray[0]     = MeshVelo_FEvect;

    int tmp_mesh_DOF    = fespace_mesh->GetN_DegreesOfFreedom();

    // Mesh - temp Velocity
    MeshTemp_FEvectArray[0]     = MeshVelo_FEvect;
    


    double start_int = GetTime();
    int n_lvel = 1;   // THIVIN - Multigrid Levels

    double** solutionArray, **rhsArray;
    solutionArray           = new double*[1];
    rhsArray                = new double*[1];

    solutionArray[0]        = solution ;
    rhsArray[0]             = rhs;

    // Calculate Volume of the initial Mesh COnfirgiuration 
    double refVolume = coll->GetVolume();
    cout << " VOLUME ( REFERENCE CONFIGURATION ) : " << refVolume <<endl;

    // End Setting up array values for the Current Parmoon system Matrix

    #ifdef _SMPI
    if(rank == 0) {
    #endif

    SystemMatrix_ALE =  new TSystemTNSE3D_ALE(n_lvel,Velocity_FeSpaceArray,Pressure_FeSpaceArray,Velocity_FEvectArray
                                            ,Pressure_FEvectArray,solutionArray,rhsArray, TDatabase::ParamDB->DISCTYPE,
                                            NSEType, TDatabase::ParamDB->SOLVER_TYPE,gridFeSpaceArray,MeshTemp_FEvectArray,FALSE );

    cout << " Constructior Completed " <<endl;
    SystemMatrix_ALE->Init(LinCoeffs, BoundCondition, U1BoundValue, U2BoundValue, U3BoundValue,
                           Grid_BoundCondition, Grid_BoundValue_X,Grid_BoundValue_Y,Grid_BoundValue_Z, GridCoeffs);


    printf("SystemMatrix_ALE constructed\n");

    feFunction_u1->Interpolate(InitialU1);
    feFunction_u2->Interpolate(InitialU2);
    feFunction_u3->Interpolate(InitialU3);


    

    // Call the Mesh Movement Module to update the mesh Velocity, and copy the appropriate Dataset to the SystemMatrix 3d Class
    // 

    // ------ Calculate Mesh Velocity


	mesh_deform->get_surface_normals(MeshVelo_FEvect,Velocity_FEvect);    

	mesh_deform->get_surface_normals_slipBoundary(MeshVelo_FEvect);

	mesh_deform->get_surface_normals_slipBoundary_EdgeNodes(MeshVelo_FEvect);

    // memcpy(meshVelocityOld,meshVelocity,sizeof(double)*3*N_DOF_mesh);

    mesh_deform->Solve_mesh_velocity_interior(fespace_mesh,MeshVelo_FEvect);

    // Function which migrates the FE3D_ALE datastructure  of grid to the existing System TNSE3D Ale Approach
    SystemMatrix_ALE->AssignMeshVelocity_parameters(meshVelocity,0);

    SystemMatrix_ALE->pickDOFsOfFreeSlipBoundaries(coll,Velocity_FeSpace,freeSlipBoundIds,BoundIds);

    SystemMatrix_ALE->Assemble();  // "Thivin - need to verify " 

    // PRINT VTK 
    os.seekp(std::ios::beg);
    os <<  "VTK/"<<VtkBaseName<<".0000" <<output_write_counter <<".vtk" << ends;
	Output->WriteVtk(os.str().c_str());
	output_write_counter =  output_write_counter + 1;

    // exit(0);
    #ifdef _SMPI
    }
    #endif

    //Set the Initial Iteration Parameters
    int m = 0,checker;
    int N_SubSteps = GetN_SubSteps();
    double oldtau = 1.;
    double impuls_residual = 0,residual = 0;
    double end_time = TDatabase::TimeDB->ENDTIME; 
    double limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE;
    double Max_It = TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE;       
    double AllErrors [7];
    memset(AllErrors, 0, 7.*SizeOfDouble);


    #ifdef _SMPI
    if(rank !=0)
    { 
    solmumps = new TMumpsSolver(0,0,NULL,NULL,0); // 
    }
    if(rank==0)
    {
     stime = TDatabase::TimeDB->CURRENTTIME ;
    }


    MPI_Bcast(&stime,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

    TDatabase::TimeDB->CURRENTTIME = stime;
    MPI_Bcast(&end_time,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&N_SubSteps,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&Max_It,1,MPI_INT,0,MPI_COMM_WORLD);
    #endif

    while(TDatabase::TimeDB->CURRENTTIME< end_time)   // time cycle
    {  
        m++;
        TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;

        for(int l=0;l<N_SubSteps;l++)
        {
            #ifdef _SMPI 
                if(rank==0) { 
            #endif
            
            SetTimeDiscParameters(1);

            if(m==1)
            {
                OutPut( " -------   For Time Step : " << TDatabase::TimeDB->CURRENTTIME <<endl );
                OutPut("Theta1: " << TDatabase::TimeDB->THETA1<< "\t");
                OutPut("Theta2: " << TDatabase::TimeDB->THETA2<< "\t");
                OutPut("Theta3: " << TDatabase::TimeDB->THETA3<< "\t");
                OutPut("Theta4: " << TDatabase::TimeDB->THETA4<< endl);
                OutPut(endl);
            }
            double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
            TDatabase::TimeDB->CURRENTTIME += tau;

            OutPut(endl << "CURRENT Sub TIME: ");
            OutPut(TDatabase::TimeDB->CURRENTTIME << endl);  

            memcpy(oldrhs, rhs, N_TotalDOF*SizeOfDouble);  

            // SystemMatrix_ALE->imposeExternalBoundaryCondition(externalBoundaryParameters, External_Velocity_FEvect, Velocity_FEvect);
            
            // assemble  matrix for NSE time-dependent domain
            SystemMatrix_ALE->Assemble(); 


            //scale B matices and assemble NSE-rhs based on the \theta time stepping scheme 
            SystemMatrix_ALE->AssembleSystMat(tau, rhs, rhs, solution); 
            oldtau = tau;

            // calculate the residual
            SystemMatrix_ALE->GetResidual(solution, impuls_residual,residual);
   	

            OutPut(" nonlinear iteration step   0");
            OutPut(setw(14) << impuls_residual);
            OutPut(setw(14) << residual-impuls_residual);
            OutPut(setw(14) << sqrt(residual) << endl);  

            #ifdef _SMPI 
                }
            #endif

            double stime;    
            #ifdef _SMPI
            if(rank==0)
            {
                stime = TDatabase::TimeDB->CURRENTTIME ;
            }
            MPI_Bcast(&stime,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
            TDatabase::TimeDB->CURRENTTIME = stime;
            #endif
            
            for(int j=1;j<=Max_It;j++)
            {
                checker = 0;
                // SystemMatrix_ALE->imposeExternalBoundaryCondition( externalBoundaryParameters,  External_Velocity_FEvect, Velocity_FEvect);
                #ifdef _SMPI 
                if(rank==0) { 
                #endif
                SystemMatrix_ALE->get_surface_normals_slipBoundary(coll,Velocity_FEvect);   
                SystemMatrix_ALE->get_surface_normals_slipBoundary_EdgeNodes(coll,Velocity_FEvect );
                SystemMatrix_ALE->impose_FreeSlip_BoundaryCondition(rhs,Velocity_FEvect->GetLength(), N_Active);
                
                #ifdef _SMPI 
                } 
                #endif


                #ifdef _SMPI
                if(rank==0)
                #endif
                {
                 SystemMatrix_ALE->Solve(solution);
                }

                #ifdef _SMPI
                else
                {  
                    solmumps->FactorizeAndSolve(NULL,NULL);
                }
                #endif

                #ifdef _SMPI 
                if(rank==0) { 
                #endif

                if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
                    IntoL20FEFunction3D(solution+3*N_U, N_P, Pressure_FeSpace); 
                
                if(TDatabase::ParamDB->FLOW_PROBLEM_TYPE==STOKES)
                    break;
                
                //restore the mass matrix for the next nonlinear iteration      
                SystemMatrix_ALE->RestoreMassMatNonLinear();  

                //assemble the system matrix with given aux, sol and rhs 
                SystemMatrix_ALE->AssembleNonLinear();

                //assemble system mat, S = M + dt\theta_1*A
                SystemMatrix_ALE->AssembleSystMatNonLinear();  

                SystemMatrix_ALE->GetResidual(solution, impuls_residual,residual);


                mesh_deform->get_surface_normals(MeshVelo_FEvect,Velocity_FEvect);    
                mesh_deform->get_surface_normals_slipBoundary(MeshVelo_FEvect);
	            mesh_deform->get_surface_normals_slipBoundary_EdgeNodes(MeshVelo_FEvect);

                mesh_deform->Solve_mesh_velocity_interior(fespace_mesh,MeshVelo_FEvect);
                memcpy(meshVelocityOld,meshVelocity,sizeof(double)*3*N_DOF_mesh);

                OutPut(" nonlinear iteration step " << setw(3) << j<< setw(14) << impuls_residual);
                OutPut(setw(14) << residual-impuls_residual << setw(14) << sqrt(residual) << endl);

                #ifdef _SMPI 
                } 
                #endif

                if(sqrt(residual)<=limit)
                {
                    #ifdef _SMPI
                        checker = -1;
                    #else
                    break;
                    #endif
                }

                #ifdef _SMPI
                MPI_Bcast(&checker,1,MPI_INT,0,MPI_COMM_WORLD);
                if(checker ==-1)
                break;
                #endif	

            }

            #ifdef _SMPI 
            if(rank==0) { 
            #endif

            double timeStep = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH ;
            cout << "Current Time Step Length : " << timeStep <<endl;
            
            mesh_deform->move_mesh_ale(MeshVelo_FEvect,meshVelocityOld,timeStep);

            

            // all mat will be assembled
            SystemMatrix_ALE->SetSystMatAssembled(FALSE);

            #ifdef _SMPI 
                }
            #endif
        }

        #ifdef _SMPI 
            if(rank==0) { 
        #endif

        if(TDatabase::ParamDB->MEASURE_ERRORS)
        {
            SystemMatrix_ALE->MeasureTNSEErrors(ExactU1, ExactU2, ExactU3, ExactP, AllErrors);
            OutPut("L2(u): " <<   AllErrors[0] << "\t");
	        OutPut("H1-semi(u): " <<  AllErrors[1] << "\t");
	        OutPut("L2(p): " <<  AllErrors[2] << "\t");
	        OutPut("H1-semi(p): " <<  AllErrors[3]   << "\t"); 
	        OutPut(AllErrors[4] <<  " l_infty(L2(u)) " <<AllErrors[5] << endl);
	        OutPut(TDatabase::TimeDB->CURRENTTIME << " L2(0,t,L2)(u) : " <<   sqrt(AllErrors[6]) << endl); 
        }

        if(m==1 || m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)  
        {
            os.seekp(std::ios::beg);
            os <<  "VTK/"<<VtkBaseName<<".0000" <<output_write_counter <<".vtk" << ends;
            Output->WriteVtk(os.str().c_str());
            output_write_counter =  output_write_counter + 1;
        }

        // Calculate the relative Mass Conservation Value 
        double currentVol = coll->GetVolume();
        cout << " MASS CONSERVATION : " <<  ( 1  -  ( (refVolume   - currentVol ) / refVolume )  ) * 100 <<endl;

        #ifdef _SMPI 
        }
        #endif
    }

  



    return 0;
}
