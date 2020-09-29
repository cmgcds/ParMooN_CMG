// =======================================================================
//
// Purpose:     main program for solving a stationary TNSE equation in ParMooN
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 17.12.2015

// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <SystemTCD3D.h>
#include <FEDatabase3D.h>
#include <LinAlg.h>
#include <FESpace3D.h>
#include <SystemTNSE3D.h>
#include <SquareStructure3D.h>
#include <Structure3D.h>
#include <Output3D.h>
#include <LinAlg.h>
#include <MainUtilities.h>
#include <TimeDiscRout.h>
// #include <TCD3D.h>

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>

#ifdef _SMPI
#include "mpi.h"
 #include "MumpsSolver.h"
#endif

#ifdef _MPI
#include "mpi.h"
#include <MeshPartition.h>
#endif
#include <fstream> 

double bound = 0;
double timeC = 0;

#define AMG 0
#define GMG 1
#define DIRECT 2
// =======================================================================
// include current example
// =======================================================================
// #include "../Examples/TNSE_3D/ChannelObstacle3D_slip_volker.h"
#include "CD_NSE_3D.h"

void printall_array ( double * Arr1,double * Arr2, int SizeOfArr )
{
	std::ofstream myfile;
	myfile.open("entries/Arr1.txt");
	for(int ii =0; ii< SizeOfArr ; ii++ )
	myfile << " " << Arr1[ii] << endl;
	myfile.close();
	
	
	myfile.open("entries/Arr2.txt");
	for(int ii =0; ii< SizeOfArr ; ii++ )
	myfile << " " << Arr2[ii] << endl;
	myfile.close();
}
     
void norm_array ( double * Arr1,double * Arr2, int SizeOfArr )
{	
	double sum_Arr1=0;
	for(int ii =0; ii< SizeOfArr ; ii++ )
	{sum_Arr1= sum_Arr1+ (Arr1[ii]*Arr1[ii]);}
	cout<<"sum_Arr1::"<<sum_Arr1<<endl;	
	
        double sum_Arr2=0;
	for(int ii =0; ii< SizeOfArr ; ii++ )
	{sum_Arr2= sum_Arr2+ (Arr2[ii]*Arr2[ii]);}
	cout<<"sum_Arr2::"<<sum_Arr2<<endl;
}
// =======================================================================
// main program
// =======================================================================
int main(int argc, char* argv[])
{
  // ======================================================================
  //  declaration of variables
  // ======================================================================
    int i, j, l, m, N_Cells, N_U, N_P, N_TotalDOF, img=1, pressure_space_code;
    int Max_It, NSEType, velocity_space_code;
    int LEVELS, mg_level, mg_type;
    int N_SubSteps, N_L;
    int ORDER, N_T;
  
    double *sol_NSE, *rhs_NSE, *sol_CD, *rhs_CD, *oldrhs_NSE, *oldrhs_CD, t1, t2, errors[4], residual, impuls_residual;
    double *oldsol;
    double **Sol_array_NSE, **Rhs_array_NSE, **Sol_array_CD, **Rhs_array_CD;
    double limit, AllErrors[7];
    double tau, oldtau, end_time;
  
    double start_time, stop_time, start_assembling_solving=0, end_assembling_solving=0,
            total_assembling_solving=0, start_int=0, end_int=0, total_int=0;
	 
    double Cd, Cl;	 
    int checker = 1;
    TDomain *Domain;
    TDatabase *Database = new TDatabase();
  
    int profiling;
    
    bool UpdateStiffnessMat, UpdateRhs,  ConvectionFirstTime;    
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
  
    TFEDatabase3D *FEDatabase = new TFEDatabase3D(); 
    TCollection *coll, *mortarcoll = NULL;

    TFESpace3D *velocity_space, *pressure_space, **Velocity_FeSpace, **Pressure_FeSpace, *fesp[2];
    TFEVectFunct3D **Velocity, *u;
    TFEFunction3D *p, *u1, *u2, *u3, **Pressure, *fefct[3];
    TSystemTNSE3D *SystemMatrix_NSE;
  
    TFESpace3D **Temperature_FeSpace;
    TFEFunction3D *Theta, **Temperature;
    TSystemTCD3D *SystemMatrix_CD;
    
    TOutput3D *Output;
    TAuxParam3D *aux_CD, *aux_NSE;
    MultiIndex3D AllDerivatives[4] = { D000, D100, D010, D001 };
  
    TFESpace3D *projection_space;
    TFESpace3D **Projection_FeSpace;

    const char vtkdir[] = "VTK"; 
    char *PsBaseName, *VtkBaseName;
 
    char Name[] = "name";
    char Description[] = "description";
    char NameString[] = "VMS";
    char UString[] = "U";
    char PString[] = "P";
    char TString[] = "T";
    char SubID[] = "";
    std::ostringstream os;    
    double stime;
#ifdef _SMPI
    TMumpsSolver * solmumps;
#endif
  
#ifdef _SMPI    
    if(rank==0)  
#endif      
    {
        os << " ";   
        mkdir(vtkdir, 0777);
    
        // ======================================================================
        // set the database values and generate mesh
        // ======================================================================    
        /** set variables' value in TDatabase using argv[1] (*.dat file), and generate the MESH based */
        Domain = new TDomain(argv[1]);
        profiling = TDatabase::ParamDB->timeprofiling;

//         if(profiling)
//         {
// #ifdef _MPI
//             start_time = MPI_Wtime();
// #else
//             start_time = GetTime();
// #endif
//         }

        OpenFiles();
        OutFile.setf(std::ios::scientific);
   
        Database->CheckParameterConsistencyNSE();

#ifdef _MPI
        out_rank=TDatabase::ParamDB->Par_P0;
        //out_rank = 0;
        if(rank == out_rank)
#endif
        {
            Database->WriteParamDB(argv[0]);
            //Database->WriteTimeDB();
            ExampleFile();
        }
   
        /** needed in the new implementation */
        if(TDatabase::ParamDB->FLOW_PROBLEM_TYPE == STOKES)
        {  TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE = 1; }
        /* meshgenerator */
  
        Domain->GmshGen(TDatabase::ParamDB->GEOFILE);//gmsh mesh 
        LEVELS = TDatabase::ParamDB->LEVELS;
  
        if(TDatabase::ParamDB->SOLVER_TYPE==DIRECT)
        {
            TDatabase::ParamDB->UNIFORM_STEPS += (LEVELS-1);
            LEVELS = 1;
        } 
        // refine grid up to the coarsest level
        for(i=0;i<TDatabase::ParamDB->UNIFORM_STEPS;i++)
            Domain->RegRefineAll();
        
#ifdef _MPI
        Domain->GenerateEdgeInfo();
  
//         if(profiling)  t1 = MPI_Wtime();
  
        if(rank == 0)
        {
            printf("\n----------------------------------------------------------------------------------------\n");
            printf("metis type set to %d\n",TDatabase::ParamDB->Par_P2);
            printf("----------------------------------------------------------------------------------------\n\n");
        }
        //this loop checks if number of cells are sufficient in the coarsest level, such that each 
        //rank get some own cells to work on
        //it does so by changing the metis type first, if not possible then refine and start again
        do
        {
            metisType[TDatabase::ParamDB->Par_P2] = 1;
            Refine = Partition_Mesh3D(Comm, Domain, MaxCpV);	//MaxCpV=maximum cell per vertex
    
            if(metisType[0]*metisType[1] == 1 && Refine)
            {
                metisType[0] = 0;      metisType[1] = 0;
                TDatabase::ParamDB->Par_P2 = 0;
                if(rank == 0)
                {
                    printf("\n----------------------------------------------------------------------------------------\n");
                    printf("Warning :: both metisType used. Now refining the mesh by one step \n");
                    printf("metis type set to 0\n");
                    printf("----------------------------------------------------------------------------------------\n\n");
                }
                Domain->RegRefineAll();
                Domain->GenerateEdgeInfo();
                TDatabase::ParamDB->UNIFORM_STEPS +=1;
            }
        }while(Refine);
  
//         if(profiling)  t2 = MPI_Wtime(); 
//   
//         if(profiling)
//         {
//             time2 = t2-t1;
//             MPI_Reduce(&time2, &time1, 1, MPI_DOUBLE, MPI_MAX, out_rank, Comm);
//             if(rank == out_rank)
//                 printf("Time taken for Domain Decomposition is %e\n", time1);
//         }

        Domain->GenerateEdgeInfo();
        MaxSubDomainPerDof = MIN(MaxCpV, size);
        TDatabase::ParamDB->WRITE_PS = 0;
#endif 
    
        if(TDatabase::ParamDB->WRITE_PS)
        {
            // write grid into an Postscript file
            os.seekp(std::ios::beg);
            os << "Domain" << ".ps" << ends;
            Domain->PS(os.str().c_str(),It_Finest,0);
        }

//=========================================================================
// set data for multigrid
//=========================================================================  

        // set type of multilevel
        mg_type = TDatabase::ParamDB->SC_MG_TYPE_SADDLE;
 
        if(TDatabase::ParamDB->SOLVER_TYPE==AMG_SOLVE|| TDatabase::ParamDB->SOLVER_TYPE==DIRECT)
        { 
            mg_type=0;
            TDatabase::ParamDB->SC_MG_TYPE_SADDLE = mg_type;
            cout << mg_type << "-------------------------------" << endl;
        }
  
        if(mg_type==1)
        { 
            mg_level =  LEVELS + 1;
            ORDER = -1;
        }
        else
        { 
            mg_level = LEVELS; 
            ORDER = TDatabase::ParamDB->ANSATZ_ORDER;      
        }
   
        if(TDatabase::ParamDB->SOLVER_TYPE==GMG)
#ifdef _MPI  
        if(rank == out_rank)
#endif   
        {   
            OutPut("=======================================================" << endl);
            OutPut("======           GEOMETRY  LEVEL ");
            OutPut(LEVELS << "              ======" << endl);
            OutPut("======           MULTIGRID LEVEL ");
            OutPut(mg_level << "              ======" << endl);
            OutPut("=======================================================" << endl);   
        }   

        Velocity_FeSpace = new TFESpace3D*[mg_level];  
        Pressure_FeSpace = new TFESpace3D*[mg_level];  
        Temperature_FeSpace   = new TFESpace3D*[mg_level];
  
#ifdef __PRIVATE__ 
        Projection_FeSpace = new TFESpace3D*[mg_level];
#endif
  
        Velocity = new TFEVectFunct3D*[mg_level];
        Pressure = new TFEFunction3D*[mg_level];
        Temperature = new TFEFunction3D*[mg_level];
  
        Sol_array_NSE = new double*[mg_level];
        Rhs_array_NSE = new double*[mg_level];
        Sol_array_CD = new double*[mg_level];
        Rhs_array_CD = new double*[mg_level];    

//=========================================================================
// loop over all levels (not a multigrid level but for convergence study)  
//=========================================================================
        for(i=0;i<LEVELS;i++)
        {   
            if(i)
            { Domain->RegRefineAll(); }
#ifdef _MPI
            if(rank == out_rank)
                printf("Level :: %d\n\n",i);
            if(i)
            {
                Domain->GenerateEdgeInfo();
                Domain_Crop(Comm, Domain);       // removing unwanted cells in the hallo after refinement 
            }
#endif
            coll=Domain->GetCollection(It_Finest, 0);
//=========================================================================
// construct all NSE finite element spaces
//=========================================================================  
            if(mg_type==1) // lower order FE on coarse grids 
            {
                Velocity_FeSpace[i] = new TFESpace3D(coll, Name, UString, BoundCondition_NSE, Non_USpace,1);
                Pressure_FeSpace[i] = new TFESpace3D(coll, Name, PString, BoundCondition_NSE, DiscP_PSpace,0);    
       
                if(i==LEVELS-1) // higher order on fine level
                {
                    GetVelocityAndPressureSpace3D(coll, BoundCondition_NSE, velocity_space,
                                                pressure_space, &pressure_space_code,
                                                TDatabase::ParamDB->VELOCITY_SPACE,
                                                TDatabase::ParamDB->PRESSURE_SPACE); 
                    Velocity_FeSpace[i+1] = velocity_space;
                    Pressure_FeSpace[i+1] = pressure_space;

                    // defaulty inf-sup pressure space will be selected based on the velocity space, so update it in database
                    TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE = pressure_space_code;
                    velocity_space_code = TDatabase::ParamDB->VELOCITY_SPACE;
#ifdef _MPI
                    Velocity_FeSpace[LEVELS]->SetMaxSubDomainPerDof(MaxSubDomainPerDof);
                    Pressure_FeSpace[LEVELS]->SetMaxSubDomainPerDof(MaxSubDomainPerDof);
#endif
                }
            }
            else
            {
                GetVelocityAndPressureSpace3D(coll, BoundCondition_NSE, velocity_space,
                                            pressure_space, &pressure_space_code,
                                            TDatabase::ParamDB->VELOCITY_SPACE,
                                            TDatabase::ParamDB->PRESSURE_SPACE); 

                Velocity_FeSpace[i] = velocity_space;
                Pressure_FeSpace[i] = pressure_space;

                TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE = pressure_space_code;
                velocity_space_code = TDatabase::ParamDB->VELOCITY_SPACE;
            }
#ifdef __PRIVATE__ 
            if(TDatabase::ParamDB->DISCTYPE == VMS_PROJECTION)
            {
                if (TDatabase::ParamDB->VMS_LARGE_VELOCITY_SPACE==0)
                    projection_space = new TFESpace3D(coll,NameString, UString, BoundCondition_NSE, DiscP_PSpace, 0);
                else
                    projection_space = new TFESpace3D(coll,NameString, UString, BoundCondition_NSE, DiscP_PSpace, 1);
          
                Projection_FeSpace[i] = projection_space;
                N_L = Projection_FeSpace[i]->GetN_DegreesOfFreedom();
                OutPut("Dof Projection : " << setw(10) << N_L << endl);     
            }
#endif
#ifdef _MPI
            Velocity_FeSpace[i]->SetMaxSubDomainPerDof(MaxSubDomainPerDof);
            Pressure_FeSpace[i]->SetMaxSubDomainPerDof(MaxSubDomainPerDof);
#endif
//=========================================================================
// construct all CD finite element spaces
//========================================================================= 
            Temperature_FeSpace[i] =  new TFESpace3D(coll, Name, Description, BoundCondition_CD, ORDER);     

#ifdef _MPI
            Temperature_FeSpace[i]->SetMaxSubDomainPerDof(MaxSubDomainPerDof);
#endif
     
        //multilevel multigrid disc
            if(i==LEVELS-1 && mg_type==1) 
            {
                ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
                Temperature_FeSpace[mg_level-1] =  new TFESpace3D(coll, Name, Description, 
                                                            BoundCondition_CD, ORDER);
#ifdef _MPI
                Temperature_FeSpace[mg_level-1]->SetMaxSubDomainPerDof(MaxSubDomainPerDof);
#endif
            } //  if(i==LEVELS-1 && i!=mg_level-1)
//======================================================================
// construct all NSE finite element functions
//======================================================================   
            N_U = Velocity_FeSpace[i]->GetN_DegreesOfFreedom();
            N_P = Pressure_FeSpace[i]->GetN_DegreesOfFreedom();
            N_TotalDOF = 3*N_U + N_P;    

            sol_NSE = new double[N_TotalDOF];
            memset(sol_NSE, 0, N_TotalDOF*SizeOfDouble);     
            Sol_array_NSE[i] = sol_NSE;
     
            rhs_NSE = new double[N_TotalDOF];
            memset(rhs_NSE, 0, N_TotalDOF*SizeOfDouble);     
            Rhs_array_NSE[i] = rhs_NSE;

            u = new TFEVectFunct3D(Velocity_FeSpace[i], UString,  UString,  sol_NSE, N_U, 3);
            Velocity[i] = u;
            p = new TFEFunction3D(Pressure_FeSpace[i], PString,  PString,  sol_NSE+3*N_U, N_P);    
            Pressure[i] = p;
 
            if(i==LEVELS-1 && mg_type==1)  
            {
                N_U = Velocity_FeSpace[i+1]->GetN_DegreesOfFreedom();
                N_P = Pressure_FeSpace[i+1]->GetN_DegreesOfFreedom();    
                N_TotalDOF = 3*N_U + N_P;    
                
                sol_NSE = new double[N_TotalDOF];
                memset(sol_NSE, 0, N_TotalDOF*SizeOfDouble);
                Sol_array_NSE[i+1] = sol_NSE;
        
                rhs_NSE = new double[N_TotalDOF];
                memset(rhs_NSE, 0, N_TotalDOF*SizeOfDouble);
                Rhs_array_NSE[i+1] = rhs_NSE;
                
                u = new TFEVectFunct3D(Velocity_FeSpace[i+1], UString,  UString,  sol_NSE, N_U, 3);
                Velocity[i+1] = u;
                p = new TFEFunction3D(Pressure_FeSpace[i+1], PString,  PString,  sol_NSE+3*N_U, N_P);
                Pressure[i+1] = p;    
            }// if(i==LEVELS-1 && mg_type==1)
            
//======================================================================
// construct all CD finite element functions
//======================================================================
            N_T = Temperature_FeSpace[i]->GetN_DegreesOfFreedom();
            sol_CD = new double[N_T];
            rhs_CD = new double[N_T];
            Sol_array_CD[i] = sol_CD;
            Rhs_array_CD[i] = rhs_CD;
        
            Theta  = new TFEFunction3D(Temperature_FeSpace[i], TString, TString, sol_CD, N_T);
            Temperature[i] = Theta;
     
            if(i==LEVELS-1 && mg_type==1) 
            {  
                N_T = Temperature_FeSpace[mg_level-1]->GetN_DegreesOfFreedom();
                sol_CD = new double[N_T];
                rhs_CD = new double[N_T];
                Sol_array_CD[mg_level-1] = sol_CD;
                Rhs_array_CD[mg_level-1] = rhs_CD;
        
                Theta = new TFEFunction3D(Temperature_FeSpace[mg_level-1], TString, TString, sol_CD, N_T);
                Temperature[mg_level-1] = Theta;
            }//   if(i==LEVELS-1 && mg_type==1)     

#ifdef _MPI
            N_Cells = coll->GetN_Cells();
            printf("rank=%d\t N_Cells   : %d\t Dof all   :%d\t Dof Velocity :%d\t Dof Pressure: %d\n",rank,N_Cells,N_TotalDOF,3*N_U,N_P);
#endif
        } //  for(i=0;i<LEVELS;i++)
   
        u1 = Velocity[mg_level-1]->GetComponent(0);
        u2 = Velocity[mg_level-1]->GetComponent(1);
        u3 = Velocity[mg_level-1]->GetComponent(2);  

        oldrhs_NSE = new double[N_TotalDOF];
        oldrhs_CD = new double[N_T];
        oldsol = new double[N_T];
    
#ifndef _MPI       
        N_Cells = coll->GetN_Cells();
        OutPut("N_Cells      : "<< setw(10) << N_Cells <<endl);
        OutPut("Dof velocity : "<< setw(10) << 3*N_U << endl);
        OutPut("Dof pressure : "<< setw(10) << N_P << endl);
        OutPut("Dof all      : "<< setw(10) << N_TotalDOF  << endl);  
#endif 
//======================================================================
// SystemMatrix CD construction
//======================================================================
//         if (TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD)
//         {
//             fesp[0] = NULL;//Velocity_FeSpace;
//             fesp[1] = NULL;//Temperature_FeSpace;
//             
//             fefct[0] = u1;
//             fefct[1] = u2;
//             fefct[2] = u3;
//             
//             aux_CD =  new TAuxParam3D(TimeCDParamsVeloFieldN_FESpaces,
//                                       TimeCDParamsVeloFieldN_Fct,
//                                       TimeCDParamsVeloFieldN_ParamFct,
//                                       TimeCDParamsVeloFieldN_FEValues,
//                                       fesp+1, fefct,
//                                       TimeCDParamsVeloFieldFct,
//                                       TimeCDParamsVeloFieldFEFctIndex,
//                                       TimeCDParamsVeloFieldFEMultiIndex,
//                                       TimeCDParamsVeloFieldN_Params,
//                                       TimeCDParamsVeloFieldBeginParam);
//         }
//         else
//         {
//             aux_CD =  new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
//         }
        /** interpolate the initial value */
        Theta->Interpolate(InitialCondition);   
        
        /** Disc type: GALERKIN (or) SDFEM  (or) UPWIND (or) SUPG (or) LOCAL_PROJECTION
         *  Solver: AMG_SOLVE (or) GMG  (or) DIRECT */
        //         if(profiling){
        // #ifdef _MPI
        //             start_int = MPI_Wtime();
        // #else
        //             start_int = GetTime();
        // #endif
        //         }
        
        SystemMatrix_CD = new TSystemTCD3D(mg_level, Temperature_FeSpace, Sol_array_CD, Rhs_array_CD,
                                           TDatabase::ParamDB->DISCTYPE, TDatabase::ParamDB->SOLVER_TYPE);
        
        #ifdef _MPI
        if(rank==0)
            #endif
            printf("SystemMatrix constructed\n");
        /** initilize the system matrix with the functions defined in Example file */
        // last argument aux is used to pass  addition fe functions (eg. mesh velocity) that is nedded for assembling,
        // otherwise, just pass it with NULL 
        SystemMatrix_CD->Init(BilinearCoeffs, BoundCondition_CD, BoundValue, NULL);
        
        //         if(profiling){
        // #ifdef _MPI
        //             end_int = MPI_Wtime();
        // #else
        //             end_int = GetTime();
        // #endif
        //             total_int = end_int-start_int;
        //         }
        
        /** assemble the system matrix with given aux, sol and rhs 
         *  aux is used to pass  addition fe functions (eg. mesh velocity) that is nedded for assembling,
         *  otherwise, just pass it with NULL  */
        
        //         if(profiling){
        // #ifdef _MPI
        //             start_assembling = MPI_Wtime();
        // #else
        //             start_assembling = GetTime();
        // #endif
        //         }
        SystemMatrix_CD->AssembleMRhs(); 
        
        //         if(profiling){
        // #ifdef _MPI
        //             end_assembling = MPI_Wtime();
        // #else
        //             end_assembling = GetTime();
        // #endif
        //             total_assembling += (end_assembling-start_assembling);
        //         }
        /** copy rhs to oldrhs before calling the solver, as rhs will change in multigrid solver */
        memcpy(oldrhs_CD, rhs_CD, N_T*SizeOfDouble); 
//======================================================================
// SystemMatrix_NSE construction
//======================================================================  
        NSEType = TDatabase::ParamDB->NSTYPE;
    
//         if(profiling)
//         {
// #ifdef _MPI
//             start_int = MPI_Wtime();
// #else
//             start_int = GetTime();
// #endif
//         }
        // get a  TNSE3D system
        SystemMatrix_NSE = new TSystemTNSE3D(mg_level, Velocity_FeSpace, Pressure_FeSpace, Velocity, Pressure, 
                                       Sol_array_NSE, Rhs_array_NSE, TDatabase::ParamDB->DISCTYPE, NSEType, TDatabase::ParamDB->SOLVER_TYPE
#ifdef __PRIVATE__ 
                                       , Projection_FeSpace
#endif   
                                       );   
      
 
        // initilize the system matrix with the functions defined in Example file
        SystemMatrix_NSE->Init(LinCoeffs, BoundCondition_NSE, U1BoundValue, U2BoundValue, U3BoundValue);
      
        // cout<<"test main " << endl;
 
#ifdef _MPI
        if(rank==0)
#endif
        printf("SystemMatrix_NSE constructed\n");
    
//         if(profiling){
// #ifdef _MPI
//             end_int = MPI_Wtime();
// #else
//             end_int = GetTime();
// #endif
//         total_int = end_int-start_int;
//         }

/*        if(profiling){
#ifdef _MPI
        start_assembling_solving = MPI_Wtime();
#else
        start_assembling_solving = GetTime();
#endif
        }        */  
      
        //  interpolate the initial solution
        u1->Interpolate(InitialU1);
        u2->Interpolate(InitialU2);
        u3->Interpolate(InitialU3);

        // assemble M, A, rhs matrices at time t=0 with initia sol and rhs 
        SystemMatrix_NSE->Assemble();
       
        //  SystemMatrix_NSE->CheckAllMat();
        //  SystemMatrix_NSE->RHS_stats();
        //  norm_array(sol,rhs,N_TotalDOF);
        //  SystemMatrix_NSE->All_levels_check();
               
//======================================================================
// prepare for outout
//======================================================================
        VtkBaseName = TDatabase::ParamDB->VTKBASENAME;    
        Output = new TOutput3D(2, 2, 1, 1, Domain);
        Output->AddFEVectFunct(u);
        Output->AddFEFunction(p);
        Output->AddFEFunction(Theta);
   
#ifdef _MPI
        if(TDatabase::ParamDB->WRITE_VTK)
            Output->Write_ParVTK(MPI_COMM_WORLD, img, SubID);
            img++;       
#else  
        if(TDatabase::ParamDB->WRITE_VTK)
        {
            os.seekp(std::ios::beg);
            if(img<10) os <<  "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
            else if(img<100) os <<  "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
            else if(img<1000) os <<  "VTK/"<<VtkBaseName<<".00"<<img<<".vtk" << ends;
            else if(img<10000) os <<  "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
            else  os <<  "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
            Output->WriteVtk(os.str().c_str());
            img++;
        }   
#endif    
        //  SystemMatrix_NSE->GetMesh(); 
        //  
//======================================================================
// time disc loop
//======================================================================    
        // parameters for time stepping scheme
        m = 0;
        N_SubSteps = GetN_SubSteps();
        oldtau = 1.;
        end_time = TDatabase::TimeDB->ENDTIME; 
        limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE;
        Max_It = TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE;        
        memset(AllErrors, 0, 7.*SizeOfDouble);
        
        UpdateStiffnessMat = FALSE; //check BilinearCoeffs in example file
        UpdateRhs = TRUE; //check BilinearCoeffs in example file
        ConvectionFirstTime=TRUE;        
    }

#ifdef _SMPI
    if(rank !=0)
    { 
        solmumps = new TMumpsSolver(0,0,NULL,NULL,0);
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

    // time loop starts
    while(TDatabase::TimeDB->CURRENTTIME< end_time)   // time cycle
    {                                             
#ifdef _SMPI
        if(rank==0)
#endif
        {
            m++;
            TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;
        }

        for(l=0;l<N_SubSteps;l++) // sub steps of fractional step theta
        {
#ifdef _SMPI
            if(rank==0)
#endif
            {
                SetTimeDiscParameters(1);
#ifdef _MPI
                //MPI_Comm_rank(MPI_COMM_WORLD,&rank);
                if(rank == 0)
#endif
                {
                    if(m==1)
                    {
                        OutPut("Theta1: " << TDatabase::TimeDB->THETA1<< endl);
                        OutPut("Theta2: " << TDatabase::TimeDB->THETA2<< endl);
                        OutPut("Theta3: " << TDatabase::TimeDB->THETA3<< endl);
                        OutPut("Theta4: " << TDatabase::TimeDB->THETA4<< endl);
                    }
                }      
                tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
                TDatabase::TimeDB->CURRENTTIME += tau;
#ifdef _MPI
                if(rank==0)
#endif
                {   
                    OutPut(endl << "CURRENT TIME: ");
                    OutPut(TDatabase::TimeDB->CURRENTTIME << endl);          
                }
//====================================================================== 
//Solve the NSE System Matrix
//======================================================================                
                // copy rhs to oldrhs
                memcpy(oldrhs_NSE, rhs_NSE, N_TotalDOF*SizeOfDouble);        
  
                // assemble only rhs, nonlinear matrix for NSE will be assemble in fixed point iteration
                // not needed if rhs is not time-dependent
                if(m!=1)
                { SystemMatrix_NSE->AssembleRhs();}
                else
                { SystemMatrix_NSE->Assemble();}
       
                // cout<<"start printing"<<endl;
                // printall_array(sol,rhs,N_TotalDOF);
                // cout<<"end printing"<<endl;
                // exit(0);
                
                //scale B matices and assemble NSE-rhs based on the \theta time stepping scheme 
                SystemMatrix_NSE->AssembleSystMat(tau/oldtau, oldrhs_NSE, rhs_NSE, sol_NSE); 
                oldtau = tau;
                
                // calculate the residual
                SystemMatrix_NSE->GetResidual(sol_NSE, impuls_residual,residual);

#ifdef _MPI
                if(rank==0)
#endif
                {
                    OutPut(" nonlinear iteration step   0");
                    OutPut(setw(14) << impuls_residual);
                    OutPut(setw(14) << residual-impuls_residual);
                    OutPut(setw(14) << sqrt(residual) << endl);  
                }

            }
      
#ifdef _SMPI

            if(rank==0)
            {
                stime = TDatabase::TimeDB->CURRENTTIME ;
            }

            MPI_Bcast(&stime,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
            TDatabase::TimeDB->CURRENTTIME = stime;
#endif
      
//====================================================================== 
//Solve the system
//Nonlinear iteration of fixed point type
//======================================================================
            for(j=1;j<=Max_It;j++)
            { 
                checker = 0;
#ifdef _SMPI
                if(rank==0)
#endif
                {
                    SystemMatrix_NSE->Solve(sol_NSE);
                }
#ifdef _SMPI
                else
                {  
                    solmumps->FactorizeAndSolve(NULL,NULL);
                }
#endif

#ifdef _SMPI
                if(rank==0)
#endif
                {
                    if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
                    IntoL20FEFunction3D(sol_NSE+3*N_U, N_P, Pressure_FeSpace[mg_level-1]);  
       
                    //no nonlinear iteration for Stokes problem  
                    if(TDatabase::ParamDB->FLOW_PROBLEM_TYPE==STOKES)
                        break;
                    // restore the mass matrix for the next nonlinear iteration      
                    SystemMatrix_NSE->RestoreMassMatNonLinear();    
                    
                    // assemble the system matrix with given aux, sol and rhs 
                    SystemMatrix_NSE->AssembleNonLinear();    
       
                    // if(j==1){
                    // SystemMatrix_NSE->CheckAllMat();
                    // SystemMatrix_NSE->RHS_stats();
                    //     	
                    // exit(0);
                    // }

                    // assemble system mat, S = M + dt\theta_1*A
                    SystemMatrix_NSE->AssembleSystMatNonLinear();    
                    // get the residual
                    SystemMatrix_NSE->GetResidual(sol_NSE, impuls_residual,residual);
#ifdef _MPI
                    if(rank==0)
#endif
                    {
                        OutPut(" nonlinear iteration step " << setw(3) << j);
                        OutPut(setw(14) << impuls_residual);
                        OutPut(setw(14) << residual-impuls_residual);
                        OutPut(setw(14) << sqrt(residual) << endl);
                    }
                    if(sqrt(residual)<=limit)
                    {
#ifdef _SMPI
                        checker = -1;
#else
                        break;
#endif
                    }
                }
      
#ifdef _SMPI
                MPI_Bcast(&checker,1,MPI_INT,0,MPI_COMM_WORLD);
                if(checker ==-1)
                    break;
#endif	
            } // for(j=1;j<=Max_It;j++)     
      
#ifdef _SMPI
            if(rank==0)
#endif
            {
                // restore the mass matrix for the next time step
                SystemMatrix_NSE->RestoreMassMat();        
            }
            
//====================================================================== 
//Solve the CD System Matrix
//======================================================================             
#ifdef _SMPI
            if(rank==0)
#endif      
            {      
                memcpy(oldsol, sol_CD, N_T*SizeOfDouble);
/*                if(profiling)
                { 
#ifdef _MPI
                    start_assembling = MPI_Wtime();
#else
                    start_assembling = GetTime();
#endif
                }    */        

                if(UpdateStiffnessMat || UpdateRhs ||  ConvectionFirstTime)
                {  
                    if(UpdateRhs)
                    { SystemMatrix_CD->AssembleARhs(); }
                    else
                    { SystemMatrix_CD->AssembleARhs(); }
         
                    SystemMatrix_CD->AssembleSystMat(oldrhs_CD, oldsol, rhs_CD, sol_CD
#ifdef _MPI
                                                                , Rhs_array_CD
#endif
                                              );
                    /** copy rhs to oldrhs before calling the solver, as rhs will change in multigrid solver */
                    memcpy(oldrhs_CD, rhs_CD, N_T*SizeOfDouble); 
                    ConvectionFirstTime = FALSE;
                }

/*                if(profiling)
                {
#ifdef _MPI
                    end_assembling = MPI_Wtime();
#else
                    end_assembling = GetTime();
#endif
                    total_assembling += (end_assembling-start_assembling);
                }
                // solve the system matrix 
                if(profiling)
                {
#ifdef _MPI
                    start_solve = MPI_Wtime();
#else
                    start_solve = GetTime();
#endif
                }   */    

                SystemMatrix_CD->Solve(sol_CD);
            
//                 if(profiling)
//                 {
// #ifdef _MPI
//                     end_solve = MPI_Wtime();
// #else
//                     end_solve = GetTime();
// #endif
//                     total_solve += (end_solve-start_solve);
//                 }
            
                if(UpdateStiffnessMat || UpdateRhs)
                {         
                    SystemMatrix_CD->RestoreMassMat();
                }
            }
        }
//======================================================================
// produce outout
//======================================================================           
#ifdef _SMPI
        if(rank==0)
#endif
        {
            if(m==1 || m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)  
       
 #ifdef _MPI
            if(TDatabase::ParamDB->WRITE_VTK)
                Output->Write_ParVTK(MPI_COMM_WORLD, img, SubID);
            img++;       
#else     
            if(TDatabase::ParamDB->WRITE_VTK)
            {
                os.seekp(std::ios::beg);
                if(img<10) os <<  "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
                else if(img<100) os <<  "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
                else if(img<1000) os <<  "VTK/"<<VtkBaseName<<".00"<<img<<".vtk" << ends;
                else if(img<10000) os <<  "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
                else  os <<  "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
                    Output->WriteVtk(os.str().c_str());
                img++;
            }
 #endif   
        }
    } // while(TDatabase::TimeDB->CURRENTTIME< e
      
    CloseFiles();
#ifdef _MPI
    MPI_Finalize();
#endif
    return 0;
} // end main
