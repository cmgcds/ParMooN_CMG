#ifndef _MPI
  #error ParConvDiff.C: program must be compiled with MPI option enabled 
#else

#include <MooNMD.h>

#include <stdlib.h>
#include <string.h>
#include <sstream>

#include "../Examples/CD_3D/Plane.h"

#define EXITMAIN(x) MPI_Barrier(comm); MPI_Finalize();  return (x)

#ifdef BEGIN_SEQ
  #undef BEGIN_SEQ
#endif

#ifdef END_SEQ
  #undef END_SEQ
#endif

#define BEGIN_SEQ \
  MPI_Barrier(comm); \
  for(int _i=0;_i<size;++_i) \
  { \
    if ( g_rank == _i ) \
    { 
    
#define END_SEQ \
    } \
    MPI_Barrier(comm); \
  }
  
#define BEGIN_SEQ_WORKER \
  MPI_Barrier(comm); \
  for(int _i=1;_i<size;++_i) \
  { \
    if ( g_rank == _i ) \
    { 
    
#define END_SEQ_WORKER \
    } \
    MPI_Barrier(comm_worker); \
  }



int g_rank;
double bound = 0;

void WriteSolution(const char *basename, int number, TOutput3D *Output);
void WriteGrid(const char *name, TDomain *Domain);

int main(int argc, char **argv)
{  
  int size;
  MPI_Comm comm=MPI_COMM_WORLD;
  MPI_Comm comm_worker;
  MPI_Group world_group, worker_group;
  int root[1] = { 0 };
  
  MPI_Init(&argc, &argv);
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &g_rank);  
  
  TDomain Domain;
  TDatabase Database;
  TFEDatabase3D FEDatabase; 
  TCollection *Coll;
  TAuxParam3D *aux_dummy;
  TFESpace3D *FESP[1], *RHSFESP[3];
  TFEFunction3D *u;
  TOutput3D *Output;
  
  TSquareMatrix3D *SQMATRICES[1];
  TSquareMatrix3D *A = NULL;
  
  BoundCondFunct3D *BOUNDCONDITIONS[3];
  BoundValueFunct3D *BOUNDVALUES[3];
  
  int MaxCpV, Order, N_Unknowns, N_DOF;
  int N_MAT, N_SQMAT, N_RHS, N_FESP;
  double *sol, *rhs;
  double *RHS[3];
  
  char *SMESH;
  
  char name[] = "name";
  char desc[] = "desc";
  
  if ( argc > 1 )
    Domain.ReadParam(argv[1]);
  else
  {
    OutPut("No readin file given!"<< endl);
    EXITMAIN(0);
  }
  
//   if ( TDatabase::ParamDB->Par_P1 == 0 )
//   {
//     MPI_Comm_group(comm, &world_group);
//     MPI_Group_excl(world_group, 1, root, &worker_group);
//     MPI_Comm_create(comm, worker_group, &comm_worker);
//     MPI_Group_free(&worker_group);
//     MPI_Group_free(&world_group);
//   }
//   else
//     comm_worker = MPI_COMM_WORLD;
  
  BOUNDCONDITIONS[0] = BoundCondition;
  BOUNDCONDITIONS[1] = BoundCondition;
  BOUNDCONDITIONS[2] = BoundCondition;
  BOUNDVALUES[0] = BoundValue;
  BOUNDVALUES[1] = BoundValue;
  BOUNDVALUES[2] = BoundValue;
  
  SMESH = TDatabase::ParamDB->SMESHFILE;
  
  Domain.Tetgen(SMESH);
  
  for (int i=0;i<TDatabase::ParamDB->UNIFORM_STEPS;++i)
    Domain.RegRefineAll();
  
  Partition_Mesh(comm, &Domain, MaxCpV);
  BEGIN_SEQ
  OutPut("MaxCpV: " << MaxCpV << endl);
  END_SEQ
  
  WriteGrid("grid", &Domain);

  // discreteform
  TDiscreteForm3D DiscreteForm(name, desc, N_Terms, Derivatives, SpacesNumbers,
			       N_Matrices, N_Rhs, RowSpace, ColumnSpace, RhsSpace,
			       BilinearAssemble, BilinearCoeffs, NULL);

  // collection
  Coll = Domain.GetCollection(It_Finest, 0);

MPI_Barrier(comm);
  // fe space
  Order = TDatabase::ParamDB->VELOCITY_SPACE;
  TFESpace3D fesp (Coll, name, desc, BoundCondition, Order);
  TParFECommunicator3D FEComm (&comm, &fesp);
  FEComm.Test();
// EXITMAIN(0);


  // structure
  TSquareStructure3D sqstructure (&fesp);
  sqstructure.Sort();
//   FEComm.WriteStructure(&sqstructure);
// EXITMAIN(0);
MPI_Barrier(comm);
  TMumpsSolver *MumpsSolver = new TMumpsSolver(comm, &sqstructure, &FEComm, 1);
MPI_Barrier(comm);

  N_DOF = fesp.GetN_DegreesOfFreedom();
  N_Unknowns = N_DOF;
  
  BEGIN_SEQ
  OutPut("N_DOF: " << N_DOF << endl);
  OutPut("N_Unknowns: " << N_Unknowns << endl);
  END_SEQ
  
  sol = new double [N_Unknowns];
  rhs = new double [N_Unknowns];
  memset(sol, 0, N_Unknowns*sizeof(sol[0]));
  memset(rhs, 0, N_Unknowns*sizeof(rhs[0]));
  
  // fe function and output (only on root)
  if ( g_rank == 0 )
  {
    u = new TFEFunction3D (&fesp, name, desc, sol, N_DOF);
    Output = new TOutput3D (0, 1, 0, 0, &Domain);
    
    Output->AddFEFunction(u);
  }
  
  TParVector3D ParSol (&FEComm, N_Unknowns, sol);
  TParVector3D ParRhs (&FEComm, N_Unknowns, rhs);
    
MPI_Barrier(comm);
//   BEGIN_SEQ
  if ( g_rank != 0 || TDatabase::ParamDB->Par_P1 == 1)
  {   
    // matrices
    A = new TSquareMatrix3D(&sqstructure);
 
    N_SQMAT = 1;
    N_MAT   = 0;
    SQMATRICES[0] = A;
    SQMATRICES[0]->Reset();
    
    // rhs
    RHS[0] = rhs;
    N_RHS = 1;
    RHSFESP[0] = &fesp;
    
    FESP[0] = &fesp;
    N_FESP = 1;
    
    // aux object
    aux_dummy = new TAuxParam3D (0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, 0, NULL);
    
    Assemble3D(N_FESP, FESP,
	       N_SQMAT, SQMATRICES,
	       N_MAT, NULL,
	       N_RHS, RHS, RHSFESP,
	       &DiscreteForm,
	       BOUNDCONDITIONS, BOUNDVALUES, aux_dummy);
  }
//   END_SEQ

MPI_Barrier(comm);
  ParRhs.AssembleAtRoot();
MPI_Barrier(comm);
  MumpsSolver->Factorize(A, &FEComm);
MPI_Barrier(comm);
  MumpsSolver->Solve(sol, rhs);
MPI_Barrier(comm); 
  if ( g_rank == 0 )
  {
    WriteSolution("sol", 0, Output);
  }

  delete Coll;
  delete MumpsSolver;
  
  delete [] sol;
  delete [] rhs;
  
  if ( g_rank != 0 )
  {
    delete A;
  }
  
  if ( g_rank == 0 )
  {
    delete Output;
    delete u;
  }
  
//   if ( TDatabase::ParamDB->Par_P1 == 0 ) 
//      MPI_Comm_free(&comm_worker);
  MPI_Finalize();  
  return 0;
}

void WriteSolution(const char *basename, int number, TOutput3D *Output)
{
  std::ostringstream os;
  
  os << basename << ".";
  os.width(5);
  os.fill('0');
  os << number << "." << g_rank << ".vtk" << ends;
  
  Output->WriteVtk(os.str().c_str());
}

void WriteGrid(const char *name, TDomain *Domain)
{
  std::string str;
  std::ostringstream os;
  
  os << g_rank;
  
  str += name;
  str += ".";
  str += os.str();
  str += ".vtk";
  
  TCollection *Coll = Domain->GetCollection(It_Finest, 0);
  TOutput3D Output(0,0,0,0,Domain, Coll);
  
  Output.WriteVtk(str.c_str());
  
//   PrintMesh(Coll);
  
  delete Coll;
}

#endif
