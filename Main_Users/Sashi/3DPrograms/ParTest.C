#include <MooNMD.h>

#ifdef _MPI

#include <ConvDiff3D.h>

#include <stdlib.h>

#include "../Examples/CD_3D/Plane.h"

#define EXITMAIN(x) MPI_Barrier(comm); MPI_Finalize();  return (x)

int g_rank = 0;
double bound = 0;

void WriteGrid(const char *name, TDomain *Domain);

int main(int argc, char **argv)
{
  int size, rank;
  MPI_Comm comm = MPI_COMM_WORLD;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);
  g_rank = rank;
  
  TDomain Domain;
  TDatabase Database;
  TFEDatabase3D FEDatabase;
  
  BoundCondFunct3D *BOUNDCONDITIONS[3];
  BoundValueFunct3D *BOUNDVALUES[3];
  
  int MaxCpV;
  
  char name[] = "name";
  char desc[] = "desc";
  
  if ( argc > 1 )
    Domain.ReadParam(argv[1]);
  else
  {
    OutPut("No readin file given!"<< endl);
    EXITMAIN(0);
  }
  
  char *SMESH = TDatabase::ParamDB->SMESHFILE;
  
  Domain.Tetgen(SMESH);
  
  for (int i=0;i<TDatabase::ParamDB->UNIFORM_STEPS;++i)
    Domain.RegRefineAll();
  
  Partition_Mesh(comm, &Domain, MaxCpV);
  
  WriteGrid("grid", &Domain);
  
  BOUNDCONDITIONS[0] = BoundCondition;
  BOUNDCONDITIONS[1] = BoundCondition;
  BOUNDCONDITIONS[2] = BoundCondition;
  BOUNDVALUES[0] = BoundValue;
  BOUNDVALUES[1] = BoundValue;
  BOUNDVALUES[2] = BoundValue;
  
  TDiscreteForm3D DiscreteForm(name, desc, N_Terms, Derivatives, SpacesNumbers,
			       N_Matrices, N_Rhs, RowSpace, ColumnSpace, RhsSpace,
			       BilinearAssemble, BilinearCoeffs, NULL);
  
  TCollection *Coll = Domain.GetCollection(It_Finest, 0);
  
  int Order = TDatabase::ParamDB->VELOCITY_SPACE;
  TFESpace3D fesp (Coll, name, desc, BoundCondition, Order);
  TParFECommunicator3D FEComm (&comm, &fesp);
  FEComm.Test();
  
  TSquareStructure3D sqstructure (&fesp);
  sqstructure.Sort();
  
  int N_Unknowns = fesp.GetN_DegreesOfFreedom();
  
  double *sol = new double [N_Unknowns];
  double *rhs = new double [N_Unknowns];
  double *rhs_root = new double [N_Unknowns];
  memset(sol, 0, N_Unknowns*sizeof(sol[0]));
  memset(rhs, 0, N_Unknowns*sizeof(rhs[0]));
  memset(rhs_root, 0 , N_Unknowns*sizeof(rhs_root[0]));
  TParVector3D ParSol (&FEComm, N_Unknowns, sol);
  TParVector3D ParRhs (&FEComm, N_Unknowns, rhs);
 
  if ( rank == 0 ) // root
  {
    for (int i=0;i<N_Unknowns;++i)
      sol[i] = i;
  }
  
  ParSol.ScatterFromRoot();
  
  TAuxParam3D aux_dummy (0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, 0, NULL);
  
  TSquareMatrix3D A (&sqstructure);
  TSquareMatrix3D *SQMATRICES[10];
  
  int N_SQMAT = 1;
  int N_MAT   = 0;
  SQMATRICES[0] = &A;
  SQMATRICES[0]->Reset();
  
  double *RHS[10];
  TFESpace3D *RHSFESP[10], *FESP[10];
  int N_RHS = 1;
  int N_FESP = 1;
  RHS[0] = rhs;
  FESP[0] = &fesp;
  RHSFESP[0] = &fesp;
  
  Assemble3D(N_FESP, FESP,
	       N_SQMAT, SQMATRICES,
	       N_MAT, NULL,
	       N_RHS, RHS, RHSFESP,
	       &DiscreteForm,
	       BOUNDCONDITIONS, BOUNDVALUES, &aux_dummy);
	  
  memset(rhs, 0, N_Unknowns*sizeof(rhs[0]));
  memset(rhs_root, 0 , N_Unknowns*sizeof(rhs_root[0]));
  
  if ( rank != 0 ) // worker
  {
    MatVectActive(&A, sol, rhs);
  }
  else // root
  {
    MatVectActive(&A, sol, rhs_root);
  }
  
  ParRhs.AssembleAtRoot();
  
  OutPut("blaaaaa" << endl);
  
  if ( rank == 0 ) // root
  {
    for (int i=0;i<N_Unknowns;++i)
    {
      OutPut("root: " << rhs_root[i] << "  |  " << rhs[i] << "  |  " << rhs_root[i] / rhs[i] << endl);
      
      rhs_root[i] = rhs_root[i] - rhs[i];
    }
    
    OutPut(Dnorm(N_Unknowns, rhs_root) << endl);
  }
  
  MPI_Finalize();  
  return 0;
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