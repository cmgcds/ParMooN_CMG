// =======================================================================
//
// Purpose:     Assignment 1
//
// Author:      Shiladitya Banerjee
//
// History:     Implementation started on 14.07.2012
// =======================================================================

 /** INCLUDES*/
#include <mpi.h>
#include <Domain.h>
#include <Database.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <MainUtilities.h>
#include <MeshPartition2D.h>
 
using namespace std;

int main(int argc, char* argv[])
{
  TDatabase *Database = new TDatabase();

  const int root = 0;
  int rank, size, len;
  double t_par1, t_par2;
  char ReadinDat [] = "readin.dat";

  MPI_Init(&argc, &argv);
  
  MPI_Comm comm;
  comm = MPI_COMM_WORLD;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
   
  /** Variable List */
  
  TDomain *Domain = new TDomain();
  TCollection *coll;
  TBaseCell **trial;

  int i,ret, ORDER, N_DOF, N_Active, N_NonActive,MaxCpV;
  int N_Unknowns, img=1, N_RDOF[0],N_cells;


  std::ostringstream os;
  
  char *PRM, *GEO;
  
/** *********************************************************************/
/** READ PARAMETER FILE readin.dat
/** *********************************************************************/
  if(argc>=2)
    ret=Domain->ReadParam(argv[1]);
  else
    ret=Domain->ReadParam(ReadinDat);

  if(ret==-1)
    exit(-1);

  OpenFiles();
  OutFile.setf(std::ios::scientific);

  Database->WriteParamDB(argv[0]);

/** *********************************************************************/
/** COPY READ PARAMETERS INTO LOCAL VARIABLE
/** *********************************************************************/
  PRM = TDatabase::ParamDB->BNDFILE;
  GEO = TDatabase::ParamDB->GEOFILE;
  Domain->Init(PRM, GEO);
  
  
  for(i=0;i<TDatabase::ParamDB->UNIFORM_STEPS;i++)
     Domain->RegRefineAll();

  /*if(rank == 0)
  {
    os.seekp(std::ios::beg);
    os << "Domain" << ".ps" << ends;
    Domain->PS(os.str().c_str(),It_Finest,0);
  }*/
 //Domain->GetTreeInfo(trial,N_cells);
 //printf(" NUMBER OF CELLS IN SUBDOMAIN %d : \n  OWN CELLS : %d\n ",rank,N_cells); 
 Partition_Mesh2D(comm,Domain,MaxCpV);//4
//   MPI_Barrier(comm);
 /** ************************************************************************/
 /** PRINTING SUBDOMAINS INTO POSTSCRIPT*/
 /** ************************************************************************/ 
//  coll = Domain->GetOwnCollection(It_Finest, 0, rank); 
//  printf(" NUMBER OF CELLS IN SUBDOMAIN %d : \n  OWN CELLS : %d\n ",rank,coll->GetN_Cells()); 
 
   os.seekp(std::ios::beg);
  os << "Domain"<<rank << ".ps" << ends;
    Domain->PS(os.str().c_str(),It_Finest,0);
  
  MPI_Finalize();
  return 0;
}

