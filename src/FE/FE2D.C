// =======================================================================
// @(#)FE2D.C        1.1 10/30/98
//
// Class:       TFE2D
// Purpose:     store all information for one finite element class
//
// Author:      Gunar Matthies  09.07.98
//
// =======================================================================
#ifdef _MPI
#  include "mpi.h"
# endif

#include <FE2D.h>
#include <FEDatabase2D.h>
#include <NodalFunctional2D.h>
#include <Database.h>


/** constructor */
TFE2D::TFE2D()
{
}

/** constructor with data */
TFE2D::TFE2D(BaseFunct2D basefunct_id, NodalFunctional2D nodalfunctional_id,
         RefTrans2D reftransid, FEDesc2D fedesc_id, int n_info)
{
  BaseFunct_ID = basefunct_id;
  BaseFunct = TFEDatabase2D::GetBaseFunct2D(BaseFunct_ID);

  NodalFunctional_ID = nodalfunctional_id;
  NodalFunctional  = TFEDatabase2D::GetNodalFunctional2D(NodalFunctional_ID);

  RefTransID = reftransid;

  FEDesc_ID = fedesc_id;
  FEDesc = TFEDatabase2D::GetFEDesc2D(FEDesc_ID);

  N_Info = n_info;
  N_DOF = FEDesc->GetN_DOF();

  Size = N_Info + N_DOF;
}

/** check N[i](b[j]) = delta[ij] */
void TFE2D::CheckNFandBF()
{
  int i,j,k,l, N_Points;
  double *xi, *eta;
  double PointValues[MaxN_PointsForNodal2D];
  double FunctionalValues[MaxN_BaseFunctions2D];
  double AllPointValues[MaxN_PointsForNodal2D][MaxN_BaseFunctions2D];

  NodalFunctional->GetPointsForAll(N_Points, xi, eta);

  for(k=0;k<N_Points;k++)
    BaseFunct->GetDerivatives(D00, xi[k], eta[k], 
                              AllPointValues[k]);
#ifdef _MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(rank==TDatabase::ParamDB->Par_P0 && TDatabase::ParamDB->SC_VERBOSE>0)
#endif 
  {
   cout << "CheckNFandBF: " << "BaseFunct_ID: " << BaseFunct_ID << " ";
   cout << "NodalFunctional_ID: " << NodalFunctional_ID << endl;
  }

  for(k=0;k<N_DOF;k++)
  {
    for(l=0;l<N_Points;l++)
      PointValues[l] = AllPointValues[l][k];

    NodalFunctional->GetAllFunctionals(NULL, NULL, PointValues,
                          FunctionalValues);

    for(i=0;i<N_DOF;i++)
    {
      if(fabs(FunctionalValues[i])<1e-10)
        FunctionalValues[i] = 0;
      // cout << k << " " << i << " " << FunctionalValues[i] << endl;
      if( i == k )
        if( fabs(FunctionalValues[i]-1) > 1e-8 )
          cout << "BF: " << k << " NF: " << i << " " << FunctionalValues[i] << endl;
      if( i != k )
        if( fabs(FunctionalValues[i]-0) > 1e-8 )
          cout << "BF: " << k << " NF: " << i << " " << FunctionalValues[i] << endl;
    }
  }
#ifdef _MPI
  if(rank==TDatabase::ParamDB->Par_P0 && TDatabase::ParamDB->SC_VERBOSE>0)
#endif 
  cout << endl;
}
