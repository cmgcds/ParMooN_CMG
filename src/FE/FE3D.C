// =======================================================================
// %W% %G%
//
// Class:       TFE3D
// Purpose:     store all information for one finite element class
//
// Author:      Gunar Matthies  19.11.99
//
// =======================================================================
#ifdef _MPI
#  include "mpi.h"
# endif

#include <FE3D.h>
#include <FEDatabase3D.h>
#include <Database.h>

/** constructor */
TFE3D::TFE3D()
{
}

/** constructor with data */
TFE3D::TFE3D(BaseFunct3D basefunct_id, NodalFunctional3D nodalfunctional_id,
         RefTrans3D reftransid, FEDesc3D fedesc_id, int n_info)
{
  BaseFunct_ID = basefunct_id;
  BaseFunct = TFEDatabase3D::GetBaseFunct3D(BaseFunct_ID);

  NodalFunctional_ID = nodalfunctional_id;
  NodalFunctional  = TFEDatabase3D::GetNodalFunctional3D(NodalFunctional_ID);

  RefTransID = reftransid;

  FEDesc_ID = fedesc_id;
  FEDesc = TFEDatabase3D::GetFEDesc3D(FEDesc_ID);

  N_Info = n_info;
  N_DOF = BaseFunct->GetDimension();

  Size = N_Info + N_DOF;
}

/** check N[i](b[j]) = delta[ij] */
void TFE3D::CheckNFandBF()
{
  int i,j,k,l, N_Points;
  double *xi, *eta, *zeta;
  
  NodalFunctional->GetPointsForAll(N_Points, xi, eta, zeta);
  
  int n_basis_functions = BaseFunct->GetDimension();
  int baseVectDim = BaseFunct->GetBaseVectDim();
  double PointValues[N_Points*baseVectDim];
  double FunctionalValues[n_basis_functions];
  double AllPointValues[N_Points][n_basis_functions*baseVectDim];

  for(k=0;k<N_Points;k++)
    BaseFunct->GetDerivatives(D000, xi[k], eta[k], zeta[k],
                              AllPointValues[k]);
#ifdef _MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(rank==TDatabase::ParamDB->Par_P0)
#endif
  if(TDatabase::ParamDB->SC_VERBOSE > 0)
  {
   cout << "CheckNFandBF: " << "BaseFunct_ID: " << BaseFunct_ID << " ";
   cout << "NodalFunctional_ID: " << NodalFunctional_ID << endl;
  }
  for(k=0;k<N_DOF;k++)
  {
    for(l=0;l<N_Points;l++)
      for(i = 0; i< baseVectDim; i++)
        PointValues[l+N_Points*i] = AllPointValues[l][k+n_basis_functions*i];

    NodalFunctional->GetAllFunctionals(NULL, NULL, PointValues,
                          FunctionalValues);

    for(i=0;i<N_DOF;i++)
    {
      if(fabs(FunctionalValues[i])<1e-10)
        FunctionalValues[i] = 0;
      if( i == k )
        if( fabs(FunctionalValues[i]-1) > 1e-7 )
          cout << "BF: " << k << " NF: " << i << " " << FunctionalValues[i] << endl;
      if( i != k )
        if( fabs(FunctionalValues[i]-0) > 1e-7 )
          cout << "BF: " << k << " NF: " << i << " " << FunctionalValues[i] << endl;
    }
  }

#ifdef _MPI
  if(rank==TDatabase::ParamDB->Par_P0)
#endif 
  if(TDatabase::ParamDB->SC_VERBOSE > 0)
  cout << endl;
}
