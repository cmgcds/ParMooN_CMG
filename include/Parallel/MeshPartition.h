// =======================================================================
// @(#)MeshPartition.h
// 
// Purpose:     partition the domain into "npart" parts for parallel computing
// 
// Author:      Sashikumaar Ganesan
// History:      start of implementation  07/09/09 (Sashikumaar Ganesan)
// =======================================================================
#include <Domain.h>

void Sort(TVertex **Array, int length);
int GetIndex(TVertex **Array, int Length, TVertex *Element);
#ifdef  _MPI 
void Partition_Mesh(MPI_Comm comm, TDomain *Domain, int &MaxCpV);

int Partition_Mesh3D(MPI_Comm comm, TDomain *Domain, int &MaxCpV);

void Domain_Crop(MPI_Comm comm, TDomain *Domain);

#endif

