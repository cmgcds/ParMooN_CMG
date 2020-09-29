#include <Database.h>
#include <Domain.h>
#include <Collection.h>
#include <MainUtilities.h>
#include <LinAlg.h>

#include <PardisoSolver.h>

#ifdef __3D__
  #include <FEDatabase3D.h>
  #include <Output3D.h>
  #include <MeshPartition.h>
  #include <DiscreteForm3D.h>
  #include <FESpace3D.h>
  #include <SquareMatrix3D.h>
  #include <Matrix3D.h>
  #include <AuxParam3D.h>
  #include <Assemble3D.h>
  #include <FreeSurface3D.h>
  #include <RefTrans3D.h>
  #include <MovingTNSE3D.h>
  #include <TetGenMeshLoader.h>
#else

#endif

#ifdef _MPI
  #include "mpi.h"
  #include <MumpsSolver.h>

  #ifdef __3D__
    #include <ParFECommunicator3D.h>
    #include <ParVector3D.h>
  #else
  
  #endif
#endif