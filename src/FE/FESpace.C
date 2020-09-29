// =======================================================================
// @(#)FESpace.C        1.5 09/15/99
// 
// Class:       TFESpace
// Purpose:     general super class for all finite element spaces
//              special spaces are implemented in subclasses
//
// Author:      Gunar Matthies (04.11.97)
//
// History:     start of implementation 04.11.97 (Gunar Matthies)
//
//              split TFESpace into TFESpacexD (15.04.1998) Volker Behns
//
// =======================================================================

#include <Constants.h>
#include <FESpace.h>
#include <MooNMD_Io.h>
#include <string.h>

#include "Vector.h"

/** copying given parameters into inner storage places */
int TFESpace::InitData(TCollection *coll, char *name, char *description)
{
  Name=strdup(name);
  Description=strdup(description);

  // cout << "Name: " << name << endl;
  // cout << "Desc: " << Description << endl;

  Collection=coll;
  N_Cells=Collection->GetN_Cells();
  N_DegreesOfFreedom=0;
  GlobalNumbers=NULL;
  BeginIndex=NULL;
  N_UsedElements=0;

  N_Dirichlet=0;
  DGSpace = 0; // use 'void SetAsDGSpace()' to change this

 # ifdef _MPI
  MaxSubDomainPerDof = -1;
  // cout << "number of cells: " << N_Cells << endl;
 # endif
  return 0;
}

/** Constructor */
TFESpace::TFESpace(TCollection *coll, char *name, char *description)
{
  InitData(coll, name, description);
}

/** destructor */
TFESpace::~TFESpace()
{
  delete Name;
  delete Description;
  delete [] GlobalNumbers;
  delete [] BeginIndex;
  delete [] BoundaryNodeTypes;
  delete [] N_BoundaryNodes;
}

/** write info on fespace into file */
int TFESpace::Write(const char *filename)
{
  int header[4];
  int N_LocalDOF;

  std::ofstream dat(filename);
  if(!dat)
  {
    cerr << "cannot open file '" << filename << "' for output" << endl;
    return -1;
  }

  N_LocalDOF = BeginIndex[N_Cells];
  header[0] = N_Cells;
  header[1] = N_DegreesOfFreedom;
  header[2] = ActiveBound;
  header[3] = N_LocalDOF;

  dat.write((char *)header, sizeof(int)*4);
  
  dat.write((char *)BeginIndex, sizeof(int)*(N_Cells+1));
  dat.write((char *)GlobalNumbers, sizeof(int)*N_LocalDOF);
  
  dat.close();

  return 0;
}
