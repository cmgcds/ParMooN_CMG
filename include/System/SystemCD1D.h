// =======================================================================
// @(#)SystemCD1D.h        4.1 13.04.20
// 
// Class:       TSYSTEMCD1D
// Purpose:     general super class for all 1D System
//
// Author:      Sashikumaar Ganesan (12.12.20)
//
// History:     start of implementation 12.12.20 (Sashikumaar Ganesan)
//
// =======================================================================

#ifndef __SYSTEMCD1D__
#define __SYSTEMCD1D__

#include <System1D.h>
// #include <Domain.h>
// #include <BaseCell.h>
// #include <Collection.h>
// #include <FESpace1D.h>
// #include <FEFunction1D.h>
// #include <SquareMatrix1D.h>

/** general super class for 1D */
class TSystemCD1D : public TSystem1D
{
  protected:


  public:
  
  /** constructor */  
  TSystemCD1D(int N_L, double start, double end, BoundCond1D *boundConLminLMax, DoubleFunctND *BdValues, char *ParamFile);

  void Init(CoeffFctND *bilinear);

  /** interpolate */
  void Interpolate(DoubleFunctND *Exact);

  /** assemble and solve the system */
  void Solve();

  /** assemble the stiffness mat */
  void AssembleARhs();
  void AssembleARhs_SUPG();
  void AssembleARhs_FD();
  void AssembleARhs_DG();

  /** set BC */
  void SetDirichletBc();

  //** Print the Solution **/
  void printSolution();

  //** Generate the Solution to VTK **/
  void generateVTK();


  // ** Generate Output in GNUPLOT ** //
  void plotGNU();



    /** destrcutor */
    ~TSystemCD1D();

};

#endif
