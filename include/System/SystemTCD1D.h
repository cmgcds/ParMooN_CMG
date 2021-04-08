// =======================================================================
// @(#)SystemCD1D.h        4.1 13.04.20
// 
// Class:       TSYSTEMTCD1D
// Purpose:     general super class for all TCD 1D System
//
// Author:      Sashikumaar Ganesan (12.12.20)
//
// History:     start of implementation 12.12.20 (Sashikumaar Ganesan)
//
// =======================================================================

#ifndef __SYSTEMTCD1D__
#define __SYSTEMTCD1D__

#include <SystemCD1D.h>
#include <Domain.h>
#include <BaseCell.h>
#include <Collection.h>
#include <FESpace1D.h>
#include <FEFunction1D.h>
#include <SquareMatrix1D.h>

/** general super class for 1D */
class TSystemTCD1D : public TSystemCD1D
{
  protected:

  public:
    /** constructor */  
    TSystemTCD1D(int N_L, double start, double end, BoundCond1D *boundConLminLMax, DoubleFunctND *BdValues, char *ParamFile);

    /** destrcutor */
    ~TSystemTCD1D();

};

#endif
