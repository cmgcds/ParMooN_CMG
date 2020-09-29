// =======================================================================
// @(#)DefineParams.h        1.3 10/18/99
// 
// Purpose:     defines for mortar
//
// Author:      Volker Behns 13.09.99
//
// =======================================================================

#ifdef __MORTAR__
  // additional link term for SDFEM
  //#define __ADD_LINK_SDFEM__
  // additional link term for upwind
  //#define __ADD_LINK_UPW__

  // continuity in cross points
  #define __CONNECT_CROSSPOINTS__
#endif // __MORTAR__

// correct new midpoints (regular refinement of quadrangles)
#define __CORRECT_MP__

// -------------- end of changable parameters ------------------

#ifdef __ADD_LINK_SDFEM__
#ifndef __ADD_LINK__
// extension of matrix structure
#define __ADD_LINK__
#endif
#endif

#ifdef __ADD_LINK_UPW__
#ifndef __ADD_LINK__
// extension of matrix structure
#define __ADD_LINK__
#endif
#endif
