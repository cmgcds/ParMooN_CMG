// =======================================================================
// Purpose:     main program for 1D with new kernels of ParMooN
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 12.12.2020
// =======================================================================

#include <Domain.h>
#include <Database.h>
#include <SystemCD1D.h>
#include <FEDatabase2D.h>

// =======================================================================
// include current example
// =======================================================================

#include "../Examples/CD_1D/cd1d.h"

// main program

int main(int argc, char* argv[])
{
  TDatabase *Database = new TDatabase();
  TSystemCD1D *SystemCD;
  TFEDatabase2D *FEDatabase = new TFEDatabase2D(); 
 
  std::ostringstream os;
  os << " ";   
  //======================================================================
  // CD1D System construction 
  //======================================================================  
  SystemCD = new TSystemCD1D(10, 0, 1, BoundCondition_LminLMax, BoundVales, argv[1]);
 
  // initiaize system
  SystemCD->Init(BilinearCoeffs);

//   // interpolate (if needed)
//   SystemCD->Interpolate(Exact);
 
  // interpolate (if needed)
  SystemCD->Solve();


CloseFiles();

return 0;
} // end main
