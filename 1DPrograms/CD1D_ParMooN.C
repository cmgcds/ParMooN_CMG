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

  // Read the Param Values with Empty constructor
  // The actual TDomain Constructor is constructed inside the System1D
  TDomain* Domain = new TDomain();
  Domain->ReadParam(argv[1]);

  // Obtain FE Mesh Details from the Param File
  int N_Elements    =  TDatabase::ParamDB->N_ELEMENTS_1D;
  double Start_X    =  TDatabase::ParamDB->START_X;
  double End_X    =  TDatabase::ParamDB->END_X;

  int FE_Order      =  TDatabase::ParamDB->ANSATZ_ORDER ;
  int N_DOF         =  N_Elements*(FE_Order) + 1;

  // cout << " NELEM : " << N_Elements<<endl;
  // cout << " Start : " << Start_X<<endl;
  // cout << " End : " << End_X<<endl;
  // cout << " FE ORDER : " << FE_Order <<endl;
  // cout << " N_DOF : " << N_DOF <<endl;
  //======================================================================
  // CD1D System construction 
  //======================================================================  
  SystemCD = new TSystemCD1D(N_Elements, Start_X, End_X, BoundCondition_LminLMax, BoundVales, argv[1]);

  // initiaize system
  SystemCD->Init(BilinearCoeffs);

//   // interpolate (if needed)
//   SystemCD->Interpolate(Exact);
 
  // interpolate (if needed)
  SystemCD->Solve();

  // Get the Mean Value of Derivatives
  SystemCD->getMeanValueDerivatives();

  // Print the Solution to terminal
  SystemCD->printSolution();

  // Generate VTK File for the problem
  SystemCD->generateVTK();

  // Generate plot.dat file to be used for using in Tecplot or GNUplot.
  SystemCD->plotGNU();


CloseFiles();

return 0;
} // end main
