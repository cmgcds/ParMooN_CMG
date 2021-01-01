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

#include <ANN.h>
// main program

int main(int argc, char* argv[])
{

  TANNParamReader paramReader(argv[1]);

  TANN ann(&paramReader);

  for (int i=0; i<paramReader.nHL+2; i++){
    std::cout << ann.layers[i].rank << "   " << ann.layers[i].dim << "   " << ann.layers[i].type << std::endl; 
  };
return 0;
} // end main
