/** =======================================================================
* @class     TSystemRTE
* @brief     stores the information of ADI system 
* @author    Sashikumaar Ganesan
* @date      13.04.2020
* @History 
* ======================================================================= */
#include <SystemADI.h>

TSystemADI::TSystemADI()
{
  
}

TSystemADI::TSystemADI(DoubleFunctND *growthAndNuc)
{
  GetGrowthAndNuc = growthAndNuc;
}



TSystemADI::~TSystemADI()
{

}