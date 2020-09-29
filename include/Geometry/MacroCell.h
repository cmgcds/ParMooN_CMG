/** ************************************************************************ 
*
* @class    TMacroCell
* @brief    represent a unit of the macro grid
* @author   Volker Behns  
* @date     15.03.98
************************************************************************  */

#ifndef __MACROCELL__
#define __MACROCELL__

#include <GridCell.h>

 /** @brief represent a unit of the macro grid */
class TMacroCell : public TGridCell
{
  protected:
     /** @brief ID of subgrid this macro belongs to */
    int SubGridID;

  public:
    // Constructor
     /** @brief initialize the refinement descriptor to refdesc */
    TMacroCell(TRefDesc *refdesc, int reflevel);

    // Methods
     /** @brief set subgrid ID */
    void SetSubGridID(int subgridid)
    { SubGridID = subgridid; }

     /** @brief return subgrid ID */
    virtual int GetSubGridID()
    { return SubGridID; }
};

#endif
