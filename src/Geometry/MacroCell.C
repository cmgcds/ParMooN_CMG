
#include <MacroCell.h>

// Constructors
TMacroCell::TMacroCell(TRefDesc *refdesc, int reflevel) :
              TGridCell(refdesc, reflevel)
{
  SubGridID = 0;
}

// Methods
