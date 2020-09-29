// =======================================================================
// @(#)It_Mortar.h        1.2 02/08/99
// 
// Class:       TIt_Mortar
// Purpose:     iterator which gets all cells on a given mortar face
//
// Author:      Volker Behns  20.03.98
//
// =======================================================================

#ifndef __IT_MORTAR__
#define __IT_MORTAR__

#include <It_Search.h>

/** iterator which gets all cells on a given mortar face */
class TIt_Mortar : public TIt_Search
{
  protected:
    /** local number of mortar face on each level */
    int LocMortarFace[MAX_ItLevel];

    /** indicator for mortar side */
    bool MortarSide;

  public:
    // Constructors

    // Methods
    /** initialize the iterator */
    virtual int Init(int level);

    /** return the maximum level */
    virtual int GetMaxLevel();

    /** return the next cell */
    virtual TBaseCell *Next(int &info);

    /** return the previous cell */
    virtual TBaseCell *Prev()
    { return NULL; }

    /** return number of moratr faces */
    int GetN_MortarFace()
    { return Domain->GetN_MortarFace(); }

    /** return coordinates of 'start'-point of mortar edge */
    void GetPoint(double &X, double &Y)
    {
      TVertex *Vert = ActiveCell->GetVertex(LocMortarFace[ActiveLevel]);
      X = Vert->GetX();
      Y = Vert->GetY();
    }
};

#endif
