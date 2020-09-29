// =======================================================================
// 
// Class:       TVector
//
// Purpose:     emulates a dynamic data structure for storing data
//              of arbitrary type (template)
//
// Author:      Gunar Matthies
// Version:     1.0
//
// =======================================================================

#ifndef __VECTOR__
#define __VECTOR__

#ifdef __AIX__
#pragma implementation("../../src/General/Vector.C")
#endif

#include <Constants.h>
#include <AllClasses.h>

template <class Type>
class TVector
{
  protected:
    /** number of all elements */
    int N_Elements;

    /** initial size */
    int init_size;
    
    /** increment */
    int init_incr;

    /** number of lists used */
    int NumberOfLists;

    /** number of free lists now */
    int NumberOfFreeLists;

    /** number of free entries in current list */
    int NumberOfFreeEntries;

    /** list for next entry */
    int ListForNext;
    
    /** index for next entry */
    int IndexForNext;

    /** array of arrays containing the data */
    Type **Lists;

    /** if there is no space left in the array Lists, create a new array
        which is IncrForListNumbers entries longer than the old */
    int IncrForListNumbers;

  public:
    /** constructor */
    TVector(int i_size, int i_incr);

    /** destructor = free all allocated memory */
    ~TVector();

    /** return the number of elements strored */
    int GetN_Elements()
    { return N_Elements; }

    /** return the element i */
    Type GetElement(int i);

    /** set the value of the already existing element i to value */
    void SetElement(int i, Type value);

    /** add a new element at the end */
    int AddElement(Type value);

};


#endif
