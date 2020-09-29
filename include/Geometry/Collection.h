/** ************************************************************************ 
*
* @class TCollection 
* @date  14.10.97
* @brief store cells in an array
* @author Gunar Matthies & Sashikumaar Ganesan
* @History: MPI methods (Sashikumaar Ganesan, 08.08.14)
   
****************************************************************************/

#ifndef __COLLECTION__
#define __COLLECTION__

#include <BaseCell.h>
#include <JointCollection.h>

/** @brief store cells in an array, used by cell iterators */
class TCollection
{
  protected:
    /** @brief number of cells stored */
    int N_Cells;

    /** @brief array containing the pointers to the cells */
    TBaseCell **Cells;

    /** @brief array with all cells sorted by pointer */
    TBaseCell **SortedCells;

    /** @brief array with index of SortedCells in Cells */
    int *Index;

    #ifdef  _MPI
    /** @brief Number of own cells (excluding Halo cells) */
    int N_OwnCells;

    /** @brief array for Globalcell number in Cells */
    int *GlobalIndex;
    #endif

  public:
    /** @brief constructor */
    TCollection(int n_cells, TBaseCell **cells);

    /** @brief return number of cells */
    int GetN_Cells() const
    { return N_Cells; }

    /** @brief return Cell with index i in Cells-array */
    TBaseCell *GetCell(int i) const
    { return Cells[i]; }


  // THIVIN
  // Get volume of the Current Cell COllection
  double GetVolume()
  {
    int N_cells = this->GetN_Cells();
    double volume = 0.;
    for ( int cellNo = 0 ; cellNo < N_cells ; cellNo++)
    {
      TBaseCell* currentCell;
      currentCell = this->GetCell(cellNo);
      volume += currentCell->GetMeasure();

    }

    return volume;

  }

    // THIVIN
  // Get volume of the Current Cell COllection
  double GetminVolumeCollection()
  {
    int N_cells = this->GetN_Cells();
    double minvolume = 99999999999999999;
    double volume = 0.;
    for ( int cellNo = 0 ; cellNo < N_cells ; cellNo++)
    {
      TBaseCell* currentCell;
      currentCell = this->GetCell(cellNo);
      volume = currentCell->GetMeasure();
      if (volume < minvolume)
        minvolume = volume ;

    }

    return minvolume;

  }



  //THIVIN 
  // Get the shortest Edge length in the Collection
   double GetShortestEdgeinCollection()
  {
    int N_cells = this->GetN_Cells();
    double minshortedge = 99999999999999999;
    double shortedge = 0.;
    for ( int cellNo = 0 ; cellNo < N_cells ; cellNo++)
    {
      TBaseCell* currentCell;
      currentCell = this->GetCell(cellNo);
      shortedge = currentCell->GetDiameter();
      if(shortedge < minshortedge)   minshortedge = shortedge;

    }
    return minshortedge;

  }

    /** @brief destructor: delete arrays */
    ~TCollection();

    /** @brief get maximal and minimal diameter */
    int GetHminHmax(double *hmin, double *hmax);

    /** @brief return Index of cell in Cells-array */
    int GetIndex(TBaseCell *cell);

    /** @brief mark the vertices that are on the boundary */
    int MarkBoundaryVertices();

    /** @brief return Index of joints in Cells-array */
    TJointCollection  *GetJointCollection();

    /** @brief Generate   Vertex Neibs for all cells in the collection */
    void GenerateCellVertNeibs();

    /** @brief return the Index of the vertex in the sorted array */
    int GetIndex(TVertex **Array, int Length, TVertex *Element);
    
#ifdef  _MPI
    void SetN_OwnCells(int n_OwnCells)
     { N_OwnCells = n_OwnCells; }

    int GetN_OwnCells()
     { return N_OwnCells; }

    int GetN_HaloCells()
     { return (N_Cells - N_OwnCells); }

    int *GetGlobalIndex()
     {
      return GlobalIndex;
     }
#endif

   void Replace_Coll(int n_cells, TBaseCell **cells)
     {
      N_Cells = n_cells;
      Cells = cells;
     }

  private:
    /** @brief provide additional arrays */
    void GenerateSortedArrays();

    /** @brief return Index of cell in SortedCells-array */
    int GetSortedIndex(TBaseCell *cell);

    // THIVIN 
    // Get Volume of the collection 



};

#endif
