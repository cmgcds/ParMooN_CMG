#ifndef __OCTREE__
#define __OCTREE__

#include <Collection.h>
#include <BaseCell.h>

class TOctTree
{
  public:
    typedef struct {
      double StartX;
      double StartY;
      double StartZ;
      double BoundX;
      double BoundY;
      double BoundZ;
    } TBoundBox;
  
  protected:   
    class TNode
    {
      protected:
      	/// parameters describing box
	TBoundBox Box;
	
	/// parent
	TNode *Parent;
	
	/// sub boxes
	TNode *Childs[8];
	
	///
	int Depth;
	
	int N_Cells;
	TBaseCell **Cells;

      protected:
	void MakeSubBoxes();
	bool Intersect(TBaseCell *Cell);
	bool PointInBox(double x, double y, double z);
		
      public:
	TNode (TNode *Parent, TBoundBox *box, int n_cells, TBaseCell **cells);	
	~TNode ();
	
	TNode *GetLeaf(double x, double y, double z);
	void GetCells(int &n_cells, TBaseCell **&cells);
    };
 
  protected:
    
    /// Collection from which the octtree is build
    TCollection *Coll;
    
    /// head of the tree
    TNode *Head;
    
    /// bounding box for whole domain
    TBoundBox Bound;
    
  protected:
    void BuildTree();
    
  public:
    TOctTree(TCollection *coll, TBoundBox *bounds); 
    ~TOctTree();
    
    void GetCells(double x, double y, double z, int &N_Cells, TBaseCell **&Cells);
};

#endif
