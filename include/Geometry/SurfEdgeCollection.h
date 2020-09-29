#ifndef __SURFEDGECOLLECTION__
#define __SURFEDGECOLLECTION__

#include <Collection.h>

// #define iFLIP(x,y) if ( x > y ) { int ivar; ivar = a; a = b; b = ivar; }

class TSurfEdgeCollection
{
  protected:
    /** number of vertices */
    int mN_Vertices;
    
    /** number of edges */
    int mN_Edges;
    
    /** hash table for edges */
    int **mEdgeHash;
  
  protected:
    void Init(TCollection *Coll, int N_SurfaceJoints,
	      int *CellNumbers, int *JointNumbers);
	      
    inline int Flip(int &a, int &b)
    {
      if ( a > b )
      { 
	int ivar;
	ivar = a; a = b; b = ivar;
      }
      
      return a+b;
    }
    
    void IncreaseBucket(int a, int b, int &counter);
    int FindEdge(int *Bucket, int a, int b);
    
    
  public:
    TSurfEdgeCollection (TCollection *Coll, int N_SurfaceJoints,
			 int *CellNumbers, int *JointNumbers);
			 
    ~TSurfEdgeCollection ();
    
    int GetN_Edges()
    { return mN_Edges; }
    
    int GetN_Vertices()
    { return mN_Vertices; }
    
    int GetEdge(int a, int b)
    { 
      if ( a < 0 || b < 0 ) return -1;
      else return FindEdge(mEdgeHash[Flip(a,b)], a, b);
    }
};
#endif
