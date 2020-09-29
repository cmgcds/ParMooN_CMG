#ifndef __ISOJOINT__
#define __ISOJOINT__

#include <Vertex.h>

class TIsoJoint
{
  public:
    
    virtual int GetN_Vertices() = 0;
    virtual TVertex **GetVertices() = 0;
    
    virtual void SetVertices(int n_vertices, TVertex **vertices) = 0; 
};
#endif
