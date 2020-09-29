/** ************************************************************************ 
*
* @class     TMesh
* @brief     stores the information of a mesh 
* @author    Sashikumaar Ganesan, 
* @date      08.08.14
* @History 
 ************************************************************************  */

#ifndef __MESH__
#define __MESH__

#include <BaseCell.h>
#include <Vertex.h>
#include <Joint.h>

/**  @brief mesh details */
class TMesh : 
{
  protected:
    /**  @brief number of vertices in the mesh */    
    int N_RootVertices;
    
    /**  @brief number of joints (2D:edge, 3D:Face) in the mesh */    
    int N_Joints;    
    
    /**  @brief number of cells in the mesh */    
    int N_Cells;   

     /**  @brief cell-vetrtices index */    
    int *CellVertices;
    
    /**  @brief cell-joints (2D:edge, 3D:Face) in the mesh */    
    int *CellJoints;      
    
    /**  @brief array of pointers to vertices in the mesh  */
    TVertex  **Vertices;
 
    /**  @brief array of pointers to joints in the mesh  */
    TJoint **Joints;
    
    /**  @brief array of pointers to cells in the mesh */
    TBaseCell **CellTree;
 
    /**  @brief grid level on with this cell was generated */
    int RefLevel;

  public:
    // Constructor
    TMesh();

    TMesh(int N_RootVertices, int N_Joints, int N_Cells, int *CellVertices, int *CellJoints, TVertex  **Vertices, TJoint **Joints, TBaseCell **CellTree);
	
    // Destructor
    ~TMesh();

    // Methods
 

#ifdef __2D__
 
#else
 
 
#endif 
 
};

#endif
