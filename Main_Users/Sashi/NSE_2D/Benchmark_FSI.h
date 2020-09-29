// Navier-Stokes problem, Benchmark problem
// 
// u(x,y) = unknown
// p(x,y) = unknown

void ExampleFile()
{
  OutPut("Example: Benchmark_FSI.h" << endl) ;
}

#define __BENCH__

#include <Joint.h>
#include <BoundEdge.h>
#include <BoundComp.h>
#include <FE2D.h>
#include <FEDesc2D.h>

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void ExactU2(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void ExactP(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int i, double t, BoundCond &cond)
{
  if (i==1)
  {
    cond = NEUMANN;
    TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
  }
  else
    cond = DIRICHLET;
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0: value = 0;
            break;
    case 1: value= 0;
            break;
    case 2: value = 0;
            break;
    case 3: value=1.2*Param*(1-Param); // 4*0.3
            break;
    case 4: value=0;
            break;
    default: cout << "wrong boundary part number: " << BdComp << endl;
  }  
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  value = 0;
  if(BdComp>4) cout << "wrong boundary part number: " << BdComp << endl;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *x, double *y,
               double **parameters, double **coeffs)
{
  double eps = 1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = 0; // f1
    coeff[2] = 0; // f2
  }
}


void DeleteDomain(TDomain *&Domain)
{
  int i, j, k, N_Cells, N_RootCells;
  int CurrVertex, N_Faces, N_Vertices, ID;
  
  TBaseCell **CellTree,  **SurfCellTree, *cell;
  TGridCell **DelCell;
  TVertex **VertexDel;  
  TCollection *coll;
  
  
    Domain->GetTreeInfo(CellTree,N_RootCells);
    coll = Domain->GetCollection(It_Finest, 0);
    N_Cells = coll->GetN_Cells();  

    
 //     cout<<"N_RootCells: "<<N_RootCells<<endl;
   // remove all existing vertices and joints
    VertexDel = new TVertex*[8*N_RootCells];
    CurrVertex = 0;   
    
   for(i=0;i<N_Cells;i++)
     {
       cell = coll->GetCell(i);
       N_Faces = cell->GetN_Joints();
       N_Vertices = cell->GetN_Vertices();
       for(j=0;j<N_Faces;j++)
         {
          if(CurrVertex==0)
           {
               VertexDel[CurrVertex++] = cell->GetVertex(j);
            }
           else
            {
             ID = 0;
             for(k=0;k<CurrVertex;k++)
             if(VertexDel[k]==cell->GetVertex(j))
              {
               ID = 1; break;
              }
             if(ID!=1)
              {
               VertexDel[CurrVertex++] = cell->GetVertex(j);
              }
          } // else if(CurrVertex==0)

           ID = 0;
           for(k=0;k<CurrVertex;k++)
           if(VertexDel[k]==cell->GetVertex((j+1)%N_Vertices))
            {
             ID = 1; break;
            }
             if(ID!=1)
            {
             VertexDel[CurrVertex++] = cell->GetVertex((j+1)%N_Vertices);
            }
 
           ID = 0;
           for(k=0;k<CurrVertex;k++)
           if(VertexDel[k]==cell->GetVertex((j+2)%N_Vertices))
            {
             ID = 1; break;
            }
             if(ID!=1)
            {
             VertexDel[CurrVertex++] = cell->GetVertex((j+2)%N_Vertices);
            }
         if(N_Faces==6) // If old cell is hexahedrol
          {  
           ID = 0;
           for(k=0;k<CurrVertex;k++)
           if(VertexDel[k]==cell->GetVertex((j+4)%N_Vertices))
            {
             ID = 1; break;
            }
             if(ID!=1)
            {
             VertexDel[CurrVertex++] = cell->GetVertex((j+4)%N_Vertices);
            }        
        }
      } // for j
    } // for i

    for(i=0;i<CurrVertex;i++)
     delete VertexDel[i];
     delete [] VertexDel;
     OutPut(CurrVertex<<" vertices were deleted"<<endl);

  // remove all existing cells and joints

   for(i=0;i<N_RootCells;i++)
    delete (TGridCell*)CellTree[i];
    delete [] CellTree;
    OutPut(N_RootCells<<" cells were deleted"<<endl);    
    
    
/*    
    cout << "DeleteDomain" <<endl;
    exit(0);
    */
}


void MEditMeshCreate(TDomain *&Domain)
{
  int i, j, k, l, dimension, N_Vertices;
  int N_FVert, N_Faces, *Facemarkerlist, *Facelist, N_RootCells;
  int v1, v2, v3, v4, CellMarker, RefLevel=0, *Triangles, *PointNeighb, maxEpV, *Triangles_loc;
  int a, b, c, Neib[2], Neighb_tmp, CurrNeib, len1, len2, len3, CurrComp, N_Points ;
  int N_Cells, MaxLen, jj, N_SourcePts;  
  const int *EdgeVertex, *TmpFV, *TmpLen;
  
  double X, Y, Z, DispX, DispY, DispZ;
  double Xmin = 1e10, Xmax = -1e10, Ymin = 1e10, Ymax = -1e10;  
  double Zmin = 1e10, Zmax = -1e10, T[4]={0,0,0,0}, S[4]={0,0,0,0}; 
  double StartX, StartY, StartZ, BoundX, BoundY, BoundZ; 

  TVertex **NewVertices;
  TBaseCell **CellTree, *cell;  
  TJoint *Joint;  
  TBoundComp3D *bdcomp; 
  TShapeDesc *ShapeDesc;
  TCollection *coll;  
  
  char *MEDITMESH, line[100];
  MEDITMESH = TDatabase::ParamDB->SMESHFILE;
 
  
  
  /** delete all cells, vertices in the old domain */
  DeleteDomain(Domain);
  
    cout << " SMESHFILE is " << MEDITMESH << endl;
  /** load the tet mesh file */  
  std::ifstream dat(MEDITMESH);

  if (!dat)
  {
    cerr << "cannot open '" << MEDITMESH << "' for input" << endl;
    exit(-1);
  }     
    
  // check the dimension
  while (!dat.eof())
  {
    dat >> line;
    
    if ( (!strcmp(line, "Dimension"))  ||  (!strcmp(line, "dimension")) ||  (!strcmp(line, "DIMENSION")))
    {
     dat.getline (line, 99);
     dat >> dimension;
     break;
    }    
    // read until end of line
    dat.getline (line, 99);   
  }

  dimension--;
  if(dimension!=2)
   {
    cerr << "dimension: " << dimension << endl;
    exit(-1);
   }
  
  // find the N_Vertices
  while (!dat.eof())
  {
    dat >> line;
    
    if ( (!strcmp(line, "Vertices")) ||  (!strcmp(line, "vertices"))   ||  (!strcmp(line, "VERTICES"))   ) 
    {
     dat.getline (line, 99);
     dat >> N_Vertices;
     break;
    }    
    // read until end of line
    dat.getline (line, 99);   
  } 
  
  
  // generate new vertices
  NewVertices = new TVertex*[N_Vertices];
     
  
  //read from file
   for(i=0;i<N_Vertices; i++)
    {
     dat.getline (line, 99);
     dat >> X >> Y >> Z;   
   
     if(Z!=0)
      {
       cerr << "Z: " << Z << endl;
       exit(-1);
      }
        
     NewVertices[i] = new TVertex(X, Y);

//       // Pt source is at (291.977 -1.87333 159.89 )
//      if(i==98 || i==123 || i==122 )      
//       cout<< i << " vert X: " <<X << " vert Y: " <<Y <<" vert Z: " <<Z <<endl; 
//           
     
      // set bounding box
      if (X > Xmax) Xmax = X;
      if (X < Xmin) Xmin = X;
      if (Y > Ymax) Ymax = Y;
      if (Y < Ymin) Ymin = Y;
//       if (Z > Zmax) Zmax = Z;
//       if (Z < Zmin) Zmin = Z;
   
     cout<< i << " vert X: " <<X << " vert Y: " <<Y <<" vert Z: " <<Z <<endl;   
    } // for(i=0;i<N_Vertices; i++)
    
   // set bounding box
    StartX = Xmin;
    StartY = Ymin;

    BoundX = Xmax - Xmin;
    BoundY = Ymax - Ymin;
    
    Domain->SetBoundBox(BoundX,BoundY);
    Domain->SetBoundBoxstart(StartX,StartY);   
    N_RootCells=0;
     
  i=0;
   // find the N_RootCells
   while (!dat.eof())
   {
    dat >> line;
      cout << "N_RootCells  " << i++ << " " << line <<endl;
    if ( (!strcmp(line, "Edges")) ||  (!strcmp(line, "triangles"))   ||  (!strcmp(line, "TRIANGLES"))   ) 
    {
     dat.getline (line, 99);
     dat >> N_RootCells;
     
//      cout << "N_RootCells  " << N_RootCells <<endl;
     break; 
    }    
    // read until end of line
    dat.getline (line, 99);   
    
    exit(0);
    
   }   
     
    
  // generate new cells
//    CellTree = new TBaseCell*[N_RootCells];
//    Triangles = new int[3*N_RootCells];
       
    
    
    
    
  cout <<"test MEditMeshCreate "<< endl;
  exit(0);
  
  
} // void MEditMeshCreate(TD
 
 
 
 
 
 