// ======================================================================
// Sine problem 3D
// ======================================================================
// #include <ConvDiff3D.h>

void ExampleFile()
{
  OutPut("Example: Cylinder.h" << endl) ;
}

// exact solution
void Exact(double x, double y, double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(double x, double y, double z, BoundCond &cond)
{

  if ( (fabs(z-10)<1e-8))
   {
    cond = NEUMANN;
   } 
   else
   {
    cond = DIRICHLET;
   }
}

// value of boundary condition
void BoundValue(double x, double y, double z, double &value)
{

   if (fabs(z)<1e-10)
   {
    value =  1. ;
   } 
   else if (fabs(z -10.)<1e-8)
   {
    value = 0;
   }
  else
   {
    value = 0.5;
   }  
}

void BilinearCoeffs(int n_points, double *x, double *y, double *z,
        double **parameters, double **coeffs)
{
  static double eps=1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, *param;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    coeff[0] = eps;
    coeff[1] = 0;
    coeff[2] = 0;
    coeff[3] = 0.25*(1. - (x[i]*x[i] +y[i]*y[i] + z[i]*z[i]));
    coeff[4] = 0;
    coeff[5] = 0;
  }
}


void ReadMeditMesh(char *SMESH, tetgenio &In)
{  
  int i, j, k, dimension, N_FVert, N_Faces, N_Vertices;
  int BDComp_Min=10000000;
  
  char line[100];
  tetgenio::facet *F;
  tetgenio::polygon *P;    

 std::ifstream dat(SMESH);

  if (!dat)
  {
    cerr << "cannot open '" << SMESH << "' for input" << endl;
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
  
  if(dimension!=3)
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
   
  In.numberofpoints = N_Vertices;
  In.pointlist = new double[3*N_Vertices];  

  //read from file
   for(i=0;i<N_Vertices; i++)
    {
     dat.getline (line, 99);
     dat >> In.pointlist[3*i] >> In.pointlist[3*i+1] >> In.pointlist[3*i+2];
//      cout<< i << " vert X: " <<In.pointlist[3*i] << " vert Y: " <<In.pointlist[3*i+1] <<endl;   
    }

  // find the N_Triangles
  // find the N_Vertices
  while (!dat.eof())
  {
    dat >> line;
    
    if ( (!strcmp(line, "Triangles")) ||  (!strcmp(line, "triangles"))   ||  (!strcmp(line, "TRIANGLES"))   ) 
    {
     N_FVert = 3;
     dat.getline (line, 99);
     dat >> N_Faces;
     break;
    }    
    else if ( (!strcmp(line, "Quadrilaterals")) ||  (!strcmp(line, "quadrilaterals"))   ||  (!strcmp(line, "QUADRILATERALS"))   ) 
    {
     N_FVert = 4;
     dat.getline (line, 99);
     dat >> N_Faces;
     break;
    }    
   
    // read until end of line
    dat.getline (line, 99);   
  } 
  cout << " NFaces " << N_Faces << endl;
     
    In.numberoffacets = N_Faces;
    In.facetlist = new tetgenio::facet[In.numberoffacets];
    In.facetmarkerlist = new int[In.numberoffacets];
    
    for(i=0;i<N_Faces; i++)
     {
      dat.getline (line, 99);       
      
      F = &In.facetlist[i];
      tetgenio::init(F);      
      F->numberofpolygons = 1;
      F->polygonlist = new tetgenio::polygon[F->numberofpolygons];
      F->numberofholes = 0;
      F->holelist = NULL;
      P = &F->polygonlist[0];
      tetgenio::init(P);
      P->numberofvertices = N_FVert;
      P->vertexlist = new int[P->numberofvertices];
      
      for(j=0;j<N_FVert;j++)
      {
        dat >> k;
        P->vertexlist[j] = k-1;  // c numbering 
      }
      
      dat >>  In.facetmarkerlist[i];        
      
      if(BDComp_Min > In.facetmarkerlist[i])
        BDComp_Min=In.facetmarkerlist[i];
      
      // cout << i<<  " P->vertexlist[j]:  " <<  In.facetmarkerlist[i]  << endl;  
     } //   for(i=0;i<N_Faces; i++)
   
    BDComp_Min--;
   
    for(i=0;i<N_Faces; i++)
      In.facetmarkerlist[i] -= BDComp_Min;
         
   cout <<    " ReadMeditMesh Done !!! "  << endl;     
     
     
       dat.close();
//    exit(0);    
} // ReadMeditMesh




void TetrameshGen(TDomain *Domain)
{
 //======================================================================
 // Tetgen for grid generation begin
 //======================================================================
  int i, j, k, l, N_Coord, N_FVert, *N_FVerts, *Facets;
  double *Vertices;
  int N, N_RootCells, N_Cells, CurrVertex, N_Vertices, ID, N_Faces, N_G, RefLevel=0;
  int CurrNeib, len1, len2, len3, maxEpV = 0, a, b, c, Neib[2], Neighb_tmp, CurrComp;
  int *Tetrahedrals, *PartMarker, *PointNeighb;
  double *Coordinates, N_x, N_y, N_z; 
  tetgenio In, Out;
  double Xmin = 1e10, Xmax = -1e10, Ymin = 1e10, Ymax = -1e10;  tetgenio::facet *F;
  double Zmin = 1e10, Zmax = -1e10, T[4]={0,0,0,0}, S[4]={0,0,0,0};  tetgenio::polygon *P;  
  double StartX, StartY, StartZ, BoundX, BoundY, BoundZ;  char *SMESH, line[100];
  TBaseCell **CellTree,  **SurfCellTree;  std::ostringstream opts;
  TGridCell **DelCell;
  TVertex **VertexDel, **NewVertices, **NewSurfVertices;  opts << " ";
  TBoundPart *BoundPart;
  TBdPlane **UpdateFaceParams;  SMESH = TDatabase::ParamDB->SMESHFILE;
  TJoint *Joint;
  TBoundComp3D *bdcomp; 
  TCollection *coll, *SurfColl;
  TBaseCell *cell;
  TBdSphere *UpdateParam;


 std::ifstream dat(SMESH);

//   if (!dat)
//   {
//     cerr << "cannot open '" << SMESH << "' for input" << endl;
//     exit(-1);
//   }
// 
//   dat.getline (line, 99);
// 
//   // determine the number of vettices and alocate memory
//   dat >> N_Vertices >> N_Coord;
//   cerr << "N_Vertices: " << N_Vertices << endl;
//   if(N_Coord!=3)  
//    {
//     cerr << "N_Vertices: " << N_Vertices << endl;
//     cerr << "N_Coord must be 3 but it has: " << N_Coord << endl;
//     exit(-1);
//   }
// 
//   Vertices = new double[3*N_Vertices];
// 
//   for(i=0;i<N_Vertices; i++)
//    {
//     dat.getline (line, 99);
//     dat >> j >>Vertices[3*i] >> Vertices[3*i+1] >> Vertices[3*i+2];
// //     cout<< i << " vert X: " <<Vertices[3*i] <<endl;
//    }
// 
//   dat.getline (line, 99);
//   // determine the number of vettices and alocate memory
//   dat >> N_Faces >> N_FVert;
//   cerr << "N_Faces: " << N_Faces << endl;
// //   N_FVert =4;
//   Facets = new int[N_Faces*N_FVert];
//   N_FVerts = new int[N_Faces];
// 
// 
//   
//  if(N_FVert==3)
//   for(i=0;i<N_Faces; i++)
//    {
//     dat.getline (line, 99);
// //     dat >> N_FVerts[i] >>Facets[3*i] >> Facets[3*i+1] >> Facets[3*i+2];
//     dat >> N_FVerts[i] >>Facets[3*i] >> Facets[3*i+1] >> Facets[3*i+2];
// //     cout<< i << " Facets X: " <<Facets[3*i] <<endl;
//    }
//   else if(N_FVert==4)
//    for(i=0;i<N_Faces; i++)
//     {
//      dat.getline (line, 99);
// //     dat >> N_FVerts[i] >>Facets[3*i] >> Facets[3*i+1] >> Facets[3*i+2];
//      dat >> N_FVerts[i] >>Facets[4*i] >> Facets[4*i+1] >> Facets[4*i+2]>> Facets[4*i+3];
// //     cout<< i << " Facets X: " <<Facets[3*i] <<endl;
//     }
//   else
//    {
//     cout<< "Number of face vertices should be 3 or 4, see TetrameshGen(TDomain *Domain) !!!!!!!!! "<< endl;
//     exit(0); 
//    }
// 
//   dat.close();

   /** load the medit mesh file into Tetgen */  
  ReadMeditMesh(SMESH, In);
  

  
//     BoundPart = Domain->GetBdPart(0);

//     N = BoundPart->GetN_BdComps();
//     UpdateFaceParams = new TBdPlane *[N];
//     for(i=0; i<N; i++)
//      {
//       UpdateFaceParams[i] = (TBdPlane*)BoundPart->GetBdComp(i);
//      }
//    cout<<"N"<<N <<endl;

//     UpdateParam = (TBdSphere*)BoundPart->GetBdComp(0);  // sphere domain
//     UpdateParam->SetParams(0.0, 0.0, 0.0, 1.0);

    opts.seekp(std::ios::beg);
    opts<<'p'; // Tetrahedralize the PLC. Switches are chosen to read a PLC (p)
    opts<<"q"<<1.25; // quality mesh generation(q) with a specified quality bound
//     opts<<"a"<<0.1; // maximum volume constraint
//     opts<<'i'; // Inserts a list of additional points into mesh.
    opts<<'z'; // numbers all output items starting from zero
//     opts<<'d'; // Detect intersections of PLC facets.
    opts<<'f'; // Outputs all  faces (including non-boundary) 
    opts<<'e'; // Outputs a list of edges of the triangulation
//     opts<<'I'; // Suppresses mesh iteration numbers.
    opts<<'C'; // Checks the consistency of the final mesh.
//     opts<<'Q'; // Quiet: No terminal output except errors.
//     opts<<'g'; // Outputs mesh to .mesh file for viewing by Medit
    opts<<'Y'; // Suppresses boundary facets/segments splitting
    opts<<'V';  //verbose mode
    opts<<ends;

/*
    In.numberofpoints = N_Vertices;
    In.pointlist = new double[3*In.numberofpoints];
    for(i=0;i<3*N_Vertices; i++)
     In.pointlist[i]= Vertices[i];

    In.numberoffacets = N_Faces;
    In.facetlist = new tetgenio::facet[In.numberoffacets];
    In.facetmarkerlist = new int[In.numberoffacets];
// cout<< " test main 1 " <<endl;
    for(i=0;i<N_Faces; i++)
     {
      F = &In.facetlist[i];
      F->numberofpolygons = 1;
      F->polygonlist = new tetgenio::polygon[F->numberofpolygons];
      F->numberofholes = 0;
      F->holelist = NULL;
      P = &F->polygonlist[0];
      P->numberofvertices = N_FVert;
      P->vertexlist = new int[P->numberofvertices];
      for(j=0;j<N_FVert;j++)
        P->vertexlist[j] = Facets[N_FVert*i + j];
      In.facetmarkerlist[i] = 1;
     }*/
// cout<< " test main " <<endl;

//     for(i=0;i<In.numberofpoints;i++)
//       OutPut(i<<" (x, y, z) =  "<<
//        In.pointlist[3*i]<<' '<<In.pointlist[3*i+1]<<' '<<In.pointlist[3*i+2]<<endl);

 // Calling  tetrahedralize function of 3dtetgen mesh generator
    tetrahedralize((char*)opts.str().c_str(), &In, &Out);

//       cout<<"Test meditmesh" <<endl;
//    exit(0);
   
 //   output: coordinates of all vertices
//  for(i=0;i<Out.numberofpoints;i++)
//   OutPut(" (x, y, z) =  "<< Out.pointlist[3*i]<<' '<<Out.pointlist[3*i+1]<<' '<<Out.pointlist[3*i+2]<<endl);

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
       N_Faces = cell->GetN_Faces();
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


   N_RootCells = Out.numberoftetrahedra;


  // allocate auxillary fields
   Coordinates = Out.pointlist;
   Tetrahedrals = Out.tetrahedronlist;

  // generate new vertices
   N_G = Out.numberofpoints;
   NewVertices = new TVertex*[N_G];

   for (i=0;i<N_G;i++)
    {
      NewVertices[i] = new TVertex(Coordinates[3*i], Coordinates[3*i+1], Coordinates[3*i+2]);

      // set bounding box
      if (Coordinates[3*i] > Xmax) Xmax = Coordinates[3*i];
      if (Coordinates[3*i] < Xmin) Xmin = Coordinates[3*i];
      if (Coordinates[3*i+1] > Ymax) Ymax = Coordinates[3*i+1];
      if (Coordinates[3*i+1] < Ymin) Ymin = Coordinates[3*i+1];
      if (Coordinates[3*i+2] > Zmax) Zmax = Coordinates[3*i+2];
      if (Coordinates[3*i+2] < Zmin) Zmin = Coordinates[3*i+2];
   }

   // set bounding box
    StartX = Xmin;
    StartY = Ymin;
    StartZ = Zmin;
    BoundX = Xmax - Xmin;
    BoundY = Ymax - Ymin;
    BoundZ = Zmax - Zmin;


   Domain->SetBoundBox(StartX, StartY, StartZ, BoundX, BoundY, BoundZ);
//        cout<<Xmin <<"  "<<Ymin <<"  "<<Zmin<<endl;
//        cout<<Xmax <<"  "<<Ymax <<"  "<<Zmax<<endl;

   CellTree = new TBaseCell*[N_RootCells];
 //  output of each tetraheron vertex indices (four vertices for each)
//   for (i=0;i<N_RootCells;i++)
//        cout<< Tetrahedrals[4*i]<<"  "<<Tetrahedrals[4*i + 1]<<"  "
//          <<Tetrahedrals[4*i + 2]<<"  "<<Tetrahedrals[4*i + 3]<<endl;

   for (i=0;i<N_RootCells;i++)
   {
     CellTree[i] = new TMacroCell(TDatabase::RefDescDB[Tetrahedron],
                                    RefLevel);

     CellTree[i]->SetVertex(0, NewVertices[Tetrahedrals[4*i    ]]);
     CellTree[i]->SetVertex(1, NewVertices[Tetrahedrals[4*i + 1]]);
     CellTree[i]->SetVertex(2, NewVertices[Tetrahedrals[4*i + 2]]);
     CellTree[i]->SetVertex(3, NewVertices[Tetrahedrals[4*i + 3]]);

    CellTree[i]->SetClipBoard(i);
     ((TMacroCell *) CellTree[i])->SetSubGridID(0);
   }

   Domain->SetTreeInfo(CellTree, N_RootCells);


   // initialize iterators
   TDatabase::IteratorDB[It_EQ]->SetParam(Domain);
   TDatabase::IteratorDB[It_LE]->SetParam(Domain);
   TDatabase::IteratorDB[It_Finest]->SetParam(Domain);
   TDatabase::IteratorDB[It_Between]->SetParam(Domain);
   TDatabase::IteratorDB[It_OCAF]->SetParam(Domain);


//     search neighbours
   N_G = Out.numberofpoints;
   PointNeighb = new int[N_G];
   cout<<"numberofpoints "<<N_G<<endl;
   memset(PointNeighb, 0, N_G *SizeOfInt);

     for (i=0;i<4*N_RootCells;i++)
     PointNeighb[Tetrahedrals[i]]++;

   for (i=0;i<N_G;i++)
     if (PointNeighb[i] > maxEpV) maxEpV = PointNeighb[i];
   delete [] PointNeighb;

   cout<<"maxEpV "<< maxEpV<<endl;

   PointNeighb = new int[++maxEpV * N_G];

   memset(PointNeighb, 0, maxEpV*N_G*SizeOfInt);

    // every vertex contains "maxEpV" columns
    // for every vertex at first colomn contains the number of cells containing this vertex
    // at further columns we set the index of corresponding cells containing this vertex
    // cout<<"maxEpV*N_G "<<maxEpV*N_G<<endl;

   for(i=0;i<4*N_RootCells;i++)
    {
     j = Tetrahedrals[i]*maxEpV;
     PointNeighb[j]++;
     //cout<<"j + PointNeighb[j] " << j <<endl;
     PointNeighb[j + PointNeighb[j]] = i / 4;
    }
 //  output of PointNeighb columns for each point
//   for (i=0;i<N_G;i++)
//    {
//     for (j=0;j<maxEpV;j++)
//     cout<<"  "<< PointNeighb[i*maxEpV+j];
//     cout<<endl;
//    }

   // generate new faces 
   N_G =  Out.numberoftrifaces;
   cout<<"numberoftrifaces "<<N_G<<endl;
   for (i=0;i<N_G;i++)
   {
     a = Out.trifacelist[3*i];
     b = Out.trifacelist[3*i+1];
     c = Out.trifacelist[3*i+2];

//      cout<<"  "<< a<<"  "<< b<<"  "<< c<<endl;

     Neib[0] = -1;
     Neib[1] = -1;
     CurrNeib = 0;

     len1 = PointNeighb[a*maxEpV];
     len2 = PointNeighb[b*maxEpV];
     len3 = PointNeighb[c*maxEpV];

   // find the index of the cells containing current face with point indices a,b,c
    for (j=1;j<=len1;j++)
     {
       Neighb_tmp = PointNeighb[a*maxEpV + j];
        for (k=1;k<=len2;k++)
         {
          if (Neighb_tmp == PointNeighb[b*maxEpV + k])
           {
            for (l=1;l<=len3;l++)
             if (Neighb_tmp == PointNeighb[c*maxEpV + l])
             {
              Neib[CurrNeib++] = Neighb_tmp;
              break;
             }
           }
          }
       if (CurrNeib == 2) break;
     }
 //   cout<<"CurrNeib " << CurrNeib <<endl;
// cout<<"Out.trifacemarkerlist[i] : "<<Out.trifacemarkerlist[i]<<endl;
     if (Out.trifacemarkerlist[i]) // 0 for inner edges and Boundcomp+1 for Boundedge respect
      {

       CurrComp = Out.trifacemarkerlist[i] - 1;
//        cout<<"Boundary face CurrComp: "<<CurrComp<<endl;

       bdcomp = Domain->GetBdPart(0)->GetBdComp(CurrComp);

       if(bdcomp->GetTSofXYZ(NewVertices[a]->GetX(), NewVertices[a]->GetY(),
                             NewVertices[a]->GetY(), T[1], S[1]) ||
          bdcomp->GetTSofXYZ(NewVertices[b]->GetX(), NewVertices[b]->GetY(),
                             NewVertices[b]->GetY(), T[2], S[2]) ||
          bdcomp->GetTSofXYZ(NewVertices[c]->GetX(), NewVertices[c]->GetY(),
                             NewVertices[c]->GetY(), T[3], S[3])    )
         {
          cerr<<"Error: could not set parameter values"<<endl;
          OutPut(NewVertices[a]<<endl);
          OutPut(NewVertices[b]<<endl);
          OutPut(NewVertices[c]<<endl);
          exit(0);
         }

       if (CurrNeib == 2)
        {
         if(bdcomp->IsFreeBoundary())
           Joint = new TIsoInterfaceJoint3D(bdcomp, T, S,
                       CellTree[Neib[0]], CellTree[Neib[1]] );
          else
           Joint = new TInterfaceJoint3D(bdcomp, T, S, 
                       CellTree[Neib[0]], CellTree[Neib[1]] );
        }
       else
        {
         if(bdcomp->IsFreeBoundary())
           Joint = new TIsoBoundFace(bdcomp, T, S);
         else
           Joint = new TBoundFace(bdcomp, T, S);
        }
      }
     else
      {
// //       cout<<"Inner face"<<endl;
       if (CurrNeib != 2)
       cerr << "Error !!!!!!!! not enough neighbours!" << endl;

       Joint = new TJointEqN(CellTree[Neib[0]], CellTree[Neib[1]]);

      }

    // First element containing the current face
    // find the local index for the point 'a' on the cell
    for (j=0;j<4;j++)
      if (Tetrahedrals[4*Neib[0]+j] == a) break;

    // find the local index for the point 'b' on the cell
    for (k=0;k<4;k++)
      if (Tetrahedrals[4*Neib[0]+k] == b) break;

       // find the local index for the point 'c' on the cell
    for (l=0;l<4;l++)
      if (Tetrahedrals[4*Neib[0]+l] == c) break;   

     l = l*100 + k*10 + j;  

//      cout<<""<< l <<endl;

     switch (l) // j will contain the local index for the current face
      {
        case 210: case 21: case 102:
        case 120: case 12: case 201:
          j = 0;
          break;  
        case 310: case 31: case 103:
        case 130: case 13: case 301:
          j = 1;
          break;  
        case 321: case 132: case 213:
        case 231: case 123: case 312:
          j = 2;
          break;  
        case 230: case 23: case 302:
        case 320: case 32: case 203:
          j = 3;
          break; 

      default:
       Error("Unable to set the face !!!!!!!!!!!!" << endl);
       exit(0);
     }
      CellTree[Neib[0]]->SetJoint(j, Joint);

   if (Neib[1] != -1) // second element containing the current face
    {
          // find the local index for the point 'a' on the cell
    for (j=0;j<4;j++)
      if (Tetrahedrals[4*Neib[1]+j] == a) break;

    // find the local index for the point 'b' on the cell
    for (k=0;k<4;k++)
      if (Tetrahedrals[4*Neib[1]+k] == b) break;

       // find the local index for the point 'c' on the cell
    for (l=0;l<4;l++)
      if (Tetrahedrals[4*Neib[1]+l] == c) break;   

     l = l*100 + k*10 + j;  

//      cout<<""<< l <<endl;

     switch (l) // j will contain the local index for the current face
      {
        case 210: case 21: case 102:
        case 120: case 12: case 201:
          j = 0;
          break;  
        case 310: case 31: case 103:
        case 130: case 13: case 301:
          j = 1;
          break;  
        case 321: case 132: case 213:
        case 231: case 123: case 312:
          j = 2;
          break;  
        case 230: case 23: case 302:
        case 320: case 32: case 203:
          j = 3;
          break; 

      default:
       Error("Unable to set the face !!!!!!!!!!!!" << endl);
       exit(0);
      }
      CellTree[Neib[1]]->SetJoint(j, Joint);
     }

  if (Joint->GetType() == InterfaceJoint3D ||
      Joint->GetType() == IsoInterfaceJoint3D)
      {
        ((TInterfaceJoint3D*)Joint)->SetMapType();
        ((TInterfaceJoint3D*)(Joint))->CheckOrientation();
      }
      else 
       if (Joint->GetType() == JointEqN)
           Joint->SetMapType();

  }

  delete [] PointNeighb;


cout<<"Out.numberofpoints "<< Out.numberofpoints <<endl;
// cout<<"cout main tetra "<< N_Faces <<endl;
// exit(0);
 //======================================================================
 // Tetgen for grid generation end
 //======================================================================


delete [] Vertices;

}



