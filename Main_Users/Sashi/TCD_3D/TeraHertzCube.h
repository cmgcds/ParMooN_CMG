// ==========================================================================
// instationary problem
// ==========================================================================

//===========================================================================
// example file
// =========================================================================
// exact solution in unit cube
void ExampleFile()
{
#define __ROBINBC__ 

  OutPut("Example: TeraHertzCube.h" << endl);
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

// initial conditon
void InitialCondition(double x, double y, double z, double *values)
{
 values[0] = TDatabase::ParamDB->P0;

}

// kind of boundary condition
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
//   cond = ROBIN;
  double beam_r = TDatabase::ParamDB->P11;
  double scale = 2.5e-4/beam_r;
  scale = 1.;
  
if(fabs(z+(60*scale))<1e-12 )
 { cond = DIRICHLET; }
 else
 { cond =  NEUMANN;  }

  
}

// value of boundary condition
void BoundValue(double x, double y, double z, double &value)
{
  value = 0;
}

/*
void BilinearCoeffs(int n_points, double *X, double *Y, double *Z,
        double **parameters, double **coeffs)
{
  double eps=0;
  int i, Region_ID;
  double *coeff;
  double x, y, z;
 
  
  double rhsfact, alpha,  char_L, beam_r, DimlessBeam_r;
  double xp=0., yp=0., zp=0., SourceCoord, Sourceradius;
    
//    if(TDatabase::ParamDB->P13==0)
//     { //Lside    //   3.606596999999999830777142	93.27791100000000312775228	78.19347199999999986630428
//       yp=93.27791100000000312775228;
//       zp=78.19347199999999986630428; 
//     } 
//    else if(TDatabase::ParamDB->P13==1)
//     { //Back side  83.11056343315070193966676	12.29866318748098663604651	97.409675520607805765394
//       xp=83.11056343315070193966676;
//       zp=97.409675520607805765394; 
//     }   
//    else if(TDatabase::ParamDB->P13==2) //  Front side
//     { //Front side  91.21127599999999802093953 214.3263880000000028758222 85.10373400000000287946023 
//       xp=91.21127599999999802093953;
//       zp=85.10373400000000287946023; 
//     }  
//    else // default is right side    
//      { //Rside    // 177.8944370000000105846993	 95.85110799999999642295734	81.84575200000000450017978
//       yp=95.85110799999999642295734;
//       zp=81.84575200000000450017978; 
//     }  
    
 
  beam_r = TDatabase::ParamDB->P11;
  char_L = TDatabase::ParamDB->P12;

  Region_ID = coeffs[0][1];
  DimlessBeam_r =  beam_r/char_L;
 
  // cout << "Region_ID = " << Region_ID << endl;
   
   if(Region_ID==1) //1-scalp and skull 
    {
     eps = TDatabase::ParamDB->P1/(1600.0*2000*char_L*char_L); ;   
     alpha = TDatabase::ParamDB->P7;     
     rhsfact = (alpha*TDatabase::ParamDB->P4)/(1600.0*2000*(22./7.)*beam_r*beam_r);   
    }
   else if(Region_ID==2)  
    {
      //2-CSF  
     eps = TDatabase::ParamDB->P2/(1030.0*3710*char_L*char_L); 
     alpha = TDatabase::ParamDB->P8;        
     rhsfact = alpha*TDatabase::ParamDB->P4/(1030.0*3710*(22./7.)*beam_r*beam_r);
     
//      eps = TDatabase::ParamDB->P1/(1600.0*2000*char_L*char_L); ;   
//      alpha = TDatabase::ParamDB->P7;     
//      rhsfact = (alpha*TDatabase::ParamDB->P4)/(1600.0*2000*(22./7.)*beam_r*beam_r);        
//      
    }
    else if(Region_ID==3) // 3-gray matter,  
    {
     eps = TDatabase::ParamDB->P3/(1030.0*3854*char_L*char_L);
     alpha = TDatabase::ParamDB->P9;           
     rhsfact = (alpha*TDatabase::ParamDB->P4)/(1030.0*3854*(22./7.)*beam_r*beam_r);     
    }
   else if(Region_ID==4) // 4-white matter
    {
     eps = TDatabase::ParamDB->P3/(1030.0*3854*char_L*char_L);;   
     alpha = TDatabase::ParamDB->P10;           
     rhsfact = (alpha*TDatabase::ParamDB->P4)/(1030.0*3854*(22./7.)*beam_r*beam_r);    
    }   
    else // by default it will take as skull
    {
      cout<< " No region info, see example file " << Region_ID << endl;
      exit(0);  
    } 
    
//       cout<< " No region info, see example file " << rhsfact  << endl;
 cout<< eps << " No "  << " region " << char_L  << " info, see example file " << rhsfact  << endl;    
    
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    x = X[i];
    y = Y[i];
    z = Z[i];

    // diffusion
    coeff[0] = eps;
    // convection in x direction
    coeff[1] = 0;
    // convection in y direction
    coeff[2] = 0;
    // convection in z direction
    coeff[3] = 0;
    // reaction
    coeff[4] = 0;
     // rhs

     
   Sourceradius = sqrt(  (x-xp)*(x-xp) + (y-yp)*(y-yp) );
   SourceCoord = z;
    
    if(Sourceradius<=DimlessBeam_r)
     {           
       
       if((SourceCoord)>0)
       {
        cout <<     " SourceCoord " <<SourceCoord  << endl;
        SourceCoord = 0.;
       }
       
      coeff[5] = rhsfact*exp(alpha*SourceCoord*char_L);// f   
//             if(TDatabase::ParamDB->P13<coeff[5]) TDatabase::ParamDB->P13=coeff[5];
      

     }
    else
     {coeff[5] = 0; }// f
    coeff[6] = 0;
    
  }
}
*/

// /*
void BilinearCoeffs(int n_points, double *X, double *Y, double *Z,
        double **parameters, double **coeffs)
{
  double eps;
  int i, Region_ID;
  double *coeff;
  double x, y, z, beam_r, P, k, DimlessBeam_r;
  
//   double t = TDatabase::TimeDB->CURRENTTIME;
  
  double rhsfact, alpha,  char_L;
  
  alpha = TDatabase::ParamDB->P1;
  P = TDatabase::ParamDB->P2;
  beam_r = TDatabase::ParamDB->P11;
  char_L = TDatabase::ParamDB->P12;
  k = TDatabase::ParamDB->P5;
  
//   if(TDatabase::ParamDB->RE_NR!=0)
//    eps = 1.0/TDatabase::ParamDB->RE_NR;
//   else
//    eps = 0;
  
   rhsfact = alpha*P/(1030.0*3710*(22./7.)*beam_r*beam_r);
   eps = k/(1030.0*3710*char_L*char_L);
     
     
//   cout << "eps  eps = " << eps << endl;
//    Region_ID = coeffs[0][1];
  // cout << "Region_ID = " << Region_ID << endl;
   
//        cout<< eps << " No " << k << " region " << char_L  << " info, see example file " << rhsfact  << endl;
   
  DimlessBeam_r =  beam_r/char_L;
   
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    x = X[i];
    y = Y[i];
    z = Z[i];

    // diffusion
    coeff[0] = eps;
    // convection in x direction
    coeff[1] = 0;
    // convection in y direction
    coeff[2] = 0;
    // convection in z direction
    coeff[3] = 0;
    // reaction
    coeff[4] = 0;
     // rhs
    if(sqrt(x*x + y*y)<=DimlessBeam_r)
     {   
      coeff[5] = rhsfact*exp(alpha*z*char_L);// f
      if(TDatabase::ParamDB->P14<coeff[5]) TDatabase::ParamDB->P14=coeff[5];
//       if(z>-1)
//       cout << z << " pow " << alpha*z*char_L << " f " << coeff[5] << endl;
     }
    else
     {coeff[5] = 0; }// f
    coeff[6] = 0;
    
  }
}
 
//  */


 
void ReadMeditMesh(char *SMESH, tetgenio &In)
{  
  int i, j, k, dimension, N_FVert, N_Faces, N_Vertices;
  int BDComp_Min=10000000;

  double beam_r = TDatabase::ParamDB->P11;
  double scale = 2.5e-4/beam_r, x, y, z;
  
  scale = 1.;
    
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
     dat >> x >> y >> z;
     In.pointlist[3*i]   = x*scale;
     In.pointlist[3*i+1] = y*scale;
     In.pointlist[3*i+2] = z*scale;
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
         
//         cout << i<<  " P->vertexlist[j]:  " <<  BDComp_Min << endl;     
     
     
       dat.close();
//    exit(0);    
} // ReadMeditMesh


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
    
    
/*    
    cout << "DeleteDomain" <<endl;
    exit(0);
    */
}


void TetrameshCreate(TDomain *&Domain)
{
  int i, j, k, l, dimension, N_Vertices;
  int N_FVert, N_Faces, *Facemarkerlist, *Facelist, N_RootCells;
  int v1, v2, v3, v4, CellMarker, RefLevel=0, *Tetrahedrals, *PointNeighb, maxEpV, *Tetrahedrals_loc;
  int a, b, c, Neib[2], Neighb_tmp, CurrNeib, len1, len2, len3, CurrComp, N_Points ;
  int N_Cells, MaxLen, jj;  
  const int *EdgeVertex, *TmpFV, *TmpLen;
  
  double X, Y, Z;
  double Xmin = 1e10, Xmax = -1e10, Ymin = 1e10, Ymax = -1e10;  
  double Zmin = 1e10, Zmax = -1e10, T[4]={0,0,0,0}, S[4]={0,0,0,0}; 
  double StartX, StartY, StartZ, BoundX, BoundY, BoundZ; 
  
  char *MESH, line[100];
  MESH = TDatabase::ParamDB->SMESHFILE;

  TVertex **NewVertices;
  TBaseCell **CellTree, *cell;  
  TJoint *Joint;  
  TBoundComp3D *bdcomp; 
  TShapeDesc *ShapeDesc;
  TCollection *coll;
  
#ifdef _MPI  
  int rank, size;  
  MPI_Comm Comm = TDatabase::ParamDB->Comm;
    
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);
 
#endif     
  
  /** delete all cells, vertices in the old domain */
  DeleteDomain(Domain);
  
  
  //   cout << " SMESHFILE is " << SMESH << endl;
  /** load the tet mesh file */  
  std::ifstream dat(MESH);

  if (!dat)
  {
    cerr << "cannot open '" << MESH << "' for input" << endl;
    exit(-1);
  }     
  
    // #ifdef _MPI    
//      if(rank==1)
//     printf(" TetrameshGen  complete \n" );  
//     MPI_Finalize(); 
// #endif    
//   exit(0);

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
  
   // generate new vertices
   NewVertices = new TVertex*[N_Vertices];
  
  //read from file
   for(i=0;i<N_Vertices; i++)
    {
     dat.getline (line, 99);
     dat >> X >> Y >> Z;
     
     NewVertices[i] = new TVertex(X, Y, Z);

      // set bounding box
      if (X > Xmax) Xmax = X;
      if (X < Xmin) Xmin = X;
      if (Y > Ymax) Ymax = Y;
      if (Y < Ymin) Ymin = Y;
      if (Z > Zmax) Zmax = Z;
      if (Z < Zmin) Zmin = Z;
   
     //cout<< i << " vert X: " <<X << " vert Y: " <<Y <<" vert Z: " <<Z <<endl;   
    } // for(i=0;i<N_Vertices; i++)
    
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
// 0.032813  3.78927  -0.883686 // brain cells range
// 181.509  216.944  174.768 // brain cells range
 
 
    
  // find the N_Triangles
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
   
 
   if(N_FVert!=3)
   {
    cerr << "Only Tria surfaces implemented N_FVert: " << N_FVert << endl;
    exit(-1);
   }
  
    Facelist = new int[N_FVert*N_Faces];
    Facemarkerlist = new int[N_Faces];
    
    for(i=0;i<N_Faces; i++)
     {
      dat.getline (line, 99);      

      for(j=0;j<N_FVert;j++)
      {
        dat >> k;
        Facelist[i*N_FVert + j]  = k-1;  // c numbering 
      }
      
      dat >>  Facemarkerlist[i];        

     } //   for(i=0;i<N_Faces; i++)
   

//     for(i=0;i<N_Faces; i++)
//      {
//        if(Facelist[3*i]==0 ||Facelist[3*i+1]==0 ||Facelist[3*i+2]==0 )
//      cout<< i << " vert X: " <<Facelist[3*i] << " vert Y: " <<Facelist[3*i+1] <<" vert Z: " <<Facelist[3*i+2] <<endl;
//    
//     
//      }
   // find the N_RootCells
  while (!dat.eof())
  {
    dat >> line;
    
    if ( (!strcmp(line, "Tetrahedron")) ||  (!strcmp(line, "tetrahedron"))   ||  (!strcmp(line, "TETRAHEDRON"))   ) 
    {
     dat.getline (line, 99);
     dat >> N_RootCells;
     break;
    }    
    // read until end of line
    dat.getline (line, 99);   
  }   
 
  // generate new cells
   CellTree = new TBaseCell*[N_RootCells];
   Tetrahedrals = new int[4*N_RootCells];
   
   for (i=0;i<N_RootCells;i++)
   {
     dat.getline (line, 99);
     dat >> v1 >> v2 >> v3 >> v4 >> CellMarker;  
     Tetrahedrals[4*i    ] = v1 -1;
     Tetrahedrals[4*i + 1] = v2 -1;
     Tetrahedrals[4*i + 2] = v3 -1;
     Tetrahedrals[4*i + 3] = v4 -1;
     
     
     CellTree[i] = new TMacroCell(TDatabase::RefDescDB[Tetrahedron], RefLevel);
     CellTree[i]->SetRegionID(CellMarker);
      
     CellTree[i]->SetVertex(0, NewVertices[Tetrahedrals[4*i    ]]);
     CellTree[i]->SetVertex(1, NewVertices[Tetrahedrals[4*i + 1]]);
     CellTree[i]->SetVertex(2, NewVertices[Tetrahedrals[4*i + 2]]);
     CellTree[i]->SetVertex(3, NewVertices[Tetrahedrals[4*i + 3]]);

     CellTree[i]->SetClipBoard(i);
     ((TMacroCell *) CellTree[i])->SetSubGridID(0);
   }   

   dat.close();
  
  
   Domain->SetTreeInfo(CellTree, N_RootCells);

   // initialize iterators
   TDatabase::IteratorDB[It_EQ]->SetParam(Domain);
   TDatabase::IteratorDB[It_LE]->SetParam(Domain);
   TDatabase::IteratorDB[It_Finest]->SetParam(Domain);
   TDatabase::IteratorDB[It_Between]->SetParam(Domain);
   TDatabase::IteratorDB[It_OCAF]->SetParam(Domain);

   // search neighbours
   PointNeighb = new int[N_Vertices];
#ifdef _MPI    
     if(rank==0)   
#endif       
   cout<<"Numberofpoints "<<N_Vertices<<endl;
   memset(PointNeighb, 0, N_Vertices*SizeOfInt);     
  
     for (i=0;i<4*N_RootCells;i++)
     PointNeighb[Tetrahedrals[i]]++;

   maxEpV = 0;
   for (i=0;i<N_Vertices;i++)
     if (PointNeighb[i] > maxEpV) maxEpV = PointNeighb[i];
   delete [] PointNeighb;
 
#ifdef _MPI    
     if(rank==0)   
#endif   
   cout<<"maxEpV "<< maxEpV<<endl;   
   
   PointNeighb = new int[++maxEpV * N_Vertices];

   memset(PointNeighb, 0, maxEpV*N_Vertices*SizeOfInt);

   // every vertex contains "maxEpV" columns
   // for every vertex at first colomn contains the number of cells containing this vertex
   // at further columns we set the index of corresponding cells containing this vertex
   // cout<<"maxEpV*N_Vertices "<<maxEpV*N_Vertices<<endl;
   for(i=0;i<4*N_RootCells;i++)
    {
     j = Tetrahedrals[i]*maxEpV;
     PointNeighb[j]++;
     //cout<<"j + PointNeighb[j] " << j <<endl;
     PointNeighb[j + PointNeighb[j]] = i / 4;
    }
 
   
   // first generate surface faces   
#ifdef _MPI     
     if(rank==0) 
#endif  
   cout<<"Surface faces "<<N_Faces<<endl;
 
   for (i=0;i<N_Faces;i++)
   {
     a = Facelist[3*i];
     b = Facelist[3*i+1];
     c = Facelist[3*i+2];

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
     
      if(CurrNeib>1)
        continue;
      
//   if(CurrNeib>1)
//    {
//     // cout<<"Face " << i <<  "No. CurrNeib " << CurrNeib << " Neib 1" << Neib[0] <<  " Neib 2" << Neib[1] <<endl;
//     //     exit(0);
//    }

     CurrComp = Facemarkerlist[i] - 1;
 
     if (CurrNeib == 1) // boundary face, scalp and skull layer
      {

       if(CurrComp!=0)
        cout<<"Not a scalp and skull layer: "<<CurrComp<<endl;
  
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
 
//        if (CurrNeib == 2)
//         {
//          if(bdcomp->IsFreeBoundary())
//            Joint = new TIsoInterfaceJoint3D(bdcomp, T, S,
//                        CellTree[Neib[0]], CellTree[Neib[1]] );
//           else
//            Joint = new TInterfaceJoint3D(bdcomp, T, S, 
//                        CellTree[Neib[0]], CellTree[Neib[1]] );
//         }
//        else
//         {
//          if(bdcomp->IsFreeBoundary())
//            Joint = new TIsoBoundFace(bdcomp, T, S);
//          else
           Joint = new TBoundFace(bdcomp, T, S);
//         }
      }
     else
      {
       // cout<<"Inner face"<<endl;
       if (CurrNeib != 2)
       cerr << "Error !!!!!!!! not enough or more neighbours!" << endl;

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

     //cout<<""<< l <<endl;

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
     } // switch (l)
     
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
     } // if (Neib[1] != -1) 

  if (Joint->GetType() == InterfaceJoint3D ||
      Joint->GetType() == IsoInterfaceJoint3D)
      {
        ((TInterfaceJoint3D*)Joint)->SetMapType();
        ((TInterfaceJoint3D*)(Joint))->CheckOrientation();
      }
      else 
       if (Joint->GetType() == JointEqN)
           Joint->SetMapType();

  }   //   for (i=0;i<N_Faces;i++)   
  
 
 /** now generate inner faces (excluding region surface faces) */  
   coll=Domain->GetCollection(It_Finest, 0);
   N_Cells = coll->GetN_Cells();       
     
    for(i=0;i<N_Cells;i++)
     {
      cell = coll->GetCell(i);
      ShapeDesc= cell->GetShapeDesc();   
      ShapeDesc->GetFaceVertex(TmpFV, TmpLen, MaxLen);      
      N_Faces = cell->GetN_Faces();
      Tetrahedrals_loc = Tetrahedrals+4*i;
      
      
       for(jj=0;jj<N_Faces;jj++)
         if(cell->GetJoint(jj) == NULL)
          {
           N_Points = TmpLen[jj];
   
           if(N_Points!=3)
            {     
             cerr << "Only Tria faces are allowed!!! N_FVert: " << N_Points << endl;
             exit(-1);     
            }
            
           //printf(" TetrameshGen  Null joint  \n" );             
           a = Tetrahedrals_loc[TmpFV[jj*MaxLen]];
           b = Tetrahedrals_loc[TmpFV[jj*MaxLen + 1]];
           c = Tetrahedrals_loc[TmpFV[jj*MaxLen + 2]];

           //cout<<"  "<< a<<"  "<< b<<"  "<< c<<endl;

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
               if(Neighb_tmp == PointNeighb[b*maxEpV + k])
                {
                 for(l=1;l<=len3;l++)
                  if(Neighb_tmp == PointNeighb[c*maxEpV + l])
                   {
                    Neib[CurrNeib++] = Neighb_tmp;
                    break;
                   }
                 }
               }
             if (CurrNeib == 2) break;
            } // for (j=1;j<=len1;j++)
     
        if(CurrNeib!=2)
         {
          cerr<<"Face " << i <<  "No. CurrNeib " << CurrNeib << " Neib 1" << Neib[0] <<  " Neib 2" << Neib[1] <<endl;
          exit(0);
         }           
           
        // inner face
        Joint = new TJointEqN(CellTree[Neib[0]], CellTree[Neib[1]]);

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

        //cout<<" l "<< l <<endl;           

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
      } // switch (l)
     
     CellTree[Neib[0]]->SetJoint(j, Joint);

     // second cell containing the current face 
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

     //cout<<" l "<< l <<endl;

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
 
      if (Joint->GetType() == InterfaceJoint3D ||
       Joint->GetType() == IsoInterfaceJoint3D)
        {
         ((TInterfaceJoint3D*)Joint)->SetMapType();
         ((TInterfaceJoint3D*)(Joint))->CheckOrientation();
        }
       else 
        if (Joint->GetType() == JointEqN)
           Joint->SetMapType();      
       
       } //  if(cell->GetJoint(j)
      }  // for(i=0;i<N_Cells;i++)
     
delete [] PointNeighb;

// #ifdef _MPI    
//      if(rank==1)
//     printf(" TetrameshGen  complete \n" );  
//     MPI_Finalize(); 
// #endif    
//   exit(0);
//       
  
}


void TetrameshGen(TDomain *&Domain)
{
 //======================================================================
 // Tetgen for grid generation begin
 //======================================================================
  int i, j, k, l, N_Coord, *N_FVerts, N_Faces, *Facets;
  int N, N_RootCells, N_Cells, CurrVertex, N_Vertices, ID, N_G, RefLevel=0;
  int CurrNeib, len1, len2, len3, maxEpV = 0, a, b, c, Neib[2], Neighb_tmp, CurrComp;
  int *Tetrahedrals, *PointNeighb, dimension;
  int *Facelist, *Facemarkerlist;
 
  double *Coordinates, N_x, N_y, N_z;   
  double *Vertices;
  double Xmin = 1e10, Xmax = -1e10, Ymin = 1e10, Ymax = -1e10;  
  double Zmin = 1e10, Zmax = -1e10, T[4]={0,0,0,0}, S[4]={0,0,0,0}; 
  double StartX, StartY, StartZ, BoundX, BoundY, BoundZ; 
  
  tetgenio In, Out;
   
  TBaseCell **CellTree,  **SurfCellTree;
  TGridCell **DelCell;
  TVertex **VertexDel, **NewVertices, **NewSurfVertices; 
  TBoundPart *BoundPart;
  TBdPlane **UpdateFaceParams;
  TJoint *Joint;
  TBoundComp3D *bdcomp; 
  TCollection *coll, *SurfColl;
  TBaseCell *cell;
  TBdSphere *UpdateParam;  

  char *SMESH, line[100];
  SMESH = TDatabase::ParamDB->SMESHFILE;

#ifdef _MPI  
  int rank, size;  
  MPI_Comm Comm = TDatabase::ParamDB->Comm;
    
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);
 
#endif    

  

  /** only the root call the mesh generator and broadcast the mesh info to all */
#ifdef _MPI    
  if(rank==0)
#endif  
   {
   
//     In.initialize();
//     Out.initialize();
   
    std::ostringstream opts;
    opts << " "; 
  
    opts.seekp(std::ios::beg);
    opts<<'p'; // Tetrahedralize the PLC. Switches are chosen to read a PLC (p)
//     opts<<'r'; // -r  Reconstructs a previously generated mesh     
    opts<<"q"<<1.20; // quality mesh generation(q) with a specified quality bound
//     opts<<"a"<<0.1; // maximum volume constraint
//     opts<<'i'; // Inserts a list of additional points into mesh.
    opts<<'z'; // numbers all output items starting from zero
//     opts<<'d'; // Detect intersections of PLC facets.
    opts<<'f'; // Outputs all  faces (including non-boundary) 
    opts<<'e'; // Outputs a list of edges of the triangulation
// //     opts<<'I'; // Suppresses mesh iteration numbers.
    opts<<'C'; // Checks the consistency of the final mesh.
//     opts<<'Q'; // Quiet: No terminal output except errors.
//     opts<<'g'; // Outputs mesh to .mesh file for viewing by Medit
    opts<<'Y'; // Suppresses boundary facets/segments splitting
//     opts<<'V';  //verbose mode
    opts<<ends;
    
  
//   cout << " SMESHFILE is " << SMESH << endl;
  /** load the medit mesh file into Tetgen */  
  ReadMeditMesh(SMESH, In);
  
  
//     for(i=0;i<In.numberofpoints;i++)
//       OutPut(i<<" (x, y, z) =  "<<
//        In.pointlist[3*i]<<' '<<In.pointlist[3*i+1]<<' '<<In.pointlist[3*i+2]<<endl);

 // Calling  tetrahedralize function of 3dtetgen mesh generator
    tetrahedralize((char*)opts.str().c_str(), &In, &Out);
    
   } // if(rank==0)
     
    
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

   
#ifdef _MPI     
     if(rank==0) 
#endif  
     {
       N_RootCells = Out.numberoftetrahedra;
       N_G = Out.numberofpoints;   
       
       if(N_G==0 || N_RootCells==0)
       {
	OutPut(N_G<<" Out.numberofpoints"<<endl);
	OutPut(N_RootCells<<" Out.N_RootCells"<<endl);
	exit(0);
       }
       
       
       
       // allocate auxillary fields
       Coordinates = Out.pointlist;  
       Tetrahedrals = Out.tetrahedronlist;       
     }

     
     
#ifdef _MPI
   MPI_Bcast(&N_RootCells, 1, MPI_INT, 0, Comm);
   MPI_Bcast(&N_G, 1, MPI_INT, 0, Comm);  
   
   
  if(rank!=0)
   {
    Coordinates = new double [3*N_G];
    Tetrahedrals = new int[4*N_RootCells];
   }
   
   MPI_Bcast(Coordinates, 3*N_G, MPI_DOUBLE, 0, Comm);  
   MPI_Bcast(Tetrahedrals, 4*N_RootCells, MPI_INT, 0, Comm);     
#endif      
   
  // generate new vertices
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
     
    CellTree[i]->SetRegionID(2); // default is CSF
        
   }

 
     

   Domain->SetTreeInfo(CellTree, N_RootCells);


   // initialize iterators
   TDatabase::IteratorDB[It_EQ]->SetParam(Domain);
   TDatabase::IteratorDB[It_LE]->SetParam(Domain);
   TDatabase::IteratorDB[It_Finest]->SetParam(Domain);
   TDatabase::IteratorDB[It_Between]->SetParam(Domain);
   TDatabase::IteratorDB[It_OCAF]->SetParam(Domain);
   
  
   // search neighbours
   PointNeighb = new int[N_G];
#ifdef _MPI    
     if(rank==0)   
#endif       
   cout<<"numberofpoints "<<N_G<<endl;
   memset(PointNeighb, 0, N_G *SizeOfInt);

     for (i=0;i<4*N_RootCells;i++)
     PointNeighb[Tetrahedrals[i]]++;

   for (i=0;i<N_G;i++)
     if (PointNeighb[i] > maxEpV) maxEpV = PointNeighb[i];
   delete [] PointNeighb;

#ifdef _MPI    
     if(rank==0)   
#endif   
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
  
#ifdef _MPI     
     if(rank==0) 
#endif  
     {    
      N_G =  Out.numberoftrifaces;
      Facelist = Out.trifacelist;
      Facemarkerlist = Out.trifacemarkerlist;
     }
      
#ifdef _MPI
   MPI_Bcast(&N_G, 1, MPI_INT, 0, Comm);  

  if(rank!=0)
   {
    Facelist = new int [3*N_G];
    Facemarkerlist = new int[N_G];
   }   
   
   
   MPI_Bcast(Facelist, 3*N_G, MPI_INT, 0, Comm); 
   MPI_Bcast(Facemarkerlist, N_G, MPI_INT, 0, Comm); 
  
  if(rank==0) 
#endif      
   cout<<"numberoftrifaces "<<N_G<<endl;
    
   for (i=0;i<N_G;i++)
   {
     a = Facelist[3*i];
     b = Facelist[3*i+1];
     c = Facelist[3*i+2];

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
// cout<<"Facemarkerlist[i] : "<<Facemarkerlist[i]<<endl;
     if ( CurrNeib != 2 ) // 0 for inner edges and Boundcomp+1 for Boundedge respect
      {

       CurrComp = Facemarkerlist[i] - 1;
//        cout<<"Boundary face CurrComp: "<<CurrComp<<endl;
// exit(0);
       CurrComp = 0; // not yet implemented fully
       
       
       
       bdcomp = Domain->GetBdPart(0)->GetBdComp(CurrComp);
       

//        if(bdcomp->GetTSofXYZ(NewVertices[a]->GetX(), NewVertices[a]->GetY(),
//                              NewVertices[a]->GetY(), T[1], S[1]) ||
//           bdcomp->GetTSofXYZ(NewVertices[b]->GetX(), NewVertices[b]->GetY(),
//                              NewVertices[b]->GetY(), T[2], S[2]) ||
//           bdcomp->GetTSofXYZ(NewVertices[c]->GetX(), NewVertices[c]->GetY(),
//                              NewVertices[c]->GetY(), T[3], S[3])    )
//          {
//           cerr<<"Error: could not set parameter values"<<endl;
//           OutPut(NewVertices[a]<<endl);
//           OutPut(NewVertices[b]<<endl);
//           OutPut(NewVertices[c]<<endl);
//           exit(0);
//          }
 
       if (CurrNeib != 1)
        {
         if(bdcomp->IsFreeBoundary())
           Joint = new TIsoInterfaceJoint3D(bdcomp, T, S,
                       CellTree[Neib[0]], CellTree[Neib[1]] );
          else
           Joint = new TInterfaceJoint3D(bdcomp, T, S, 
                       CellTree[Neib[0]], CellTree[Neib[1]] );
          cerr<<"Error: could not be interface"<<endl;   
   
                exit(0);                 
//                        
        }
       else
        {
//          if(bdcomp->IsFreeBoundary())
//            Joint = new TIsoBoundFace(bdcomp, T, S);
//          else
           Joint = new TBoundFace(bdcomp, T, S);
        }
      }
     else
      {
// //       cout<<"Inner face"<<endl;
//        if (CurrNeib != 2)
//        cerr << "Error !!!!!!!! not enough neighbours!" << endl;

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

#ifdef _MPI    
  if(rank==0)    
#endif
  {
//    In.deinitialize();
//    Out.deinitialize();
  }
  

 
 
#ifdef _MPI    
  if(rank!=0) 
   {
    delete [] Tetrahedrals ;  
    delete [] Coordinates;      
    delete [] Facelist;
    delete [] Facemarkerlist;
   }
#endif  
  
  delete [] PointNeighb;

// #ifdef _MPI    
//      if(rank==3)
//     printf(" TetrameshGen  %d \n",   N_G );  
//     MPI_Finalize(); 
// #endif    
//   exit(0);
//      

}


/** same as TetrameshGen but without using the face list info from tetgen */
void TetraGen(TDomain *&Domain)
{
 //======================================================================
 // Tetgen for grid generation begin
 //======================================================================
  int i, j, k, l, N_Coord, *N_FVerts, N_Faces, *Facets;
  int N, N_RootCells, N_Cells, CurrVertex, N_Vertices, ID, N_G, RefLevel=0;
  int CurrNeib, len1, len2, len3, maxEpV = 0, a, b, c, Neib[2], Neighb_tmp, CurrComp;
  int *Tetrahedrals, *PointNeighb, dimension;
  int jj, N_Points, MaxLen, *Tetrahedrals_loc, v[4], N_Layercells=0;
  const int *EdgeVertex, *TmpFV, *TmpLen;

  double *Coordinates, N_x, N_y, N_z, x_mean, y_mean, z_mean;   
  double *Vertices;
  double Xmin = 1e10, Xmax = -1e10, Ymin = 1e10, Ymax = -1e10;  
  double Zmin = 1e10, Zmax = -1e10, T[4]={0,0,0,0}, S[4]={0,0,0,0}; 
  double StartX, StartY, StartZ, BoundX, BoundY, BoundZ; 
  
  tetgenio In, Out;
   
  TBaseCell **CellTree,  **SurfCellTree;
  TGridCell **DelCell;
  TVertex **VertexDel, **NewVertices, **NewSurfVertices; 
  TBoundPart *BoundPart;
  TBdPlane **UpdateFaceParams;
  TJoint *Joint;
  TBoundComp3D *bdcomp; 
  TCollection *coll, *SurfColl;
  TBaseCell *cell;
  TBdSphere *UpdateParam;  
  TShapeDesc *ShapeDesc;
   
  char *SMESH, line[100];
  SMESH = TDatabase::ParamDB->SMESHFILE;

#ifdef _MPI  
  int rank, size;  
  MPI_Comm Comm = TDatabase::ParamDB->Comm;
    
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);
 
#endif    

  /** only the root call the mesh generator and broadcast the mesh info to all */
#ifdef _MPI    
  if(rank==0)
#endif  
   {
   
//     In.initialize();
//     Out.initialize();
   
    std::ostringstream opts;
    opts << " "; 
  
    opts.seekp(std::ios::beg);
    opts<<'p'; // Tetrahedralize the PLC. Switches are chosen to read a PLC (p)
//     opts<<'r'; // -r  Reconstructs a previously generated mesh     
    opts<<"q"<<1.20; // quality mesh generation(q) with a specified quality bound
//     opts<<"a"<<10; // maximum volume constraint
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
//     opts<<'V';  //verbose mode
    opts<<ends;
    
  
//   cout << " SMESHFILE is " << SMESH << endl;
  /** load the medit mesh file into Tetgen */  
  ReadMeditMesh(SMESH, In);
  
  
//     for(i=0;i<In.numberofpoints;i++)
//       OutPut(i<<" (x, y, z) =  "<<
//        In.pointlist[3*i]<<' '<<In.pointlist[3*i+1]<<' '<<In.pointlist[3*i+2]<<endl);

 // Calling  tetrahedralize function of 3dtetgen mesh generator
    tetrahedralize((char*)opts.str().c_str(), &In, &Out);

 //   output: coordinates of all vertices
//  for(i=0;i<Out.numberofpoints;i++)
//   OutPut(" (x, y, z) =  "<< Out.pointlist[3*i]<<' '<<Out.pointlist[3*i+1]<<' '<<Out.pointlist[3*i+2]<<endl);    
    
   } // if(rank==0)
     
  /** delete all cells, vertices in the old domain */
  DeleteDomain(Domain);

   
#ifdef _MPI     
     if(rank==0) 
#endif  
     {
       N_RootCells = Out.numberoftetrahedra;
       N_G = Out.numberofpoints;   

       // allocate auxillary fields
       Coordinates = Out.pointlist;  
       Tetrahedrals = Out.tetrahedronlist;       
     }
 
#ifdef _MPI
   MPI_Bcast(&N_RootCells, 1, MPI_INT, 0, Comm);
   MPI_Bcast(&N_G, 1, MPI_INT, 0, Comm);  
   
   
  if(rank!=0)
   {
    Coordinates = new double [3*N_G];
    Tetrahedrals = new int[4*N_RootCells];
   }
   
   MPI_Bcast(Coordinates, 3*N_G, MPI_DOUBLE, 0, Comm);  
   MPI_Bcast(Tetrahedrals, 4*N_RootCells, MPI_INT, 0, Comm);     
#endif      
   
  // generate new vertices
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

     v[0] = Tetrahedrals[4*i];
     v[1] = Tetrahedrals[4*i+1];
     v[2] = Tetrahedrals[4*i+2];
     v[3] = Tetrahedrals[4*i+3];
     
     CellTree[i]->SetVertex(0, NewVertices[v[0]]);
     CellTree[i]->SetVertex(1, NewVertices[v[1]]);
     CellTree[i]->SetVertex(2, NewVertices[v[2]]);
     CellTree[i]->SetVertex(3, NewVertices[v[3]]);

     CellTree[i]->SetClipBoard(i);
     ((TMacroCell *) CellTree[i])->SetSubGridID(0);

     
    CellTree[i]->SetRegionID(2); // default is CSF
    CellTree[i]->SetAsLayerCell(1);    
//            
//     if(TDatabase::ParamDB->P15)
//      {
//       x_mean=0; y_mean=0; z_mean=0;
//       for (j=0;j<4;j++)
//       {
//        x_mean += Coordinates[3*v[j]];
//        y_mean += Coordinates[3*v[j] + 1];       
//        z_mean += Coordinates[3*v[j] + 2];           
//       }
//       
//        x_mean /=4.;
//        y_mean /=4.;  
//        z_mean /=4.;  
//        
//        if((sqrt(x_mean*x_mean + y_mean*y_mean)<=2) && z_mean>=-10.)
//         {
// 
//          CellTree[i]->SetAsLayerCell(1);
// //          cout << N_Layercells  << " x_mean " << x_mean << " y_mean " << y_mean << " z_mean " << z_mean << endl;
// //          N_Layercells++;
// 
//         }
//         
//       } //if(TDatabase::ParamDB->P15)
    } //  for (i=0;i<N_RootCells;i++)

   Domain->SetTreeInfo(CellTree, N_RootCells);

// exit(0);

   // initialize iterators
   TDatabase::IteratorDB[It_EQ]->SetParam(Domain);
   TDatabase::IteratorDB[It_LE]->SetParam(Domain);
   TDatabase::IteratorDB[It_Finest]->SetParam(Domain);
   TDatabase::IteratorDB[It_Between]->SetParam(Domain);
   TDatabase::IteratorDB[It_OCAF]->SetParam(Domain);
   
  
   // search neighbours
   PointNeighb = new int[N_G];
#ifdef _MPI    
     if(rank==0)   
#endif       
   cout<<"numberofpoints "<<N_G<<endl;
   memset(PointNeighb, 0, N_G *SizeOfInt);

     for (i=0;i<4*N_RootCells;i++)
     PointNeighb[Tetrahedrals[i]]++;

   for (i=0;i<N_G;i++)
     if (PointNeighb[i] > maxEpV) maxEpV = PointNeighb[i];
   delete [] PointNeighb;

#ifdef _MPI    
     if(rank==0)   
#endif   
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


   // generate new faces    
   coll=Domain->GetCollection(It_Finest, 0);
   N_Cells = coll->GetN_Cells();       
   N_G = 0;
   for (i=0;i<N_Cells;i++)
   {
    cell = coll->GetCell(i);
    ShapeDesc= cell->GetShapeDesc();   
    ShapeDesc->GetFaceVertex(TmpFV, TmpLen, MaxLen);      
    N_Faces = cell->GetN_Faces();
    Tetrahedrals_loc = Tetrahedrals+4*i;
           
    for(jj=0;jj<N_Faces;jj++)
     if(cell->GetJoint(jj) == NULL)
      {
       N_G++;    
       N_Points = TmpLen[jj];
   
       if(N_Points!=3)
        {     
         cerr << "Only Tria faces are allowed!!! N_FVert: " << N_Points << endl;
         exit(-1);     
        }     
     
     
       //printf(" TetrameshGen  Null joint  \n" );             
       a = Tetrahedrals_loc[TmpFV[jj*MaxLen]];
       b = Tetrahedrals_loc[TmpFV[jj*MaxLen + 1]];
       c = Tetrahedrals_loc[TmpFV[jj*MaxLen + 2]];

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
        }// for (j=1;j<=len1;j++)
        
      //   cout<<"CurrNeib " << CurrNeib <<endl;
      // cout<<"Facemarkerlist[i] : "<<Facemarkerlist[i]<<endl;
      if (CurrNeib == 1) // 0 for inner edges and Boundcomp+1 for Boundedge respect
       {
        CurrComp = 0; // not yet implemented fully
       
       
       
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
 
//        if (CurrNeib == 2)
//         {
//          if(bdcomp->IsFreeBoundary())
//            Joint = new TIsoInterfaceJoint3D(bdcomp, T, S,
//                        CellTree[Neib[0]], CellTree[Neib[1]] );
//           else
//            Joint = new TInterfaceJoint3D(bdcomp, T, S, 
//                        CellTree[Neib[0]], CellTree[Neib[1]] );
//         }
//        else
        {
         if(bdcomp->IsFreeBoundary())
           Joint = new TIsoBoundFace(bdcomp, T, S);
         else
           Joint = new TBoundFace(bdcomp, T, S);
        }
      }
     else
      {
       // cout<<"Inner face"<<endl;
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
     } // if (Neib[1] != -1)

  if (Joint->GetType() == InterfaceJoint3D ||
      Joint->GetType() == IsoInterfaceJoint3D)
      {
        ((TInterfaceJoint3D*)Joint)->SetMapType();
        ((TInterfaceJoint3D*)(Joint))->CheckOrientation();
      }
      else 
       if (Joint->GetType() == JointEqN)
           Joint->SetMapType();
   } // if(cell->GetJoint(jj) == NULL)
  } //    for (i=0;i<N_Cells;i++)

  
#ifdef _MPI    
    if(rank==0) 
#endif      
   cout<<"numberoftrifaces after "<<N_G<<endl;
  
#ifdef _MPI    
  if(rank==0)    
#endif
  {
//    In.deinitialize();
//    Out.deinitialize();
  }
  

 
 
#ifdef _MPI    
  if(rank!=0) 
   {
    delete [] Tetrahedrals ;  
    delete [] Coordinates;      
//     delete [] Facelist;
//     delete [] Facemarkerlist;
   }
#endif  
  
  delete [] PointNeighb;

// #ifdef _MPI    
//      if(rank==3)
//     printf(" TetrameshGen  %d \n",   N_G );  
//     MPI_Finalize(); 
// #endif    
//   exit(0);
//      

}

