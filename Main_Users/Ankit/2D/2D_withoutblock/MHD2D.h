// Navier-Stokes problem, Poiseuille-Problem
// 
// u(x,y) = (4*y*(1-y), 0)
// p(x,y) = x-1/2

void ExampleFile()
{
  OutPut("Example: MHD2D.h" << endl) ;
  
  TDatabase::ParamDB->INTERNAL_STEADY_STATE_MATRICES_OR_RHS = 0;
  TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD = 1;
}

#include <Joint.h>
#include <BoundEdge.h>
#include <BoundComp.h>
#include <FE2D.h>
#include <FEDesc2D.h>

#include <MainUtilities.h>


#include <MacroCell.h>
#include <BoundEdge.h>
#include <IsoBoundEdge.h>

#include <IsoInterfaceJoint.h>
#include <BdLine.h>
#include <BdCircle.h>
#include <GridCell.h>

#include <QuadAffin.h>
#include <QuadBilinear.h>
#include <QuadIsoparametric.h>
#include <TriaAffin.h>
#include <TriaIsoparametric.h>

extern "C"
{
#include <gridgen.h>
  void triangulate(char*, struct triangulateio*,
		   struct triangulateio*, struct triangulateio*);
}

// ========================================================================
// initial solution
// ========================================================================

void InitialU1(double x, double y, double *values)
{
  values[0] = 0;
}

void InitialU2(double x, double y, double *values)
{
  values[0] = 0;
}

void InitialP(double x, double y, double *values)
{
  values[0] = 0;
}

// ========================================================================
// boundary conditions for Convection-Diffusion
// ========================================================================
void BoundCondition_CD(int BdComp, double t, BoundCond &cond)
{    
  if(BdComp == 0 || BdComp == 2 || BdComp == 3)
    cond = DIRICHLET;
  else
    cond = NEUMANN;
}

// value of boundary condition
void BoundValue_CD(int BdComp, double Param, double &value)
{

  switch(BdComp)
  {
    case 0: value = 1;
            break;
    case 1: value = 0;
            break;
    case 2: value = 0;
            break;
    case 3: value = 0;
            break;
    case 4: value = 0;
            break;      
    default: cout << "CD wrong boundary part number" << endl;
            break;
  }
}


void BilinearCoeffs(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  int i;
  double *coeff, *param;
  double Pe = TDatabase::ParamDB->RE_NR*Prandtl_Number;
  double eps=1/Pe;
  
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];
    
    coeff[0] = eps;
    if(TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD)
     {
      coeff[1] = param[0];  // u1
      coeff[2] = param[1];  // u2
//      cout << "Param " << param[0] << "  " << param[1] << endl;
     }
    else
     {
      coeff[1] =  0.;  // u1
      coeff[2] =  0.;  // u2
     }    
    coeff[3] = 0;
    coeff[4] = 0;
    // rhs from previous time step
  }
    
}

// ========================================================================
// boundary conditions for Navier-Stokes
// ========================================================================
void BoundCondition(int i, double t, BoundCond &cond)
{
  if(i == 1)
    cond = NEUMANN;
  else
    cond = DIRICHLET;
  
  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0; 
}

void U1BoundValue(int BdComp, double Param, double &value)
{
//   double y = (2*Param-1)/2;
//   double k = sqrt(2*Ha);
//   double h_ = a_duct/d_duct;
  switch(BdComp)
  {
    case 0: value= 0;
            break;
    case 1: value= 0;
            break;
    case 2: value=0;
            break;
    case 3: value = 4*(Param)*(1-Param);//( cosh(y*k/h_) - cosh(.5*k/h_) )/( 1 - cosh(.5*k/h_) ) ; //4*(Param)*(1-Param);//0.5*Param*(1-Param);//(TDatabase::ParamDB->RE_NR*Kinematic_Viscosity/d_duct)
            //cout << "Param ------------ " << Param << endl;
            break;
    case 4: value=0;
            break;          
    default: cout << "U1 wrong boundary part number" << endl;
            break;
  }
}


void U2BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0: value=0;
            break;
    case 1: value=0;
            break;
    case 2: value=0;
            break;
    case 3: value=0;
            break;
    case 4: value=0;
            break;        
    default: cout << "U2 wrong boundary part number" << endl;
            break;
      
  }
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *x, double *y,
               double **parameters, double **coeffs)
{
  static double eps = 1./TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = 0;// f1
    coeff[2] = 0; // f2
  }
}

// **************************************************************************
// Triangular Mesh Generation
// **************************************************************************
void  TriaReMeshGen(TDomain *&Domain)
{
  int j, ID, k, N_G, *PartMarker, *PointNeighb, maxEpV=0;
  int a, b, len1, len2, Neighb_tmp, BDpart;
  int i, temp, N_Cells, N_RootCells, CurrVertex, N_Joints, N_Vertices;
  int N_Interface_Vert, N_Verti, N_Hori, N_SlipBound_Vert, N_BDVertices;
  int CurrComp, In_Index, *Triangles, Neib[2], CurrNeib;
  
  double deviation, hi, x0, y0, x, y, phi1, phi2;
  double T_a, T_b, C_x, C_y, s, theta;  
  double dt, area = TDatabase::ParamDB->Area;
  double *Coordinates, *Hole_List, *S_BX, *S_BY;  
  
  
  double d = .01;
  double w = 4;      // w = w/d  
  double Lu = 12;   // Lu = 12d/d
  double Ld = 42;   // Lu = 42d/d
  double Xi[4] = {-Lu, Ld, Ld, -Lu};
  double Yi[4] = {-w/2, -w/2, w/2, w/2};
  

  TBaseCell **CellTree, *cell;
  TBoundPart *BoundPart;
  TJoint *Joint;
  TCollection *coll;
  TVertex **VertexDel, **NewVertices;  
  TBdLine *UpdateBound[12];
  TBdCircle *UpdateIntface;

  
  boolean AllowEdgeRef = (boolean) TDatabase::ParamDB->MESHGEN_ALLOW_EDGE_REF;
  
  struct triangulateio In, Out;
  std::ostringstream opts;
  opts << " ";
  
  Out.pointlist = NULL;
  Out.pointattributelist = NULL;
  Out.pointmarkerlist = NULL;
  Out.trianglelist = NULL;
  Out.triangleattributelist = NULL;
  Out.trianglearealist = NULL;
  Out.neighborlist = NULL;
  Out.segmentlist = NULL;
  Out.segmentmarkerlist = NULL;
  Out.holelist = NULL;
  Out.regionlist = NULL;
  Out.edgelist = NULL;
  Out.edgemarkerlist = NULL;
  Out.normlist = NULL;

  opts.seekp(std::ios::beg); 

  

//   char *SMESH, line[100];
     char *SMESH, line[309];
  
  SMESH = TDatabase::ParamDB->SMESHFILE; 
  
  std::ifstream dat(SMESH);

  if (!dat)
  {
    cerr << "cannot open '" << SMESH << "' for input" << endl;
    exit(-1);
  }      
  
   // check the dimension
  N_Interface_Vert = -1;
  while (!dat.eof())
  {
    dat >> line;
    N_Interface_Vert++; 
    
       // read until end of line
    dat.getline (line, 309); 
//        dat.getline (line, 99);
  }   
  
   dat.close();
  
//   cout<< " TriaReMeshGen " << N_G << endl;

  
  S_BX = new double[N_Interface_Vert];
  S_BY = new double[N_Interface_Vert];
  
  std::ifstream dat1(SMESH);
  
  // check the dimension
 

  //read from file
   for(i=0;i<N_Interface_Vert; i++)
    {
     dat1 >> S_BX[i] >> S_BY[i]; 
     dat1.getline (line, 309);
//      dat1.getline (line, 99);
      S_BX[i]= S_BX[i];
      S_BY[i]= S_BY[i];
//      cout<< i << " vert X: " <<S_BX[i] << " vert Y: " << S_BY[i] <<endl;   
    }
    
  dat1.close();  
// exit(0);

  BoundPart = Domain->GetBdPart(0);
  UpdateBound[0]  = (TBdLine*)BoundPart->GetBdComp(0);
  UpdateBound[1]  = (TBdLine*)BoundPart->GetBdComp(1);
  UpdateBound[2]  = (TBdLine*)BoundPart->GetBdComp(2);
  UpdateBound[3]  = (TBdLine*)BoundPart->GetBdComp(3);

  BoundPart = Domain->GetBdPart(1);
  UpdateIntface = (TBdCircle*)BoundPart->GetBdComp(0); 


  cout<<area<<endl;

  
// OutPut("MESHGEN_REF_QUALIT " << TDatabase::ParamDB->MESHGEN_REF_QUALITY << endl);

  opts<<'p'; // Constrained Delaunay Triangulation:
           // initial values - only points defined on the boundary of the domain;
           // triangulation near boundary may variate from Delaunay criterion
  opts<<"q"<<  TDatabase::ParamDB->MESHGEN_REF_QUALITY;
              // Quality mesh generation with no angles smaller than 20 degrees;
  opts<<"a"<< area; // Imposes a maximum triangle area.
  opts<<'e'; // Outputs a list of edges of the triangulation
  opts<<'z'; // Numbers if items starting from 0
  //opts<<"VVVV"; // Gives detailed information about what Triangle is doing
  opts<<'Q'; // Supress all explanation of what Triangle is doing, unless an error occurs
//   opts<<'Y'; // Supress adding vertices on boundary edges
  opts<<ends;
  
//   N_Interface_Vert = int (TDatabase::ParamDB->P6);    //Freesurf except end point
  N_Hori  = 400;      // number of horrizontal BD vertices
  N_Verti = 100;       // number of vertical BD vertices
  N_SlipBound_Vert = 2*N_Hori + 2*N_Verti;  

  N_BDVertices = N_Interface_Vert+N_SlipBound_Vert;
  In.numberofpoints = N_BDVertices;
  In.pointlist = new double[2*In.numberofpoints];
  In.pointmarkerlist = new int[In.numberofpoints];
  In.numberofpointattributes = 0;

  In.numberofsegments = In.numberofpoints;
  In.segmentlist = new int[2*In.numberofsegments];
  In.segmentmarkerlist = new int[In.numberofsegments]; 
  In.numberofregions = 0;
  In.regionlist = NULL;
  
  In.numberofholes = 0;
  In.holelist = NULL;

  Hole_List = new double[2* In.numberofholes];
//   Hole_List[0] = 0.;
//   Hole_List[1] = 0.;
  Hole_List[0] = 0;
  Hole_List[1] = 0;
  In.holelist = Hole_List;
    
  In_Index = 0;
  CurrComp = 1;  
  
  hi = (Xi[1] - Xi[0])/(double)N_Hori;
  x0 = Xi[0];
  y0 = Yi[0];
  x  = Xi[0];
  
  // points and segments on the horizontal boundary (marker=1)
  for(i=0;i<N_Hori;i++) // without last point
   {
    x = x0 + (double)i*hi;
    In.pointlist[2*In_Index] = x;
    In.pointlist[2*In_Index+1] = y0;
    // cout<<" x : "<< In.pointlist[2*In_Index] << " y : "<< In.pointlist[2*In_Index+1] <<endl;
    In.pointmarkerlist[In_Index] = CurrComp;
    In.segmentlist[2*In_Index] = In_Index;
    In.segmentlist[2*In_Index+1] = In_Index+1;
    In.segmentmarkerlist[In_Index] = CurrComp;
    In_Index++;
   }
   
  CurrComp++;  
  
  hi = (Yi[2] - Yi[1])/(double)N_Verti;
  x0 = Xi[1];
  y0 = Yi[1];
  y  = Yi[1];
  // points and segments on the horizontal boundary (marker=1)
  for(i=0;i<N_Verti;i++) // without last point
   {
    y = y0 + (double)i*hi; 
    In.pointlist[2*In_Index] = x0;
    In.pointlist[2*In_Index+1] = y;
    //cout<<" x : "<< In.pointlist[2*In_Index] << " y : "<< In.pointlist[2*In_Index+1] <<endl;
    In.pointmarkerlist[In_Index] = CurrComp;
    In.segmentlist[2*In_Index] = In_Index;
    In.segmentlist[2*In_Index+1] = In_Index+1;
    In.segmentmarkerlist[In_Index] = CurrComp;
    In_Index++;
   }
   
  CurrComp++;  
//   cout<<endl;

  hi = (Xi[3] - Xi[2])/(double)N_Hori;
  x0 = Xi[2];
  y0 = Yi[2];
  x  = Xi[2];
  // points and segments on the horizontal boundary (marker=1)
 for(i=0;i<N_Hori;i++) // without last point
   {
    x = x0 + (double)i*hi;
    In.pointlist[2*In_Index] = x;
    In.pointlist[2*In_Index+1] = y0;
    // cout<<" x : "<< In.pointlist[2*In_Index] << " y : "<< In.pointlist[2*In_Index+1] <<endl;
    In.pointmarkerlist[In_Index] = CurrComp;
    In.segmentlist[2*In_Index] = In_Index;
    In.segmentlist[2*In_Index+1] = In_Index+1;
    In.segmentmarkerlist[In_Index] = CurrComp;
    In_Index++;

   }
   
  CurrComp++;
//   cout<<endl;

  hi = (Yi[0] - Yi[3])/(double)N_Verti;
  x0 = Xi[3];
  y0 = Yi[3];
  y  = Yi[3];
  // points and segments on the horizontal boundary (marker=1)
 for(i=0;i<N_Verti;i++) // without last point
   {
    y = y0 + (double)i*hi;     
    In.pointlist[2*In_Index] = x0;
    In.pointlist[2*In_Index+1] = y;
    // cout<<" x : "<< In.pointlist[2*In_Index] << " y : "<< In.pointlist[2*In_Index+1] <<endl;
    In.pointmarkerlist[In_Index] = CurrComp;
    In.segmentlist[2*In_Index] = In_Index;
    In.segmentlist[2*In_Index+1] = In_Index+1;
    In.segmentmarkerlist[In_Index] = CurrComp;
    In_Index++;
 
   }
  CurrComp++;  
  
  In.segmentlist[2*(In_Index-1)+1] = 0;
//   temp=In_Index;
//  
//   for(i=0;i<N_Interface_Vert;i++) // without last point
//    { 
//       In.pointlist[2*In_Index] = S_BX[i];
//       In.pointlist[2*In_Index+1] = S_BY[i];
//       //cout<<" x : "<< S_BX[i] << " y : "<< S_BY[i]<<endl;
//       In.pointmarkerlist[In_Index] = CurrComp;
//       In.segmentlist[2*In_Index] = In_Index;
//       In.segmentlist[2*In_Index+1] = In_Index+1;
// //       if (AllowEdgeRef)
//        { In.segmentmarkerlist[In_Index] = CurrComp; }
// //       else
// //        { In.segmentmarkerlist[In_Index] = 100000 + CurrComp; }
// 
//       In_Index++;  
//    }
// 
//   In.segmentlist[2*(In_Index-1)+1] = temp;  
//   
//   exit(0);
 
  if(Out.pointlist!=NULL) {
    free(Out.pointlist); Out.pointlist = NULL;}
  if(Out.pointattributelist!=NULL) {
    free(Out.pointattributelist); Out.pointattributelist = NULL;}
  if(Out.pointmarkerlist!=NULL) {
    free(Out.pointmarkerlist); Out.pointmarkerlist = NULL;}
  if(Out.trianglelist!=NULL) {
    free(Out.trianglelist); Out.trianglelist = NULL;}
  if(Out.triangleattributelist!=NULL) {
    free(Out.triangleattributelist); Out.triangleattributelist = NULL;}
  if(Out.trianglearealist!=NULL) {
    free(Out.trianglearealist); Out.trianglearealist = NULL;}
  if(Out.neighborlist!=NULL) {
    free(Out.neighborlist); Out.neighborlist = NULL;}
  if(Out.segmentlist!=NULL) {
    free(Out.segmentlist); Out.segmentlist = NULL;}
  if(Out.segmentmarkerlist!=NULL) {
    free(Out.segmentmarkerlist); Out.segmentmarkerlist = NULL;}
  if(Out.holelist!=NULL) {
    free(Out.holelist); Out.holelist = NULL;}
  if(Out.regionlist!=NULL) {
    free(Out.regionlist); Out.regionlist = NULL;}
  if(Out.edgelist!=NULL) {
    free(Out.edgelist); Out.edgelist = NULL;}
  if(Out.edgemarkerlist!=NULL) {
    free(Out.edgemarkerlist); Out.edgemarkerlist = NULL;}
  if(Out.normlist!=NULL) {
    free(Out.normlist); Out.normlist = NULL;}
//    cout<<"asd\n";
  // call triangle
  triangulate((char*)opts.str().c_str(), &In, &Out, (struct triangulateio *)NULL);
//    cout<<"hi\n";
//      exit(0);
  
 
  Domain->GetTreeInfo(CellTree, N_RootCells);
  coll = Domain->GetCollection(It_Finest, 0);
  N_Cells = coll->GetN_Cells();
  
  // remove all existing vertices and joints
  VertexDel = new TVertex*[3*N_RootCells];
  
   CurrVertex = 0;

  for(i=0;i<N_Cells;i++)
    {
      cell = coll->GetCell(i);
      N_Joints = cell->GetN_Joints();
      N_Vertices = cell->GetN_Vertices();
      for(j=0;j<N_Joints;j++)
        {
         if(CurrVertex==0)
          {
              VertexDel[CurrVertex] = cell->GetVertex(j);
              CurrVertex++;
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
              VertexDel[CurrVertex] = cell->GetVertex(j);
              CurrVertex++;
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
            VertexDel[CurrVertex] = cell->GetVertex((j+1)%N_Vertices);
            CurrVertex++;
           }
        } // for j
    } // for i
  for(i=0;i<CurrVertex;i++)
    delete VertexDel[i];

  delete []VertexDel;
  OutPut(CurrVertex<<" vertices were deleted"<<endl);

 // remove all existing cells and joints
  for(i=0;i<N_RootCells;i++)
    delete (TGridCell*)CellTree[i];
  OutPut(N_RootCells<<" cells were deleted"<<endl);
   delete CellTree;
   delete coll;

   
   
 // Solid Bound startx, starty, x length and y length
  UpdateBound[0]->SetParams(Xi[0], Yi[0], Xi[1]-Xi[0],Yi[1]-Yi[0]);
  UpdateBound[1]->SetParams(Xi[1], Yi[1], Xi[2]-Xi[1],Yi[2]-Yi[1]);
  UpdateBound[2]->SetParams(Xi[2], Yi[2], Xi[3]-Xi[2],Yi[3]-Yi[2]);
  UpdateBound[3]->SetParams(Xi[3], Yi[3], Xi[0]-Xi[3],Yi[0]-Yi[3]);
  

// Free boundary xmid, ymid, radius_a, radius_b, start angle, end angle
  UpdateIntface->SetParams(C_x, C_y, T_a, T_b, phi1, phi2);   
   
   
  N_RootCells = Out.numberoftriangles; 

  // allocate auxillary fields
  Coordinates = Out.pointlist;
  Triangles = Out.trianglelist;
  PartMarker = new int[Out.numberofpoints];

  // generate new vertices
  N_G = Out.numberofpoints;
  NewVertices = new TVertex*[N_G];

  for (i=0;i<N_G;i++)
     NewVertices[i] = new TVertex(Coordinates[2*i], Coordinates[2*i+1]);
 
 // generate cells
  CellTree = new TBaseCell*[N_RootCells];

  for (i=0;i<N_RootCells;i++)
  {
    CellTree[i] = new TMacroCell(TDatabase::RefDescDB[Triangle], 0);

    CellTree[i]->SetVertex(0, NewVertices[Out.trianglelist[3*i    ]]);
    CellTree[i]->SetVertex(1, NewVertices[Out.trianglelist[3*i + 1]]);
    CellTree[i]->SetVertex(2, NewVertices[Out.trianglelist[3*i + 2]]);

      ((TMacroCell *) CellTree[i])->SetSubGridID(0);
  }

  Domain->SetTreeInfo(CellTree, N_RootCells);

  // initialize iterators
  TDatabase::IteratorDB[It_EQ]->SetParam(Domain);
  TDatabase::IteratorDB[It_LE]->SetParam(Domain);
  TDatabase::IteratorDB[It_Finest]->SetParam(Domain);
  TDatabase::IteratorDB[It_Between]->SetParam(Domain);
  TDatabase::IteratorDB[It_OCAF]->SetParam(Domain);   
  
  // search neighbours
  N_G = Out.numberofpoints;
  PointNeighb = new int[N_G];

  memset(PointNeighb, 0, N_G *SizeOfInt);

  for (i=0;i<3*N_RootCells;i++)
    PointNeighb[Triangles[i]]++;

  maxEpV=0;

  for (i=0;i<N_G;i++)
    if (PointNeighb[i] > maxEpV) maxEpV = PointNeighb[i];
  delete [] PointNeighb;

  PointNeighb = new int[++maxEpV * N_G];

  memset(PointNeighb, 0, maxEpV * N_G *SizeOfInt);
      
   // first colomn contains the number of following elements
   // for every point at first column we set the number of neighbour points
   // at further columns we set the index of corresponding cells
  for(i=0;i<3*N_RootCells;i++)
  {
    j = Triangles[i]*maxEpV;
    PointNeighb[j]++;
    PointNeighb[j + PointNeighb[j]] = i / 3;
  }
 
  // generate new edges
  N_G = Out.numberofedges;
  for (i=0;i<N_G;i++)
  {
    a = Out.edgelist[2*i];
    b = Out.edgelist[2*i+1];
    Neib[0] = -1;
    Neib[1] = -1;
    CurrNeib = 0;

    len1 = PointNeighb[a*maxEpV];
    len2 = PointNeighb[b*maxEpV];

  // find indexes of cells containing the current edge
   for (j=1;j<=len1;j++)
    {
      Neighb_tmp = PointNeighb[a*maxEpV + j];
       for (k=1;k<=len2;k++)
        if (Neighb_tmp == PointNeighb[b*maxEpV + k])
        {
          Neib[CurrNeib++] = Neighb_tmp;
          break;
        }
      if (CurrNeib == 2) break;
    }
 

    
    if (Out.edgemarkerlist[i]) // 0 for inner edges and Boundcomp+1 for Boundedge respect
    {
      CurrComp = Out.edgemarkerlist[i] - 1;
      if (CurrComp >= 100000) CurrComp -= 100000;
      
      BDpart=Domain->GetBdPartID(CurrComp);
      CurrComp= Domain->GetLocalBdCompID(CurrComp);

      if(Domain->GetBdPart(BDpart)->GetBdComp(CurrComp)->GetTofXY(
            NewVertices[a]->GetX(), NewVertices[a]->GetY(), T_a) ||
          Domain->GetBdPart(BDpart)->GetBdComp(CurrComp)->GetTofXY(
            NewVertices[b]->GetX(), NewVertices[b]->GetY(), T_b))
       {
          cerr<<"Error: could not set parameter values"<<endl;
          OutPut(NewVertices[a]<<endl);
          OutPut(NewVertices[b]<<endl);
        //  exit(0);
       }
       
// ===============================================================================       
 //  due to the hole the orientation of the circle is colck-wise  from out side
// the parameter of the starting edge (in 0, 2pi) is for e.g (0, 0.9) (colck-wise) 
// while refining to get the mid point we should get  the mid point parameter as 0.95 
// but we will get only (0+0.9)/2 =0.45 (wrong mid point). So we change.
//  Note only if the orientation is colck-wise  !!!!!!!
// ===============================================================================  

     if(BDpart==1 && CurrComp==0  && fabs(T_a)==0 ) T_a=1;      
           
      if (CurrNeib == 2)    // 2 cells contain the current edge
       {
        if(Domain->GetBdPart(BDpart)->GetBdComp(CurrComp)->IsFreeBoundary())
	 {
          Joint = new TIsoInterfaceJoint(Domain->GetBdPart(BDpart)->GetBdComp(CurrComp),
                  T_a, T_b, CellTree[Neib[0]], CellTree[Neib[1]]);
	 }
        else
	 {
          Joint = new TInterfaceJoint(Domain->GetBdPart(BDpart)->GetBdComp(CurrComp),
                  T_a, T_b, CellTree[Neib[0]], CellTree[Neib[1]]);
	 }
       }
      else
       {
        if(Domain->GetBdPart(BDpart)->GetBdComp(CurrComp)->IsFreeBoundary())
	 {
          Joint = new TIsoBoundEdge(Domain->GetBdPart(BDpart)->GetBdComp(CurrComp), T_a, T_b);
// 	   cout<<" FreeBoundary"<<endl;
	 }
        else
         { 
          Joint = new TBoundEdge(Domain->GetBdPart(BDpart)->GetBdComp(CurrComp), T_a, T_b);
         }
       }
    }
    else // inner edge
    {
    if (CurrNeib != 2)
        cerr << "Error!!!!!!!! not enough neighbours!" << endl;

    Joint = new TJointEqN(CellTree[Neib[0]], CellTree[Neib[1]]);
    }

    // find the local index for the point 'a' on the cell
    for (j=0;j<3;j++)
      if (Triangles[3*Neib[0]+j] == a) break;

    // find the local index for the point 'b' on the cell
    for (k=0;k<3;k++)
      if (Triangles[3*Neib[0]+k] == b) break;

   k = k*10 + j;

    switch (k)
    {
      case  1:
      case 10:
        j = 0;
        break;
      case 12:
      case 21:
        j = 1;
        break;
      case  2:
      case 20:
        j = 2;
        break;
    }

  CellTree[Neib[0]]->SetJoint(j, Joint);

   if (Neib[1] != -1)
    {
      // find the local index for the point 'a' on the cell
      for (j=0;j<3;j++)
        if (Triangles[3*Neib[1]+j] == a) break;

      // find the local index for the point 'b' on the cell
      for (k=0;k<3;k++)
        if (Triangles[3*Neib[1]+k] == b) break;

      k = k*10 + j;

      switch (k) // j will contain the local index for the current
      {
        case  1:
        case 10:
          j = 0;
          break;
        case 12:
        case 21:
          j = 1;
          break;
        case  2:
        case 20:
          j = 2;
          break;
      }

      CellTree[Neib[1]]->SetJoint(j, Joint);
    }

  if (Joint->GetType() == InterfaceJoint ||
        Joint->GetType() == IsoInterfaceJoint)
      ((TInterfaceJoint *) Joint)->CheckOrientation();
  }
 

  delete [] NewVertices;
  delete [] PointNeighb;
  delete [] In.pointlist;
  delete [] In.pointmarkerlist;
  delete [] In.segmentlist;
  delete [] In.segmentmarkerlist;
  
  if(Out.pointlist!=NULL) {
    free(Out.pointlist); Out.pointlist = NULL;}
  if(Out.pointattributelist!=NULL) { 
    free(Out.pointattributelist); Out.pointattributelist = NULL;}
  if(Out.pointmarkerlist!=NULL) {
    free(Out.pointmarkerlist); Out.pointmarkerlist = NULL;}
  if(Out.trianglelist!=NULL) {
    free(Out.trianglelist); Out.trianglelist = NULL;}
  if(Out.triangleattributelist!=NULL) {
    free(Out.triangleattributelist); Out.triangleattributelist = NULL;}
  if(Out.trianglearealist!=NULL) {
    free(Out.trianglearealist); Out.trianglearealist = NULL;}
  if(Out.neighborlist!=NULL) {
    free(Out.neighborlist); Out.neighborlist = NULL;}
  if(Out.segmentlist!=NULL) {
    free(Out.segmentlist); Out.segmentlist = NULL;}
  if(Out.segmentmarkerlist!=NULL) {
    free(Out.segmentmarkerlist); Out.segmentmarkerlist = NULL;}
  if(Out.holelist!=NULL) {
    free(Out.holelist); Out.holelist = NULL;}
  if(Out.regionlist!=NULL) {
    free(Out.regionlist); Out.regionlist = NULL;}
  if(Out.edgelist!=NULL) {
    free(Out.edgelist); Out.edgelist = NULL;}
  if(Out.edgemarkerlist!=NULL) {
    free(Out.edgemarkerlist); Out.edgemarkerlist = NULL;}
  if(Out.normlist!=NULL) {
    free(Out.normlist); Out.normlist = NULL;} 
  
//======================================================================
// Triangular for grid generation --end
//======================================================================


//   cout<< "tetgen" <<endl;
//   exit(0);
  
 
} // TriaReMeshGen
