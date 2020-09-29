 
// ==========================================================================
// instationary problem
// ==========================================================================
#include <MacroCell.h>
#include <BoundEdge.h>
#include <IsoBoundEdge.h>
// #include <gridgen.h>
#include <IsoInterfaceJoint.h>
#include <BdLine.h>
#include <BdCircle.h>
#include <GridCell.h>
#include <TimeConvDiff2D.h>

// #include <TimeConvDiff2D.h>

#define __HEMKER__
#define  __MOVINGMESH__
// #define __CONSERVATIVEALE__

extern "C"
{
  #include <gridgen.h>    
  void triangulate(char*, struct triangulateio*,
                   struct triangulateio*, struct triangulateio*);
}


void ExampleFile()
{
  OutPut("Example: Hemker.h" << endl) ;
  
#ifdef  __CONSERVATIVEALE__  
  TDatabase::ParamDB->P6 = 1;// con-ALE
  OutPut("Conservative ALE form with, that is, with - div w" << endl);   
#else
  TDatabase::ParamDB->P6=0;// non-conservative ALE
  OutPut("Non-Conservative ALE form, that is, witout - div w term" << endl);  
#endif 
}

// exact solution
void Exact(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void BoundCondition(int BdComp, double t, BoundCond &cond)
{
    switch(BdComp)
    {
	case 0:
	case 1:
	case 2:
	    cond = NEUMANN;
	    break;
	default:
	    cond = DIRICHLET;
    }
}


// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 1:
      value = 0;
      break;
    case 4:
      value = 1;
      break;
    default:
      value = 0;
  }
}

// initial conditon
void InitialCondition(double x,  double y, double *values)
{
  values[0] = 0;
}

void BoundConditionAdjoint(int BdComp, double t, BoundCond &cond)
{
    switch(BdComp)
    {
	case 0:
	case 2:
	case 3:
	    cond = NEUMANN;
	    break;
	default:
	    cond = DIRICHLET;
    }
}


// value of boundary condition
void BoundValueAdjoint(int BdComp, double Param, double &value)
{
    value = 0;
}

void BilinearCoeffs(int n_points, double *x, double *y,
        double **parameters, double **coeffs)
{
  double eps=1./TDatabase::ParamDB->PE_NR;
  double angle = 0, v1, v2;
  int i;
  double *coeff, *param;

  v1 = cos(angle);
  v2 = sin(angle);

//   v1 = 0.;
//   v2 = 0.;
  
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = v1;
    coeff[2] = v2;
    coeff[3] = 0;

    coeff[4] = 0;
  }
}

// kind of boundary condition (for FE space needed)
void GridBoundCondition(int BdComp, double t, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void GridBoundValue(int BdComp, double Param, double &value)
{
  value = 0;
}

void GridCoeffs(int n_points, double *x, double *y,
        double **parameters, double **coeffs)
{
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    coeff[0] = 1;

    coeff[1] = 0;
    coeff[2] = 0;
    coeff[3] = 0;

    coeff[4] = 0;
  }
}

 
void  ModifyBdCoords(int *N_MovVert, TVertex **MovBoundVert, TIsoBoundEdge **Free_Joint,  double * Iso_refX, double t)
{
 int i, j, k, l, m;
   
 double x, y;//, t=TDatabase::TimeDB->CURRENTTIME;
 double disp, RefMaxy=-1e20;
 
 TVertex **Vertices;
 TBoundPart *BoundPart; 
 TBdCircle *UpdateIntface;
 TDomain *Domain = TDatabase::ParamDB->Domain;
 
  BoundPart = Domain->GetBdPart(1);
  UpdateIntface = (TBdCircle*)BoundPart->GetBdComp(0); 
  
  
  disp = 0.5*sin((2./5.)*Pi*t);
//    disp = 0;
   for(i=0;i<N_MovVert[0];i++)
    {
     MovBoundVert[i]->GetCoords(x, y);
     if(RefMaxy<y) RefMaxy=y;
//      if(  fabs(x - 0.970031)<1e-3    && fabs(y - 0.24298)<1e-3)
//         cout<<i<<  "x " << x << " y : "<< y<<   " y_new : "<< y + disp<<endl;
     y += disp;
     MovBoundVert[i]->SetCoords(x, y);    
 
    }
   
 // update iso points
   m = 0;     
   for(i=0;i<N_MovVert[0];i++)
    { 
     Vertices = Free_Joint[i]->GetVertices();
     k = Free_Joint[i]->GetN_Vertices();
   
     for(l=0;l<k;l++)
      { 
       x= Iso_refX[2*m];
       y = Iso_refX[2*m+1] + disp;
 
       Vertices[l]->SetCoords(x, y);    
       m++;
      }      
    }// for(i=0;i<N_MovVert[0];i++) 
 
//    cout<<  N_MovVert[0] <<  " RefMaxy : "<< RefMaxy<<endl;
//    exit(0);
} // ModifyBdCoords



void  GetMovingBoundData(TCollection *coll, int *N_MovVert, TVertex ** &MovBoundVert,
                         TIsoBoundEdge ** &Free_Joint, TBaseCell ** &Free_Cells, int ** &IsoCellEdgeNos, double * &Iso_refX)
{
 int i, j, k, l, m0, m1, m2, N_Cells, comp;
 int ORDER, VSP, temp;
 
 double  TX[2], TY[2], x, y, x1, y1;
  
 TBaseCell *Me, *temp_cell;
 TJoint *Joint;
 TBoundComp *BoundComp;  
 TVertex *temp_vert, **Vertices;

 TIsoBoundEdge *temp_isoJoint;
 
  ORDER = 0;
  VSP = TDatabase::ParamDB->ANSATZ_ORDER;

  if (abs(VSP) > 20)
   {ORDER = abs(VSP) - 20;}
  else if ( abs(VSP) > 10)
    {ORDER = abs(VSP) - 10;}
  else ORDER = abs(VSP);

  N_Cells = coll->GetN_Cells();
  N_MovVert[0] = 0;

    for(j=0;j<N_Cells;j++)
      {
        Me = coll->GetCell(j);
        k = Me->GetN_Edges();
        for(l=0;l<k;l++)
        {
         Joint = Me->GetJoint(l);
         if(Joint->GetType() == IsoBoundEdge)
          N_MovVert[0]++; 
            
          }// endfor l
        }// endfor j

//       cout<<   " N_MovVert[0] " << N_MovVert[0] << endl; 
// exit(0);
     // free bound
     Free_Joint = new TIsoBoundEdge*[N_MovVert[0]];
     MovBoundVert = new TVertex*[N_MovVert[0]];
     Free_Cells = new TBaseCell*[N_MovVert[0]];
     IsoCellEdgeNos[0]  = new int [N_MovVert[0]];
     IsoCellEdgeNos[1]  = new int [N_MovVert[0]];
     Iso_refX = new double [2*(ORDER-1)*N_MovVert[0]];

     m2 = 0;  
     for(j=0;j<N_Cells;j++)
      {
        Me = coll->GetCell(j);
        k = Me->GetN_Edges();
        for(l=0;l<k;l++)
        {
          Joint = Me->GetJoint(l);
          if(Joint->GetType() == IsoBoundEdge)
           {
            Free_Joint[m2] = (TIsoBoundEdge *)Joint;
            MovBoundVert[m2] = Me->GetVertex(l);
            Free_Cells[m2] = Me;
            IsoCellEdgeNos[0][m2] = j;
            IsoCellEdgeNos[1][m2] = l;
            Free_Joint[m2]->GenerateVertices(ORDER-1);
            m2++;
           }
          }// endfor l
         }// endfor j    
         
         
  // sort         
//    Sort_Imping(Free_Cells, MovBoundVert[2], IsoCellEdgeNos[0], IsoCellEdgeNos[1], N_MovVert[2], x, y);

  for(i=0;i<N_MovVert[0]-1;i++)
   {
    m2 = 0; 
    k = Free_Cells[i]->GetN_Edges();
    Free_Cells[i]->GetVertex((IsoCellEdgeNos[1][i]+1) % k)->GetCoords(x, y);

     for(j=i+1;j<N_MovVert[0];j++)
      {
       MovBoundVert[j]->GetCoords(x1, y1);
       if((x==x1) && (y==y1))
        {
         temp_vert = MovBoundVert[j];
         MovBoundVert[j] = MovBoundVert[i+1];
         MovBoundVert[i+1] = temp_vert;
       
         temp_cell = Free_Cells[j];
         Free_Cells[j] = Free_Cells[i+1];
         Free_Cells[i+1] = temp_cell;
       
         temp = IsoCellEdgeNos[0][j];
         IsoCellEdgeNos[0][j] = IsoCellEdgeNos[0][i+1];
         IsoCellEdgeNos[0][i+1] = temp; 
         
         temp = IsoCellEdgeNos[1][j];
         IsoCellEdgeNos[1][j] = IsoCellEdgeNos[1][i+1];
         IsoCellEdgeNos[1][i+1] = temp;        
         
         temp_isoJoint = Free_Joint[j];
         Free_Joint[j] =  Free_Joint[i+1];
         Free_Joint[i+1]= temp_isoJoint;
 
         m2++;
        }
      if(m2) break;  
     }
   }

   //Iso refVertices
   m2 = 0;     
   for(i=0;i<N_MovVert[0];i++)
    { 
     Vertices = Free_Joint[i]->GetVertices();
     k = Free_Joint[i]->GetN_Vertices();
   
     for(l=0;l<k;l++)
      { 
       Vertices[l]->GetCoords(Iso_refX[2*m2], Iso_refX[2*m2+1]);
       m2++;
      }      
    }// for(i=0;i<N_MovVert[0];i++)
   
   
//    print   
//   for(i=0;i<N_MovVert;i++)
//    {
//      MovBoundVert[i]->GetCoords(x, y);
//     cout<<i<<  " x : "<< x << " y : "<< y<<  " Angle of free Vertices "<<(180/Pi)*atan2(y,x)<<endl;
//    }
//    
//    for(i=0;i<N_MovVert;i++)
//    {
//      Vertices = Free_Joint[i]->GetVertices();
//      k = Free_Joint[i]->GetN_Vertices();
//    
//      for(l=0;l<k;l++)
//       { 
//        Vertices[l]->GetCoords(x, y);
//        cout<<i<<  " x : "<< x << " y : "<< y<<  " Angle of free Vertices "<<(180/Pi)*atan2(y,x)<<endl;
//       }
//    }
   
//   exit(0); 


}// GetMovingBoundData


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
  double *Coordinates, *Hole_List;  
  
  double Xi[4] = {-3., 9., 9.,-3.};
  double Yi[4] = {-3.,-3., 3., 3.}; 

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
  

  BoundPart = Domain->GetBdPart(0);
  UpdateBound[0]  = (TBdLine*)BoundPart->GetBdComp(0);
  UpdateBound[1]  = (TBdLine*)BoundPart->GetBdComp(1);
  UpdateBound[2]  = (TBdLine*)BoundPart->GetBdComp(2);
  UpdateBound[3]  = (TBdLine*)BoundPart->GetBdComp(3);
  BoundPart = Domain->GetBdPart(1);
  UpdateIntface = (TBdCircle*)BoundPart->GetBdComp(0);
  
  
//OutPut("MESHGEN_REF_QUALIT " << TDatabase::ParamDB->MESHGEN_REF_QUALITY << endl);

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
  
  N_Interface_Vert = 200;    //Freesurf except end point
  N_Hori  = 60;      // number of horrizontal BD vertices
  N_Verti = 20;       // number of vertical BD vertices
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
  
  In.numberofholes = 1;
  In.holelist = NULL;

  Hole_List = new double[2* In.numberofholes];
  Hole_List[0] = 0.;
  Hole_List[1] = 0.;
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
//     cout<<" x : "<< In.pointlist[2*In_Index] << " y : "<< In.pointlist[2*In_Index+1] <<endl;
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
//     cout<<" x : "<< In.pointlist[2*In_Index] << " y : "<< In.pointlist[2*In_Index+1] <<endl;
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
//     cout<<" x : "<< In.pointlist[2*In_Index] << " y : "<< In.pointlist[2*In_Index+1] <<endl;
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
//     cout<<" x : "<< In.pointlist[2*In_Index] << " y : "<< In.pointlist[2*In_Index+1] <<endl;
    In.pointmarkerlist[In_Index] = CurrComp;
    In.segmentlist[2*In_Index] = In_Index;
    In.segmentlist[2*In_Index+1] = In_Index+1;
    In.segmentmarkerlist[In_Index] = CurrComp;
    In_Index++;
 
   }
  CurrComp++;  
  
  In.segmentlist[2*(In_Index-1)+1] = 0;
  temp=In_Index;
 
  
//   T_a = TDatabase::ParamDB->P1; // x axis value in ellipse
//   T_b = TDatabase::ParamDB->P2; // y axis value in ellipse !!! set in modifyGausscoord funct also
  T_a = 1.;
  T_b = 1.;
//   deviation = fabs( T_a - T_b);  


   C_x = 0.0;  // center of the inner phase
   C_y = 0.0, // center of the inner phase
   phi1 = 0.000000E+0000; // end degree value of interface
   phi2 = 2.*Pi;
 
//   C_x = 0.0;  // center of the inner phase
//   C_y = 0.0, // center of the inner phase
  s = (phi2- phi1)/(double)N_Interface_Vert; 
 
  // points and segments on the interface (marker=2)
//   theta = 0.;
//   t0 = theta;
//   cout << " s " << s << endl;
   for(i=0;i<N_Interface_Vert;i++)
    {
      theta = phi1 + (double)i*s;       
//      cout<<" theta : "<< theta <<endl;
      
//       if(fabs(I_FaceX[i]) < 1e-10 ) I_FaceX[i] = 0.0;
      In.pointlist[2*In_Index] =   T_a*cos(theta);;
      In.pointlist[2*In_Index+1] =  T_b*sin(theta);

//       if(i==0) 
//        {
//         FreeBD_X = I_FaceX[i];   FreeBD_Y = I_FaceY[i]; 
//        cout << " sorting " << FreeBD_X << ' ' << FreeBD_Y<<endl;
//        }

//       if(i==1)
//        {
//         h_interface = sqrt((I_FaceX[i-1]-I_FaceX[i])*(I_FaceX[i-1]-I_FaceX[i]) +
//                             (I_FaceY[i-1]-I_FaceY[i])*(I_FaceY[i-1]-I_FaceY[i]) );
//         OutPut("h_interface " <<h_interface << endl);
//        }
//       cout<<(180./Pi)*theta<< " x : "<< In.pointlist[2*In_Index] << " y : "<< In.pointlist[2*In_Index+1] <<endl;

      In.pointmarkerlist[In_Index] = CurrComp;
      In.segmentlist[2*In_Index] = In_Index;
      In.segmentlist[2*In_Index+1] = In_Index+1;
      In.segmentmarkerlist[In_Index] = CurrComp;
      In_Index++;
    }

  In.segmentlist[2*(In_Index-1)+1] = temp;  
  
//  exit(0);
 
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
    
    
  // call triangle
  triangulate((char*)opts.str().c_str(), &In, &Out, (struct triangulateio *)NULL);
  
  
 
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

//       // set bounding box
//   left = bottom = 1e8;
//   right = top = -1e8;
// 
//    for(i=0;i<In.numberofpoints;i++)
//     {
//       if(left>In.pointlist[2*i]) left = In.pointlist[2*i];
//       if(right<In.pointlist[2*i]) right = In.pointlist[2*i];
//       if(top<In.pointlist[2*i+1]) top = In.pointlist[2*i+1];
//       if(bottom>In.pointlist[2*i+1]) bottom = In.pointlist[2*i+1];
//     }  
 
/* // Solid Bound startx, starty, x length and y length
  UpdateSlipBound[0]->SetParams(Xi[0], Yi[0], Xi[1]-Xi[0],Yi[1]-Yi[0]);
  UpdateSlipBound[1]->SetParams(Xi[1], Yi[1], Xi[2]-Xi[1],Yi[2]-Yi[1]);
  UpdateSlipBound[2]->SetParams(Xi[2], Yi[2], Xi[3]-Xi[2],Yi[3]-Yi[2]);
  UpdateSlipBound[3]->SetParams(Xi[3], Yi[3], Xi[0]-Xi[3],Yi[0]-Yi[3]);

// Free boundary xmid, ymid, radius_a, radius_b, start angle, end angle
  UpdateIntface->SetParams(C_x, C_y, T_a, T_b, phi1, phi2);*/  
  
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

//          cout<<BDpart << " BDpart CurrComp "<< CurrComp <<endl;
 
	 
    
    if (Out.edgemarkerlist[i]) // 0 for inner edges and Boundcomp+1 for Boundedge respect
    {
      CurrComp = Out.edgemarkerlist[i] - 1;
      if (CurrComp >= 100000) CurrComp -= 100000;
      
      BDpart=Domain->GetBdPartID(CurrComp);
      CurrComp= Domain->GetLocalBdCompID(CurrComp);

//      cout<<BDpart << " BDpart CurrComp "<< CurrComp <<endl;
 
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
       
       
//       cout<<BDpart << " BDpart CurrComp "<< CurrComp <<endl;
     
     
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

void ComputeExtremalValues(int N, double *sol, double  *values)
{
   int i;
   double max, min;

   min = 1e10;
   max = -1e10;
   
   for(i=0;i<N;i++)
   {
      if(sol[i]-1 > max)
         max = sol[i]-1;
      if(sol[i] < min)
         min = sol[i];
   }

   values[0] = min;
   values[1] = max;
}
 
/** compute curve of the outflow boundary */
void ComputeOutflowBoundary(int level, TFEFunction2D *ufct)
{
  double h, x=4,values[3],y;
  int i, bound_points = 401;
  h = 6.0/(bound_points-1);
  for (i=0;i<bound_points; i++)
  {
      y = -3+i*h;
      ufct->FindGradient(x,y,values);
      OutPut("cutline " << x << " " <<  y << 
            " " <<  values[0] << endl);
  }
}

// computation of some global errors, only for P1 or Q1 !!!
//
// values[0] : absolute value of largest negative undershoot in 
//             a circle around the cylinder
// values[1] : difference of largest positive value in a circle
//             around the cylinder and 1
// values[2] : absolute value of largest negative undershoot 
//             for x > 2
// values[3] : difference of largest positive value for x>2 and 1
// 
void ComputeLocalExtrema(TFEFunction2D *ufct, double *values)
{
  TBaseCell *cell;
  TCollection *Coll;
  TFESpace2D *FESpace2D;
  TJoint *joint;
  RefTrans2D RefTrans;
  TBaseFunct2D *bf;
  FE2D FE_ID;
  TFE2D *FE_Obj;
  TBoundEdge *boundedge;
  TBoundComp *BoundComp;
  int N_BaseFunct, comp;
  double xi, eta, eps = 1e-6, edgeval[4];
  double *uorig, *uxorig, *uyorig, *uref, *uxiref, *uetaref, u;
  double linfty = 0, l2 = 0, l2interior = 0;
  double *Values, fevalue[4], area, diam, local, x, y, val;
  int *GlobalNumbers, *BeginIndex, index, N_Cells, N_Edges;
  int i, j, k, boundary_cell, *Numbers;
  double extr[4];

  extr[0] = -1;
  extr[1] = -1;
  extr[2] = -1;
  extr[3] = 0;

  FESpace2D = ufct->GetFESpace2D();
  BeginIndex = FESpace2D->GetBeginIndex();
  GlobalNumbers = FESpace2D->GetGlobalNumbers();
  Values = ufct->GetValues();  

  Coll = FESpace2D->GetCollection();
  N_Cells = Coll->GetN_Cells();

  // loop over all edges
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_Edges=cell->GetN_Edges();
    diam = cell->GetDiameter();
    if (N_Edges==3)
	area = diam*diam / 4.0;
    else
	area = diam*diam/2.0;
    for(j=0;j<N_Edges;j++)              // loop over all edges of cell
    {                                  
      fevalue[j] = 0;
    }

    FE_ID = FESpace2D->GetFE2D(i, cell);
    FE_Obj = TFEDatabase2D::GetFE2D(FE_ID);
    RefTrans = FE_Obj->GetRefTransID();

    // get base function object
    bf = FE_Obj->GetBaseFunct2D();
    N_BaseFunct = bf->GetDimension();
    
    uorig = new double[N_BaseFunct];
    uxorig = new double[N_BaseFunct];
    uyorig = new double[N_BaseFunct];
    
    uref = new double[N_BaseFunct];
    uxiref = new double[N_BaseFunct];
    uetaref = new double[N_BaseFunct];
    
    // set cell for reference transformation
    TFEDatabase2D::SetCellForRefTrans(cell, RefTrans);
    for (j=0;j<N_Edges;j++)
    {
      // compute coordinates
      x = cell->GetVertex(j)->GetX();
      y = cell->GetVertex(j)->GetY();
      if (x<-1.5)
	  continue;
      // find local coordinates of the given point
      //cout << " x: " << x << endl;
      //cout << " y: " << y << endl;
      TFEDatabase2D::GetRefFromOrig(RefTrans, x, y, xi, eta);
      //cout << " xi: " << xi << endl;
      //cout << "eta: " << eta << endl;

      bf->GetDerivatives(D00, xi, eta, uref);
      bf->GetDerivatives(D10, xi, eta, uxiref);
      bf->GetDerivatives(D01, xi, eta, uetaref);
      
      TFEDatabase2D::GetOrigValues(RefTrans, xi, eta, bf, Coll, (TGridCell *)cell,
                uref, uxiref, uetaref, uorig, uxorig, uyorig);
      // compute value of fe function at (x,y)
      u = 0;
      Numbers = GlobalNumbers + BeginIndex[i];
      for(k=0;k<N_BaseFunct;k++)
      {
        val = Values[Numbers[k]];
        u += uorig[k]*val;
      }
      //OutPut(x << " " << y << " " << u << endl);
      // strip with the circle
      if ((x>-1.5)&&(x<1.5))
      {
	  if ((u<=0)&&(fabs(u)>extr[0]))
	      extr[0] = fabs(u);
	  if ((u>=1)&&(fabs(u-1)>extr[1]))
	      extr[1] = fabs(u-1);
      }
      if (x>2)
      {
	  if ((u<=0)&&(fabs(u)>extr[2]))
	      extr[2] = fabs(u);
	  if ((u>=1)&&(fabs(u-1)>extr[3]))
	      extr[3] = fabs(u-1);
      }
    } // endfor (j) N_Edges
  
    delete uorig;
    delete uxorig;
    delete uyorig;
    delete uref;
    delete uxiref;
    delete uetaref;
  } // endfor
    
  values[0] = extr[0];
  values[1] = extr[1];
  values[2] = extr[2];
  values[3] = extr[3];
 }

/****************************************************************/
//
// for FEM_TVD
//
/****************************************************************/

void CheckWrongNeumannNodes(TCollection *Coll, TFESpace2D *fespace,
			    int &N_neum_to_diri, int* &neum_to_diri,
			    int* &neum_to_diri_bdry, 
			    double* &neum_to_diri_param)
{
    const int max_entries = 50000;  
    int i, j, N_, min_val, type;
  int N_Cells, N_V, diri_counter = 0, found, diri_counter_1 = 0;
  int *global_numbers, *begin_index, *dof;
  int boundary_vertices[4], tmp_diri[max_entries], tmp_bdry[max_entries];
  double x[4], y[4], eps = 1e-6, tmp_param[max_entries];
  TBaseCell *cell;
  TVertex *vertex;
  FE2D CurrentElement;

  // number of mesh cells
  N_Cells = Coll->GetN_Cells();
  // array with global numbers of d.o.f.
  global_numbers = fespace->GetGlobalNumbers();
  // array which points to the beginning of the global numbers in
  // global_numbers for each mesh cell
  begin_index = fespace->GetBeginIndex();

  diri_counter = 0;
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
     N_V = cell->GetN_Vertices();
    found = 0;
    for (j=0;j<N_V;j++)
    {
      // read coordinates of the mesh cell
      boundary_vertices[j] = 0;
      vertex = cell->GetVertex(j);
      vertex->GetCoords(x[j], y[j]);
      if ((fabs(x[j]+3)<eps)//||(fabs(y[j]+3)<eps)||(fabs(y[j]-3)<eps)
	  || (fabs(sqrt(x[j]*x[j]+y[j]*y[j])-1)<eps))
      {
	   boundary_vertices[j] = 1;
	   found++;
      }
    }
    // no cell with edge with vertex on the boundary
    if (found<2) 
	continue;
       // finite element on the mesh cell
    CurrentElement = fespace->GetFE2D(i, cell);
    // number of basis functions (= number of d.o.f.)
    N_ = TFEDatabase2D::GetN_BaseFunctFromFE2D(CurrentElement);
    // the array which gives the mapping of the local to the global d.o.f.
    dof = global_numbers+begin_index[i];
    switch(CurrentElement)
    {
	// P_1, Q_1
	case C_P1_2D_T_A:
	case C_Q1_2D_Q_A:
	case C_Q1_2D_Q_M:
	    for (j=0;j<N_V;j++)
	    {
		// vertex on the boundary
		if (boundary_vertices[j])
		{
		    if (CurrentElement==C_P1_2D_T_A)
			tmp_diri[diri_counter] = dof[j];
		    else
		    {
			if (j<2){
			    tmp_diri[diri_counter] = dof[j];
			}
			else
			{
			    if (j==2)
				tmp_diri[diri_counter] = dof[3];
			    else
				tmp_diri[diri_counter] = dof[2];
			}
		    }
		    if (diri_counter > max_entries)
		    {
			OutPut("tmp_diri too short !!!"<<endl);
			exit(4711);
		    }
		    // inflow x = -3
		    if (fabs(x[j]+3)<eps) 
		    {
			tmp_bdry[diri_counter] = 3;
			tmp_param[diri_counter] = (-y[j]+3)/6.0;
		    }
		    // circle
		    if (fabs(sqrt(x[j]*x[j]+y[j]*y[j])-1)<eps) 
		    {
			tmp_bdry[diri_counter] = 4;
			// parameter does not matter, since b.c. equal to 1
			tmp_param[diri_counter] =  0;
		    }
		    diri_counter++;
		}
	    }
	    break;
	// P_2, Q_2
	case C_P2_2D_T_A:
	case C_Q2_2D_Q_A:
	case C_Q2_2D_Q_M:
            // loop over the edges
 	    for (j=0;j<N_V;j++)
	    {
              // check of edge j is on boundary  
              if (boundary_vertices[j] && boundary_vertices[(j+1)%N_V])
              {
		// check if this is a boundary edge
		type = cell->GetJoint(j)->GetType();
		if (!((type == BoundaryEdge)||(type == IsoBoundEdge)))
		  continue;
	        switch(j)
                {
                   case 0:
                     tmp_diri[diri_counter] = dof[0];
                     tmp_diri[diri_counter+1] = dof[1];
                     tmp_diri[diri_counter+2] = dof[2];
                   break;
                  case 1:
                     if (N_V==3)
                     {
                       tmp_diri[diri_counter] = dof[2];
                       tmp_diri[diri_counter+1] = dof[4];
                       tmp_diri[diri_counter+2] = dof[5];
                     }
                     else
                     {
                       tmp_diri[diri_counter] = dof[2];
                       tmp_diri[diri_counter+1] = dof[5];
                       tmp_diri[diri_counter+2] = dof[8];
                     }
                   break;
                  case 2:
                     if (N_V==3)
                     {
                       tmp_diri[diri_counter] = dof[5];
                       tmp_diri[diri_counter+1] = dof[3];
                       tmp_diri[diri_counter+2] = dof[0];
                     }
                     else
                     {
                       tmp_diri[diri_counter] = dof[8];
                       tmp_diri[diri_counter+1] = dof[7];
                       tmp_diri[diri_counter+2] = dof[6];
                     }
                   break;
                   case 3:
                     tmp_diri[diri_counter] = dof[6];
                     tmp_diri[diri_counter+1] = dof[3];
                     tmp_diri[diri_counter+2] = dof[0];
                   break;

                }
              
		if (diri_counter+2 > max_entries)
		{
			OutPut("tmp_diri too short !!!"<<endl);
			exit(4711);
		}

		// inflow x = -3
		if ((fabs(x[j]+3)<eps)&&(fabs(x[(j+1)%N_V]+3)<eps)) 
		{
		    tmp_bdry[diri_counter] = 3;
		    tmp_bdry[diri_counter+1] = 3;
		    tmp_bdry[diri_counter+2] = 3;
		    tmp_param[diri_counter] = (-y[j]+3)/6.0;
		    tmp_param[diri_counter+2] = (-y[(j+1)%N_V]+3)/6.0;
		    tmp_param[diri_counter+1] = (tmp_param[diri_counter] +  tmp_param[diri_counter+2])/2.0;
		}
		// circle
		if ((fabs(sqrt(x[j]*x[j]+y[j]*y[j])-1)<eps) && 
		    (fabs(sqrt(x[(j+1)%N_V]*x[(j+1)%N_V]+y[(j+1)%N_V]*y[(j+1)%N_V])-1)<eps))
		{
		    tmp_bdry[diri_counter] = 4;
		    tmp_bdry[diri_counter+1] = 4;
		    tmp_bdry[diri_counter+2] = 4;
		    // parameter does not matter, since b.c. equal to 1
		    tmp_param[diri_counter] = 0;
		    tmp_param[diri_counter+1] = 0;
		    tmp_param[diri_counter+2] = 0;
		}
		diri_counter +=3;
	      }
	    }
	    break;
	// P_3, Q_3
	case C_P3_2D_T_A:
	case C_Q3_2D_Q_A:
	case C_Q3_2D_Q_M:
            // loop over the edges
 	    for (j=0;j<N_V;j++)
	    {
              // check of edge j is on boundary  
              if (boundary_vertices[j] && boundary_vertices[(j+1)%N_V])
              {
		// check if this is a boundary edge
		type = cell->GetJoint(j)->GetType();
		if (!((type == BoundaryEdge)||(type == IsoBoundEdge)))
		  continue;

               // P3: local dof 0, 1, 2, 3 are on the boundary
               // Q3: local dof 0, 1, 2, 3 are on the boundary
	        switch(j)
                {
                   case 0:
                     tmp_diri[diri_counter] = dof[0];
                     tmp_diri[diri_counter+1] = dof[1];
                     tmp_diri[diri_counter+2] = dof[2];
		     tmp_diri[diri_counter+3] = dof[3];
                   break;
                  case 1:
                     if (N_V==3)
                     {
                       tmp_diri[diri_counter] = dof[3];
                       tmp_diri[diri_counter+1] = dof[6];
                       tmp_diri[diri_counter+2] = dof[8];
		       tmp_diri[diri_counter+3] = dof[9];
                     }
                     else
                     {
                       tmp_diri[diri_counter] = dof[3];
                       tmp_diri[diri_counter+1] = dof[7];
                       tmp_diri[diri_counter+2] = dof[11];
		       tmp_diri[diri_counter+3] = dof[15];
                     }
                   break;
                  case 2:
                     if (N_V==3)
                     {
                       tmp_diri[diri_counter] = dof[9];
                       tmp_diri[diri_counter+1] = dof[7];
                       tmp_diri[diri_counter+2] = dof[4];
                       tmp_diri[diri_counter+3] = dof[0];
		     }
                     else
                     {
                       tmp_diri[diri_counter] = dof[15];
                       tmp_diri[diri_counter+1] = dof[14];
                       tmp_diri[diri_counter+2] = dof[13];
			tmp_diri[diri_counter+3] = dof[12];
                     }
                   break;
                   case 3:
                     tmp_diri[diri_counter] = dof[12];
                     tmp_diri[diri_counter+1] = dof[8];
                     tmp_diri[diri_counter+2] = dof[4];
		     tmp_diri[diri_counter+3] = dof[0];
                   break;
                }
              
		if (diri_counter+3 > max_entries)
		{
			OutPut("tmp_diri too short !!!"<<endl);
			exit(4711);
		}

		// inflow x = -3
		if ((fabs(x[j]+3)<eps)&&(fabs(x[(j+1)%N_V]+3)<eps)) 
		{
		    tmp_bdry[diri_counter] = 3;
		    tmp_bdry[diri_counter+1] = 3;
		    tmp_bdry[diri_counter+2] = 3;
		    tmp_bdry[diri_counter+3] = 3;
		    tmp_param[diri_counter] = (-y[j]+3)/6.0;
		    tmp_param[diri_counter+3] = (-y[(j+1)%N_V]+3)/6.0;
		    tmp_param[diri_counter+1] = (2*tmp_param[diri_counter] +  tmp_param[diri_counter+3])/3.0;
		    tmp_param[diri_counter+2] = (tmp_param[diri_counter] +  2*tmp_param[diri_counter+3])/2.0;
		}
		// circle
		if ((fabs(sqrt(x[j]*x[j]+y[j]*y[j])-1)<eps) && 
		    (fabs(sqrt(x[(j+1)%N_V]*x[(j+1)%N_V]+y[(j+1)%N_V]*y[(j+1)%N_V])-1)<eps))
		{
		    tmp_bdry[diri_counter] = 4;
		    tmp_bdry[diri_counter+1] = 4;
		    tmp_bdry[diri_counter+2] = 4;
		    tmp_bdry[diri_counter+3] = 4;
		    // parameter does not matter, since b.c. equal to 1
		    tmp_param[diri_counter] = 0;
		    tmp_param[diri_counter+1] = 0;
		    tmp_param[diri_counter+2] = 0;
		    tmp_param[diri_counter+3] = 0;
		}
		diri_counter +=4;
	      }
	    }
	    break;
	default:
	    OutPut("CheckNeumannNodesForVelocity not implemented for element "
		   << CurrentElement << endl);
	    OutPut("code can be run without CheckNeumannNodesForVelocity, just delete the exit" << endl);
	    exit(4711);
    }	    
  }
 
  // condense
  for (i=0;i<diri_counter;i++)
  {
      if (tmp_diri[i] == -1)
	  continue;
      diri_counter_1++;
      for (j=i+1;j<diri_counter;j++)
      {
	  if (tmp_diri[i] == tmp_diri[j])
	  {
	      tmp_diri[j] = -1;
	  }
      }
  }

  //OutPut("CheckNeumannNodesForVelocity: N_neum_to_diri " << diri_counter_1 << endl);
  N_neum_to_diri = diri_counter_1;
  // allocate array for the indices
  neum_to_diri = new int[diri_counter_1];
  // allocate array for the corresponding boundary numbers
  neum_to_diri_bdry = new int[diri_counter_1];
  // allocate array for the corresponding boundary parameters
  neum_to_diri_param = new double[diri_counter_1];
  // fill array and sort
  for (i=0;i<diri_counter_1;i++)
  {
      min_val = tmp_diri[0];
      found = 0;
      for (j=1;j<diri_counter;j++)
      {
	  if ((tmp_diri[j]>0) && ((tmp_diri[j] < min_val) || 
				  (min_val == -1)))
	  {
	       min_val =  tmp_diri[j];
	       found = j;
	  }
      }
      neum_to_diri[i] = tmp_diri[found];
      neum_to_diri_bdry[i] = tmp_bdry[found];
      neum_to_diri_param[i] = tmp_param[found];
      tmp_diri[found] = -1;
  }
/*
  for (i=0;i<diri_counter_1;i++)
  {
      OutPut(i << " " << neum_to_diri[i] << " " << neum_to_diri_bdry[i] <<
	     " " << neum_to_diri_param[i] <<  endl);
  }
*/
}

void SetDirichletNodesFromNeumannNodes(TSquareMatrix2D **SQMATRICES, 
				       double *rhs,
				       int N_neum_to_diri,
				       int *neum_to_diri,
				       int *neum_to_diri_bdry,
				       double *neum_to_diri_param)
{
    TSquareMatrix2D *MatrixA;
    double *Entries_A, value;
    int i, j, l, l0, l1, index, *RowPtr_A, *KCol_A;
    
    MatrixA = SQMATRICES[0];
    RowPtr_A      = MatrixA->GetRowPtr();
    KCol_A        = MatrixA->GetKCol();
    Entries_A     = MatrixA->GetEntries();
    // loop over dof to change
    for (i=0;i<N_neum_to_diri;i++)
    {
	index = neum_to_diri[i];
	l0 = RowPtr_A[index];
	l1 = RowPtr_A[index+1];
	for (l=l0;l<l1;l++)
	{
	    // diagonal entry
	    if (KCol_A[l]==index)
		Entries_A[l] = 1;  
	    else
		Entries_A[l] = 0;
	}
	// set boundary condition
	BoundValue(neum_to_diri_bdry[i], neum_to_diri_param[i], rhs[index]);
    }
}

void ComputeDifferenceToCoarseLevel(TCollection *Coll_fine,
				    TCollection *Coll_coarse,
				    TFEFunction2D *u_fine, 
				    TFEFunction2D *u_coarse)
{
    int i, j, k, N_Cells, N_Edges, coarse_no;
    double x, y, x_c, y_c, val_fine[4], val_coarse[4], c1err = -1, c1err_coarse = -1;
    double x_err, y_err, x_err_c, y_err_c;
    TBaseCell *cell, *parent;
    
    // number of cells
    N_Cells = Coll_fine->GetN_Cells();
    
    // loop over all edges
    for(i=0;i<N_Cells;i++)
    {
	// cell
	cell = Coll_fine->GetCell(i);
	// parent cell
	parent = cell->GetParent();
	coarse_no = Coll_coarse->GetIndex(parent);
	//OutPut(coarse_no << " ");
	// number of edges
	N_Edges=cell->GetN_Edges();
	for (j=0;j<N_Edges;j++)
	{
	    cell->GetVertex(j)->GetCoords(x, y);
	    u_fine->FindGradientLocal(cell, i, x, y, val_fine);
	    u_coarse->FindGradientLocal(parent, coarse_no, x, y, val_coarse);
	    if (fabs(val_fine[0] - val_coarse[0]) > c1err)
	    {
		c1err = fabs(val_fine[0] - val_coarse[0]);
		x_err = x;
		y_err = y;
	    }
	    for (k=0;k<N_Edges;k++)
	    {
		parent->GetVertex(k)->GetCoords(x_c, y_c);
		if ((fabs(x_c -x ) < 1e-6) && (fabs(y_c -y ) < 1e-6))
		{
		    if (fabs(val_fine[0] - val_coarse[0]) > c1err_coarse)
		    {
			c1err_coarse = fabs(val_fine[0] - val_coarse[0]);
			x_err_c = x;
			y_err_c = y;
		    }
		}
	    }
	}
    } 

    //OutPut("C1 error f " << c1err  << " \\& " << x_err <<"," << y_err << endl);
    //OutPut(" C1 error c " << c1err_coarse  << " \\& " << x_err_c <<"," << y_err_c << endl);

    OutPut("C1 error f "<<" & " << c1err  << " &  ( " << x_err <<"," << y_err << ")"<< "\\\\\\hline" << endl);
    OutPut("C1 error c "<< " & "<< c1err_coarse  << " &  ( " << x_err_c <<"," << y_err_c << ")" << "\\\\\\hline"<< endl);
}
void ComputeCutLines_X(TCollection *Coll, TFEFunction2D *ufct, int level)
{
    double h, val[3], x, y, *cutvalues, y0=0;
  int i, j, N_Cells, bound_points = 10001;
    TBaseCell *cell;
  
  h = 3.0/(bound_points-1);

  cutvalues = new double[6*bound_points];
  memset(cutvalues , 0 , 6*bound_points*SizeOfDouble);

  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {   
    cell = Coll->GetCell(i);
    x = -1;
    for (j=0;j<bound_points;j++)
    {
	y = y0 + h * j;
	cutvalues[j] = y;
	if(cell->PointInCell(x,y))
	{
	    ufct->FindGradientLocal(cell, i, x, y, val);
	    cutvalues[j+bound_points] = val[0];
	}
    }
    x = 0;
    for (j=0;j<bound_points;j++)
    {
	y = y0 + h * j;
	if(cell->PointInCell(x,y))
	{
	    ufct->FindGradientLocal(cell, i, x, y, val);
	    cutvalues[j+2*bound_points] = val[0];
	}
    }
    x = 1;
    for (j=0;j<bound_points;j++)
    {
	y = y0 + h * j;
	if(cell->PointInCell(x,y))
	{
	    ufct->FindGradientLocal(cell, i, x, y, val);
	    cutvalues[j+3*bound_points] = val[0];
	}
    }
    x = 2;
    for (j=0;j<bound_points;j++)
    {
	y = y0 + h * j;
	if(cell->PointInCell(x,y))
	{
	    ufct->FindGradientLocal(cell, i, x, y, val);
	    cutvalues[j+4*bound_points] = val[0];
	}
    }
    x = 4;
    for (j=0;j<bound_points;j++)
    {
	y = y0 + h * j;
	if(cell->PointInCell(x,y))
	{
	    ufct->FindGradientLocal(cell, i, x, y, val);
	    cutvalues[j+5*bound_points] = val[0];
	}
    }
  }

  for (j=0;j<bound_points;j++)
  {
      OutPut("cutx"<< level << " " << cutvalues[j] << " " << cutvalues[j+bound_points] 
	      << " " << cutvalues[j+2*bound_points]  << " " << cutvalues[j+3*bound_points]  << " " 
	     << cutvalues[j+4*bound_points]  << " " << cutvalues[j+5*bound_points] << endl);
  }
}

void ComputeCutLines_Y(TCollection *Coll, TFEFunction2D *ufct, int level)
{
    double h, val[3], x, y, *cutvalues;
  int i, j, N_Cells, bound_points = 20001;
    TBaseCell *cell;
  
  h = 10.0/(bound_points-1);

  cutvalues = new double[3*bound_points];
  memset(cutvalues , 0 , 3*bound_points*SizeOfDouble);
  
  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {   
    cell = Coll->GetCell(i);
    y = 0;
    for (j=0;j<bound_points;j++)
    {
	x = -2 + h * j;
	cutvalues[j] = x;
	if(cell->PointInCell(x,y))
	{
	    ufct->FindGradientLocal(cell, i, x, y, val);
	    cutvalues[j+bound_points] = val[0];
	}
    }
    y = 1;
    for (j=0;j<bound_points;j++)
    {
	x = -2 + h * j;
	if(cell->PointInCell(x,y))
	{
	    ufct->FindGradientLocal(cell, i, x, y, val);
	    /*if (cutvalues[j+2*bound_points]!=4711)
	    {
		OutPut("belegt"<<endl);
		exit(1);
	    }
	    else*/
	    cutvalues[j+2*bound_points] = val[0];
	}
    }
  }

  for (j=0;j<bound_points;j++)
  {
      OutPut("cuty"<< level << " " << cutvalues[j] << " " << cutvalues[j+bound_points] 
	      << " " << cutvalues[j+2*bound_points] << endl);
  }
}

void ComputeCutLines_epsY(TCollection *Coll, TFEFunction2D *ufct, int level)
{
    double h, val[3], x, y, *cutvalues;
  int i, j, N_Cells, bound_points = 20001;
    TBaseCell *cell;
    double eps = 1.0/TDatabase::ParamDB->RE_NR;
  
  h = 10.0/(bound_points-1);

  cutvalues = new double[3*bound_points];
  memset(cutvalues , 0 , 3*bound_points*SizeOfDouble);

  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {   
    cell = Coll->GetCell(i);
    y = 1-sqrt(eps);
    for (j=0;j<bound_points;j++)
    {
	x = -2 + h * j;
	cutvalues[j] = x;
	if(cell->PointInCell(x,y))
	{
	    ufct->FindGradientLocal(cell, i, x, y, val);
	    cutvalues[j+bound_points] = val[0];
	}
    }
    y = 1+sqrt(eps);
    for (j=0;j<bound_points;j++)
    {
	x = -2 + h * j;
	if(cell->PointInCell(x,y))
	{
	    ufct->FindGradientLocal(cell, i, x, y, val);
	    cutvalues[j+2*bound_points] = val[0];
	}
    }
  }

  for (j=0;j<bound_points;j++)
  {
      OutPut("cutyeps"<< level << " " << cutvalues[j] << " " << cutvalues[j+bound_points] 
	      << " " << cutvalues[j+2*bound_points] << endl);
  }
}
void ComputeCutLines_eps_radial(TCollection *Coll, TFEFunction2D *ufct, int level)
{
    double h, val[3], x, y, *cutvalues, tmp, r;
  int i, j, N_Cells, bound_points = 10001;
    TBaseCell *cell;
    double eps = 1.0/TDatabase::ParamDB->RE_NR;
  
  h = 2*Pi/(bound_points-1);

  cutvalues = new double[10*bound_points];
  memset(cutvalues , 0 , 10*bound_points*SizeOfDouble);

  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {   
    cell = Coll->GetCell(i);
    r = (1+eps)*(1+eps);
    for (j=0;j<bound_points;j++)
    {
	x = r*cos(h * j-Pi);
	y = r*sin(h*j-Pi);
	cutvalues[j] = h*j -Pi;
	cutvalues[j+bound_points] = x;
	cutvalues[j+2*bound_points] = y;
	if(cell->PointInCell(x,y))
	{
	    ufct->FindGradientLocal(cell, i, x, y, val);
	    cutvalues[j+3*bound_points] = val[0];
	}
    }
    tmp=pow(eps,2.0/3.0);
    r = (1+tmp)*(1+tmp);
    for (j=0;j<bound_points;j++)
    {
	x = r*cos(h * j-Pi);
	y = r*sin(h*j-Pi);
	cutvalues[j+4*bound_points] = x;
	cutvalues[j+5*bound_points] = y;
	if(cell->PointInCell(x,y))
	{
	    ufct->FindGradientLocal(cell, i, x, y, val);
	    cutvalues[j+6*bound_points] = val[0];
	}
    }
    tmp=sqrt(eps);
    r = (1+tmp)*(1+tmp);
    for (j=0;j<bound_points;j++)
    {
	x = r*cos(h * j-Pi);
	y = r*sin(h*j-Pi);
	cutvalues[j+7*bound_points] = x;
	cutvalues[j+8*bound_points] = y;
	if(cell->PointInCell(x,y))
	{
	    ufct->FindGradientLocal(cell, i, x, y, val);
	    cutvalues[j+9*bound_points] = val[0];
	}
    }

  }

  for (j=0;j<bound_points;j++)
  {
      OutPut("cutradeps"<< level << " " << cutvalues[j] << " " << cutvalues[j+bound_points] 
	      << " " << cutvalues[j+2*bound_points] << " " << cutvalues[j+3*bound_points] 
	     << " " << cutvalues[j+4*bound_points] << " " << cutvalues[j+5*bound_points] 
	     << " " << cutvalues[j+6*bound_points] << " " << cutvalues[j+7*bound_points] 
	     << " " << cutvalues[j+8*bound_points] << " " << cutvalues[j+9*bound_points]  << endl);
  }
}
