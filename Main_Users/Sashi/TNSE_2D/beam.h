// Navier-Stokes problem, Benchmark problem
// 
// u(x,y) = unknown
// p(x,y) = unknown

void ExampleFile()
{
  OutPut("Example: beam_masudslip.h" << endl) ;
}

// #define __BENCH__
#define __BEAM__
// #define __HEMKER__

#include <Joint.h>
#include <BoundEdge.h>
#include <BoundComp.h>
#include <FE2D.h>
#include <FEDesc2D.h>


extern "C"
{
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
//   values[0] = values[0];
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int i, double t, BoundCond &cond)
{
  switch(i)
  {
    case 1:
      cond = NEUMANN;
      TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
      break;
    default:
      cond = DIRICHLET;
      break;
  }   
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0: value = 1;
            break;
    case 1: value=0; // 6
            break;
    case 2: value = 1;
            break;
    case 3: value=1; // 6
            break;
    case 4: value=0;
            break;
    default: cout << "wrong boundary part number: " << BdComp << endl;
  }
//   double t=TDatabase::TimeDB->CURRENTTIME;

//   switch(BdComp)
//   {
//     case 3: 
// //        value=sin(Pi*t/20)*10*Param*(1-Param);
//          value= 2*Param*(1-Param); // 6
//        break;
//     default: 
//        value = 0;
//        break;
//   }
  
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  value = 0;
  if(BdComp>4) cout << "wrong boundary part number: " << BdComp << endl;
}

// ========================================================================
// Description for moving grid
// ========================================================================
void GridU1(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void GridU2(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}
 
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

void ModifyCoords(double x, double y, double &X, double &Y, double t)
{
// //  double t = TDatabase::TimeDB->CURRENTTIME; 
// //  double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH; 
//  double r,d, norm, r1, n1, n2;
//  
//   d = 2. - cos(20.*Pi*t);
//  
//  X = x*d;
//  Y = y*d;
   
 X = x;
 Y = y;
  }

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *x, double *y,
               double **parameters, double **coeffs)
{
  double eps = 1./TDatabase::ParamDB->RE_NR;
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

 
void ModifyBdCoords(int N_MovVert, TVertex **MovBoundVert, TIsoBoundEdge **Free_Joint, 
		    double * Iso_refX, double t)
{
 int i, j, k, l, m;
   
 double x, y, Fact=0.1;
 double disp, theta;
 
 
 TVertex **Vertices;

   for(i=0;i<N_MovVert;i++)
    {
     MovBoundVert[i]->GetCoords(x, y);
     
     if( fabs(0.545-x)<1e-8 )
       continue;

	     disp = Fact*(x-0.545)*(x-0.545)*sin(2.*Pi*(1/7.5)*t);
     theta = atan2(y,(x-0.545));
     x += (-0.25*disp*tan(theta) -y*sin(theta));
     y += disp; 
     
     MovBoundVert[i]->SetCoords(x, y);    
    }
   
 // update iso points
   m = 0;     
   for(i=0;i<N_MovVert;i++)
    { 
     Vertices = Free_Joint[i]->GetVertices();
     k = Free_Joint[i]->GetN_Vertices();
   
     for(l=0;l<k;l++)
      { 
       x = Iso_refX[2*m];
       y = Iso_refX[2*m+1];
             disp = Fact*(x-0.545)*(x-0.545)*sin(2.*Pi*((2.*Pi)/7.5)*t);
	     
       theta = atan2(y,(x-0.545));
       x += (-0.25*disp*tan(theta) -y*sin(theta));
       y += disp;          
 
       Vertices[l]->SetCoords(x, y);    
       m++;
      }      
    }// for(i=0;i<N_MovVert;i++) 
 
//    cout<<   " RefMaxy : "<< RefMaxy<<endl;
//    exit(0);
} // ModifyBdCoords

// sorting curved surfrace vertices - general case
void Sort_Imping(TBaseCell **cell, TVertex **Vertex, int *CellNo, int *EdgeNo, int N,
                     double X0, double Y0 )
{
 int i, j, k, temp, test=0;
 
 double x, y, x1, y1;
 
 TVertex *temp_vert;
 TBaseCell *temp_cell;

 // finding the right(starting) vertex
  for(i=0;i<N;i++)
   {
    Vertex[i]->GetCoords(x, y);
    if(  sqrt((x-X0)*(x-X0) +(y-Y0)* (y-Y0))<1e-5  )
      {
//        cout << " sorting " << x << ' ' << y<<endl;
       temp_vert = Vertex[0];
       Vertex[0] = Vertex[i];
       Vertex[i] = temp_vert;

       temp_cell = cell[0];
       cell[0] = cell[i];
       cell[i] = temp_cell;
       
       temp = EdgeNo[0];
       EdgeNo[0] = EdgeNo[i];
       EdgeNo[i] = temp;
       
       temp = CellNo[0];
       CellNo[0] = CellNo[i];
       CellNo[i] = temp;        

       test++;
      }
    if(test) break;
   }

  if(i==N)
  {
   cout<< "Error in finding start vert " << X0 << endl;   
   exit(0);
  }

  for(i=0;i<N-1;i++)
   {
    test = 0; 
    k = cell[i]->GetN_Edges();
    cell[i]->GetVertex((EdgeNo[i]+1) % k)->GetCoords(x, y);

     for(j=i+1;j<N;j++)
      {
       Vertex[j]->GetCoords(x1, y1);
       if((x==x1) && (y==y1))
        {
         temp_vert = Vertex[j];
         Vertex[j] = Vertex[i+1];
         Vertex[i+1] = temp_vert;
       
         temp_cell = cell[j];
         cell[j] = cell[i+1];
         cell[i+1] = temp_cell;
       
         temp = EdgeNo[j];
         EdgeNo[j] = EdgeNo[i+1];
         EdgeNo[i+1] = temp; 
         
         temp = CellNo[j];
         CellNo[j] = CellNo[i+1];
         CellNo[i+1] = temp;          
         
         test++;
        }
      if(test) break;  
     }
   }

//    print   
//   for(i=0;i<N;i++)
//    {
//     Vertex[i]->GetCoords(x, y);
//     cout<<i<<  " x : "<< x << " y : "<< y <<endl;
// //     cout<<i<<  " x : "<< x << " y : "<< y<<  " Angle of free Vertices "<<(180/Pi)*atan2(y,x)<<endl;
//    }
//   exit(0); 
}

  
void  GetMovingBoundData(TCollection *coll, int &N_MovVert, TVertex ** &MovBoundVert,
                         TIsoBoundEdge ** &Free_Joint, TBaseCell ** &Free_Cells,
			 int ** &IsoCellEdgeNos, double * &Iso_refX)
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
  N_MovVert = 0;

    for(j=0;j<N_Cells;j++)
      {
        Me = coll->GetCell(j);
        k = Me->GetN_Edges();
        for(l=0;l<k;l++)
        {
         Joint = Me->GetJoint(l);
	 //loop over all the cells and edges and consider only those as moving vertices , which are of type 
	 //IsoBoundEdge, and and increment the number of moving vertices
         if(Joint->GetType() == IsoBoundEdge)
          {
           (Me->GetVertex(l))->GetCoords(x, y);
           (Me->GetVertex((l+1)%k))->GetCoords(x1, y1);   
	   
           if(((x+x1)/2.)>0.545)
	    {
//              cout << "x " << x << " y " << y <<endl;
             N_MovVert++; 
	    }
           } 
          }// endfor l
        }// endfor j


     // free bound
     Free_Joint = new TIsoBoundEdge*[N_MovVert];
     MovBoundVert = new TVertex*[N_MovVert];
     Free_Cells = new TBaseCell*[N_MovVert];
     IsoCellEdgeNos[0]  = new int [N_MovVert];
     IsoCellEdgeNos[1]  = new int [N_MovVert];
     Iso_refX = new double [2*(ORDER-1)*N_MovVert];
     
 //Once the information about the moving vertices are obtained,pass the info to the function Sort_Imping
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
           (Me->GetVertex(l))->GetCoords(x, y);
           (Me->GetVertex((l+1)%k))->GetCoords(x1, y1);       

           if(((x+x1)/2.)>0.545)     
           {
//              cout << "x " << x << " y " << y <<endl;	     
            Free_Joint[m2] = (TIsoBoundEdge *)Joint;
            MovBoundVert[m2] = Me->GetVertex(l);
            Free_Cells[m2] = Me;
            IsoCellEdgeNos[0][m2] = j;
            IsoCellEdgeNos[1][m2] = l;
            Free_Joint[m2]->GenerateVertices(ORDER-1);
            m2++;
            }
           } //
          }// endfor l
         }// endfor j    
         
//         cout<<  m2 << " N_MovVert " << N_MovVert << endl; 
   
  // sort   
   x=0.545; y=0.1;
   Sort_Imping(Free_Cells, MovBoundVert, IsoCellEdgeNos[0], IsoCellEdgeNos[1], N_MovVert, x, y);
   
// exit(0);    
//   for(i=0;i<N_MovVert-1;i++)
//    {
//     m2 = 0; 
//     k = Free_Cells[i]->GetN_Edges();
//     Free_Cells[i]->GetVertex((IsoCellEdgeNos[1][i]+1) % k)->GetCoords(x, y);
// 
//      for(j=i+1;j<N_MovVert;j++)
//       {
//        MovBoundVert[j]->GetCoords(x1, y1);
//        if((x==x1) && (y==y1))
//         {
//          temp_vert = MovBoundVert[j];
//          MovBoundVert[j] = MovBoundVert[i+1];
//          MovBoundVert[i+1] = temp_vert;
//        
//          temp_cell = Free_Cells[j];
//          Free_Cells[j] = Free_Cells[i+1];
//          Free_Cells[i+1] = temp_cell;
//        
//          temp = IsoCellEdgeNos[0][j];
//          IsoCellEdgeNos[0][j] = IsoCellEdgeNos[0][i+1];
//          IsoCellEdgeNos[0][i+1] = temp; 
//          
//          temp = IsoCellEdgeNos[1][j];
//          IsoCellEdgeNos[1][j] = IsoCellEdgeNos[1][i+1];
//          IsoCellEdgeNos[1][i+1] = temp;        
//          
//          temp_isoJoint = Free_Joint[j];
//          Free_Joint[j] =  Free_Joint[i+1];
//          Free_Joint[i+1]= temp_isoJoint;
//  
//          m2++;
//         }
//       if(m2) break;  
//      }
//    }

   //Iso refVertices
   m2 = 0;     
   for(i=0;i<N_MovVert;i++)
    { 
     Vertices = Free_Joint[i]->GetVertices();
     k = Free_Joint[i]->GetN_Vertices();
   
     for(l=0;l<k;l++)
      { 
       Vertices[l]->GetCoords(Iso_refX[2*m2], Iso_refX[2*m2+1]);
       m2++;
      }      
    }// for(i=0;i<N_MovVert;i++)
   
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
  double *Coordinates, *Hole_List, *S_BX, *S_BY;  
  
  
  
//   double Xi[4] = {-3., 9., 9.,-3.};
//   double Yi[4] = {-3.,-3., 3., 3.}; 
//   double Xi[4] = {0., 10, 10,0};
//   double Yi[4] = {0.,0., 4.1, 4.1}; 
  double Xi[4] = {-10., 36., 36.,-10.};
  double Yi[4] = {-10.,-10., 10., 10.};
  

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
  opts<<'Y'; // Supress adding vertices on boundary edges
  opts<<ends;
  
//   N_Interface_Vert = int (TDatabase::ParamDB->P6);    //Freesurf except end point
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
//   Hole_List[0] = 0.;
//   Hole_List[1] = 0.;
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
     cout<<" x : "<< In.pointlist[2*In_Index] << " y : "<< In.pointlist[2*In_Index+1] <<endl;
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
    cout<<" x : "<< In.pointlist[2*In_Index] << " y : "<< In.pointlist[2*In_Index+1] <<endl;
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
     cout<<" x : "<< In.pointlist[2*In_Index] << " y : "<< In.pointlist[2*In_Index+1] <<endl;
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
     cout<<" x : "<< In.pointlist[2*In_Index] << " y : "<< In.pointlist[2*In_Index+1] <<endl;
    In.pointmarkerlist[In_Index] = CurrComp;
    In.segmentlist[2*In_Index] = In_Index;
    In.segmentlist[2*In_Index+1] = In_Index+1;
    In.segmentmarkerlist[In_Index] = CurrComp;
    In_Index++;
 
   }
  CurrComp++;  
  
  In.segmentlist[2*(In_Index-1)+1] = 0;
  temp=In_Index;
 
  for(i=0;i<N_Interface_Vert;i++) // without last point
   { 
      In.pointlist[2*In_Index] = S_BX[i];
      In.pointlist[2*In_Index+1] = S_BY[i];
      cout<<" x : "<< S_BX[i] << " y : "<< S_BY[i]<<endl;
      In.pointmarkerlist[In_Index] = CurrComp;
      In.segmentlist[2*In_Index] = In_Index;
      In.segmentlist[2*In_Index+1] = In_Index+1;
//       if (AllowEdgeRef)
       { In.segmentmarkerlist[In_Index] = CurrComp; }
//       else
//        { In.segmentmarkerlist[In_Index] = 100000 + CurrComp; }

      In_Index++;  
   }

  In.segmentlist[2*(In_Index-1)+1] = temp;  
  
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

 



 

/** calculate characteristic values */
 void GetCdCl(TFEFunction2D *u1fct, TFEFunction2D *u2fct,TFEFunction2D *pfct, TFEFunction2D *u1old,TFEFunction2D *u2old,double &cd, double &cl)

//void GetCdCl(TFEFunction2D *u1fct, TFEFunction2D *u2fct,TFEFunction2D *pfct, double &cd, double &cl)

{  
  cout << "test21 " <<endl;
  
  int i,j,k,l, N_;
  int N_Points,N_Edges,comp,ed_nr;
  double *weights, *xi, *eta;
  double X[MaxN_QuadPoints_2D];
  double Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  int N_LocalUsedElements;
  FE2D LocalUsedElements[2], CurrentElement;
  int *DOF;
  double **OrigFEValues, *Orig;
  bool SecondDer[2] = { FALSE, FALSE };
  double *u1, *u2, *p;
  TFESpace2D *USpace, *PSpace;
  int *UGlobalNumbers, *UBeginIndex;
  int *PGlobalNumbers, *PBeginIndex;
  int *N_BaseFunct, N_Cells;
  BaseFunct2D BaseFunct, *BaseFuncts;
  TCollection *Coll;
  TBaseCell *cell;
  double value, value1, value2, value3;
  double FEFunctValues[MaxN_BaseFunctions2D];
  double FEFunctValues1[MaxN_BaseFunctions2D];
  double FEFunctValues2[MaxN_BaseFunctions2D];
  double FEFunctValues3[MaxN_BaseFunctions2D];
  int N_DerivativesU = 3;
  double *Derivatives[MaxN_BaseFunctions2D];
  MultiIndex2D NeededDerivatives[3] = { D00, D10, D01 };
  TFEFunction2D *vfct;
  double *v, nu = 1./TDatabase::ParamDB->RE_NR;
  double *Der, *aux;
  TJoint *joint;
  TBoundEdge *boundedge;
  TBoundComp *BoundComp;
  TFE2D *eleCell;
  FE2D FEEle;
  TFEDesc2D *FEDesc;
  int N_DOF_Circ, *DOF_Circ;
  char VString[] = "v";

  u1 = u1fct->GetValues();
  u2 = u2fct->GetValues();
  p = pfct->GetValues();

  USpace = u1fct->GetFESpace2D();
  PSpace = pfct->GetFESpace2D();

  UGlobalNumbers = USpace->GetGlobalNumbers();
  UBeginIndex = USpace->GetBeginIndex();

  PGlobalNumbers = PSpace->GetGlobalNumbers();
  PBeginIndex = PSpace->GetBeginIndex();

  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  aux = new double [MaxN_QuadPoints_2D*10];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    Derivatives[j] = aux + j*10;

  N_ = u1fct->GetLength();
  v = new double[N_];
  memset(v,0,N_*SizeOfDouble);
  vfct = new TFEFunction2D(USpace, VString, VString, v, N_);

// ########################################################################
// loop over all cells
// ########################################################################
  Coll = USpace->GetCollection(); // all spaces use same Coll
  N_Cells = Coll->GetN_Cells();
 
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_Edges=cell->GetN_Edges();
    for(j=0;j<N_Edges;j++)              // loop over all edges of cell
    {
      joint=cell->GetJoint(j);
      if ((joint->GetType() == BoundaryEdge)||
          (joint->GetType() == IsoBoundEdge)) // boundary edge 
      {
        
        boundedge=(TBoundEdge *)joint;  
        BoundComp=boundedge->GetBoundComp();  // get boundary component
        comp=BoundComp->GetID();              // boundary id 
        if (comp==4) 
          {
            FEEle = USpace->GetFE2D(i,cell);   // finite element of cell
            eleCell =  TFEDatabase2D::GetFE2D(FEEle); 
            FEDesc = eleCell->GetFEDesc2D();   // fe descriptor
            N_DOF_Circ = FEDesc->GetN_JointDOF(); // # local dofs on joints
            DOF_Circ = FEDesc->GetJointDOF(j); // local dofs on joint j
            DOF = UGlobalNumbers + UBeginIndex[i]; // pointer to global dofs
            for (k=0;k<N_DOF_Circ;k++)         // set fe on circle to 1 
              v[DOF[DOF_Circ[k]]] = 1;
          }
      }      
    }
  }
  
  cd = 0;
  cl = 0;

// ########################################################################
// loop over all cells
// ########################################################################
  Coll = USpace->GetCollection(); // all spaces use same Coll
  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);

    // ####################################################################
    // find local used elements on this cell
    // ####################################################################
    N_LocalUsedElements = 2;
    LocalUsedElements[0] = USpace->GetFE2D(i, cell);
    LocalUsedElements[1] = PSpace->GetFE2D(i, cell);

    // ####################################################################
    // calculate values on original element
    // ####################################################################
    TFEDatabase2D::GetOrig(N_LocalUsedElements, LocalUsedElements, 
                         Coll, cell, SecondDer,
                         N_Points, xi, eta, weights, X, Y, AbsDetjk);

    // calculate all needed values of p 
    CurrentElement = LocalUsedElements[1];
    BaseFunct = BaseFuncts[CurrentElement];
    N_ = N_BaseFunct[CurrentElement];

    DOF = PGlobalNumbers + PBeginIndex[i];
    for(l=0;l<N_;l++)
      FEFunctValues[l] = p[DOF[l]];

    OrigFEValues = TFEDatabase2D::GetOrigElementValues(BaseFunct, D00);

    for(j=0;j<N_Points;j++)
    {
      Orig = OrigFEValues[j];
      value = 0;
      for(l=0;l<N_;l++)
        value += FEFunctValues[l] * Orig[l];

      Derivatives[j][0] = value;
    }

    // calculate all needed values of u1, u2 
    CurrentElement = LocalUsedElements[0];
    BaseFunct = BaseFuncts[CurrentElement];
    N_ = N_BaseFunct[CurrentElement];

    DOF = UGlobalNumbers + UBeginIndex[i];
    for(l=0;l<N_;l++)
    {
      FEFunctValues1[l] = u1[DOF[l]];
      FEFunctValues2[l] = u2[DOF[l]];
      FEFunctValues3[l] = v[DOF[l]];
    }

    for(k=0;k<N_DerivativesU;k++)
    {
      OrigFEValues = TFEDatabase2D::GetOrigElementValues(BaseFunct,
                                      NeededDerivatives[k]);
      for(j=0;j<N_Points;j++)
      {
        Orig = OrigFEValues[j];
        value1 = 0;
        value2 = 0;
        value3 = 0;
        for(l=0;l<N_;l++)
        {
          value1 += FEFunctValues1[l] * Orig[l];
          value2 += FEFunctValues2[l] * Orig[l];
          value3 += FEFunctValues3[l] * Orig[l];
        } // endfor l
        Derivatives[j][k+1] = value1;
        Derivatives[j][k+4] = value2;
        Derivatives[j][k+7] = value3;
      } // endfor j
    } // endfor k

    // calculation
    for(j=0;j<N_Points;j++)
    {
      Der = Derivatives[j];
      
      // nu * (u1_x*v_x, u1_y*v_y), v= (v,0)
      value1  = nu*(Der[2]*Der[8]+Der[3]*Der[9]);
      // (u1 * u1_x + u2* u1_y) * (1,0)
      value1 += (Der[1]*Der[2]+Der[4]*Der[3])*Der[7];
      // pressure times divergence of test function (1,0)
      value1 -= Der[0]*Der[8];

      value2  = nu*(Der[5]*Der[8]+Der[6]*Der[9]);
      value2 += (Der[1]*Der[5]+Der[4]*Der[6])*Der[7];
      value2 -= Der[0]*Der[9];

      cd += AbsDetjk[j]*weights[j] * value1;
      cl += AbsDetjk[j]*weights[j] * value2;
    }

  } // endfor i

  cd *= -500;
  cl *= -500;

  delete Derivatives[0];
  delete vfct;
  delete v;
}
