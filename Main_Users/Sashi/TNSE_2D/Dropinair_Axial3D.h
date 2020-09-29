


void ExampleFile()
{
  OutPut("Dropinair_axial3D.h" << endl) ;

  #define __AXIAL3D__

}


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
//   const double theta=TDatabase::ParamDB->IMPACT_ANGLE;

//   values[0] = cos((Pi/180)*theta);
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void InitialU2(double x, double y, double *values)
{
//   const double theta=TDatabase::ParamDB->IMPACT_ANGLE;

//   if(y==0) values[0] = 0;
//   else values[0] = -(1-cos((Pi/180)*theta));
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
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
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}


// ========================================================================
// Grid velocity
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


// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int i, double t, BoundCond &cond)
{

//   TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
//    cond =  FREESURF; // DIRICHLET; //
//  
  switch(i)
      {
        case 0:
      cond =  SLIP_FRICTION_PENETRATION_RESISTANCE;
       TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION = 1;
       TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE=0;
//       cond =  DIRICHLET;
        break;

        case 1:
             cond = FREESURF;
             TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
        break;
  default:
            Error("Unknown Boundary component ref example file" << endl);
            exit(0);
       }
       
 }

void U1BoundValue(int BdComp, double Param, double &value)
{
  value = 0;
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  value = 0;
}

// ========================================================================
// Modify coordinate, starting from reference geometry
// ========================================================================
void ModifyCoords(double &x, double &y, double &theta)
{

  const double a=TDatabase::ParamDB->P1;
  const double b=TDatabase::ParamDB->P2;
  const double k=TDatabase::ParamDB->P4;
  const  int temp1 = int(TDatabase::ParamDB->P3 + 0.1);
  double midr, dr;
  double phi, r;
  int temp2;

  switch(temp1)
  {
  case 0:
     dr=(a*a*b*b)/(b*b*cos(theta)*cos(theta)+a*a*sin(theta)*sin(theta));

     x = sqrt(dr)*cos(theta);
     if (fabs(x)< 1e-9) x = 0;
     y = sqrt(dr)*sin(theta);
     if (fabs(y)< 1e-9) y = 0;
      break;
     
 case 1:
     midr=0.5*(a+b);
     dr=0.5*(a-b);

     phi = atan2(y, x);
     r = midr + dr * cos(k*phi);

     x = r*cos(phi);
     if (fabs(x)< 1e-6) x = 0;
     y = r*sin(phi);
     if (fabs(y)< 1e-6) y = 0;
  //    cout << "x = " << x << "  y = " <<  y << endl;
      
   break;

  default:
        OutPut("Specify P3 value " << endl);
        exit(1);
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
// inclined angle of the plane on which the droplet deforms
//   static double theta = TDatabase::ParamDB->P1;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
//     if(TDatabase::ParamDB->FR_NR == 0)
//        coeff[1] = sin((Pi/180)*theta);
//     else
//     coeff[1] = sin((Pi/180)*theta)/TDatabase::ParamDB->FR_NR; // f1
//     if(TDatabase::ParamDB->FR_NR == 0)
//        coeff[2] = 0.;
//     else
//      coeff[2] = -cos((Pi/180)*theta)/TDatabase::ParamDB->FR_NR; // f2 - => gravity in opposite direction

     coeff[1]= 0.;
     coeff[2]= 0.;

  }
}


void InitialS(double x, double y, double *values)
{
 double t = TDatabase::TimeDB->CURRENTTIME;
 double Re = TDatabase::ParamDB->RE_NR;
 double Pr = TDatabase::ParamDB->PR_NR;

  static double eps = 1.0/(Re*Pr);

  double phi = atan2(y, x);
  phi = Pi/2. - phi;

//   values[0] = 1. -  cos(phi)*exp(-2.*t*eps); // diffusion test
//   values[0] = 0.5 + 0.5*cos(phi); // expanding static sphere
  values[0] = 0.9*fabs(cos(phi)); // Oscllating sphere
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void ExactS(double x, double y, double *values)
{
 double L, t = TDatabase::TimeDB->CURRENTTIME;
 double Re = TDatabase::ParamDB->RE_NR;
 double Pr = TDatabase::ParamDB->PR_NR;

  static double eps = 1.0/(Re*Pr);

  double phi = atan2(y, x);
//   phi = Pi/2. - phi;

//   L = ( 1./(16.*Pi*t + 1.) ); // expanding static sphere
//   values[0] = 0.5 + 0.5*cos(phi)*exp(-2.*t*eps); // diffusion test
//   values[0] = (0.5 + 0.5*cos(phi)) * (pow(L, 1./3.)); // expanding static sphere
  values[0] = cos(phi); // Oscllating sphere
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
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
  double r2;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    // if( (r2 = (x[i]-0.2)*(x[i]-0.2) + (y[i]-0.2)*(y[i]-0.2)) < 0.01)
    //   coeff[0] = 10*sqrt(r2);
    // else
    //   coeff[0] = 10*0.1;
    coeff[0] = 1;

    coeff[1] = 0;
    coeff[2] = 0;
    coeff[3] = 0;

    coeff[4] = 0;
  }
}
void HeatCoeffs(int n_points, double *x, double *y,
        double **parameters, double **coeffs)
{
  int i;
  double *coeff;
  double r2;
  double Re = TDatabase::ParamDB->RE_NR;
  double Pr = TDatabase::ParamDB->PR_NR;

  static double eps = 1.0/(Re*Pr);

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = 0;
    coeff[2] = 0;
    coeff[3] = 0;
    coeff[4] = 0;
  }
}
 
// **************************************************************************
// Triangular Mesh Generation
// **************************************************************************

void  TriaReMeshGen(TDomain *&Domain, int &N_LineBDComp, double &FreeBD_X,  double &FreeBD_Y)
{
  TBoundComp *BoundComp;
  TBdLine *UpdateSlipBound;
  TBdCircle *UpdateFreeBound;
  TBaseCell **CellTree, *cell;
  TBoundPart *BoundPart;
  TJoint *Joint;
  TCollection *coll;
  TVertex **VertexDel, **NewVertices;

  double dt, area = TDatabase::ParamDB->Area;
  double mode_diff = TDatabase::ParamDB->P4;
  double r, phi, t, T_a, T_b, x, y, theta;
  double *Coordinates, temp;
  double left, right, bottom, top;
  double y_begin, y_end, dx, deviation;
  
  int i, j, k, ID, In_Index, N_Cells, N_G, *PointNeighb, maxEpV=0;
  int N_FreeBound_Vert, a, b, len1, len2, Neighb_tmp, mode = int (TDatabase::ParamDB->P4);
  int N_Interf_Vertices, CurrComp, CurrVertex, N_Joints, N_Vertices;
  int N_RootCells, *PartMarker, *Triangles, Neib[2], CurrNeib;
  int N_SlipBound_Vert;

  boolean AllowEdgeRef = (boolean) TDatabase::ParamDB->MESHGEN_ALLOW_EDGE_REF;
 // cout<< " MESHGEN_ALLOW_EDGE_REF" <<AllowEdgeRef<<endl;

  struct triangulateio In, Out;
  std::ostringstream opts;
  opts << " ";

  deviation = fabs(TDatabase::ParamDB->P1 -  TDatabase::ParamDB->P2);
  BoundPart = Domain->GetBdPart(0);
  UpdateFreeBound = (TBdCircle*)BoundPart->GetBdComp(0);
  UpdateSlipBound = (TBdLine*)BoundPart->GetBdComp(1);


  N_LineBDComp = 1;

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

  N_FreeBound_Vert = int (TDatabase::ParamDB->P6);    //Freesurf except end point
  N_SlipBound_Vert = 50;     // Initially only three points on solid bound (except end point)
  N_Interf_Vertices = N_FreeBound_Vert+N_SlipBound_Vert;
  In.numberofpoints = N_Interf_Vertices;
  In.pointlist = new double[2*In.numberofpoints];
  In.pointmarkerlist = new int[In.numberofpoints];
  In.numberofpointattributes = 0;

  In.numberofsegments = In.numberofpoints;
  In.segmentlist = new int[2*In.numberofsegments];
  In.segmentmarkerlist = new int[In.numberofsegments];
  In.numberofholes = 0;
  In.holelist = NULL;
  In.numberofregions = 0;
  In.regionlist = NULL;

  In_Index = 0;
  CurrComp = 1;

    T_a = 1.; // x axis value in ellipse
    T_b = 1.; // y axis value in ellipse !!! set in modifyGausscoord funct also
    x = 0.0;     // Axial bound

    theta = -Pi/2.; // end degree value of freesurface
    dt = Pi/N_FreeBound_Vert;

 // points and segments on the interface (marker=2)
    t = theta;

   for(i=0;i<N_FreeBound_Vert;i++) // without last point
    {
    //  cout<<" theta : "<< theta <<endl;
      x = T_a*cos(theta);
      y = T_b*sin(theta);

// spherical harmonic of order 2
// spherical harmonic of order 4
    phi = atan2(y, x);
    if(mode==2)
      r = 1.0 + deviation*(sqrt(5./Pi)*(3.*cos(phi+Pi/2.)*cos(phi+Pi/2.) -1. )/4.);
    else if(mode==4) 
     {
      temp = cos(phi+Pi/2.);
      r = 1.0 + deviation*(sqrt(1./Pi)*(35.*temp*temp*temp*temp - 30.*temp*temp +3. )*3./16.);
     }
    else
      r=1.;

      x = r*cos(phi);
      y = r*sin(phi);
//       ModifyCoords( x,  y);
      if (fabs(y)<1e-10) y = 0.;
      if (fabs(x)<1e-10) x = 0.;

      if(i==0) 
       {
        FreeBD_X = x;   FreeBD_Y = y; 
       }


      In.pointlist[2*In_Index] = x;
      In.pointlist[2*In_Index+1] = y;
//       cout<<" x : "<< x << " y : "<< y<<endl;
      //In.pointmarkerlist[In_Index] = CurrComp;
      In.segmentlist[2*In_Index] = In_Index;
      In.segmentlist[2*In_Index+1] = In_Index+1;
      if (AllowEdgeRef)
       {
    //    if (i)
         In.segmentmarkerlist[In_Index] = CurrComp;
       }
      else
       {
    //    if (i)
          In.segmentmarkerlist[In_Index] = 100000 + CurrComp;
       }

      In_Index++;
      theta = t + double(i+1)*dt;
//       if(theta>360) theta -= 360.;
    }
   CurrComp++;


    T_a = 1.; // x axis value in ellipse
    T_b = 1.; // y axis value in ellipse !!! set in modifyGausscoord funct also
    x = 0.0;     // Axial bound

    y = T_b;
// spherical harmonic of order 2 
// spherical harmonic of order 4
    phi = atan2(y, x);
    if(mode==2)
      r = 1.0 + deviation*(sqrt(5./Pi)*(3.*cos(phi+Pi/2.)*cos(phi+Pi/2.) -1. )/4.);
    else if(mode==4) 
     {
      temp = cos(phi+Pi/2.);
      r = 1.0 + deviation*(sqrt(1./Pi)*(35.*temp*temp*temp*temp - 30.*temp*temp +3. )*3./16.);
     }
    else 
      OutPut("No. of mode is taken as 0 check main programme"<<endl);

    y = r*sin(phi);

    y_begin = y;



    y = -T_b;
    phi = atan2(y, x);
// spherical harmonic of order 2
// spherical harmonic of order 4
    phi = atan2(y, x);
    if(mode==2)
      r = 1.0 + deviation*(sqrt(5./Pi)*(3.*cos(phi+Pi/2.)*cos(phi+Pi/2.) -1. )/4.);
    else if(mode==4) 
     {
      temp = cos(phi+Pi/2.);
      r = 1.0 + deviation*(sqrt(1./Pi)*(35.*temp*temp*temp*temp - 30.*temp*temp +3. )*3./16.);
     }

    y = r*sin(phi);
    y_end = y;

    dx = (y_end - y_begin)/N_SlipBound_Vert;
    y = y_begin;
  // points and segments on the horizontal boundary (marker=1)
  for(i=0;i<N_SlipBound_Vert;i++) // without last point
   {
    if (fabs(y)<1e-10) y = 0.;
    if (fabs(x)<1e-10) x = 0.;
    In.pointlist[2*In_Index] = 0.0;
    In.pointlist[2*In_Index+1] = y;
//     cout<<" x : "<< x << " y : "<< y<<endl;
    In.pointmarkerlist[In_Index] = CurrComp;
    In.segmentlist[2*In_Index] = In_Index;
    In.segmentlist[2*In_Index+1] = In_Index+1;
    In.segmentmarkerlist[In_Index] = CurrComp;
    In_Index++;
    y = y_begin + double(i+1)*dx;
   }
   CurrComp++;




  In.segmentlist[2*(In.numberofsegments-1)+1] = 0;

// exit(0);
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

/*
  for(i=0;i<In.numberofpoints;i++)
    OutPut(i<<' '<<In.pointmarkerlist[i]<<' '<<
	   In.pointlist[2*i]<<' '<<In.pointlist[2*i+1]<<endl);
cout<<endl;
//exit(0);
*/
  triangulate((char*)opts.str().c_str(), &In, &Out, (struct triangulateio *)NULL);

/*
  for(i=0;i<Out.numberofpoints;i++)
     OutPut(i<<' '<<Out.pointmarkerlist[i]<<' '<<
      Out.pointlist[2*i]<<' '<<Out.pointlist[2*i+1]<<endl);
  */


  Domain->GetTreeInfo(CellTree,N_RootCells);
  coll = Domain->GetCollection(It_Finest, 0);
  N_Cells = coll->GetN_Cells();

  // remove all existing vertices and joints
  VertexDel = new TVertex*[3*N_RootCells];
  // DelCell = new TGridCell*[N_Cells];

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

      // set bounding box
  left = bottom = 1e8;
  right = top = -1e8;

   for(i=0;i<In.numberofpoints;i++)
    {
      if(left>In.pointlist[2*i]) left = In.pointlist[2*i];
      if(right<In.pointlist[2*i]) right = In.pointlist[2*i];
      if(top<In.pointlist[2*i+1]) top = In.pointlist[2*i+1];
      if(bottom>In.pointlist[2*i+1]) bottom = In.pointlist[2*i+1];
    }

  OutPut("left: "<<left<<" right: "<<right<<" top: "<<top<<" bottom: "<<bottom<<endl);

  Domain->SetBoundBox(right-left,top-bottom);
  Domain->SetBoundBoxstart(left,bottom);

 // Solid Bound startx, starty, x length and y length
 UpdateSlipBound->SetParams(0., y_begin, 0, y_end - y_begin);

// Free boundary xmid, ymid, radius_a, radius_b, start angle, end angle
//  UpdateFreeBound ->SetParams( 0.0, 0.0, T_a, T_b, -Pi/2., Pi/2.);
 
//  cout << " In.pointlist[N_SlipBound_Vert] : " << In.pointlist[2*N_SlipBound_Vert]<< endl;
 // cout << " In.pointlist[N_SlipBound_Vert] : " << In.pointlist[0]<< endl;
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


      if(Domain->GetBdPart(0)->GetBdComp(CurrComp)->GetTofXY(
            NewVertices[a]->GetX(), NewVertices[a]->GetY(), T_a) ||
          Domain->GetBdPart(0)->GetBdComp(CurrComp)->GetTofXY(
            NewVertices[b]->GetX(), NewVertices[b]->GetY(), T_b))
       {
          cerr<<"Error: could not set parameter values"<<endl;
          OutPut(NewVertices[a]<<endl);
          OutPut(NewVertices[b]<<endl);
        //  exit(0);
       }

      if (CurrNeib == 2)    // 2 cells contain the current edge
        if(Domain->GetBdPart(0)->GetBdComp(CurrComp)->IsFreeBoundary())
          Joint = new TIsoInterfaceJoint(Domain->GetBdPart(0)->GetBdComp(CurrComp),
                  T_a, T_b, CellTree[Neib[0]], CellTree[Neib[1]]);
        else
          Joint = new TInterfaceJoint(Domain->GetBdPart(0)->GetBdComp(CurrComp),
                  T_a, T_b, CellTree[Neib[0]], CellTree[Neib[1]]);
      else
        if(Domain->GetBdPart(0)->GetBdComp(CurrComp)->IsFreeBoundary())
          Joint = new TIsoBoundEdge(Domain->GetBdPart(0)->GetBdComp(CurrComp), T_a, T_b);
        else
          Joint = new TBoundEdge(Domain->GetBdPart(0)->GetBdComp(CurrComp), T_a, T_b);
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
// Triangular for grid generation end
//======================================================================

}