

void ExampleFile()
{  
  #define __IMPINGINGDROPLET__
  #define __AXIAL3D__
//   #define __EXPERIMENTAL__
 #define __SURFACT__
//  #define __SURFACTDIFFTEST__
 
 #define __SURFACTCONVTEST__
 
  
  
 TDatabase::ParamDB->REACTOR_P9=1;
 
  TDatabase::ParamDB->Axial3D=1;
  OutPut("Drop_Imping_Axial3D.h" << endl) ;
  OutPut("TDatabase::ParamDB->Axial3D = " << TDatabase::ParamDB->Axial3D <<endl) ;
}


// ========================================================================
// initial solution
// ========================================================================

void InitialU1(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void InitialU2(double x, double y, double *values)
{
 
  if(y==0)
  { values[0] = 0.; }
  else 
  { values[0] = TDatabase::ParamDB->P15;}
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void InitialP(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
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

void InitialS(double x, double y, double *values)
{
//  double t=0;
//  double Re = TDatabase::ParamDB->RE_NR;
//  double Pr = TDatabase::ParamDB->Pr_NR;

//   static double eps = 1./TDatabase::ParamDB->REACTOR_P17;

  double phi = atan2(y, x);
  phi += Pi/2.;
  
//   if(fabs(sqrt(x*x + y*y))
//    cout << t << " Phi " << (180./Pi)*phi << " X " << x << " Y " << y << endl;
//   values[0] = 0.5 +  0.5*cos(phi)*exp(-2.*t*eps); // diffusion test
  
//   if(phi<1e-8 && x<1.)
//     values[0] = 0.;
//     
  values[0] = 0.5 + 0.5*cos(2.*phi - Pi); // expanding static sphere
//   values[0] = 0.9*fabs(cos(phi)); // Oscllating sphere
//   values[0] = sqrt((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5));

  values[0] = TDatabase::ParamDB->REACTOR_P12;  // static bubble
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void ExactS(double x, double y, double *values)
{
//  double L, t = TDatabase::TimeDB->CURRENTTIME;
// //  double Re = TDatabase::ParamDB->RE_NR;
// //  double Pr = TDatabase::ParamDB->PR_NR;
// 
//   static double eps = 1./TDatabase::ParamDB->REACTOR_P17;
// 
//   double phi = atan2(y-1, x);
//   phi += Pi/2.;
// //    cout << "Phi " << phi << " X " << x << " Y " << y << endl;
// //   L = ( 1./(16.*Pi*t + 1.) ); // expanding static sphere
//   values[0] = 0.5 + 0.5*cos(phi)*exp(-2.*t*eps); // diffusion test
//   
// //   if(phi<1e-8 && x<1.)
// //     values[0] = 0.;
//   
// //   values[0] = (0.5 + 0.5*cos(phi)) * (pow(L, 1./3.)); // expanding static sphere
// //   values[0] = 0.9*fabs(cos(phi)); // Oscllating sphere
// //   values[0] = 0.;
// 
// //   values[0] = TDatabase::ParamDB->REACTOR_P12;  // static bubble  
//   values[1] = 0;
//   values[2] = 0;
//   values[3] = 0;
  cout << " No exact case in Convetion test " <<endl;
  exit(0);
    
}

void InitialSuract(double x, double y, double *values)
{
//   values[0] = TDatabase::ParamDB->REACTOR_P19;
  values[0] = 0.;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void ExactSuract(double x, double y, double *values)
{
 double Pe_c = TDatabase::ParamDB->REACTOR_P18;
 static double eps = 1.0/(Pe_c);

//   double phi = atan2(y, x);
//   phi = Pi/2. - phi;
//   L = ( 1./(16.*Pi*t + 1.) ); // expanding static sphere
//   values[0] = 0.5 + 0.5*cos(phi)*exp(-2.*t*eps); // diffusion test
//   values[0] = (0.5 + 0.5*cos(phi)) * (pow(L, 1./3.)); // expanding static sphere
//   values[0] = cos(phi); // Oscllating sphere
  values[0] = 0.;
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
        case 0:
        case 2:
             cond = SLIP_FRICTION_PENETRATION_RESISTANCE;
             TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
             TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION = 1;
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


void SurfactBoundCondition(int i, double t, BoundCond &cond)
{
 cond = NEUMANN;
}

void SurfactBoundValue(int BdComp, double Param, double &value)
{
 value = 0.0;
}

void SurfactCoeffs(int n_points, double *x, double *y,
        double **parameters, double **coeffs)
{
  int i;
  double *coeff;
  static double eps = 1./ TDatabase::ParamDB->REACTOR_P18; // 1/Pe_c

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
    coeff[1] = 0.; // f1
    
    if(TDatabase::ParamDB->FR_NR == 0)
    {  coeff[2] = 0.; }
    else
    { coeff[2] = -1./TDatabase::ParamDB->FR_NR;} // f2 - => gravity in opposite direction
  }
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
 
void Solver_3dia(int N_Splines, double *a, double *b, double *c, double *rhs, double *sol)
{
  double *alpha, *beta, *y;
  int i, N;

  N = N_Splines+1;
  alpha = new double[N]; beta = new double[N]; y = new double[N];

  alpha[0] = a[0]; y[0] = rhs[0];
  for(i=1;i<N;i++)
  {
    beta[i] = b[i]/alpha[i-1];
    alpha[i] = a[i]-beta[i]*c[i-1];
    y[i] = rhs[i]-beta[i]*y[i-1];
  }

  sol[N-1] = y[N-1]/alpha[N-1];
  for(i=N-2;i>=0;i--)
    sol[i] = (y[i]-c[i]*sol[i+1])/alpha[i];

  delete [] alpha; delete [] beta; delete [] y;
}




  
void ReParam_axial3D_U(int N_E, TBaseCell **cell, int *EdgeNo,  int *CellNo, 
                       TFEVectFunct2D *Velocity,  TFEFunction2D *Surfactant, 
                       TFEFunction2D *SurfSurfactant, bool UpdateU)
{
  int i, j, ISpline, N_Splines, N_V, ORDER, VSP, m, m1, k, i3, USpline, FeDof;
  int *VeloBeginIndex, *VeloGlobalNumbers, *JointDOF, *DOF, N_DOF_Joint, *U_DOF;
  int *SurfactGlobalNumbers, *SurfactBeginIndex, *SJointDOF, *SDOF, SN_DOF_Joint, *Surf_DOF;
  int *SurfSurfactGlobalNumbers, *SurfSurfactBeginIndex, *SurfSJointDOF, *SurfSDOF;
  int *SurfSurf_DOF, SurfSN_DOF_Joint;
  
  double *h, *t, u0, u1, u2;
  double *a, *b, *c, *x, *y, teta;
  double *rhs, *u1rhs, *u2rhs, *Mx, *My,*Mu1, *Mu2, *Params, *Param9, *FEParams;
  double phi1, phi2, phi3, phi4, X, Y, T;
  double dx0, dy0, dx1, dy1, *u1_spl, *u2_spl;
  double *ValuesUX, *ValuesUY, tx, ty, *Surfact, *surf_spl, *srhs, *Msurf, surf;
  double *SurfSurfact, *Surfsurf_spl, *Surfsrhs, *SurfMsurf, Surfsurf;
   
  TIsoBoundEdge *isojoint;
  TVertex **IsoVertices;
  TJoint *Joint;
  TBaseCell *Me;
  TFESpace2D *VelocitySpace, *SurfactSpace, *SurfSurfactSpace;
  FE2D FEId, SFEId, SurfSFEId;
  TFE2D *ele, *Sele, *SurfSele;
  TFEDesc2D *FeDesc, *SFeDesc, *SurfSFeDesc;
  TCollection *coll;

  ORDER = 0;
  VSP = TDatabase::ParamDB->VELOCITY_SPACE;

  if (abs(VSP) > 20)
   {ORDER = abs(VSP) - 20;}
  else if ( abs(VSP) > 10)
   {ORDER = abs(VSP) - 10;}
  else ORDER = abs(VSP);

  N_V = N_E+1 + N_E*(ORDER-1);

  N_Splines = N_V-1;
  h = new double[N_Splines+1];
  t = new double[N_Splines+1];
  a = new double[N_Splines+1];
  b = new double[N_Splines+1];
  c = new double[N_Splines+1];
  rhs = new double[N_Splines+1];
  Mx = new double[N_Splines+1];
  My = new double[N_Splines+1];
  Params = new double [10*N_Splines];
  Param9 = new double [N_Splines+1];

  x = new double[N_V];
  y = new double[N_V];
  
  
  VelocitySpace = Velocity->GetFESpace2D();  
  coll = VelocitySpace->GetCollection();  
  
  if(UpdateU)
  {
   u1rhs = new double[N_Splines+1];
   u2rhs = new double[N_Splines+1];
   u1_spl = new double[N_Splines+1];
   u2_spl = new double[N_Splines+1];
   Mu1 = new double[N_Splines+1];
   Mu2 = new double[N_Splines+1];  
   FEParams = new double [4*2*N_Splines]; // u1, u2, surfact, surfsurfact
   U_DOF = new int[N_V];

   VeloBeginIndex = VelocitySpace->GetBeginIndex();
   VeloGlobalNumbers = VelocitySpace->GetGlobalNumbers();
   ValuesUX = Velocity->GetValues();
   ValuesUY = ValuesUX + Velocity->GetLength();
   
#ifdef __SURFACT__      
   Surf_DOF = new int[N_V];
   surf_spl = new double[N_V];  
   srhs = new double[N_Splines+1];   
   Msurf = new double[N_Splines+1];
 
   SurfSurf_DOF = new int[N_V];
   Surfsurf_spl = new double[N_V];  
   Surfsrhs = new double[N_Splines+1];   
   SurfMsurf = new double[N_Splines+1];  
   
   
   SurfactSpace = Surfactant->GetFESpace2D();
   SurfactBeginIndex = SurfactSpace->GetBeginIndex();
   SurfactGlobalNumbers = SurfactSpace->GetGlobalNumbers();
   Surfact = Surfactant->GetValues(); 
   
   SurfSurfactSpace = SurfSurfactant->GetFESpace2D();
   SurfSurfactBeginIndex = SurfSurfactSpace->GetBeginIndex();
   SurfSurfactGlobalNumbers = SurfSurfactSpace->GetGlobalNumbers();
   SurfSurfact = SurfSurfactant->GetValues();    
#endif      
   
  }
  
   m = 0;
   m1 = 0;
   for(i=0;i<N_E;i++) // i<N_E
   {
    Me = cell[i];
    Me->GetVertex(EdgeNo[i])->GetCoords(x[m], y[m]);
    m++;

    Joint = cell[i]->GetJoint(EdgeNo[i]);
    isojoint = (TIsoBoundEdge *)Joint;
    k = isojoint->GetN_Vertices();
    if(k==ORDER-1)
     {
      IsoVertices = isojoint->GetVertices();
      for(i3=0;i3<k;i3++)
       {   
        IsoVertices[i3]->GetCoords(x[m], y[m]);
//          cout<< i<<" FreeGaus " << (180/Pi)*atan2(y[m], x[m]) <<endl;
        m++;
       } 
     }
    else
     {
      // only second order conforming elements implimented
      cout<< " No match in isopoints per free edge "<<endl;
      exit(0);
     }

// for velocity
   if(UpdateU)
   {
    FEId = VelocitySpace->GetFE2D(CellNo[i], Me);
    ele = TFEDatabase2D::GetFE2D(FEId);
    FeDesc = ele->GetFEDesc2D();   // fe descriptor
    JointDOF = FeDesc->GetJointDOF(EdgeNo[i]);
    N_DOF_Joint = FeDesc->GetN_JointDOF();
    DOF = VeloGlobalNumbers + VeloBeginIndex[CellNo[i]];

#ifdef __SURFACT__
    // for surfactant
    SFEId = SurfactSpace->GetFE2D(CellNo[i], Me);
    Sele = TFEDatabase2D::GetFE2D(SFEId);
    SFeDesc = Sele->GetFEDesc2D();   // fe descriptor
    SJointDOF = SFeDesc->GetJointDOF(EdgeNo[i]);
    SN_DOF_Joint = SFeDesc->GetN_JointDOF();
    SDOF = SurfactGlobalNumbers + SurfactBeginIndex[CellNo[i]];    
    
    // for Surfsurfactant
    SurfSFEId = SurfSurfactSpace->GetFE2D(CellNo[i], Me);
    SurfSele = TFEDatabase2D::GetFE2D(SurfSFEId);
    SurfSFeDesc = SurfSele->GetFEDesc2D();   // fe descriptor
    SurfSJointDOF = SurfSFeDesc->GetJointDOF(EdgeNo[i]);
    SurfSN_DOF_Joint = SurfSFeDesc->GetN_JointDOF();
    SurfSDOF = SurfSurfactGlobalNumbers + SurfSurfactBeginIndex[CellNo[i]];       
#endif
    
    if((N_DOF_Joint-1)!=ORDER)
     {
      // only second order conforming elements implimented
      cout<< " No match in no. velo-nodal functunals and vertices on the free surface edge  "<<endl;
      cout<< " (N_DOF_Joint-1) " << N_DOF_Joint-1 << " ORDER " << ORDER <<endl;
      exit(0);
     }

    if(i !=N_E-1)// -1 due to end dof will be the start dof of the next edge except on last edge
     N_DOF_Joint--; // assumed that velocity and surfactant having same no. of dof on edge

     // //   cout << " CellNo[i] " << CellNo[i] << endl;
     for (i3=0;i3<N_DOF_Joint;i3++)
       {
         U_DOF[m1] = DOF[JointDOF[i3]]; // needed for later update
         u1_spl[m1] = ValuesUX[DOF[JointDOF[i3]]];
         u2_spl[m1] = ValuesUY[DOF[JointDOF[i3]]];

#ifdef __SURFACT__
         Surf_DOF[m1] = SDOF[SJointDOF[i3]]; // needed for later update
         surf_spl[m1] = Surfact[SDOF[SJointDOF[i3]]];
 
         SurfSurf_DOF[m1] = SurfSDOF[SurfSJointDOF[i3]]; // needed for later update
         Surfsurf_spl[m1] = SurfSurfact[SurfSDOF[SurfSJointDOF[i3]]];	 
#endif 

         m1++;
       }
     } //  if(UpdateU)       
              
    } // for(i=0;i<N_E
  
  
//   end vertex of the freeboundary
  k = cell[N_E-1]->GetN_Edges();
  cell[N_E-1]->GetVertex((EdgeNo[N_E-1]+1) % k)->GetCoords(x[m], y[m]);
//   cout << "x " << x[m] << " y " << y[m] << endl;
//   exit(0);
  
  
 
  if(m+1!=m1 && UpdateU)
    {
      // only second order conforming elements implimented
      cout<< " No match in no. velo-nodal functunals and vertices on the free surface edge  "<<endl;
      cout<< " m " << m << " m1 " << m1 <<endl;
      exit(0);
     }

//  for(i=0;i<N_V;i++)
//    OutPut("OldX: "<< i <<' '<<x[i] <<' '<< y[i] <<endl);
//  OutPut("OldX: "<< i <<' '<<x[0] <<' '<< y[0] <<endl);
//  OutPut("OldX: "<< i <<' '<<x[N_V-2] <<' '<< y[N_V-2] <<endl);
//  OutPut("OldX: "<< i <<' '<<x[N_V-1] <<' '<< y[N_V-1] <<endl);
// cout << "Surfact[Surf_DOF[0]] " <<Surfact[Surf_DOF[0]] << " Surfact[Surf_DOF[m1]] " << Surfact[Surf_DOF[m1-1]] <<endl;
  
  h[0] = 0.0; t[0] = 0.0;

 for(i=1;i<=N_Splines;i++)
  {
    h[i] = sqrt((x[i]-x[i-1])*(x[i]-x[i-1])+(y[i]-y[i-1])*(y[i]-y[i-1]));
    t[i] = t[i-1] + h[i];
  }

  dx0 = (x[1]-x[0])/h[1];
  dy0 = (y[1]-y[0])/h[1];

  dx1 = (x[N_Splines]-x[N_Splines-1])/h[N_Splines];
  dy1 = (y[N_Splines]-y[N_Splines-1])/h[N_Splines];


  a[0] = 2.; c[0] = 1.; rhs[0] = -6./h[1]*(dx0 - (x[1]-x[0])/h[1]);
  for(i=1;i<N_Splines;i++)
  {
    a[i] = 2.;  
    b[i] = h[i]/(h[i]+h[i+1]); // \mu_i in PhD thesis
    c[i] = h[i+1]/(h[i]+h[i+1]); // \lambda_i in PhD thesis
    rhs[i] = 6./(h[i]+h[i+1])*((x[i+1]-x[i])/h[i+1]-(x[i]-x[i-1])/h[i]);
  }
  b[N_Splines] = 1.; a[N_Splines] = 2.;
  rhs[N_Splines] = 6./h[N_Splines]*(dx1 - (x[N_Splines]-x[N_Splines-1])/h[N_Splines]);

  Solver_3dia(N_Splines, a, b, c, rhs, Mx);

  rhs[0] = -6./h[1]*(dy0 - (y[1]-y[0])/h[1]);
  for(i=1;i<N_Splines;i++)
  {
    rhs[i] = 6./(h[i]+h[i+1])*((y[i+1]-y[i])/h[i+1]-(y[i]-y[i-1])/h[i]);
  }
  rhs[N_Splines] = 6./h[N_Splines]*(dy1 - (y[N_Splines]-y[N_Splines-1])/h[N_Splines]);

  Solver_3dia(N_Splines, a, b, c, rhs, My);

  for(i=0;i<N_Splines;i++)
  {
    ISpline = i*10;
    Params[ISpline    ] = x[i]; 
    Params[ISpline + 1] = y[i];
    Params[ISpline + 2] = x[i+1]; 
    Params[ISpline + 3] = y[i+1];
    Params[ISpline + 4] = -Mx[i]*h[i+1]*h[i+1]/2. +
                          ((x[i+1]-x[i])/h[i+1]-h[i+1]/6.*(Mx[i+1]-Mx[i]))*h[i+1];

//     Params[ISpline + 4] = Mx[i]*h[i];
    Params[ISpline + 5] = -My[i]*h[i+1]*h[i+1]/2. +
                          ((y[i+1]-y[i])/h[i+1]-h[i+1]/6.*(My[i+1]-My[i]))*h[i+1];
    Params[ISpline + 6] = Mx[i+1]*h[i+1]*h[i+1]/2. +
                          ((x[i+1]-x[i])/h[i+1]-h[i+1]/6.*(Mx[i+1]-Mx[i]))*h[i+1];

//     Params[ISpline + 6] = -Mx[i+1];
    Params[ISpline + 7] = My[i+1]*h[i+1]*h[i+1]/2. +
                          ((y[i+1]-y[i])/h[i+1]-h[i+1]/6.*(My[i+1]-My[i]))*h[i+1];
    Params[ISpline + 8] = t[i+1]/t[N_Splines];
    Params[ISpline + 9] = 0.;

   //cout<<"  "<<Params[ISpline + 8]<<'\t'<<Params[ISpline + 9]<<endl;
  }
  
 if(UpdateU)
 {
  // ===============================================================
  // u1 component
  // ===============================================================
  for(i=1;i<N_Splines;i++)
   {
     u0 = u1_spl[i-1];
     u1 = u1_spl[i];
     u2 = u1_spl[i+1];

     u1rhs[i] = 6./(h[i]+h[i+1])*((u2-u1)/h[i+1]-(u1-u0)/h[i]);
    }

   u1rhs[0] = u1rhs[1];
   u1rhs[N_Splines] = u1rhs[N_Splines-1];


  
  Solver_3dia(N_Splines, a, b, c, u1rhs, Mu1);     
      
  // ===============================================================
  // u2 component
  // ===============================================================
  for(i=1;i<N_Splines;i++)
   {
     u0 = u2_spl[i-1];
     u1 = u2_spl[i];
     u2 = u2_spl[i+1];

     u2rhs[i] = 6./(h[i]+h[i+1])*((u2-u1)/h[i+1]-(u1-u0)/h[i]);
    }

   u2rhs[0] = u2rhs[1];
   u2rhs[N_Splines] = u2rhs[N_Splines-1];

  Solver_3dia(N_Splines, a, b, c, u2rhs, Mu2);
  
  
#ifdef __SURFACT__  
// ===============================================================
// surfactant
  for(i=1;i<N_Splines;i++)
   {
     u0 = surf_spl[i-1];
     u1 = surf_spl[i];
     u2 = surf_spl[i+1];

     srhs[i] = 6./(h[i]+h[i+1])*((u2-u1)/h[i+1]-(u1-u0)/h[i]);
    }

   srhs[0] = srhs[1];
   srhs[N_Splines] = srhs[N_Splines-1];

  Solver_3dia(N_Splines, a, b, c, srhs, Msurf);
// ===============================================================
// surfactant
  for(i=1;i<N_Splines;i++)
   {
     u0 = Surfsurf_spl[i-1];
     u1 = Surfsurf_spl[i];
     u2 = Surfsurf_spl[i+1];

     Surfsrhs[i] = 6./(h[i]+h[i+1])*((u2-u1)/h[i+1]-(u1-u0)/h[i]);
    }

   Surfsrhs[0] = Surfsrhs[1];
   Surfsrhs[N_Splines] = Surfsrhs[N_Splines-1];

  Solver_3dia(N_Splines, a, b, c, Surfsrhs, SurfMsurf);
  // ===============================================================
#endif    
  
// ===============================================================

  for(i=0;i<N_Splines;i++)
  {
    ISpline = i*8;

    FEParams[ISpline  ] = -Mu1[i]*h[i+1]*h[i+1]/2. +
                          ((u1_spl[i+1]-u1_spl[i])/h[i+1]-h[i+1]/6.*(Mu1[i+1]-Mu1[i]))*h[i+1];
    FEParams[ISpline + 1] = Mu1[i+1]*h[i+1]*h[i+1]/2. +
                          ((u1_spl[i+1]-u1_spl[i])/h[i+1]-h[i+1]/6.*(Mu1[i+1]-Mu1[i]))*h[i+1];


    FEParams[ISpline + 2  ] = -Mu2[i]*h[i+1]*h[i+1]/2. +
                          ((u2_spl[i+1]-u2_spl[i])/h[i+1]-h[i+1]/6.*(Mu2[i+1]-Mu2[i]))*h[i+1];
    FEParams[ISpline + 3] = Mu2[i+1]*h[i+1]*h[i+1]/2. +
                          ((u2_spl[i+1]-u2_spl[i])/h[i+1]-h[i+1]/6.*(Mu2[i+1]-Mu2[i]))*h[i+1];

#ifdef __SURFACT__
    FEParams[ISpline + 4  ] = -Msurf[i]*h[i+1]*h[i+1]/2. +
                          ((surf_spl[i+1]-surf_spl[i])/h[i+1]-h[i+1]/6.*(Msurf[i+1]-Msurf[i]))*h[i+1];
    FEParams[ISpline + 5] = Msurf[i+1]*h[i+1]*h[i+1]/2. +
                          ((surf_spl[i+1]-surf_spl[i])/h[i+1]-h[i+1]/6.*(Msurf[i+1]-Msurf[i]))*h[i+1];
  
    FEParams[ISpline + 6  ] = -SurfMsurf[i]*h[i+1]*h[i+1]/2. +
                          ((Surfsurf_spl[i+1]-Surfsurf_spl[i])/h[i+1]-h[i+1]/6.*(SurfMsurf[i+1]-SurfMsurf[i]))*h[i+1];
    FEParams[ISpline + 7] = SurfMsurf[i+1]*h[i+1]*h[i+1]/2. +
                          ((Surfsurf_spl[i+1]-Surfsurf_spl[i])/h[i+1]-h[i+1]/6.*(SurfMsurf[i+1]-SurfMsurf[i]))*h[i+1]; 
#endif

  }
  // ===================================================================
 } // if(UpdateU)

   teta = 1.0/N_Splines;
   T = 0;

   Param9[0] = 0;
   for(i=1;i<=N_Splines;i++) 
    Param9[i] = Params[(i-1)*10+8];

   m = 0;
   m1 = 0;
int iso = 0;
   for(j=0;j<N_E;j++)
    {
   
     T = double(m)*teta;
     for(i=1;i<=N_Splines;i++)
      {
       ISpline = (i-1)*10;
       USpline = (i-1)*8;
       FeDof   = i-1;
       if((T>=Param9[i-1]) && (T<=Param9[i]))
        {
      // further T must be from [0;1] on a subspline
         T = (T-Param9[i-1])/(Param9[i]-Param9[i-1]);
         break;
        }
      }

   phi1 = (2.*T*T - 3.*T)*T + 1.;
   phi2 = (-2.*T + 3.)*T*T;
   phi3 = (T*T - 2.*T + 1.)*T;
   phi4 = (T - 1)*T*T;

   X = Params[ISpline    ]*phi1 + Params[ISpline + 2]*phi2 +
       Params[ISpline + 4]*phi3 + Params[ISpline + 6]*phi4;
   Y = Params[ISpline + 1]*phi1 + Params[ISpline + 3]*phi2 +
       Params[ISpline + 5]*phi3 + Params[ISpline + 7]*phi4;

    if(Y < 0 || fabs(Y)<1e-8) Y = 0.0; // no penetration on solid boundary
//     if(j==0)
//     {
//       cell[j]->GetVertex(EdgeNo[j])->GetCoords(tx, ty);
//       cout << "x " <<  tx << " y " << ty << " X " << X << " Y " << Y << endl;       
//     }
    
    if(j!=0)
    cell[j]->GetVertex(EdgeNo[j])->SetCoords(X, Y);

//         if(fabs(dx0-X)>1e-4 || fabs(dy0-Y)>1e-4)
    
//     if(j==N_E -1)
//     cout << "x " << X << " y " << Y << endl;
 
//  if(UpdateU)
//        cout<<"NewX :"<<' '<< iso++ <<' '<<X<< ' '<< Y<< endl;
  
//     OutPut("NewX:"<<' '<< m <<' '<<X<<' '<< Y<<endl);
    m++;

// =========================================================================
 if(UpdateU)
 {
  // for fe values

     u0 = u1_spl[FeDof]*phi1 + u1_spl[FeDof+1]*phi2 +
              FEParams[USpline]*phi3 + FEParams[USpline + 1]*phi4;
     u1 = u2_spl[FeDof]*phi1 + u2_spl[FeDof+1]*phi2 +
              FEParams[USpline+2]*phi3 + FEParams[USpline + 3]*phi4;

#ifdef __SURFACT__
    surf = surf_spl[FeDof]*phi1 + surf_spl[FeDof+1]*phi2 +
              FEParams[USpline+4]*phi3 + FEParams[USpline + 5]*phi4;


    Surfsurf = Surfsurf_spl[FeDof]*phi1 + Surfsurf_spl[FeDof+1]*phi2 +
              FEParams[USpline+6]*phi3 + FEParams[USpline + 7]*phi4;      
#endif   
      
//    if(j==0)     
//     {
// //      cout<<"NewX :"<<' '<< ValuesUX[U_DOF[m1]] <<' '<<ValuesUY[U_DOF[m1]]<< ' '<< u0<<  ' '<< u1<< endl;
//  #ifdef __SURFACT__    
// //      cout<<"NewS :"<<' '<< surf <<' '<< Surfact[Surf_DOF[m1]]<< ' '<< Surfsurf  <<  ' '<< SurfSurfact[SurfSurf_DOF[m1]] << endl;   
// #endif       
//     }
    
    if(m1!=0) // endpoints no need to set
     {
      ValuesUX[U_DOF[m1]] = u0;
      ValuesUY[U_DOF[m1]] = u1;
      

       
#ifdef __SURFACT__      
      Surfact[Surf_DOF[m1]] = surf;
      SurfSurfact[SurfSurf_DOF[m1]] = Surfsurf;   
#endif           
     }
    m1++;
 }
// ====================================================================
// interpolation for isopoints
// no need if reparam is only for grid velo calculation
// ====================================================================
    Joint = cell[j]->GetJoint(EdgeNo[j]);
    isojoint = (TIsoBoundEdge *)Joint;
    k = isojoint->GetN_Vertices();
        
  if(UpdateU)
   {    
    if(k==ORDER-1)
     {
      IsoVertices = isojoint->GetVertices();
      for(i3=0;i3<k;i3++)
       {
       T = double(m)*teta;
       for(i=1;i<=N_Splines;i++)
        {
         ISpline = (i-1)*10;
         USpline = (i-1)*8;
         FeDof   = i-1;
         if((T>=Param9[i-1]) && (T<=Param9[i]))
          {
           // further T must be from [0;1] on a subspline
           //          cout<< ISpline << ' ' << T;
           T = (T-Param9[i-1])/(Param9[i]-Param9[i-1]);
           //          cout<< ' ' << T <<endl;
          break;
         }
       }
  
     phi1 = (2.*T*T - 3.*T)*T + 1.;
     phi2 = (-2.*T + 3.)*T*T;
     phi3 = (T*T - 2.*T + 1.)*T;
     phi4 = (T - 1)*T*T;

     X = Params[ISpline    ]*phi1 + Params[ISpline + 2]*phi2 +
         Params[ISpline + 4]*phi3 + Params[ISpline + 6]*phi4;
     Y = Params[ISpline + 1]*phi1 + Params[ISpline + 3]*phi2 +
         Params[ISpline + 5]*phi3 + Params[ISpline + 7]*phi4;

      if(j==0)
      {
//        cout<<"NewX iso:"<<' '<< iso++ <<' '<<X<<' '<< Y<<endl;
       IsoVertices[i3]->GetCoords(tx, ty);
//        OutPut("Iso x " <<  tx << " y " << ty << " X " << X << " Y " << Y << endl);            
      }
      
      
     IsoVertices[i3]->SetCoords(X, Y);     
     m++;

  // ====================================================================
  // for fe values
  // ====================================================================

    u0 = u1_spl[FeDof]*phi1 + u1_spl[FeDof+1]*phi2 +
              FEParams[USpline]*phi3 + FEParams[USpline + 1]*phi4;
    u1 = u2_spl[FeDof]*phi1 + u2_spl[FeDof+1]*phi2 +
              FEParams[USpline+2]*phi3 + FEParams[USpline + 3]*phi4;
     
#ifdef __SURFACT__
    surf = surf_spl[FeDof]*phi1 + surf_spl[FeDof+1]*phi2 +
              FEParams[USpline+4]*phi3 + FEParams[USpline + 5]*phi4;

    Surfsurf = Surfsurf_spl[FeDof]*phi1 + Surfsurf_spl[FeDof+1]*phi2 +
              FEParams[USpline+6]*phi3 + FEParams[USpline + 7]*phi4;
#endif      
//     cout<<"NewX :"<<' '<< FEParams[USpline]*phi3 <<' '<<FEParams[USpline + 1]*phi4<< ' '<< u0<<  ' '<< u1<< endl;   
     

//    if(j==0)     
//     {
//      cout<<"NewX :"<<' '<< ValuesUX[U_DOF[m1]] <<' '<<ValuesUY[U_DOF[m1]]<< ' '<< u0<<  ' '<< u1<< endl;
//  #ifdef __SURFACT__    
//      cout<<"NewS :"<<' '<< surf <<' '<< Surfact[Surf_DOF[m1]]<< ' '<< Surfsurf  <<  ' '<< SurfSurfact[SurfSurf_DOF[m1]] << endl;   
//  #endif      
//     }
    
    ValuesUX[U_DOF[m1]] = u0;
    ValuesUY[U_DOF[m1]] = u1;
#ifdef __SURFACT__
    Surfact[Surf_DOF[m1]] = surf;    
    SurfSurfact[SurfSurf_DOF[m1]] = Surfsurf;        
#endif     
    
    m1++;


// ====================================================================
    }  // for(i3=0;i3<k
    }   // if(k==ORDER-1)    

  }  // if(UpdateU) 
  else if (k==ORDER-1)
  {m += k; }
  
 }  //  for(j=0;j<N_E
//   exit(0);
   delete [] h; delete [] t; delete [] a; delete [] b;
   delete [] c; delete [] rhs; delete [] Mx; delete [] My;
   delete [] Params; delete [] Param9;  delete [] x; delete [] y;
  
 if(UpdateU)
  { 
   delete [] u1rhs;  delete [] u2rhs; delete [] u1_spl;
   delete [] u2_spl; delete [] Mu1; delete [] Mu2;
   delete [] U_DOF;  delete [] FEParams;
#ifdef __SURFACT__   
   delete [] Surf_DOF; delete [] surf_spl;  
   delete [] srhs; delete [] Msurf;
   delete [] SurfSurf_DOF; delete [] Surfsurf_spl;  
   delete [] Surfsrhs; delete [] SurfMsurf;  
#endif     
//    exit(0);
  }
  
}
 
// void ReParam_axial3D(int N_E, TBaseCell **cell, int *EdgeNo,  int *CellNo)
// {
//   int i, j, ISpline, N_Splines, N_V, ORDER, VSP, m, m1, k, i3;
//   double *h, *t;
//   double *a, *b, *c, *x, *y, teta;
//   double *rhs, *Mx, *My, *Params, *Param9;
//   double phi1, phi2, phi3, phi4, X, Y, T;
//   double dx0, dy0, dx1, dy1;
//   TIsoBoundEdge *isojoint;
//   TVertex **IsoVertices;
//   TJoint *Joint;
//   TBaseCell *Me;
// 
//   ORDER = 0;
//   VSP = TDatabase::ParamDB->VELOCITY_SPACE;
// 
//   if (abs(VSP) > 20)
//    {ORDER = abs(VSP) - 20;}
//   else if ( abs(VSP) > 10)
//    {ORDER = abs(VSP) - 10;}
//   else ORDER = abs(VSP);
// 
//   N_V = N_E+1 + N_E*(ORDER-1);
// 
//   N_Splines = N_V-1;
//   h = new double[N_Splines+1];
//   t = new double[N_Splines+1];
//   a = new double[N_Splines+1];
//   b = new double[N_Splines+1];
//   c = new double[N_Splines+1];
//   rhs = new double[N_Splines+1];
// 
//   Mx = new double[N_Splines+1];
//   My = new double[N_Splines+1];
//   Params = new double [10*N_Splines];
//   Param9 = new double [N_Splines+1];
// 
//   x = new double[N_V];
//   y = new double[N_V];
// 
//    m = 0;
//    m1 = 0;
//    for(i=0;i<N_E;i++) // i<N_E
//    {
//     Me = cell[i];
//     Me->GetVertex(EdgeNo[i])->GetCoords(x[m], y[m]);
//     m++;
// 
//     Joint = cell[i]->GetJoint(EdgeNo[i]);
//     isojoint = (TIsoBoundEdge *)Joint;
//     k = isojoint->GetN_Vertices();
//     if(k==ORDER-1)
//      {
//       IsoVertices = isojoint->GetVertices();
//       for(i3=0;i3<k;i3++)
//        {   
//         IsoVertices[i3]->GetCoords(x[m], y[m]);
//         m++;
//        } 
//      }
//     else
//      {
//       // only second order conforming elements implimented
//       cout<< " No match in isopoints per free edge "<<endl;
//       exit(0);
//      }
// 
//    } // for(i=0;i<N_E
// 
// 
// //   end vertex of the freeboundary
//   k = cell[N_E-1]->GetN_Edges();
//   cell[N_E-1]->GetVertex((EdgeNo[N_E-1]+1) % k)->GetCoords(x[m], y[m]);
// 
// //  for(i=0;i<N_V;i++)
// //    OutPut("OldX: "<< i <<' '<<x[i] <<' '<< y[i] <<endl);
// //  OutPut("OldX: "<< i <<' '<<x[0] <<' '<< y[0] <<endl);
// //  OutPut("OldX: "<< i <<' '<<x[N_V-2] <<' '<< y[N_V-2] <<endl);
// //  OutPut("OldX: "<< i <<' '<<x[N_V-1] <<' '<< y[N_V-1] <<endl);
// // cout << "Surfact[Surf_DOF[0]] " <<Surfact[Surf_DOF[0]] << " Surfact[Surf_DOF[m1]] " << Surfact[Surf_DOF[m1-1]] <<endl;
// 
//   h[0] = 0.0; t[0] = 0.0;
// 
//  for(i=1;i<=N_Splines;i++)
//   {
//     h[i] = sqrt((x[i]-x[i-1])*(x[i]-x[i-1])+(y[i]-y[i-1])*(y[i]-y[i-1]));
//     t[i] = t[i-1] + h[i];
//   }
// 
//   dx0 = (x[1]-x[0])/h[1];
//   dy0 = (y[1]-y[0])/h[1];
// 
//   dx1 = (x[N_Splines]-x[N_Splines-1])/h[N_Splines];
//   dy1 = (y[N_Splines]-y[N_Splines-1])/h[N_Splines];
// 
// 
//   a[0] = 2.; c[0] = 1.; rhs[0] = -6./h[1]*(dx0 - (x[1]-x[0])/h[1]);
//   for(i=1;i<N_Splines;i++)
//   {
//     a[i] = 2.;  
//     b[i] = h[i]/(h[i]+h[i+1]); // \mu_i in PhD thesis
//     c[i] = h[i+1]/(h[i]+h[i+1]); // \lambda_i in PhD thesis
//     rhs[i] = 6./(h[i]+h[i+1])*((x[i+1]-x[i])/h[i+1]-(x[i]-x[i-1])/h[i]);
//   }
//   b[N_Splines] = 1.; a[N_Splines] = 2.;
//   rhs[N_Splines] = 6./h[N_Splines]*(dx1 - (x[N_Splines]-x[N_Splines-1])/h[N_Splines]);
// 
//   Solver_3dia(N_Splines, a, b, c, rhs, Mx);
// 
//   rhs[0] = -6./h[1]*(dy0 - (y[1]-y[0])/h[1]);
//   for(i=1;i<N_Splines;i++)
//   {
//     rhs[i] = 6./(h[i]+h[i+1])*((y[i+1]-y[i])/h[i+1]-(y[i]-y[i-1])/h[i]);
//   }
//   rhs[N_Splines] = 6./h[N_Splines]*(dy1 - (y[N_Splines]-y[N_Splines-1])/h[N_Splines]);
// 
//   Solver_3dia(N_Splines, a, b, c, rhs, My);
// 
//   for(i=0;i<N_Splines;i++)
//   {
//     ISpline = i*10;
//     Params[ISpline    ] = x[i]; 
//     Params[ISpline + 1] = y[i];
//     Params[ISpline + 2] = x[i+1]; 
//     Params[ISpline + 3] = y[i+1];
//     Params[ISpline + 4] = -Mx[i]*h[i+1]*h[i+1]/2. +
//                           ((x[i+1]-x[i])/h[i+1]-h[i+1]/6.*(Mx[i+1]-Mx[i]))*h[i+1];
// 
// //     Params[ISpline + 4] = Mx[i]*h[i];
//     Params[ISpline + 5] = -My[i]*h[i+1]*h[i+1]/2. +
//                           ((y[i+1]-y[i])/h[i+1]-h[i+1]/6.*(My[i+1]-My[i]))*h[i+1];
//     Params[ISpline + 6] = Mx[i+1]*h[i+1]*h[i+1]/2. +
//                           ((x[i+1]-x[i])/h[i+1]-h[i+1]/6.*(Mx[i+1]-Mx[i]))*h[i+1];
// 
// //     Params[ISpline + 6] = -Mx[i+1];
//     Params[ISpline + 7] = My[i+1]*h[i+1]*h[i+1]/2. +
//                           ((y[i+1]-y[i])/h[i+1]-h[i+1]/6.*(My[i+1]-My[i]))*h[i+1];
//     Params[ISpline + 8] = t[i+1]/t[N_Splines];
//     Params[ISpline + 9] = 0.;
// 
//    //cout<<"  "<<Params[ISpline + 8]<<'\t'<<Params[ISpline + 9]<<endl;
//   }
// 
//    teta = 1.0/N_Splines;
//    T = 0;
// 
//    Param9[0] = 0;
//    for(i=1;i<=N_Splines;i++) 
//     Param9[i] = Params[(i-1)*10+8];
// 
//    m = 0;
//    for(j=0;j<N_E;j++)
//     {
//      T = double(m)*teta;
//      for(i=1;i<=N_Splines;i++)
//       {
//        ISpline = (i-1)*10;
//        if((T>=Param9[i-1]) && (T<=Param9[i]))
//         {
//       // further T must be from [0;1] on a subspline
//          T = (T-Param9[i-1])/(Param9[i]-Param9[i-1]);
//          break;
//         }
//       }
// 
//    phi1 = (2.*T*T - 3.*T)*T + 1.;
//    phi2 = (-2.*T + 3.)*T*T;
//    phi3 = (T*T - 2.*T + 1.)*T;
//    phi4 = (T - 1)*T*T;
// 
//    X = Params[ISpline    ]*phi1 + Params[ISpline + 2]*phi2 +
//        Params[ISpline + 4]*phi3 + Params[ISpline + 6]*phi4;
//    Y = Params[ISpline + 1]*phi1 + Params[ISpline + 3]*phi2 +
//        Params[ISpline + 5]*phi3 + Params[ISpline + 7]*phi4;
// 
//    if(Y < 0 || fabs(Y)<1e-8) Y = 0.0; // no penetration on solid boundary
//     cell[j]->GetVertex(EdgeNo[j])->SetCoords(X, Y);
// 
// //     OutPut("NewX:"<<' '<< m <<' '<<X<<' '<< Y<<endl);
//     m++;
// 
//     Joint = cell[j]->GetJoint(EdgeNo[j]);
//     isojoint = (TIsoBoundEdge *)Joint;
//     k = isojoint->GetN_Vertices();
//     if(k==ORDER-1)
//      {
//       IsoVertices = isojoint->GetVertices();
//       for(i3=0;i3<k;i3++)
//        {
//        T = double(m)*teta;
//        for(i=1;i<=N_Splines;i++)
//         {
//          ISpline = (i-1)*10;
//          if((T>=Param9[i-1]) && (T<=Param9[i]))
//           {
//            // further T must be from [0;1] on a subspline
//            T = (T-Param9[i-1])/(Param9[i]-Param9[i-1]);
//           break;
//          }
//        }
//   
//      phi1 = (2.*T*T - 3.*T)*T + 1.;
//      phi2 = (-2.*T + 3.)*T*T;
//      phi3 = (T*T - 2.*T + 1.)*T;
//      phi4 = (T - 1)*T*T;
// 
//      X = Params[ISpline    ]*phi1 + Params[ISpline + 2]*phi2 +
//          Params[ISpline + 4]*phi3 + Params[ISpline + 6]*phi4;
//      Y = Params[ISpline + 1]*phi1 + Params[ISpline + 3]*phi2 +
//          Params[ISpline + 5]*phi3 + Params[ISpline + 7]*phi4;
// 
//         IsoVertices[i3]->SetCoords(X, Y);
// //         OutPut("NewX:"<<' '<< m <<' '<<X<<' '<< Y<<endl);
//         m++;
//        }
//      }
//    }  //  for(j=0;j<N_E
// 
// 
//    delete [] h; delete [] t; delete [] a; delete [] b;
//    delete [] c; delete [] rhs; delete [] Mx; delete [] My;
//    delete [] Params; delete [] Param9;  delete [] x; delete [] y;
// 
// } // ReParam_axial3D
// 

void EqDist_Pts(int N_V, double *newx, double *newy)
{
  int i, j, ISpline, N_Splines,   m, m1, k, i3;
  double *h, *t;
  double *a, *b, *c,   teta;
  double *rhs, *Mx, *My, *Params, *Param9;
  double phi1, phi2, phi3, phi4, X, Y, T;
  double dx0, dy0, dx1, dy1;
  double *x, *y;

  N_Splines = N_V-1;
  h = new double[N_Splines+1];
  t = new double[N_Splines+1];
  a = new double[N_Splines+1];
  b = new double[N_Splines+1];
  c = new double[N_Splines+1];
  rhs = new double[N_Splines+1];

  Mx = new double[N_Splines+1];
  My = new double[N_Splines+1];
  Params = new double [10*N_Splines];
  Param9 = new double [N_Splines+1];

  x = new double[N_V];
  y = new double[N_V]; 
  
  memcpy(x, newx, N_V*SizeOfDouble);
  memcpy(y, newy, N_V*SizeOfDouble);  
  
  
//   cout << "N_V " << N_V <<endl;
//   exit(0);
  
  h[0] = 0.0; t[0] = 0.0;

 for(i=1;i<=N_Splines;i++)
  {
    h[i] = sqrt((x[i]-x[i-1])*(x[i]-x[i-1])+(y[i]-y[i-1])*(y[i]-y[i-1]));
    t[i] = t[i-1] + h[i];
  }

  dx0 = (x[1]-x[0])/h[1];
  dy0 = (y[1]-y[0])/h[1];

  dx1 = (x[N_Splines]-x[N_Splines-1])/h[N_Splines];
  dy1 = (y[N_Splines]-y[N_Splines-1])/h[N_Splines];


  a[0] = 2.; c[0] = 1.; rhs[0] = -6./h[1]*(dx0 - (x[1]-x[0])/h[1]);
  for(i=1;i<N_Splines;i++)
  {
    a[i] = 2.;  
    b[i] = h[i]/(h[i]+h[i+1]); // \mu_i in PhD thesis
    c[i] = h[i+1]/(h[i]+h[i+1]); // \lambda_i in PhD thesis
    rhs[i] = 6./(h[i]+h[i+1])*((x[i+1]-x[i])/h[i+1]-(x[i]-x[i-1])/h[i]);
  }
  b[N_Splines] = 1.; a[N_Splines] = 2.;
  rhs[N_Splines] = 6./h[N_Splines]*(dx1 - (x[N_Splines]-x[N_Splines-1])/h[N_Splines]);

  Solver_3dia(N_Splines, a, b, c, rhs, Mx);

  rhs[0] = -6./h[1]*(dy0 - (y[1]-y[0])/h[1]);
  for(i=1;i<N_Splines;i++)
  {
    rhs[i] = 6./(h[i]+h[i+1])*((y[i+1]-y[i])/h[i+1]-(y[i]-y[i-1])/h[i]);
  }
  rhs[N_Splines] = 6./h[N_Splines]*(dy1 - (y[N_Splines]-y[N_Splines-1])/h[N_Splines]);

  Solver_3dia(N_Splines, a, b, c, rhs, My);

  for(i=0;i<N_Splines;i++)
  {
    ISpline = i*10;
    Params[ISpline    ] = x[i]; 
    Params[ISpline + 1] = y[i];
    Params[ISpline + 2] = x[i+1]; 
    Params[ISpline + 3] = y[i+1];
    Params[ISpline + 4] = -Mx[i]*h[i+1]*h[i+1]/2. +
                          ((x[i+1]-x[i])/h[i+1]-h[i+1]/6.*(Mx[i+1]-Mx[i]))*h[i+1];

//     Params[ISpline + 4] = Mx[i]*h[i];
    Params[ISpline + 5] = -My[i]*h[i+1]*h[i+1]/2. +
                          ((y[i+1]-y[i])/h[i+1]-h[i+1]/6.*(My[i+1]-My[i]))*h[i+1];
    Params[ISpline + 6] = Mx[i+1]*h[i+1]*h[i+1]/2. +
                          ((x[i+1]-x[i])/h[i+1]-h[i+1]/6.*(Mx[i+1]-Mx[i]))*h[i+1];

//     Params[ISpline + 6] = -Mx[i+1];
    Params[ISpline + 7] = My[i+1]*h[i+1]*h[i+1]/2. +
                          ((y[i+1]-y[i])/h[i+1]-h[i+1]/6.*(My[i+1]-My[i]))*h[i+1];
    Params[ISpline + 8] = t[i+1]/t[N_Splines];
    Params[ISpline + 9] = 0.;

   //cout<<"  "<<Params[ISpline + 8]<<'\t'<<Params[ISpline + 9]<<endl;
  }

   teta = 1.0/N_Splines;
   T = 0;

   Param9[0] = 0;
   for(i=1;i<=N_Splines;i++) 
    Param9[i] = Params[(i-1)*10+8];

   
   
//    cout <<"test spline "<<endl;
//    exit(0);
   
   m = 0;
   for(j=0;j<N_V;j++)
    {
     T = double(m)*teta;
     for(i=1;i<=N_Splines;i++)
      {
       ISpline = (i-1)*10;
       if((T>=Param9[i-1]) && (T<=Param9[i]))
        {
      // further T must be from [0;1] on a subspline
         T = (T-Param9[i-1])/(Param9[i]-Param9[i-1]);
         break;
        }
      }

     phi1 = (2.*T*T - 3.*T)*T + 1.;
     phi2 = (-2.*T + 3.)*T*T;
     phi3 = (T*T - 2.*T + 1.)*T;
     phi4 = (T - 1)*T*T;

     newx[j] = Params[ISpline    ]*phi1 + Params[ISpline + 2]*phi2 +
               Params[ISpline + 4]*phi3 + Params[ISpline + 6]*phi4;
     newy[j] = Params[ISpline + 1]*phi1 + Params[ISpline + 3]*phi2 +
               Params[ISpline + 5]*phi3 + Params[ISpline + 7]*phi4;
        
//      cout<<" diff x : "<< x[j] << " " << newx[j] << " y : "<< y[j] << " " << newy[j]<<endl;     
             m++;  
    }  //  for(j=0;j<N_V

  newx[0] = x[0];
  newy[0] = y[0];
  
  newx[N_V-1] = x[N_V-1];
  newy[N_V-1] = y[N_V-1];  

   delete [] h; delete [] t; delete [] a; delete [] b;
   delete [] c; delete [] rhs; delete [] Mx; delete [] My;
   delete [] Params; delete [] Param9;  delete [] x; delete [] y; 
  
} // ReParam_axial3D



/** convert current grid to vector-values FE function */
void GridTOVelocity(TFEVectFunct2D *velo, double t, double dt, double &max_r)
{
 int N_U, i;
 TFESpace2D *Velospace;
 double *x, *y, *ValuesU, u, v, newx, newy, L, r;


 Velospace = velo->GetFESpace2D();

 ValuesU = velo->GetValues();
 N_U = velo->GetLength();


 x = new double[N_U];
 y = new double[N_U];

 Velospace->GetDOFPosition(x, y);

//  cout << " dt " << dt << " time " << t <<endl;
// exit(0);
max_r = -1e10;
 for(i=0;i<N_U;i++)
  {
   r = sqrt(x[i]*x[i] + y[i]*y[i]);
   if(max_r<r) max_r =r;
 }

OutPut("Time, Radius of the Sphere : " << TDatabase::TimeDB->CURRENTTIME<< " " << max_r<< "  "  << 2.*Pi*max_r*max_r <<endl;)

 for(i=0;i<N_U;i++)
  {
//    newx = cos(100.*t)*x[i] - sin(100.*t)*y[i];
//    newy = sin(100.*t)*x[i] + cos(100.*t)*y[i];
   u = 0.;
   v = 0.;

   r = x[i]*x[i] + y[i]*y[i] ;
   if(sqrt(r)>0.25 )
    {
     L = pow(r,3./2.);
     u = x[i]/(L);
     v = y[i]/(L);
    }


   ValuesU[i] = u;
   ValuesU[N_U+i] =  v; 
    }


  delete []x;
  delete []y;
}

void  MoveGrid(TFEVectFunct2D *GridPos, TFEVectFunct2D *AuxGridPos, 
               TFESpace2D  *Surf_space, double t, double dt)
{
 int i, j, k, l, N_G, N_Cells, N_Vertices, *GridBeginIndex, *GridGlobalNumbers, *DOF;
 int N_Edges;
 double *X, *Y, *NewX, *NewY, u, L, r, u_r, u_z, x, y, max_r=-1e10;
 TCollection *Coll;
 TBaseCell *cell;
 TFESpace2D  *GridSpace;
 TJoint *joint;
 TIsoBoundEdge *isojoint;
 TVertex **Vertices;
 int IIso, m, *VeloGlobalNumbers, *VeloBeginIndex;
 FE2D FEId;
 TFEDesc2D *FEDesc;
 int *VeloDOF, *JointDOF;

 GridSpace = GridPos->GetFESpace2D();
 Coll = GridSpace->GetCollection();
 N_Cells = Coll->GetN_Cells();
 GridBeginIndex = GridSpace->GetBeginIndex();
 GridGlobalNumbers = GridSpace->GetGlobalNumbers();

 GridPos->GridToData();
 AuxGridPos->GridToData();


 N_G=GridPos->GetLength();
 X = GridPos->GetValues();
 Y = X + N_G;

 NewX = AuxGridPos->GetValues();
 NewY = NewX + N_G;


 VeloGlobalNumbers = Surf_space->GetBeginIndex();
 VeloBeginIndex = Surf_space->GetGlobalNumbers();


//  for(i=0;i<N_G;i++)
//   {
//    r = X[i]*X[i] + Y[i]*Y[i] ;
//    if(max_r<r) max_r =r;
//  }


//  u = 10.;
 for(i=0;i<N_G;i++)
  {
   r = X[i]*X[i] + Y[i]*Y[i] ;
   if(sqrt(r)>0.25 )
    {
     L = pow(r,3./2.);
     u_r = X[i]/(L);
     u_z = Y[i]/(L);

   NewX[i] = X[i] + dt*u_r;
   NewY[i] = Y[i] + dt*u_z;

//   cout<<  " X[i] " <<X[i] << " Y[i] "<< Y[i] << " NewX[i] "  <<NewX[i] << " NewY[i] " <<NewY[i]<<endl;

   } //    if(sqrt(r)>0.5 
  } // for(i=0;i<N_G;
//   IIso = 0;
  for(i=0;i<N_Cells;i++)
   {
    cell = Coll->GetCell(i);
    DOF = GridGlobalNumbers + GridBeginIndex[i];
    N_Vertices = cell->GetN_Vertices();
    switch(N_Vertices)
    {
      case 3:
        for(j=0;j<N_Vertices;j++)
        {
          k = DOF[j];
          cell->GetVertex(j)->SetCoords(NewX[k], NewY[k]);
        }
      break;

      case 4:
        k = DOF[0];
        cell->GetVertex(0)->SetCoords(NewX[k], NewY[k]);
        k = DOF[1];
        cell->GetVertex(1)->SetCoords(NewX[k], NewY[k]);
        k = DOF[3];
        cell->GetVertex(2)->SetCoords(NewX[k], NewY[k]);
        k = DOF[2];
        cell->GetVertex(3)->SetCoords(NewX[k], NewY[k]);
      break;
    } // endswitch



    N_Edges = cell->GetN_Edges();
    for(j=0;j<N_Edges;j++)
     {
      joint = cell->GetJoint(j);
      if(joint->GetType() == IsoBoundEdge)
      {
        isojoint = (TIsoBoundEdge *)joint;
        k = isojoint->GetN_Vertices();
        Vertices = isojoint->GetVertices();
        FEId = Surf_space->GetFE2D(i, cell);
        FEDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
        m = FEDesc->GetN_JointDOF();
        if(m == k+2)
         {
          JointDOF = FEDesc->GetJointDOF(j);
          VeloDOF =  VeloGlobalNumbers+VeloBeginIndex[i];
          for(l=0;l<k;l++)
          {
            m = VeloDOF[JointDOF[l+1]];
            Vertices[l]->GetCoords(x, y);
            r = x*x + y*y;
            if(sqrt(r)>0.25)
             {
              L = pow(r,3./2.);
              u_r = x/(L);
              u_z = y/(L);
              x += dt*u_r;
              y += dt*u_z;

             Vertices[l]->SetCoords(x, y);
            } //  if(sqrt(r)
//             IIso++;
          } // endfor l

         } // if(m == k+2)
        else
        {
          // approximation order of isoparametric boundary and velocity
          // element must be the same
          Error("No match in isoparametric case" << endl);
          exit(-1);
        }
       } // endif
     } // endfor j
   }
// cout<< " r " << r <<endl;
//  AuxGridPos->DataToGrid();
}

