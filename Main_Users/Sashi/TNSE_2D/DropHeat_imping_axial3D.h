// Navier-Stokes problem, Benchmark channel
// circle moves up and down
//
// u(x,y) = unknown
// p(x,y) = unknown

// // #include <MovingNavierStokesGradDiv.h>
#define __ENERGY__
// #define __ConvergenceStudy__

void ExampleFile()
{
  #define __AXIAL3D__ 

  OutPut("Droplet Impinging Heat_axial3D.h" << endl) ;
    
}

// ========================================================================
// initial solution
// ========================================================================

void InitialU1(double x, double y, double *values)
{
  const double teta=TDatabase::ParamDB->IMPACT_ANGLE;

//   values[0] = cos((Pi/180)*teta);
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void InitialU2(double x, double y, double *values)
{

  const double teta=TDatabase::ParamDB->IMPACT_ANGLE;

  if(y==0.) values[0] = 0.;
  else values[0] = -(1-cos((Pi/180)*teta));

#ifdef __ConvergenceStudy__
  values[0] = 0.;
#endif 
  
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
  switch(i)
      {
        case 0:
        case 2:
             cond = SLIP_FRICTION_PENETRATION_RESISTANCE;
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

void BoundCondition_output(int i, double t, BoundCond &cond)
{
  switch(i)
      {
        case 0:
        case 2:
             cond = SLIP_FRICTION_PENETRATION_RESISTANCE;
             TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION = 1;
        break;

        case 1:
             cond = FREESURF;
             TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
        break;


       default:
          cond = NEUMANN;
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
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *x, double *y,
               double **parameters, double **coeffs)
{
  static double eps = 1./TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff;
// inclined angle of the plane on which the droplet deforms
  static double theta = TDatabase::ParamDB->P1;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    
    if(TDatabase::ParamDB->FR_NR == 0)
     { coeff[1] = sin((Pi/180)*theta); }
    else
     { coeff[1] = sin((Pi/180)*theta)/TDatabase::ParamDB->FR_NR; } // f1 
     
    if(TDatabase::ParamDB->FR_NR == 0)
     {  coeff[2] = 0.; }
    else
     { coeff[2] = -cos((Pi/180)*theta)/TDatabase::ParamDB->FR_NR; } // f2 - => gravity in opposite direction 
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

// ========================================================================
// boundary conditions for heat
// ========================================================================
void HeatBoundCondition(int i, double t, BoundCond &cond)
{
  switch(i)
      {
        case 0:
//              cond = NEUMANN;
//              cond = DIRICHLET;
             cond = FREESURF;
        break;

        case 1:
             cond = FREESURF;
        break;

        case 2:
             cond = NEUMANN;
        break;

        case 3:
	case 4:
        case 5:
        case 6:  
               cond = NEUMANN;
	 break;    
	 
//          case 4:
//              cond = DIRICHLET;
//  
//         break;


      default:
            Error("Unknown Boundary component ref example file" << endl);
         exit(0);
       }
     }

void InitialT(double x, double y, double *values)
{
  
double T_W = TDatabase::ParamDB->P15;
  if(y<=0)
   values[0] = T_W;
  else
   values[0] = 0;
  
  
#ifndef __ENERGY__
   values[0] = 0;
#endif
  
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void TBoundValue(int BdComp, double Param, double &value)
{
//  x = Param; ref Assemble2D.C
//   double T_W = TDatabase::ParamDB->P15;
//   double r = fabs(Param);

  switch(BdComp)
      {
       case 0:
       case 1:
       case 2:
       case 3:
       case 4:  
       case 5:
       case 6: 
           value = 0.;
        break;

        
//         case 4:  
//         
//            value = T_W;
//         break;


      default:
            Error("ref example file" << endl);
         exit(0);
     }
}

void HeatCoeffs(int n_points, double *x, double *y,
        double **parameters, double **coeffs)
{
  int i;
  double *coeff;
  double r2;
 
  double Pe = TDatabase::ParamDB->PE_NR;
  double Pe_SolidFact = TDatabase::ParamDB->HEAT_SOLID_SURFACE_FACTOR;
 
  static double eps_l = 1.0/Pe;
  static double eps_s = eps_l*Pe_SolidFact;
  
  int Phase_ID = (int)coeffs[0][0]; //see Line 226 of DiscreteForm2D.C
  
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    if(Phase_ID==0) // liquid
    { coeff[0] = eps_l; }
    else
    { coeff[0] = eps_s; }  
    
    coeff[1] = 0; // f
    coeff[2] = 0; // c
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
                       TFEVectFunct2D *Velocity, TFEFunction2D *Energy, bool UpdateU)
{
  int i, j, ISpline, N_Splines, N_V, ORDER, VSP, m, m1, k, i3, USpline, FeDof;
  int *VeloBeginIndex, *VeloGlobalNumbers, *JointDOF, *DOF, N_DOF_Joint, *U_DOF;
  int *TBeginIndex, *TGlobalNumbers, *TJointDOF, *TDOF, TN_DOF_Joint,*T_DOF;  
  
  double *h, *t, u0, u1, u2;
  double *a, *b, *c, *x, *y, teta;
  double *rhs, *u1rhs, *u2rhs, *Trhs, *Mx, *My,*Mu1, *Mu2, *MT, *Params, *Param9, *FEParams;
  double phi1, phi2, phi3, phi4, X, Y, T;
  double dx0, dy0, dx1, dy1, surf, *u1_spl, *u2_spl, *T_spl;
  double *ValuesUX, *ValuesUY, *Heat;  
  
  TIsoBoundEdge *isojoint;
  TVertex **IsoVertices;
  TJoint *Joint;
  TBaseCell *Me;
  TFESpace2D *VelocitySpace, *TSpace;
  FE2D FEId, TFEId;
  TFE2D *ele, *Tele;
  TFEDesc2D *FeDesc, *TFeDesc;
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
   FEParams = new double [3*2*N_Splines]; // 3 fe functions, u1, u2, energy
   U_DOF = new int[N_V];
  
   VeloBeginIndex = VelocitySpace->GetBeginIndex();
   VeloGlobalNumbers = VelocitySpace->GetGlobalNumbers();
   ValuesUX = Velocity->GetValues();
   ValuesUY = ValuesUX + Velocity->GetLength();

   T_DOF = new int[N_V];
   T_spl = new double[N_Splines+1];
   Trhs = new double[N_Splines+1];   
   MT = new double[N_Splines+1];
    
   TSpace = Energy->GetFESpace2D();
   TBeginIndex = TSpace->GetBeginIndex();
   TGlobalNumbers = TSpace->GetGlobalNumbers();
   Heat = Energy->GetValues();
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
  
    // for energy
    TFEId = TSpace->GetFE2D(Me->GetGlobalCellNo(), Me);
    Tele = TFEDatabase2D::GetFE2D(TFEId);
    TFeDesc = Tele->GetFEDesc2D();   // fe descriptor
    TJointDOF = TFeDesc->GetJointDOF(EdgeNo[i]);
    TN_DOF_Joint = TFeDesc->GetN_JointDOF();
    TDOF = TGlobalNumbers + TBeginIndex[Me->GetGlobalCellNo()];
    
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
         T_DOF[m1] = TDOF[TJointDOF[i3]]; // needed for later update
         T_spl[m1] = Heat[T_DOF[m1]];
         m1++;
       }

     } //  if(UpdateU)       
              
    } // for(i=0;i<N_E
 
//   end vertex of the freeboundary
  k = cell[N_E-1]->GetN_Edges();
  cell[N_E-1]->GetVertex((EdgeNo[N_E-1]+1) % k)->GetCoords(x[m], y[m]);

 
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
  
  //energy
  for(i=1;i<N_Splines;i++)
   {
     u0 = T_spl[i-1];
     u1 = T_spl[i];
     u2 = T_spl[i+1];

     Trhs[i] = 6./(h[i]+h[i+1])*((u2-u1)/h[i+1]-(u1-u0)/h[i]);
    }

   Trhs[0] = Trhs[1];
   Trhs[N_Splines] = Trhs[N_Splines-1];

  Solver_3dia(N_Splines, a, b, c, Trhs, MT);  

// ===============================================================

  for(i=0;i<N_Splines;i++)
  {
    ISpline = i*6;

    FEParams[ISpline  ] = -Mu1[i]*h[i+1]*h[i+1]/2. +
                          ((u1_spl[i+1]-u1_spl[i])/h[i+1]-h[i+1]/6.*(Mu1[i+1]-Mu1[i]))*h[i+1];
    FEParams[ISpline + 1] = Mu1[i+1]*h[i+1]*h[i+1]/2. +
                          ((u1_spl[i+1]-u1_spl[i])/h[i+1]-h[i+1]/6.*(Mu1[i+1]-Mu1[i]))*h[i+1];


    FEParams[ISpline + 2] = -Mu2[i]*h[i+1]*h[i+1]/2. +
                          ((u2_spl[i+1]-u2_spl[i])/h[i+1]-h[i+1]/6.*(Mu2[i+1]-Mu2[i]))*h[i+1];
    FEParams[ISpline + 3] = Mu2[i+1]*h[i+1]*h[i+1]/2. +
                          ((u2_spl[i+1]-u2_spl[i])/h[i+1]-h[i+1]/6.*(Mu2[i+1]-Mu2[i]))*h[i+1];
 
 
    FEParams[ISpline + 4] = -MT[i]*h[i+1]*h[i+1]/2. +
                          ((T_spl[i+1]-T_spl[i])/h[i+1]-h[i+1]/6.*(MT[i+1]-MT[i]))*h[i+1];
    FEParams[ISpline + 5] = MT[i+1]*h[i+1]*h[i+1]/2. +
                          ((T_spl[i+1]-T_spl[i])/h[i+1]-h[i+1]/6.*(MT[i+1]-MT[i]))*h[i+1]; 
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
 
   for(j=0;j<N_E;j++)
    {
   
     T = double(m)*teta;
     for(i=1;i<=N_Splines;i++)
      {
       ISpline = (i-1)*10;
       USpline = (i-1)*6;
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
    cell[j]->GetVertex(EdgeNo[j])->SetCoords(X, Y);

//         if(fabs(dx0-X)>1e-4 || fabs(dy0-Y)>1e-4)

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

     T = T_spl[FeDof]*phi1 + T_spl[FeDof+1]*phi2 +
              FEParams[USpline+4]*phi3 + FEParams[USpline + 5]*phi4;
                   
    if(j!=0) // endpoints no need to set
     {
      ValuesUX[U_DOF[m1]] = u0;
      ValuesUY[U_DOF[m1]] = u1;
      Heat[T_DOF[m1]] = T;
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
         USpline = (i-1)*6;
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

     IsoVertices[i3]->SetCoords(X, Y);
     m++;

  // ====================================================================
  // for fe values
  // ====================================================================

    u0 = u1_spl[FeDof]*phi1 + u1_spl[FeDof+1]*phi2 +
              FEParams[USpline]*phi3 + FEParams[USpline + 1]*phi4;
    u1 = u2_spl[FeDof]*phi1 + u2_spl[FeDof+1]*phi2 +
              FEParams[USpline+2]*phi3 + FEParams[USpline + 3]*phi4;

    T = T_spl[FeDof]*phi1 + T_spl[FeDof+1]*phi2 +
              FEParams[USpline+4]*phi3 + FEParams[USpline + 5]*phi4;



    ValuesUX[U_DOF[m1]] = u0;
    ValuesUY[U_DOF[m1]] = u1;
    Heat[T_DOF[m1]] = T;    
    m1++;
 
// ====================================================================
    }  // for(i3=0;i3<k
    }   // if(k==ORDER-1)    

  }  // if(UpdateU) 
  else if (k==ORDER-1)
  {m += k; }
  
 }  //  for(j=0;j<N_E


//   if(UpdateU)       
//   {
//     cout << " Energy UpdateU" << endl;
//    exit(0); 
//   }  
//   
  
   delete [] h; delete [] t; delete [] a; delete [] b;
   delete [] c; delete [] rhs; delete [] Mx; delete [] My;
   delete [] Params; delete [] Param9;  delete [] x; delete [] y;
  
 if(UpdateU)
  { 
   delete [] u1rhs;  delete [] u2rhs; delete [] u1_spl;
   delete [] u2_spl; delete [] Mu1; delete [] Mu2;
   delete [] U_DOF;  delete [] FEParams;
   
   delete []  T_spl; delete [] MT; delete [] T_DOF;
   delete []  Trhs;
  }
  
}

  