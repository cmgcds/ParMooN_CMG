

void ExampleFile()
{  
  #define __IMPINGINGDROPLET__
  #define __AXIAL3D__
//   #define __EXPERIMENTAL__
 #define __SURFACT__
 #define __SURFACTDIFFTEST__
 
 
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
 double t=0;
//  double Re = TDatabase::ParamDB->RE_NR;
//  double Pr = TDatabase::ParamDB->Pr_NR;

  static double eps = 1./TDatabase::ParamDB->REACTOR_P17;

  double phi = atan2(y-1., x);
  phi += Pi/2.;
  
//   if(fabs(sqrt(x*x + y*y))
//    cout << t << " Phi " << (180./Pi)*phi << " X " << x << " Y " << y << endl;
  values[0] = 0.5 +  0.5*cos(phi)*exp(-2.*t*eps); // diffusion test
  
//   if(phi<1e-8 && x<1.)
//     values[0] = 0.;
//     
//   values[0] = 0.5 + 0.5*cos(phi); // expanding static sphere
//   values[0] = 0.9*fabs(cos(phi)); // Oscllating sphere
//   values[0] = sqrt((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5));

//   values[0] = TDatabase::ParamDB->REACTOR_P12;  // static bubble
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void ExactS(double x, double y, double *values)
{
 double L, t = TDatabase::TimeDB->CURRENTTIME;
//  double Re = TDatabase::ParamDB->RE_NR;
//  double Pr = TDatabase::ParamDB->PR_NR;

  static double eps = 1./TDatabase::ParamDB->REACTOR_P17;

  double phi = atan2(y-1, x);
  phi += Pi/2.;
//    cout << "Phi " << phi << " X " << x << " Y " << y << endl;
//   L = ( 1./(16.*Pi*t + 1.) ); // expanding static sphere
  values[0] = 0.5 + 0.5*cos(phi)*exp(-2.*t*eps); // diffusion test
  
//   if(phi<1e-8 && x<1.)
//     values[0] = 0.;
  
//   values[0] = (0.5 + 0.5*cos(phi)) * (pow(L, 1./3.)); // expanding static sphere
//   values[0] = 0.9*fabs(cos(phi)); // Oscllating sphere
//   values[0] = 0.;

//   values[0] = TDatabase::ParamDB->REACTOR_P12;  // static bubble  
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
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

void FreeSurf_axial3D_new(TSquareMatrix2D *A11, TSquareMatrix2D *A22,
                          double *rhs1, double *rhs2,
                          BoundCondFunct2D *BoundaryCondition,
                          double dt, double *Ucl, TFEFunction2D *Surfact, double *param)
{
  int i, j, k, l, DOF_R, DOF_L, m, mm;
  int *KCol, *RowPtr, *JointDOF, N_DOF;
  int N_LinePoints;
  int *BeginIndex, *GlobalNumbers, *DOF, TestDOF, AnsatzDOF;
  int index1, index2;
  int count=0, count1=0, count2=0;
  int N_BaseFunct, *N_BaseFuncts;
  int JointNumbers[MAXN_JOINTS], IJoint, N_IsoJoints;  
  int comp, N_U, test_L=1, test_R=1;
  int N_Cells, N_Vertices, N_Edges, Semi_implicit=0;  
    
  double r2, r;  
  double  X_B[100], Y_B[100], r_axial, d1, d2, e1, e2, ngrad_test, ngrad_ansatz, tmp;
  double val, theta, factor1, factor2, angle;
  double *ValuesA11, *ValuesA12, *ValuesA21, *ValuesA22;  
  double *LineWeights, *zeta;
  double x0, y0, x1, y1,tx,ty,mod_t, x, y;
  double **uref, **uxiref, **uetaref;
  double uorig[MaxN_BaseFunctions2D], uxorig[MaxN_BaseFunctions2D];
  double uyorig[MaxN_BaseFunctions2D];
  double t0, t1, n0, n1, normn, line_wgt;
  double  Gamma, ngrad_Gamma;    
  double D = TDatabase::ParamDB->REACTOR_P12 / TDatabase::ParamDB->REACTOR_P16 ;

  
  TBaseCell *cell;
  TFEDesc2D *FeDesc;
  BaseFunct2D *BaseFuncts;
  TCollection *Coll;
  TJoint *joint;
  TIsoBoundEdge *isoboundedge;
  TBoundComp *BoundComp;
  BoundCond Cond0, Cond1;
  FE2D FEId;
  TFE2D *ele;
  RefTrans2D RefTrans;
  TRefTrans2D *F_K;
  TFESpace2D *fespace;
  BF2DRefElements RefElement;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1;  
  
#ifdef __SURFACT__ 
  TFESpace2D *surfactantspace;  
  FE2D TFEId;
  TFEDesc2D *TFeDesc;  
  
  double **Turef, **Tuxiref, **Tuetaref;  
  double Tuorig[MaxN_BaseFunctions2D], Tuxorig[MaxN_BaseFunctions2D];
  double Tuyorig[MaxN_BaseFunctions2D];
  double T_val[3], *S_Values;
  
  
  int  *TGlobalNumbers, *TBeginIndex, local_dof;
  int TN_BaseFunct, *TJointDOF, TN_DOF_Local, *TDOF;

// surfactant elasticity E
  double E = TDatabase::ParamDB->REACTOR_P10;
//Equation of state, 0 linear, 1 non-linear
  int EOS = int(TDatabase::ParamDB->REACTOR_P11);
//\Gamma_1/Gamma_\infty

//   double CHAR_L = TDatabase::ParamDB->CHAR_L0;
  
  surfactantspace = Surfact->GetFESpace2D();
  S_Values=Surfact->GetValues();
  TGlobalNumbers = surfactantspace->GetGlobalNumbers();
  TBeginIndex = surfactantspace->GetBeginIndex();  
  
//   Gamma_Max = 0;
//   N_Surf = Surfact->GetLength();

//   for(i=0;i<N_Surf;i++)
//    if(Gamma_Max<S_Values[i]) Gamma_Max=S_Values[i];
// 
//    OutPut("Gamma_Max " << Gamma_Max<<endl);  
  
#endif   
  
 
  
  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFuncts = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  fespace = A11->GetFESpace();
  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  RowPtr = A11->GetRowPtr();
  KCol = A11->GetKCol();

  ValuesA11 = A11->GetEntries();
  ValuesA22 = A22->GetEntries();

  double Re = TDatabase::ParamDB->RE_NR;
  double We = TDatabase::ParamDB->WB_NR, U;
  double Ca = We/Re, D_Angle;
  double beta = TDatabase::ParamDB->FRICTION_CONSTANT;

  double EQ_Angle = TDatabase::ParamDB->EQ_CONTACT_ANGLE;
  
  EQ_Angle = (3.141592654/180)*EQ_Angle;

 for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_Edges = cell->GetN_Edges();
    IJoint = 0;
    for(j=0;j<N_Edges;j++)
     {
      joint = cell->GetJoint(j);
      if(joint->GetType() == IsoBoundEdge)
      {
        isoboundedge = (TIsoBoundEdge *)joint;
        BoundComp = isoboundedge->GetBoundComp();
        isoboundedge->GetParameters(t0, t1);
        comp=BoundComp->GetID();
        BoundaryCondition(comp, t0, Cond0);
        BoundaryCondition(comp, t1, Cond1);

        if(Cond0 == FREESURF)
        {
          JointNumbers[IJoint] = j;
          IJoint++;
        }
      } // endif
     } // endfor j

    N_IsoJoints = IJoint;
    if(N_IsoJoints > 0)
    {
      FEId = fespace->GetFE2D(i, cell);
      
#ifdef __SURFACT__ 
      TFEId = surfactantspace->GetFE2D(i, cell);
#endif
       
      for(j=0;j<N_IsoJoints;j++)
      {
//      cout << "Cell " << i << " has free surface." << endl;
        IJoint = JointNumbers[j];
        // cout << "joint number: " << IJoint << endl;
        cell->GetVertex(IJoint)->GetCoords(x0, y0);
        cell->GetVertex((IJoint+1) % N_Edges)->GetCoords(x1, y1);
        //   if(y0==0||y1==0)
        //   cout<< " y0= " <<y0<<" y1= "<<y1<<"  x0= "<<x0<<"  x1= "<<x1<<endl;
        //   cout<< " N_LinePoints= " <<N_LinePoints<<endl;

     // entries for wetting DOF
      if(y0==0) // right wett point edge (bottom)
       {
        FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
        JointDOF = FeDesc->GetJointDOF(IJoint);
        N_DOF = FeDesc->GetN_JointDOF();
        for(m=0;m<N_DOF;m++)
         {
          DOF_R =  GlobalNumbers[BeginIndex[i]+JointDOF[m]];
          fespace->GetDOFPosition(DOF_R, x, y);
          if(y==0) // right wett point
          {
           U = Ucl[DOF_R];
          switch((int)TDatabase::ParamDB->CONTACT_ANGLE_TYPE)
          {
           case 0:
              D_Angle = EQ_Angle;
           break;

           case 1:
            if(U>0.)
	     { D_Angle = (3.141592654/180.)*TDatabase::ParamDB->AD_CONTACT_ANGLE; } // advanving angle
            else if(U<0.)
	     { D_Angle = (3.141592654/180.)*TDatabase::ParamDB->RE_CONTACT_ANGLE; }// receding angle
            else
             {
              D_Angle = EQ_Angle;
//                        + tanh(50.*U)*(TDatabase::ParamDB->AD_CONTACT_ANGLE
//                                        - TDatabase::ParamDB->RE_CONTACT_ANGLE);
             }
             
             
           break;

           case 2:
//   Hocking's expression
              D_Angle = pow(EQ_Angle, (double)3.0) 
                           + 9.0 * Ca * (fabs(U))* log(beta);

              D_Angle = pow(fabs(D_Angle), 1./3.);
           break;

           case 3:
              Ca *= fabs(U); // capillary number w.r.t contact line velocity
//   Jiang et al (1979) expression
              D_Angle = acos( cos(EQ_Angle) - (cos(EQ_Angle) + 1.)*tanh( 4.96*pow(Ca,0.702) )   );
           break;
           case 4:
              Ca *= fabs(U); // capillary number w.r.t contact line velocity
//   Bracke et al (1989) expression
              D_Angle = acos( cos(EQ_Angle) - (cos(EQ_Angle) + 1.)*( 2.*pow(Ca,0.5) )   );
           break;
           case 5:
              Ca *= fabs(U); // capillary number w.r.t contact line velocity
//   Berg et al (1992) expression
              D_Angle = acos( cos(EQ_Angle) - (cos(EQ_Angle) + 1.)*( 2.24*pow(Ca,0.54) )   );
           break;

           case 6:
//  Berg et al (1992) expression
              Ca *= fabs(U); // capillary number w.r.t contact line velocity
              D_Angle = acos( cos(EQ_Angle) - (cos(EQ_Angle) + 1.)*( 4.47*pow(Ca,0.42) )   );
           break;

          }
// OutPut("  x= "<< x <<"  y= "<< y << " U " << U<<  " D_Angle: " << (180./Pi)*D_Angle<< endl);
// exit(0);
          param[0] = x;
          param[1] = y;
          param[2] = U;
          r_axial = x;       // r value in the axial symmetric integral
          rhs1[DOF_R] +=  r_axial*((cos(D_Angle))/We);   break;
         }
        }
       }

      DOF = GlobalNumbers + BeginIndex[i];
      N_BaseFunct = N_BaseFuncts[FEId];
      ele = TFEDatabase2D::GetFE2D(FEId);
      RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);

      l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
      LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
      qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
      qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);
      TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(LineQuadFormula);
   
#ifdef __SURFACT__ 
       TFEDatabase2D::GetBaseFunct2DFromFE2D(TFEId)->MakeRefElementData(LineQuadFormula);
       TFeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(TFEId);
       TN_BaseFunct = N_BaseFuncts[TFEId];
       TJointDOF = TFeDesc->GetJointDOF(IJoint);
       TN_DOF_Local = TFeDesc->GetN_JointDOF();
       TDOF = TGlobalNumbers + TBeginIndex[i];
#endif

      switch(RefElement)
      {
        case BFUnitSquare:
          RefTrans = QuadIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TQuadIsoparametric *)F_K)->SetCell(cell);
          ((TQuadIsoparametric *)F_K)->GetOrigBoundFromRef(IJoint, N_LinePoints, zeta, X_B, Y_B);
        break;

        case BFUnitTriangle:
          RefTrans = TriaIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TTriaIsoparametric *)F_K)->SetCell(cell);
          ((TTriaIsoparametric *)F_K)->GetOrigBoundFromRef(IJoint, N_LinePoints, zeta, X_B, Y_B);

        break;
      } // endswitch


       uref = TFEDatabase2D::GetJointValues2D(BaseFuncts[FEId], LineQuadFormula, IJoint);
       uxiref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId], LineQuadFormula, IJoint, D10);
       uetaref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId], LineQuadFormula, IJoint, D01);

#ifdef __SURFACT__ 
       Turef = TFEDatabase2D::GetJointValues2D(BaseFuncts[TFEId], LineQuadFormula, IJoint);
//        Tuxiref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[TFEId], LineQuadFormula, IJoint, D10);
//        Tuetaref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[TFEId], LineQuadFormula, IJoint, D01);       
#endif
       for(k=0;k<N_LinePoints;k++)
        {
         switch(RefElement)
          {
            case BFUnitSquare:
              ((TQuadIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
// #ifdef __SURFACT__                         
//               ((TQuadIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
//                         TN_BaseFunct, Turef[k], Tuxiref[k], Tuetaref[k],
//                         Tuorig, Tuxorig, Tuyorig);
// #endif 
            break;

            case BFUnitTriangle:
              ((TTriaIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
                        
// #ifdef __SURFACT__                         
//               ((TTriaIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
//                         TN_BaseFunct, Turef[k], Tuxiref[k], Tuetaref[k],
//                         Tuorig, Tuxorig, Tuyorig);
// #endif                         
            break;
          } // endswitch

          // modify matrices
         F_K->GetTangent(IJoint, zeta[k], t0, t1);  // old line
         r_axial = fabs(X_B[k]);   // r value in the axial symmetric integral
         normn = sqrt(t0*t0+t1*t1);
         n0 =  t1/normn;
         n1 = -t0/normn;
         t1 /= normn;
         t0 /= normn;     

#ifdef __SURFACT__   
      for(mm=0;mm<TN_BaseFunct;mm++)
           Tuorig[mm] = Turef[k][mm];  
      
       T_val[0] = 0.;
//        T_val[1] = 0.; T_val[2] = 0.;
       
       for(l=0;l<TN_DOF_Local;l++)
        {
          // assumed that the velo space and 2D surfactant space are same fe space
          local_dof   = TJointDOF[l];
          m = TDOF[local_dof];

//           if(S_Values[m]<0)
// 	  {
// 	   OutPut(i<< "x : "<<X_B[k]<< " y: " << Y_B[k] <<"  Surfactant exceeds the reference value, T_val= " <<S_Values[m]<<endl);
// 	    S_Values[m] = 0.;
// 	  }
	  
	  
          val = S_Values[m];
          T_val[0] += val*Tuorig[local_dof];  // Surfactant C
        } // for(l=0;l<TN_

        if(T_val[0]<0. )
	 {
          OutPut(i<< "x : "<<X_B[k]<< " y: " << Y_B[k] <<"  Surfactant exceeds the reference value, T_val= " <<T_val[0]<<endl);
//        for(l=0;l<TN_DOF_Local;l++)
//         {
//           // assumed that the velo space and 2D surfactant space are same fe space
//           local_dof   = TJointDOF[l];
//           m = TDOF[local_dof];
//           val = S_Values[m];
// 	  
// 	   OutPut(l << " sval  "<<val<<endl);
// 	
// 	   
// 	}
          //numerical correction
          T_val[0]=0.; 
         }  

       // Marangoni effect in weak formulation
         if(EOS==0)
          {
           Gamma =(1. + E*(D - T_val[0]) ); 
           }
         else
          {
           Gamma =(1. + E*log(1. - T_val[0]) );
   
           //see SolubleSurf JCP paper
           if(Gamma<0.1)
            {
             Gamma = 0.1;
            }
          }              
#else
    Gamma = 1;
#endif       
          // Multiply with time step dt in the main program not here
          r = normn/We;
          for(l=0;l<N_BaseFunct;l++)
          {
           TestDOF = DOF[l];

           // updating rhs
            ngrad_test= n0*uxorig[l] + n1*uyorig[l];
            d1 = uxorig[l] - ngrad_test*n0;
            d2 = uyorig[l] - ngrad_test*n1;

// rhs1
//             val = r_axial*( (1.-n0*n0)*(d1 - GammaE1) - n0*n1*(d2 -GammaE2) );
	    
            val = r_axial*( (1.-n0*n0)*d1 - n0*n1*d2 );
            val += uorig[l]; // due to axialsymmetric
            val *= LineWeights[k]*r*Gamma;
            rhs1[TestDOF] -= val;

// rhs2
            val =  r_axial*( -n1*n0*d1 + (1.-n1*n1)*d2 );   
            val *= LineWeights[k]*r*Gamma;
            rhs2[TestDOF] -= val;
 
            index2 = RowPtr[TestDOF+1];
//               cout << TestDOF  << " RhsA  " << rhs1[TestDOF] << " RhsB  " << rhs2[TestDOF] << endl;
	    
            for(m=0;m<N_BaseFunct;m++)
            {
              AnsatzDOF = DOF[m];
              // cout << AnsatzDOF << " -- " << TestDOF << endl;
              index1 = RowPtr[TestDOF];
              if(index1+1 == index2) continue;
              while(KCol[index1] != AnsatzDOF) index1++;

              ngrad_ansatz= n0*uxorig[m] + n1*uyorig[m];
              e1 = uxorig[m] - ngrad_ansatz*n0;
              e2 = uyorig[m] - ngrad_ansatz*n1;

              val =d1*e1 + d2*e2 + (uorig[l]*uorig[m]/(r_axial*r_axial));
              val *= dt*LineWeights[k]*r*Gamma*r_axial;
              // cout << "A11: " << TestDOF << " ";
              // cout << AnsatzDOF << " " << val << endl;
              ValuesA11[index1] += val;

              val = d1*e1 + d2*e2;
              val *= dt*LineWeights[k]*r*Gamma*r_axial;

              // cout << "A22: " << TestDOF << " ";
              // cout << AnsatzDOF << " " << val << endl;
              ValuesA22[index1] += val;
            } // endfor m
          } // endfor l
        } // endfor k
      } // endfor j

    } // end (N_IsoJoints > 0)
  } // endfor i
 } //FreeSurf_axial3D_new(TSquareMatr

 


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
  double *ValuesUX, *ValuesUY, *Surfact, *surf_spl, *srhs, *Msurf, surf;
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
   
   Surf_DOF = new int[N_V];
   surf_spl = new double[N_V];  
   srhs = new double[N_Splines+1];   
   Msurf = new double[N_Splines+1];
   
   SurfSurf_DOF = new int[N_V];
   Surfsurf_spl = new double[N_V];  
   Surfsrhs = new double[N_Splines+1];   
   SurfMsurf = new double[N_Splines+1];  
  
   VeloBeginIndex = VelocitySpace->GetBeginIndex();
   VeloGlobalNumbers = VelocitySpace->GetGlobalNumbers();
   ValuesUX = Velocity->GetValues();
   ValuesUY = ValuesUX + Velocity->GetLength();
   
   SurfactSpace = Surfactant->GetFESpace2D();
   SurfactBeginIndex = SurfactSpace->GetBeginIndex();
   SurfactGlobalNumbers = SurfactSpace->GetGlobalNumbers();
   Surfact = Surfactant->GetValues(); 
   
   SurfSurfactSpace = SurfSurfactant->GetFESpace2D();
   SurfSurfactBeginIndex = SurfSurfactSpace->GetBeginIndex();
   SurfSurfactGlobalNumbers = SurfSurfactSpace->GetGlobalNumbers();
   SurfSurfact = SurfSurfactant->GetValues();    
// SurfSurfactant   
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

         Surf_DOF[m1] = SDOF[SJointDOF[i3]]; // needed for later update
         surf_spl[m1] = Surfact[SDOF[SJointDOF[i3]]];
 
         SurfSurf_DOF[m1] = SurfSDOF[SurfSJointDOF[i3]]; // needed for later update
         Surfsurf_spl[m1] = SurfSurfact[SurfSDOF[SurfSJointDOF[i3]]];	 
 
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

    FEParams[ISpline + 4  ] = -Msurf[i]*h[i+1]*h[i+1]/2. +
                          ((surf_spl[i+1]-surf_spl[i])/h[i+1]-h[i+1]/6.*(Msurf[i+1]-Msurf[i]))*h[i+1];
    FEParams[ISpline + 5] = Msurf[i+1]*h[i+1]*h[i+1]/2. +
                          ((surf_spl[i+1]-surf_spl[i])/h[i+1]-h[i+1]/6.*(Msurf[i+1]-Msurf[i]))*h[i+1];
  
    FEParams[ISpline + 6  ] = -SurfMsurf[i]*h[i+1]*h[i+1]/2. +
                          ((Surfsurf_spl[i+1]-Surfsurf_spl[i])/h[i+1]-h[i+1]/6.*(SurfMsurf[i+1]-SurfMsurf[i]))*h[i+1];
    FEParams[ISpline + 7] = SurfMsurf[i+1]*h[i+1]*h[i+1]/2. +
                          ((Surfsurf_spl[i+1]-Surfsurf_spl[i])/h[i+1]-h[i+1]/6.*(SurfMsurf[i+1]-SurfMsurf[i]))*h[i+1]; 

  
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
// int iso = 0;
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

    surf = surf_spl[FeDof]*phi1 + surf_spl[FeDof+1]*phi2 +
              FEParams[USpline+4]*phi3 + FEParams[USpline + 5]*phi4;


    Surfsurf = Surfsurf_spl[FeDof]*phi1 + Surfsurf_spl[FeDof+1]*phi2 +
              FEParams[USpline+6]*phi3 + FEParams[USpline + 7]*phi4;      
    
    if(j!=0) // endpoints no need to set
     {
      ValuesUX[U_DOF[m1]] = u0;
      ValuesUY[U_DOF[m1]] = u1;
      Surfact[Surf_DOF[m1]] = surf;
      SurfSurfact[SurfSurf_DOF[m1]] = Surfsurf;      
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

     IsoVertices[i3]->SetCoords(X, Y);

//       if(fabs(dx0-X)>1e-4 || fabs(dy0-Y)>1e-4)
//        cout<<"NewX iso:"<<' '<< iso++ <<' '<<X<<' '<< Y<<endl;
//  
     m++;

  // ====================================================================
  // for fe values
  // ====================================================================

    u0 = u1_spl[FeDof]*phi1 + u1_spl[FeDof+1]*phi2 +
              FEParams[USpline]*phi3 + FEParams[USpline + 1]*phi4;
    u1 = u2_spl[FeDof]*phi1 + u2_spl[FeDof+1]*phi2 +
              FEParams[USpline+2]*phi3 + FEParams[USpline + 3]*phi4;

    surf = surf_spl[FeDof]*phi1 + surf_spl[FeDof+1]*phi2 +
              FEParams[USpline+4]*phi3 + FEParams[USpline + 5]*phi4;

    Surfsurf = Surfsurf_spl[FeDof]*phi1 + Surfsurf_spl[FeDof+1]*phi2 +
              FEParams[USpline+6]*phi3 + FEParams[USpline + 7]*phi4;
  
    ValuesUX[U_DOF[m1]] = u0;
    ValuesUY[U_DOF[m1]] = u1;
    Surfact[Surf_DOF[m1]] = surf;    
    SurfSurfact[SurfSurf_DOF[m1]] = Surfsurf;        
    m1++;


// ====================================================================
    }  // for(i3=0;i3<k
    }   // if(k==ORDER-1)    

  }  // if(UpdateU) 
  else if (k==ORDER-1)
  {m += k; }
  
 }  //  for(j=0;j<N_E

   delete [] h; delete [] t; delete [] a; delete [] b;
   delete [] c; delete [] rhs; delete [] Mx; delete [] My;
   delete [] Params; delete [] Param9;  delete [] x; delete [] y;
  
 if(UpdateU)
  { 
   delete [] u1rhs;  delete [] u2rhs; delete [] u1_spl;
   delete [] u2_spl; delete [] Mu1; delete [] Mu2;
   delete [] U_DOF;  delete [] FEParams;
   delete [] Surf_DOF; delete [] surf_spl;  
   delete [] srhs; delete [] Msurf;
   delete [] SurfSurf_DOF; delete [] Surfsurf_spl;  
   delete [] Surfsrhs; delete [] SurfMsurf;  
   
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



