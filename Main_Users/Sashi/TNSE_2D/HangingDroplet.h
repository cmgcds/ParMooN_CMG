

#include <MacroCell.h>
#include <BoundEdge.h>
#include <IsoBoundEdge.h>
#include <IsoInterfaceJoint.h>
#include <BdLine.h>
#include <BdCircle.h>
#include <GridCell.h>
#include <TimeConvDiff2D.h>

extern "C"
{
  #include <gridgen.h>    
  void triangulate(char*, struct triangulateio*,
                   struct triangulateio*, struct triangulateio*);
}

void ExampleFile()
{  
// needed only if we consider the thickness of the pipe in meshing (needed for experimental comparison)
// #define __WITHCYLINDERTICKNESS__  
  
//   #define __IMPINGINGDROPLET__
  #define __AXIAL3D__
//   #define __EXPERIMENTAL__
    
//  #define __ENERGY__
 #define __SURFACT__
 #define __SOLUBLE__  
  
// #define __SURFTENSIONOLDFORM__
  
  TDatabase::ParamDB->Axial3D=1;
  OutPut("HangingDroplet.h" << endl) ;
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
  values[0] = 0.;
  values[0] = 0;  
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

void InitialS(double x, double y, double *values)
{
//  double t = TDatabase::TimeDB->CURRENTTIME;
//  double Re = TDatabase::ParamDB->RE_NR;
//  double Pr = TDatabase::ParamDB->Pr_NR;

//   static double eps = 1.0/(Re*Pr);

  double phi = atan2(y+(4.90056), x);
//   phi = Pi/2. + phi;
//      cout << " phi " << phi << endl;
//   values[0] = 0.5 +  0.5*cos(phi)*exp(-2.*t*eps); // diffusion test
//   values[0] = 0.5 + 0.5*cos(phi); // expanding static sphere
//   values[0] = 0.9*fabs(cos(phi)); // Oscllating sphere
//   values[0] = sqrt((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5));
  
//   values[0] =  TDatabase::ParamDB->REACTOR_P12*fabs(cos(phi)); 
  values[0] = TDatabase::ParamDB->REACTOR_P12;  // static bubble
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
//   values[0] = cos(phi); // Oscllating sphere
  values[0] = 0.;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void InitialSuract(double x, double y, double *values)
{
  values[0] = TDatabase::ParamDB->REACTOR_P19;
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
             cond = FREESURF;
             TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
        break;
        case 1:
           cond = DIRICHLET;
        break;
        case 2:   
#ifdef  __WITHCYLINDERTICKNESS__   
        case 3:
#endif  
           cond = DIRICHLET;
         break;       

#ifdef  __WITHCYLINDERTICKNESS__  
        case 4:   
#else  
        case 3:
#endif  
             cond = SLIP_FRICTION_PENETRATION_RESISTANCE;
             TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
             TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION = 1;
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
  
//   if(BdComp==2)
//   { 
//     double r = 1. - Param;
//     value =  1.22732*(r*r - 1.);
// //    cout << "  value " << value <<endl;
// 
//   }
//   else
  { value = 0; }
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
    { coeff[2] = -1./TDatabase::ParamDB->FR_NR;
      
//     cout<< " coeff[2] " << coeff[2] << endl; 
    
    } // f2 - => gravity in opposite direction
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


void FreeSurf_axial3D_new(TSquareMatrix2D *A11, TSquareMatrix2D *A22,
                          double *rhs1, double *rhs2,
                          BoundCondFunct2D *BoundaryCondition,
                          double dt, double *Ucl, double *param
#ifdef __SURFACT__  
                          , TFEFunction2D *Surfact
#endif  
  
)
{
  int i, j, k, l, DOF_R, DOF_L, m, N_Surf, mm;
  int *KCol, *RowPtr, *JointDOF, N_DOF;
  int N_LinePoints, N_Active;
  int *BeginIndex, *GlobalNumbers, *DOF, TestDOF, AnsatzDOF;
  int index1, index2;
  int count=0, count1=0, count2=0;
  int N_BaseFunct, *N_BaseFuncts, TN_BaseFunct, *TJointDOF, TN_DOF_Local, *TDOF;
  int JointNumbers[MAXN_JOINTS], IJoint, N_IsoJoints;  
  int comp, N_U, test_L=1, test_R=1;
  int N_Cells, N_Vertices, N_Edges, Semi_implicit=0;  
  int *TGlobalNumbers, *TBeginIndex, local_dof;
   
  double r2, r;  
  double  X_B[100], Y_B[100], r_axial, d1, d2, e1, e2, ngrad_test, ngrad_ansatz, tmp;
  double val, theta, factor1, factor2, angle;
  double *ValuesA11, *ValuesA12, *ValuesA21, *ValuesA22, *S_Values; 
  double *LineWeights, *zeta;
  double tx,ty,mod_t, x, y;
  double **uref, **uxiref, **uetaref;
  double uorig[MaxN_BaseFunctions2D], uxorig[MaxN_BaseFunctions2D];
  double uyorig[MaxN_BaseFunctions2D];
  double t0, t1, n0, n1, normn, line_wgt;    
  double Gamma_Max, Gamma_Infty, Gamma_Min;
  double **Turef, **Tuxiref, **Tuetaref;
  double Tuorig[MaxN_BaseFunctions2D], Tuxorig[MaxN_BaseFunctions2D];
  double Tuyorig[MaxN_BaseFunctions2D];
  double x0, x1, y0, y1, T_val[3], Gamma;
  double ngrad_Gamma, GammaE1, GammaE2, SurfForce1=0., SurfForce2=0., TGamm1, TGamm2;
  
  TBaseCell *cell;
  TFEDesc2D *FeDesc, *TFeDesc;
  BaseFunct2D *BaseFuncts;
  TCollection *Coll;
  TJoint *joint;
  TIsoBoundEdge *isoboundedge;
  TBoundComp *BoundComp;
  BoundCond Cond0, Cond1;
  FE2D FEId, TFEId;
  TFE2D *ele;
  RefTrans2D RefTrans;
  TRefTrans2D *F_K;
  TFESpace2D *fespace, *surfactantspace;
  BF2DRefElements RefElement;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1;  
  TIsoBoundEdge *isojoint;
  TVertex **Vertices;
     
  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFuncts = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  fespace = A11->GetFESpace();
  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();
  N_Active =  fespace->GetActiveBound();
  
  RowPtr = A11->GetRowPtr();
  KCol = A11->GetKCol();

  ValuesA11 = A11->GetEntries();
  ValuesA22 = A22->GetEntries();
  
#ifdef __SURFACT__  
  surfactantspace = Surfact->GetFESpace2D();
  S_Values=Surfact->GetValues();
  TGlobalNumbers = surfactantspace->GetGlobalNumbers();
  TBeginIndex = surfactantspace->GetBeginIndex();  
 
  Gamma_Max = 0; Gamma_Min=1;
  N_Surf = Surfact->GetLength();

  for(i=0;i<N_Surf;i++)
   {
    if(Gamma_Max<S_Values[i]) Gamma_Max=S_Values[i];
    if(Gamma_Min>S_Values[i]) Gamma_Min=S_Values[i];
   }
   
   OutPut("Gamma_Max, Gamma_Min: " <<  Gamma_Max<< " "  << Gamma_Min << endl);
  
// surfactant elasticity E
  double E = TDatabase::ParamDB->REACTOR_P10;
//Equation of state, 0 linear, 1 non-linear
  int EOS = int(TDatabase::ParamDB->REACTOR_P11);
//\Gamma_1/Gamma_\infty
  double D = TDatabase::ParamDB->REACTOR_P21 / TDatabase::ParamDB->REACTOR_P16 ;
  
#endif
  
  double Re = TDatabase::ParamDB->RE_NR;
  double We = TDatabase::ParamDB->WB_NR, U;
  double Ca = We/Re;
//   double beta = TDatabase::ParamDB->FRICTION_CONSTANT;

//   double EQ_Angle = TDatabase::ParamDB->EQ_CONTACT_ANGLE;
  
//   EQ_Angle = (3.141592654/180)*EQ_Angle;
 

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

      for(j=0;j<N_IsoJoints;j++)
      {
//      cout << "Cell " << i << " has free surface." << endl;
      IJoint = JointNumbers[j];
      joint = cell->GetJoint(IJoint);
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
       TFEId = surfactantspace->GetFE2D(i, cell);
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
#ifdef __SURFTENSIONOLDFORM__          
       Tuxiref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[TFEId], LineQuadFormula, IJoint, D10);
       Tuetaref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[TFEId], LineQuadFormula, IJoint, D01);
#endif       
#endif       

       for(k=0;k<N_LinePoints;k++)
        {
         switch(RefElement)
          {
            case BFUnitSquare:
              ((TQuadIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
#ifdef __SURFACT__
#ifdef __SURFTENSIONOLDFORM__         
              ((TQuadIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        TN_BaseFunct, Turef[k], Tuxiref[k], Tuetaref[k],
                        Tuorig, Tuxorig, Tuyorig);
            break;
#else      
      for(mm=0;mm<TN_BaseFunct;mm++)
           Tuorig[mm] = Turef[k][mm];
#endif      
#endif   
            break;

            case BFUnitTriangle:
              ((TTriaIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
#ifdef __SURFACT__     
#ifdef __SURFTENSIONOLDFORM__        
              ((TTriaIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        TN_BaseFunct, Turef[k], Tuxiref[k], Tuetaref[k],
                        Tuorig, Tuxorig, Tuyorig);
            break;
#else    
      for(mm=0;mm<TN_BaseFunct;mm++)
           Tuorig[mm] = Turef[k][mm]; 
#endif    
#endif
            break;
          } // endswitch

          // modify matrices
         F_K->GetTangent(IJoint, zeta[k], t0, t1);  // old line
         r_axial = fabs(X_B[k]);   // r value in the axial symmetric integral
         if(X_B[k]<=0)
          {
           cell->GetVertex(IJoint)->GetCoords(x0, y0);       
           cell->GetVertex((IJoint+1) % N_Edges)->GetCoords(x1, y1);    
           cout <<"X_B[k] negative in freesurface int, change Quad rule " <<  X_B[k] <<endl;
           cout << IJoint << " " << k <<" X  " << x0<<" Y  " << y0<<" X  " << x1<<" Y  " << y1 <<endl;
           cout << "If the quad formula is correct, then the time step might be large !!!! " <<endl; 
	   
	   isojoint = (TIsoBoundEdge *)joint;
           Vertices = isojoint->GetVertices();
  	   Vertices[0]->GetCoords(x0, y0);    
	   cout <<"Iso Vertex " <<" X  " << x0<<" Y  " << y0<<endl;
           exit(0);
          }
          
         normn = sqrt(t0*t0+t1*t1);
         t1 /= normn;
         t0 /= normn;    	 
         n0 =  t1;
         n1 = -t0;
         r = normn/We;
//            cout << " t1  "<< t1*t0 << endl;
#ifdef __SURFACT__   
       T_val[0] = 0.; T_val[1] = 0.; T_val[2] = 0.;
       for(l=0;l<TN_DOF_Local;l++)
        {
          // assumed that the velo space and 2D surfactant space are same fe space
          local_dof   = TJointDOF[l];
          m = TDOF[local_dof];
          val = S_Values[m];

          T_val[0] += val*Tuorig[local_dof];  // Surfactant Gamma
#ifdef __SURFTENSIONOLDFORM__             
          T_val[1] += val*Tuxorig[local_dof];  // Gamma_x
          T_val[2] += val*Tuyorig[local_dof];  // Gamma_y
#endif
        } // for(l=0;l<TN_ 
        
//           cout <<" T:  " << T_val[0]<<" Tx:  " << T_val[1] <<" Ty:  " << T_val[2] <<endl;
	  
        if(T_val[0]<0.)
	 {
          OutPut(i<< "x : "<<X_B[k]<< " y: " << Y_B[k] <<"  Surfactant exceeds the reference value, T_val= " <<T_val[0]<<endl);
          //numerical correction
          T_val[0]=0.; 
// 	  exit(0);
         }    
        
        
#ifdef __SURFTENSIONOLDFORM__  
         ngrad_Gamma = n0*T_val[1] + n1*T_val[2];
         GammaE1 = (T_val[1] - ngrad_Gamma*n0);
         GammaE2 = (T_val[2] - ngrad_Gamma*n1);
#endif  
         if(EOS==0)
          {
 #ifdef __SURFTENSIONOLDFORM__  
           GammaE1 *=E;
           GammaE2 *=E;
#endif     
           Gamma =(1. + E*(D - T_val[0]) ); 
          }
         else
          {
 #ifdef __SURFTENSIONOLDFORM__  	    
           GammaE1 *=(E/(1.- T_val[0]));
           GammaE2 *=(E/(1.- T_val[0]));    
#endif  
           Gamma =(1. + E*log(1. - T_val[0]) );

           //see SolubleSurf JCP paper
//            if(Gamma<0.1)
//             {
//              Gamma = 0.1;
//              T_Marangoni = 0.;
//             }
          }                 
#else
//     D = 0.;
    Gamma = 1;
//     T_Marangoni = 0.;          
        
#endif 
          // Multiply with time step dt in the main program not here
          for(l=0;l<N_BaseFunct;l++)
           {
           TestDOF = DOF[l];
           
           //wetting point
           if(TestDOF>=N_Active)
             continue;
             
           // updating rhs
            ngrad_test= n0*uxorig[l] + n1*uyorig[l];
            d1 = uxorig[l] - ngrad_test*n0;
            d2 = uyorig[l] - ngrad_test*n1;
#ifdef __SURFACT__   
            d1 *= Gamma;
            d2 *= Gamma;
#endif

#ifdef __SURFTENSIONOLDFORM__  
             TGamm1 = GammaE1*uorig[l];
             TGamm2 = GammaE2*uorig[l];

            val = r_axial*( (1-n0*n0)*(d1 - TGamm1) - n0*n1*(d2 -TGamm2));

#else    
            val = r_axial*( (1-n0*n0)*d1 - n0*n1*d2 );
	    
// 	    alternative form
// 	    val = r_axial*d1; ???????????  
#endif    
            val +=(Gamma*uorig[l]);
            val *= LineWeights[k]*r;
// val = r_axial*LineWeights[k]*r*1.*uorig[l]*n0; // exact curvature for bubble
            rhs1[TestDOF] -= val;
	    
#ifdef __SURFTENSIONOLDFORM__  	    
        //additional term only
//          SurfForce1 -= r*r_axial*LineWeights[k]*( (1.-n0*n0)*(- TGamm1) - n0*n1*(-TGamm2));;
//     Marangoni convection
         val  = r*r_axial*LineWeights[k]*TGamm1;
         rhs1[TestDOF] -= val;
//          SurfForce1-= val;  
#endif
   
// 	     rhs 2
#ifdef __SURFTENSIONOLDFORM__
            val =  r_axial*( -n1*n0*(d1 - TGamm1) + (1-n1*n1)*(d2 -TGamm2) );
#else    
            val =  r_axial*( -n1*n0*d1 + (1-n1*n1)*d2 );
//         alternative form
//         val = r_axial*d2;  ???????????  
#endif    
            val *= LineWeights[k]*r;
// val = r_axial*LineWeights[k]*r*1.*uorig[l]*n1; // exact curvature for bubble
            rhs2[TestDOF] -= val;

//     Marangoni convection
#ifdef __SURFTENSIONOLDFORM__ 
// 	    SurfForce2 -= r*r_axial*LineWeights[k]*( -n1*n0*( - TGamm1) + (1-n1*n1)*( -TGamm2) ); 

            val = r*r_axial*LineWeights[k]*TGamm2;
            rhs2[TestDOF] -= val;
// 	    SurfForce2 -= val;      
#endif
            index2 = RowPtr[TestDOF+1];

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

#ifdef __SURFTENSIONOLDFORM__  
             val =(d1 - TGamm1)*e1 + (d2 -TGamm2)*e2
                   + Gamma*(uorig[l]*uorig[m]/(r_axial*r_axial));
#else
              val =d1*e1 + d2*e2 + Gamma*(uorig[l]*uorig[m]/(r_axial*r_axial));
#endif      
              val *= dt*LineWeights[k]*r*r_axial;
              // cout << "A11: " << TestDOF << " ";
              // cout << AnsatzDOF << " " << val << endl;
              ValuesA11[index1] += val;
 
#ifdef __SURFTENSIONOLDFORM__  
              val = (d1 - TGamm1)*e1 + (d2 -TGamm2)*e2;
#else
              val = d1*e1 + d2*e2;
#endif      
              val *= dt*LineWeights[k]*r*r_axial;

              // cout << "A22: " << TestDOF << " ";
              // cout << AnsatzDOF << " " << val << endl;
              ValuesA22[index1] += val;
            } // endfor m
          } // endfor l
        } // endfor k
      } // endfor j

    } // end (N_IsoJoints > 0)
  } // endfor i
  
// #ifdef  __SURFTENSIONOLDFORM__  
//   OutPut("SURFTENSION OLD FORM  SurfForceA " << SurfForce1 << "  SurfForceB " << SurfForce2 << endl);
// #else
//   OutPut("SURFTENSION NEW FORM  SurfForceA " << SurfForce1 << "  SurfForceB " << SurfForce2 << endl);
// #endif  
//   exit(0);
  
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


void SolveGridEquation(double **Entries, double *sol, double *rhs,
                       int *KCol, int *RowPtr, int N_DOF)
{
  int i,j,k, col, Diognal;
  double *Entries11, *Entries12, *Entries21, *Entries22;
  double sum1, sum2, max_sum1, max_sum2;
  int start, end;

  double max_error, error=1.e-12;
  int iter;

  Entries11 = Entries[0];
  Entries12 = Entries[1];
  Entries21 = Entries[2];
  Entries22 = Entries[3];
  
  max_error = 1.; iter = 0;
  while(max_error>error)
  {
    max_error = 0.0; iter++;
    for(i=0;i<N_DOF;i++)
    {
      start = RowPtr[i];
      end = RowPtr[i+1];
      sum1 = rhs[i];
      sum2 = rhs[i+N_DOF];
      for(k=start;k<end;k++)
      {
        col = KCol[k];
        if (col==i) Diognal = k;
        sum1 -= Entries11[k] * sol[col]
              +Entries12[k] * sol[col+N_DOF];
        sum2 -= Entries21[k] * sol[col]
              +Entries22[k] * sol[col+N_DOF];
      } // endfor k
     // sol[i] += sum1/Entries11[start];
     // sol[i+N_DOF] += sum2/Entries22[start];
        sol[i] += sum1/Entries11[Diognal];
        sol[i+N_DOF] += sum2/Entries22[Diognal];
      if(max_error<fabs(sum1/Entries11[Diognal])) max_error = fabs(sum1/Entries11[Diognal]);
      if(max_error<fabs(sum2/Entries22[Diognal])) max_error = fabs(sum2/Entries22[Diognal]);
    } // endfor i
    if(iter == 1000) break;
  } // end while
//OutPut("Grid Solver: Number iteration "<<iter<<endl);
}

void ReParam_Coord(int &N_V, double *&x, double *&y)
{
  int i, j, ISpline, N_Splines, m, m1, k, i3;

  double *h, *t;
  double *a, *b, *c, teta;
  double *rhs,  *Mx, *My, *Params, *Param9;
  double phi1, phi2, phi3, phi4, X, Y, T;
  double dx0, dy0, dx1, dy1, x0, y0, xn, yn;

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
  
  
   x0 = x[0];
   y0 = y[0];
   xn = x[N_V-1]; 
   yn = y[N_V-1];    
   
   N_V = (int)TDatabase::ParamDB->P6; 

   delete [] x;
   delete [] y;
   
   x = new double[N_V];
   y = new double[N_V];  
   
   x[0] = x0;
   y[0] = y0;
   x[N_V-1] = xn;
   y[N_V-1] = yn;
   
   teta = 1.0/((double)N_V-1.);
   T = 0;

   Param9[0] = 0;
   for(i=1;i<=N_Splines;i++) 
    Param9[i] = Params[(i-1)*10+8];

   m = 0;
   m1 = 0;   
   // no need to set end vert
   for(j=0;j<N_V-1;j++)
    {
   
     T = (double)m*teta;
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

   X = Params[ISpline    ]*phi1 + Params[ISpline + 2]*phi2 +
       Params[ISpline + 4]*phi3 + Params[ISpline + 6]*phi4;
   Y = Params[ISpline + 1]*phi1 + Params[ISpline + 3]*phi2 +
       Params[ISpline + 5]*phi3 + Params[ISpline + 7]*phi4;

   if(m)
    {
     x[m] = X;
     y[m] = Y;
    }
//     OutPut("NewX:"<<' '<< m <<' '<<x[m]<<' '<< y[m]<<endl);
    m++;
 }  //  for(j=0;j<N_E
 
       
//  cout<< "ReParam_Coord " << endl;
//  exit(0);
   delete [] h; delete [] t; delete [] a; delete [] b;
   delete [] c; delete [] rhs; delete [] Mx; delete [] My;
   delete [] Params; delete [] Param9;
  
 
}
 
 

// **************************************************************************
// Triangular Mesh Generation
// **************************************************************************
void TriaReMeshGen_hori(TDomain *&Domain, int &N_FreeBound_Vert, double *&S_BX, double *&S_BY)
{
  TBoundComp *BoundComp;
  TBdLine *UpdateAxialBound, *UpdatelBoundTop, *UpdateBoundWall, *UpdateHorizBoundWall;
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
  double y_begin, y_end, dx, deviation, y_AxialEnd;
  double CHAR_L=TDatabase::ParamDB->CHAR_L0, Wall_R;
  
  int i, j, k, ID, In_Index, N_Cells, N_G, *PointNeighb, maxEpV=0;
  int a, b, len1, len2, Neighb_tmp, mode = int (TDatabase::ParamDB->P4);
  int N_Interf_Vertices, CurrComp, CurrVertex, N_Joints, N_Vertices;
  int N_RootCells, *PartMarker, *Triangles, Neib[2], CurrNeib;
  int N_SlipBound_Vert, N_TopBound_Vert, N_HorizBound_Vert, N_WallBound_Vert;

  boolean AllowEdgeRef = (boolean) TDatabase::ParamDB->MESHGEN_ALLOW_EDGE_REF;
  
  char *SMESH, line[100];
  SMESH = TDatabase::ParamDB->SMESHFILE;
  
  struct triangulateio In, Out;
  std::ostringstream opts;
  opts << " ";

  
  std::ifstream dat(SMESH);

  if (!dat)
  {
    cerr << "cannot open '" << SMESH << "' for input" << endl;
    exit(-1);
  }    
  
   // check the dimension
  N_G = -1;
  while (!dat.eof())
  {
    dat >> line;
    N_G++; 
    
       // read until end of line
    dat.getline (line, 99);    
  } 
    
  dat.close();
  
  cout<< " TriaReMeshGen " << N_G << endl;
  delete [] S_BX;
  delete [] S_BY;  
  
  S_BX = new double[N_G];
  S_BY = new double[N_G];
  
  std::ifstream dat1(SMESH);
  
  // check the dimension

  //read from file
   for(i=0;i<N_G; i++)
    {
//         dat1 >> S_BX[i] >> S_BY[i];  
     dat1 >> x >> y; 
//      S_BX[i] = x/3.5951829e-01;
//      S_BY[i] = y/3.5951829e-01;
     
     S_BX[i] = x/CHAR_L; // scaling with inner radius of the tube
     S_BY[i] = y/CHAR_L; // scaling with inner radius of the tube
     
     dat1.getline (line, 99);
//      cout<< i << " vert X: " <<S_BX[i]  << " vert Y: " << S_BY[i] <<endl;   
    }
  dat1.close();  
  
  ReParam_Coord(N_G, S_BX, S_BY);
 
  N_FreeBound_Vert = N_G-1;  //excluding last vert
  
//    for(i=0;i<N_G; i++)
//     cout<< i << " vert X: " <<S_BX[i]  << " vert Y: " << S_BY[i] <<endl;   
  
//   exit(0);
  
  BoundPart = Domain->GetBdPart(0);
  
#ifdef  __WITHCYLINDERTICKNESS__
  UpdateFreeBound = (TBdCircle*)BoundPart->GetBdComp(0);
  UpdateHorizBoundWall = (TBdLine*)BoundPart->GetBdComp(1); 
  UpdateBoundWall = (TBdLine*)BoundPart->GetBdComp(2);
  UpdatelBoundTop = (TBdLine*)BoundPart->GetBdComp(3);    
  UpdateAxialBound = (TBdLine*)BoundPart->GetBdComp(4);
#else
  UpdateFreeBound = (TBdCircle*)BoundPart->GetBdComp(0);
  UpdateBoundWall = (TBdLine*)BoundPart->GetBdComp(1);
  UpdatelBoundTop = (TBdLine*)BoundPart->GetBdComp(2);    
  UpdateAxialBound = (TBdLine*)BoundPart->GetBdComp(3); 
#endif    
  
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

//   N_FreeBound_Vert = int (TDatabase::ParamDB->P6);    //Freesurf except end point
  TDatabase::ParamDB->P6 = N_FreeBound_Vert;
  N_TopBound_Vert = 10;
#ifdef  __WITHCYLINDERTICKNESS__    
  N_HorizBound_Vert = 5;
#else
  N_HorizBound_Vert = 0;// needed only if we consider the thickness of the pipe in meshing (needed for experimental comparison)
#endif  
  N_SlipBound_Vert = 25;      
  N_WallBound_Vert = 10;  
  N_Interf_Vertices = N_FreeBound_Vert+N_WallBound_Vert+N_SlipBound_Vert + N_TopBound_Vert + N_HorizBound_Vert;
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
  
  for(i=0;i<N_FreeBound_Vert;i++) // without last point
   { 
    In.pointlist[2*In_Index] = S_BX[i];
    In.pointlist[2*In_Index+1] = S_BY[i];
    cout<<" x : "<< S_BX[i] << " y : "<< S_BY[i]<<endl;
    //In.pointmarkerlist[In_Index] = CurrComp;
    In.segmentlist[2*In_Index] = In_Index;
    In.segmentlist[2*In_Index+1] = In_Index+1;
    if (AllowEdgeRef)
     { In.segmentmarkerlist[In_Index] = CurrComp; }
    else
     { In.segmentmarkerlist[In_Index] = 100000 + CurrComp; }
    In_Index++;  
   }
  CurrComp++;   
    
  cout<< endl;
  
  /** horizontal BD (marker=1)*/
   Wall_R = 1./CHAR_L;
#ifdef  __WITHCYLINDERTICKNESS__   
   x = S_BX[N_FreeBound_Vert]; // last point from data file
   dx = (Wall_R  - x)/(double)N_HorizBound_Vert;
   y=0.;
   for(i=0;i<N_HorizBound_Vert;i++) // without last point
    {
     In.pointlist[2*In_Index] = x;
     In.pointlist[2*In_Index+1] = y;
     cout<<" x : "<< x << " y : "<< y<<endl;
     In.pointmarkerlist[In_Index] = CurrComp;
     In.segmentlist[2*In_Index] = In_Index;
     In.segmentlist[2*In_Index+1] = In_Index+1;
     In.segmentmarkerlist[In_Index] = CurrComp;
     In_Index++;
     x += dx;
    }
 
   x = S_BX[N_FreeBound_Vert]; // last point from data file
   UpdateHorizBoundWall->SetParams(x, 0, (Wall_R - x), 0);     
   CurrComp++;
   cout << endl;    
#endif    
  
    TDatabase::ParamDB->P7 =  S_BX[N_FreeBound_Vert]; // last point from data file
    x= Wall_R;
    y_begin = y = 0;
    y_end = TDatabase::ParamDB->P8/CHAR_L;  
    dx = y_end/(double)N_WallBound_Vert; 

   // points and segments on the vert Wall boundary (marker=2)
   for(i=0;i<N_WallBound_Vert;i++) // without last point
    {
     In.pointlist[2*In_Index] = x;
     In.pointlist[2*In_Index+1] = y;
     cout<<" x : "<< x << " y : "<< y<<endl;
     In.pointmarkerlist[In_Index] = CurrComp;
     In.segmentlist[2*In_Index] = In_Index;
     In.segmentlist[2*In_Index+1] = In_Index+1;
     In.segmentmarkerlist[In_Index] = CurrComp;
     In_Index++;
     y = y_begin + double(i+1)*dx;
    }
   CurrComp++;
   cout << endl;
  
    x= Wall_R;
    y= TDatabase::ParamDB->P8/CHAR_L;  
    dx = x/(double)N_TopBound_Vert; 

   UpdateBoundWall->SetParams(x, y_begin, 0, y - y_begin);   
    
   // points and segments on the horizontal boundary (marker=1)
   for(i=0;i<N_TopBound_Vert;i++) // without last point
    {
     In.pointlist[2*In_Index] = x;
     In.pointlist[2*In_Index+1] = y;
     cout<<" x : "<< x << " y : "<< y<<endl;
     In.pointmarkerlist[In_Index] = CurrComp;
     In.segmentlist[2*In_Index] = In_Index;
     In.segmentlist[2*In_Index+1] = In_Index+1;
     In.segmentmarkerlist[In_Index] = CurrComp;
     In_Index++;
     x -= dx;
    }
   CurrComp++;
  cout << endl;
  
   UpdatelBoundTop->SetParams(Wall_R, y,  -y, 0.);  
 
    x= 0.0;
    y_begin = y = TDatabase::ParamDB->P8/CHAR_L;  
    y_AxialEnd = S_BY[0];  
    dx = (y_AxialEnd-y_begin)/(double)N_SlipBound_Vert; 
//     len1 = N_SlipBound_Vert/2;
//  
//    // points and segments on the horizontal boundary (marker=1)
//    for(i=0;i<len1;i++) // without last point
//     {
//      In.pointlist[2*In_Index] = x;
//      In.pointlist[2*In_Index+1] = y;
//      cout<<" x : "<< x << " y : "<< y<<endl;
//      In.pointmarkerlist[In_Index] = CurrComp;
//      In.segmentlist[2*In_Index] = In_Index;
//      In.segmentlist[2*In_Index+1] = In_Index+1;
//      In.segmentmarkerlist[In_Index] = CurrComp;
//      In_Index++;
//      y = y_begin + double(i+1)*dx;
//     }
//    CurrComp++; 
//   cout << endl;
    
   double BD_Y = y;
   // points and segments on the horizontal boundary (marker=1)
   for(i=0;i<N_SlipBound_Vert;i++) // without last point
    {
     In.pointlist[2*In_Index] = x;
     In.pointlist[2*In_Index+1] = y;
     cout<<" x : "<< x << " y : "<< y<<endl;
     In.pointmarkerlist[In_Index] = CurrComp;
     In.segmentlist[2*In_Index] = In_Index;
     In.segmentlist[2*In_Index+1] = In_Index+1;
     In.segmentmarkerlist[In_Index] = CurrComp;
     In_Index++;
     y = y_begin + double(i+1)*dx;
    }
   CurrComp++;    
   
//   cout<< " TriaReMeshGen " << N_G << endl;
//   exit(0);
//   cout << endl;

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

  // Axial Bound startx, starty, x length and y length
  UpdateAxialBound->SetParams(0., BD_Y, 0, y_AxialEnd - BD_Y);
//   UpdateAxialBoundTop->SetParams(0., y_begin, 0, BD_Y - y_begin);
  
  
// Free boundary xmid, ymid, radius_a, radius_b, start angle, end angle
//  UpdateFreeBound ->SetParams( 0.0, 0.0, T_a, T_b, -Pi/2., Pi/2.);
  
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

//   cout<< "TriaReMeshGen: " << endl;
//   exit(0); 
} // TriaReMeshGen_hori(
  

void ReParam_axial3D_Data(int &N_E, TBaseCell **cell, int *EdgeNo,  int *CellNo, 
                          TFEVectFunct2D *Velocity,  double *&Intpol_Coord, double *&Intpol_VeloValues,
                          double h_min, double **&FreePts)
{
  int i, j, ISpline, N_Splines, N_V, ORDER, VSP, m, m1, k, i3, USpline, FeDof;
  double *h, *t, u0, u1, u2;
  double *a, *b, *c, *x, *y, teta;
  double *rhs, *u1rhs, *u2rhs, *Mx, *My,*Mu1, *Mu2, *Params, *Param9, *FEParams;
  double phi1, phi2, phi3, phi4, X, Y, T;
  double dx0, dy0, dx1, dy1, *u1_spl, *u2_spl;
  TIsoBoundEdge *isojoint;
  TVertex **IsoVertices;
  TJoint *Joint;
  TBaseCell *Me;
  TFESpace2D *VelocitySpace;
  int *VeloBeginIndex, *VeloGlobalNumbers, *JointDOF, *DOF, N_DOF_Joint, *U_DOF;
//   int *SurfactBeginIndex, *SurfactGlobalNumbers, *SJointDOF, *SDOF, SN_DOF_Joint,*Surf_DOF;
  double *ValuesUX, *ValuesUY;
  FE2D FEId;
  TFE2D *ele;
  TFEDesc2D *FeDesc;
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
  u1rhs = new double[N_Splines+1];
  u2rhs = new double[N_Splines+1];
//   srhs = new double[N_Splines+1];
  u1_spl = new double[N_Splines+1];
  u2_spl = new double[N_Splines+1];
//   surf_spl = new double[N_Splines+1];
  Mu1 = new double[N_Splines+1];
  Mu2 = new double[N_Splines+1];
//   Msurf = new double[N_Splines+1];
  Mx = new double[N_Splines+1];
  My = new double[N_Splines+1];
  Params = new double [10*N_Splines];
  Param9 = new double [N_Splines+1];
  FEParams = new double [3*2*N_Splines]; // 3 fe functions, u1, u2, surfact

  x = new double[N_V];
  y = new double[N_V];
  U_DOF = new int[N_V];
//   Surf_DOF = new int[N_V];

  VelocitySpace = Velocity->GetFESpace2D();
  VeloBeginIndex = VelocitySpace->GetBeginIndex();
  VeloGlobalNumbers = VelocitySpace->GetGlobalNumbers();
  ValuesUX = Velocity->GetValues();
  ValuesUY = ValuesUX + Velocity->GetLength();


//   SurfactSpace = Surfactant->GetFESpace2D();
//   SurfactBeginIndex = SurfactSpace->GetBeginIndex();
//   SurfactGlobalNumbers = SurfactSpace->GetGlobalNumbers();
//   Surfact = Surfactant->GetValues();

  coll = VelocitySpace->GetCollection();

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
        // cout<< i<<" FreeGaus " << (180/Pi)*atan2(y[m], x[m]) <<endl;
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
    FEId = VelocitySpace->GetFE2D(CellNo[i], Me);
    ele = TFEDatabase2D::GetFE2D(FEId);
    FeDesc = ele->GetFEDesc2D();   // fe descriptor
    JointDOF = FeDesc->GetJointDOF(EdgeNo[i]);
    N_DOF_Joint = FeDesc->GetN_JointDOF();
    DOF = VeloGlobalNumbers + VeloBeginIndex[CellNo[i]];

//     // for surfactant
//     SFEId = SurfactSpace->GetFE2D(CellNo[i], Me);
//     Sele = TFEDatabase2D::GetFE2D(SFEId);
//     SFeDesc = Sele->GetFEDesc2D();   // fe descriptor
//     SJointDOF = SFeDesc->GetJointDOF(EdgeNo[i]);
//     SN_DOF_Joint = SFeDesc->GetN_JointDOF();
//     SDOF = SurfactGlobalNumbers + SurfactBeginIndex[CellNo[i]];


    if((N_DOF_Joint-1)!=ORDER)
     {
      // only second order conforming elements implimented
      cout<< " No match in no. velo-nodal functunals and vertices on the free surface edge  "<<endl;
      cout<< " (N_DOF_Joint-1) " << N_DOF_Joint-1 << " ORDER " << ORDER <<endl;
      exit(0);
     }

    if(i != N_E-1)// -1 due to end dof will be the start dof of the next edge except on last edge
     N_DOF_Joint--; // assumed that velocity and surfactant having same no. of dof on edge

    //   cout << " CellNo[i] " << CellNo[i] << endl;
     for (i3=0;i3<N_DOF_Joint;i3++)
      {
       U_DOF[m1] = DOF[JointDOF[i3]]; // needed for later update
       u1_spl[m1] = ValuesUX[DOF[JointDOF[i3]]];
       u2_spl[m1] = ValuesUY[DOF[JointDOF[i3]]];
//        Surf_DOF[m1] = SDOF[SJointDOF[i3]]; // needed for later update
//        surf_spl[m1] = Surfact[SDOF[SJointDOF[i3]]];
       //  cout << "  SJointDOf " << JointDOF[i3] << " DOF " << DOF[JointDOF[i3]] <<  endl;
       m1++;
      }
   } // for(i=0;i<N_E

// exit(0);

//   end vertex of the freeboundary
  k = cell[N_E-1]->GetN_Edges();
  cell[N_E-1]->GetVertex((EdgeNo[N_E-1]+1) % k)->GetCoords(x[m], y[m]);

   if(m+1!=m1)
    {
     // only second order conforming elements implimented
     cout<< " No match in no. velo-nodal functunals and vertices on the free surface edge  "<<endl;
     cout<< " m " << m << " m1 " << m1 <<" N_Splines " << N_Splines <<endl;
     exit(0);
    }
//           cout<< " m " << m << " m1 " << m1 <<" N_Splines " << N_Splines <<endl;
// exit(0);

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

// ******************************************************************
// u1 component
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


// u2 component
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

// // surfactant
//   for(i=1;i<N_Splines;i++)
//    {
//      u0 = surf_spl[i-1];
//      u1 = surf_spl[i];
//      u2 = surf_spl[i+1];
// 
//      srhs[i] = 6./(h[i]+h[i+1])*((u2-u1)/h[i+1]-(u1-u0)/h[i]);
//     }

//    srhs[0] = srhs[1];
//    srhs[N_Splines] = srhs[N_Splines-1];
// 
//   Solver_3dia(N_Splines, a, b, c, srhs, Msurf);
// 


  for(i=0;i<N_Splines;i++)
  {
    ISpline = i*6;

    FEParams[ISpline  ] = -Mu1[i]*h[i+1]*h[i+1]/2. +
                          ((u1_spl[i+1]-u1_spl[i])/h[i+1]-h[i+1]/6.*(Mu1[i+1]-Mu1[i]))*h[i+1];
    FEParams[ISpline + 1] = Mu1[i+1]*h[i+1]*h[i+1]/2. +
                          ((u1_spl[i+1]-u1_spl[i])/h[i+1]-h[i+1]/6.*(Mu1[i+1]-Mu1[i]))*h[i+1];


    FEParams[ISpline + 2  ] = -Mu2[i]*h[i+1]*h[i+1]/2. +
                          ((u2_spl[i+1]-u2_spl[i])/h[i+1]-h[i+1]/6.*(Mu2[i+1]-Mu2[i]))*h[i+1];
    FEParams[ISpline + 3] = Mu2[i+1]*h[i+1]*h[i+1]/2. +
                          ((u2_spl[i+1]-u2_spl[i])/h[i+1]-h[i+1]/6.*(Mu2[i+1]-Mu2[i]))*h[i+1];


//     FEParams[ISpline + 4  ] = -Msurf[i]*h[i+1]*h[i+1]/2. +
//                           ((surf_spl[i+1]-surf_spl[i])/h[i+1]-h[i+1]/6.*(Msurf[i+1]-Msurf[i]))*h[i+1];
//     FEParams[ISpline + 5] = Msurf[i+1]*h[i+1]*h[i+1]/2. +
//                           ((surf_spl[i+1]-surf_spl[i])/h[i+1]-h[i+1]/6.*(Msurf[i+1]-Msurf[i]))*h[i+1];

  }

// =====================================================================================
// spline construction done
// =====================================================================================
   delete [] FreePts[0];
   delete [] FreePts[1];
   
   N_E = int(t[N_Splines]/h_min);
 
   FreePts[0] = new double[N_E+1];
   FreePts[1] = new double[N_E+1];

   Intpol_Coord = new double[4*(N_E+1)];  
//    Intpol_Values = new double[2*(N_E+1)];   
   Intpol_VeloValues = new double[4*(N_E+1)];
   
   FreePts[0][0] = x[0];
   FreePts[1][0] = y[0];   
   Intpol_Coord[0] = x[0];
   Intpol_Coord[1] = y[0];

   Intpol_VeloValues[0] = u1_spl[0];
   Intpol_VeloValues[1] = u2_spl[0];    
//    Intpol_Values[0] = surf_spl[0];
   
   teta = 1.0/(double)(2.*N_E); 
   
   //cout<< "X " <<FreePts[0][0] <<" Y " <<FreePts[1][0] <<endl; 
   
//    m=0;
//    for(j=0;j<N_E;j++)
//     {
//      T = (double)m*teta;   
//      cout<< "T " <<T << endl;  
//      m++;
//      
//      T = (double)m*teta;   
//      cout<< "T " <<T << endl;   
//      m++;
//     }
    
    
//    cout<< "Spline " <<endl;
//    exit(0);
   
   T = 0;
   Param9[0] = 0;
   for(i=1;i<=N_Splines;i++) 
    Param9[i] = Params[(i-1)*10+8];

   m = 0;
   m1 = 0;
 
   for(j=0;j<N_E;j++)
    {
     T = (double)m*teta;
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
      } // for(i=1;i<=N_Splines;i++)

   phi1 = (2.*T*T - 3.*T)*T + 1.;
   phi2 = (-2.*T + 3.)*T*T;
   phi3 = (T*T - 2.*T + 1.)*T;
   phi4 = (T - 1)*T*T;

   X = Params[ISpline    ]*phi1 + Params[ISpline + 2]*phi2 +
       Params[ISpline + 4]*phi3 + Params[ISpline + 6]*phi4;
   Y = Params[ISpline + 1]*phi1 + Params[ISpline + 3]*phi2 +
       Params[ISpline + 5]*phi3 + Params[ISpline + 7]*phi4;

   if(fabs(X)<1.e-12)
    X = 0.;       

   if(m!=0)
    {
     FreePts[0][j] = X;
     FreePts[1][j] = Y;
     
     Intpol_Coord[2*m] = X;
     Intpol_Coord[2*m+1] = Y; 
    }
   
    m++;
// ================================================================================
// for fe values

     u0 = u1_spl[FeDof]*phi1 + u1_spl[FeDof+1]*phi2 +
              FEParams[USpline]*phi3 + FEParams[USpline + 1]*phi4;
     u1 = u2_spl[FeDof]*phi1 + u2_spl[FeDof+1]*phi2 +
              FEParams[USpline+2]*phi3 + FEParams[USpline + 3]*phi4;

//      surf = surf_spl[FeDof]*phi1 + surf_spl[FeDof+1]*phi2 +
//               FEParams[USpline+4]*phi3 + FEParams[USpline + 5]*phi4;

    if(j!=0) // endpoints no need to set
     {      
      Intpol_VeloValues[2*m1] = u0;
      Intpol_VeloValues[2*m1+1] = u1;
//       Intpol_Values[m1] = surf;
     }
    m1++;
// ================================================================================
//     Joint = cell[j]->GetJoint(EdgeNo[j]);
//     isojoint = (TIsoBoundEdge *)Joint;
//     k = isojoint->GetN_Vertices();
       k = ORDER-1;
//     if(k==ORDER-1)
//      {
//       IsoVertices = isojoint->GetVertices();
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
           // cout<< ISpline << ' ' << T;
           T = (T-Param9[i-1])/(Param9[i]-Param9[i-1]);
           // cout<< ' ' << T <<endl;
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

     Intpol_Coord[2*m] = X;
     Intpol_Coord[2*m+1] = Y; 

     // OutPut("NewX:"<<' '<< m <<' '<<X<<' '<< Y<<endl);
     m++;

     // ====================================================================================
     // for fe values

     u0 = u1_spl[FeDof]*phi1 + u1_spl[FeDof+1]*phi2 +
              FEParams[USpline]*phi3 + FEParams[USpline + 1]*phi4;
     u1 = u2_spl[FeDof]*phi1 + u2_spl[FeDof+1]*phi2 +
              FEParams[USpline+2]*phi3 + FEParams[USpline + 3]*phi4;

//      surf = surf_spl[FeDof]*phi1 + surf_spl[FeDof+1]*phi2 +
//                FEParams[USpline+4]*phi3 + FEParams[USpline + 5]*phi4;
    
     Intpol_VeloValues[2*m1] = u0;
     Intpol_VeloValues[2*m1+1] = u1;   
//      Intpol_Values[m1] = surf;      
     m1++;
     // ====================================================================================
     }
//     }
   }  //  for(j=0;j<N_E
   
//    T = ((double)m)*teta;   
//    cout<< "T " <<T << endl;    

   FreePts[0][N_E] = x[N_Splines];
   FreePts[1][N_E] = y[N_Splines];   
   
  Intpol_Coord[2*m] = x[N_Splines];
  Intpol_Coord[2*m+1] = y[N_Splines];
  
  Intpol_VeloValues[2*m1] = u1_spl[N_Splines];;
  Intpol_VeloValues[2*m1+1] = u2_spl[N_Splines];  
//   Intpol_Values[m1] = surf_spl[N_Splines];     

//    cout<< "EndX " <<FreePts[0][N_E] <<" EndY " <<FreePts[1][N_E] <<endl;  
//    cout<< m1 << "VeloX " <<Intpol_VeloValues[2*m1] <<" VeloY " <<Intpol_VeloValues[2*m1+1] <<endl;     

   delete [] h; delete [] t; delete [] a; delete [] b;
   delete [] c; delete [] rhs; delete [] Mx; delete [] My;
   delete [] Params; delete [] Param9;  delete [] x; delete [] y;
   delete [] U_DOF; 
   delete []  u1rhs ;
   delete []  u2rhs;
//    delete []  srhs;
   delete []  u1_spl;
   delete []  u2_spl;
//    delete []  surf_spl;
   delete []  Mu1;
   delete []  Mu2;
//    delete [] Msurf;
   delete [] FEParams;   
   
//       cout<< "Spline " << Intpol_Values[m1] <<endl;
   exit(0);

} // ReParam_axial3D_Data

  
  
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
  
      // end point of the cylinder (free boundary)
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

    if(X < 0 || fabs(X)<1e-8) X = 0.0; // no penetration on axial boundary
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
       OutPut("Iso x " <<  tx << " y " << ty << " X " << X << " Y " << Y << endl);            
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
 
 
void GridVelo_imping(double **Entries, double *Sol, double *d, double *Rhs,
                     int *KCol, int *RowPtr,
                     TFEVectFunct2D *GridPos,
                     TFEVectFunct2D *AuxGridPos,
                     TFEVectFunct2D *Velocity, double dt,
                     TFEVectFunct2D *GridVelocity, 
                     TVertex ***MovBoundVert, int *N_MovVert,
                     TBaseCell **Free_Cells, int **IsoCellEdgeNos,
                     bool &reparam, TFEVectFunct2D *RefGridPos)
{
  int i,j,k,l,m,comp, N;
  int *VeloBeginIndex, *VeloGlobalNumbers;
  int *GridBeginIndex, *GridGlobalNumbers;
  int N_Cells, N_Vertices, N_Edges, N_LocalDOFs;
  int N_Levels, *DOF, *JointDOF, GridLength;
  int N_BoundaryNodes, N_LinePoints, IIso, N_Inner, N_;
  
  double xi[4], eta[4], X[4], Y[4], VX[4], VY[4];
  double FunctValues[4][MaxN_BaseFunctions2D];
  double FEValuesX[MaxN_BaseFunctions2D];
  double FEValuesY[MaxN_BaseFunctions2D];
  double *ValuesX, *ValuesY, *RefValueX, *RefValueY;
  double *ValuesVX, *ValuesVY, *NewValuesX, *NewValuesY;
  double s, t, x, y, h_tot, x0, x1, y0, y1;
  double res, oldres, *gridvelo, *Nx, *Ny, Ay;
  double *LineWeights, *zeta;
  double normalx, normaly, tangenx, tangeny, nx, ny, tx, ty;
  double un, hE, t0,t1, temp2, eps=1e-6;
  double h, hmin, hmax, hlimit, *IsoX, *IsoY;
   
  BF2DRefElements RefElement;
  TRefTrans2D *F_K;
  RefTrans2D RefTrans;
  TIsoBoundEdge *isojoint;  
  TMGLevel2D *Level;
  TFESpace2D *VelocitySpace, *GridSpace;
  TCollection *Coll;
  TBaseCell *cell;
  BoundTypes bdtype;
  TBoundEdge *BoundEdge;
  TBoundComp2D *BoundComp;
  FE2D FEId;
  TFE2D *Element;
  TFEDesc2D *FEDesc;
  BaseFunct2D BF;
  TBaseFunct2D *bf;
  boolean OnBoundary;
  TJoint *joint;
  TVertex **Vertices;
  TQuadFormula1D *qf1;
  QuadFormula1D LineQuadFormula;
  BoundCond Cond0, Cond1;
    
  VelocitySpace = Velocity->GetFESpace2D();
  VeloBeginIndex = VelocitySpace->GetBeginIndex();
  VeloGlobalNumbers = VelocitySpace->GetGlobalNumbers();
  ValuesVX = Velocity->GetValues();
  ValuesVY = ValuesVX + Velocity->GetLength();

  GridPos->GridToData();
  GridSpace = GridPos->GetFESpace2D();
  GridBeginIndex = GridSpace->GetBeginIndex();
  GridGlobalNumbers = GridSpace->GetGlobalNumbers();
  GridLength = GridPos->GetLength();
  ValuesX = GridPos->GetValues();
  ValuesY = ValuesX + GridLength;

  N_Inner = GridSpace->GetN_Inner();
  N_BoundaryNodes = GridLength - GridSpace->GetN_Inner();
  
  if(TDatabase::ParamDB->P5 > 0)
  {
   Nx = new double[N_BoundaryNodes];
   Ny = new double[N_BoundaryNodes];
   memset(Nx, 0, N_BoundaryNodes*SizeOfDouble);
   memset(Ny, 0, N_BoundaryNodes*SizeOfDouble);
  }

 // cout << "N_BoundaryNodes: " << N_BoundaryNodes << endl;
//   cout << GridLength << " " << N_DOF << endl;

  NewValuesX = AuxGridPos->GetValues();
  NewValuesY = NewValuesX + GridLength;

  memcpy(NewValuesX, ValuesX, 2*GridLength*SizeOfDouble);

  RefValueX = RefGridPos->GetValues();
  RefValueY = RefValueX + GridLength;   
  
  Coll = VelocitySpace->GetCollection();
  N_Cells = Coll->GetN_Cells();
//  cout << "N_Cells: " <<N_Cells<< endl;
  // determine outer normal vectors

  IIso = N_BoundaryNodes;
  // Outward normal no need if we move boundary with velocity
 if(TDatabase::ParamDB->P5 > 0)
  {
  // determine outer normal vectors
   for(i=0;i<N_Cells;i++)
    {
    // cout << "cell: " << i << endl;
    cell = Coll->GetCell(i);
    cell->SetClipBoard(i);
    N_Edges = cell->GetN_Edges();

    for(j=0;j<N_Edges;j++)
    {
      if( !(cell->GetJoint(j)->InnerJoint()) )
      {
        // cout << "joint: " << j << endl;
        cell->GetVertex(j)->GetCoords(x0, y0);
        cell->GetVertex((j+1)%N_Edges)->GetCoords(x1, y1);
        t = x1-x0;
        s = y1-y0;
        FEId = VelocitySpace->GetFE2D(i, cell);
        l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
        LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
        qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
        qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);

        RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
        switch(RefElement)
        {
          case BFUnitTriangle:
            RefTrans = TriaIsoparametric;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            ((TTriaIsoparametric *)F_K)->SetCell(cell);
          break;

          case BFUnitSquare:
            RefTrans = QuadIsoparametric;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            ((TQuadIsoparametric *)F_K)->SetCell(cell);
          break;

          default:
            Error("only triangles and quadrilaterals are allowes" << endl);
            Error("file: " << __FILE__ << " line " << __LINE__ << endl);
            exit(-1);
        } // endswitch

        normalx = 0;
        normaly = 0;
        hE = 0;
        for(k=0;k<N_LinePoints;k++)
        {
          F_K->GetOuterNormal(j, zeta[k], nx, ny);
          F_K->GetTangent(j, zeta[k], tx, ty);
          t = sqrt(tx*tx+ty*ty);
          normalx += t*LineWeights[k]*nx;
          normaly += t*LineWeights[k]*ny;
          hE += t*LineWeights[k];
          // cout << "k= " << k << " " << nx << " " << ny << endl;
        } // endfor k

        DOF = GridGlobalNumbers + GridBeginIndex[i];

        switch(N_Edges)
        {
          case 3:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 2:
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;
            } // endswitch j
          break;

          case 4:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[3] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 2:
                l = DOF[3] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 3:
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;
            } // endswitch j
          break;
        } // endswitch N_Edges
      } // !InnerJoint
    } // endfor j
  } // endfor i

  N_ = IIso;
  // normalize normal vector
  for(i=0;i<N_;i++)
  {
    x = Nx[i];
    y = Ny[i];
    t = sqrt(x*x+y*y);
    Nx[i] /= t;
    Ny[i] /= t;
  }
 }
  // determine new position of boundary vertices
  for(i=0;i<N_Cells;i++)
  {
    cell  = Coll->GetCell(i);
    N_Edges = cell->GetN_Edges();
    OnBoundary = FALSE;
    for(j=0;j<N_Edges;j++)
    {
      if( !(cell->GetJoint(j)->InnerJoint()) || 
           (cell->GetJoint(j)->GetType() == InterfaceJoint)  || 
           (cell->GetJoint(j)->GetType() == IsoInterfaceJoint) )
       {
        OnBoundary = TRUE;
        joint = cell->GetJoint(j);
       }
    } // endfor j

    if(OnBoundary)
    {

      FEId = VelocitySpace->GetFE2D(i, cell);
      Element = TFEDatabase2D::GetFE2D(FEId);
      bf = Element->GetBaseFunct2D();
      N_LocalDOFs = Element->GetN_DOF();

      switch(N_Edges)
      {
        case 3:
          xi[0]  = 0; xi[1]  = 1; xi[2]  = 0;
          eta[0] = 0; eta[1] = 0; eta[2] = 1;
          VX[0] = VX[1] = VX[2] = 0;
          VY[0] = VY[1] = VY[2] = 0;
          N_Vertices = N_Edges;
        break;

        case 4:
          // NOTE: sorted due to number of DOFs !!!
          xi[0]  = -1; xi[1]  =  1; xi[2]  = -1; xi[3]  =  1;
          eta[0] = -1; eta[1] = -1; eta[2] =  1; eta[3] =  1;
          VX[0] = VX[1] = VX[2] = VX[3] = 0;
          VY[0] = VY[1] = VY[2] = VY[3] = 0;
          N_Vertices = N_Edges;
        break;

        default:
          Error("only triangles and quadrilaterals are allowed!" << endl);
          exit(-1);
      } // endswitch

      for(j=0;j<N_Vertices;j++)
        bf->GetDerivatives(D00, xi[j], eta[j], FunctValues[j]);

      DOF = VeloGlobalNumbers + VeloBeginIndex[i];

      for(j=0;j<N_LocalDOFs;j++)
      {
        k = DOF[j];
        s = ValuesVX[k];
        t = ValuesVY[k];

        for(l=0;l<N_Vertices;l++)
        {
          VX[l] += FunctValues[l][j]*s;
          VY[l] += FunctValues[l][j]*t;
        } // endfor l
      } // endfor j

      FEId = GridSpace->GetFE2D(i, cell);
      Element = TFEDatabase2D::GetFE2D(FEId);
      BF = Element->GetBaseFunct2D_ID();
      if( (BF != BF_C_T_P1_2D) && (BF != BF_C_Q_Q1_2D) )
      {
        Error("Grid Space must be conforming and of first order!" << endl);
        exit(-1);
      }  // endif
      bf = Element->GetBaseFunct2D();
      N_LocalDOFs = Element->GetN_DOF();

      switch(N_Edges)
      {
        case 3:
          xi[0]  = 0; xi[1]  = 1; xi[2]  = 0;
          eta[0] = 0; eta[1] = 0; eta[2] = 1;
          X[0] = X[1] = X[2] = 0;
          Y[0] = Y[1] = Y[2] = 0;
          N_Vertices = N_Edges;
        break;

        case 4:
          // NOTE: sorted due to number of DOF !!!
          xi[0]  = -1; xi[1]  =  1; xi[2]  = -1; xi[3]  =  1;
          eta[0] = -1; eta[1] = -1; eta[2] =  1; eta[3] =  1;
          X[0] = X[1] = X[2] = X[3] = 0;
          Y[0] = Y[1] = Y[2] = Y[3] = 0;
          N_Vertices = N_Edges;
        break;

        default:
          Error("only triangles and quadrilaterals are allowed!" << endl);
          exit(-1);
      } // endswitch

      DOF = GridGlobalNumbers + GridBeginIndex[i];

      for(j=0;j<N_Vertices;j++)
      {
        l = DOF[j];
        k = l - N_Inner;
        if(k>=0)
        {
          if((TDatabase::ParamDB->P5 > 0) && (ValuesY[l] != 0) )
          {
            un = VX[j]*Nx[k] + VY[j]*Ny[k];
            NewValuesX[l] = ValuesX[l] + dt*un*Nx[k];
            NewValuesY[l] = ValuesY[l] + dt*un*Ny[k];
          }
          else
          {
// 	    if(ValuesX[l] == 0 || (ValuesX[l]==TDatabase::ParamDB->P7 && ValuesY[l]>=0) )
// 	    { NewValuesX[l] = ValuesX[l]; }
// 	    else
	    { NewValuesX[l] = ValuesX[l] + dt*VX[j]; }

//             if( (ValuesX[l]==0 && ValuesY[l] >(0.5*TDatabase::ParamDB->P8)) || ValuesY[l] == TDatabase::ParamDB->P8 
// 	        || (ValuesX[l] == TDatabase::ParamDB->P7 && ValuesY[l]>=0) ) 
//              { NewValuesY[l] = ValuesY[l];  }
//             else    
	     { NewValuesY[l] = ValuesY[l] + dt*VY[j];  }
	    
	  }
       } //  if(k>=0)
 //    Due to spline approximation solid boundary end vertices may take negative y value
//         if(NewValuesY[l]<0.0 ) NewValuesY[l] = 0.0;
        if( fabs(NewValuesX[l]) < 1e-10 ) NewValuesX[l] = 0.0;
      } // endfor j
    } // endif
  } // endfor i

// cout << " dt " << dt <<endl;
/*
   for(i=0;i<GridLength;i++)
 cout << i <<"  ---  " <<ValuesX[i] << "  ---  " << NewValuesX[i] << endl;
*/
// exit(0);
  
   MovBoundVert[1][0]->GetCoords(x, Ay);   

   AuxGridPos->DataToGrid();
   
  //======================================================================  
  //  Reparametrization of free surface - Begin
  //======================================================================  
  if(!reparam)
  {
   h_tot = 0;
   hmin = 100;
   hmax = 0.0;
   
   for(k=0;k<N_MovVert[0];k++)
    {
     MovBoundVert[0][k]->GetCoords(x1, y1);

     if(k==N_MovVert[0]-1)
      {  
       x= TDatabase::ParamDB->P7;
       y= 0.;
      }
     else
     { MovBoundVert[0][k+1]->GetCoords(x, y); }

     h = sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1));
     h_tot +=h;
     if (h < hmin) hmin = h;
     if (h > hmax) hmax = h;
    } // for(k=0;k<N_MovVert[0];k++) 
   
    h_tot /= (double)N_MovVert[0];
    hlimit =  h_tot;   
   }
 
   
   if( (hmin <  0.5*hlimit) || (hmax > 1.5*hlimit) ||  reparam )
   { 
    OutPut("Reparam b "  << hmin << " " << 0.5*hlimit << " " << hmax << " " << 1.5*hlimit << endl);
    //before reparam update iso points, since we use interpolated cubic spline
    //which pass through iso points also
    IsoX = new double[IIso];
    IsoY = new double[IIso];
    
    IIso = 0;
    for(i=0;i<N_MovVert[0];i++)
    {
     cell = Free_Cells[i];
      j = IsoCellEdgeNos[1][i];
      joint = cell->GetJoint(j);
      isojoint = (TIsoBoundEdge *)joint;
      k = isojoint->GetN_Vertices();
      Vertices = isojoint->GetVertices();
      FEId = VelocitySpace->GetFE2D(IsoCellEdgeNos[0][i], cell);
      FEDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
      m = FEDesc->GetN_JointDOF();
      if(m == k+2)
        {
          JointDOF = FEDesc->GetJointDOF(j);
          DOF =  VeloGlobalNumbers+VeloBeginIndex[IsoCellEdgeNos[0][i]];
          for(l=0;l<k;l++)
          {
            m = DOF[JointDOF[l+1]];
            Vertices[l]->GetCoords(IsoX[IIso], IsoY[IIso]);
            if(TDatabase::ParamDB->P5 > 0)
            {
              un = ValuesVX[m]*Nx[IIso+N_BoundaryNodes]
                  + ValuesVY[m]*Ny[IIso+N_BoundaryNodes];
              x  = IsoX[IIso] + dt*un*Nx[IIso+N_BoundaryNodes];
              y  = IsoY[IIso] + dt*un*Ny[IIso+N_BoundaryNodes];
            }
            else
            {
             x  = IsoX[IIso] + dt*ValuesVX[m];
             y  = IsoY[IIso] + dt*ValuesVY[m];
            }
            
//            if(y<=0) y = 1e-5;
//            if(fabs(x)<1e-12) x = 0.;
           
           Vertices[l]->SetCoords(x, y);
           IIso++;
          } // endfor l
        }
       else
        {
          // approximation order of isoparametric boundary and velocity
          // element must be the same
          Error("No match in isoparametric case" << endl);
          exit(-1);
        } //  if(m == k+2)     
    }// for(i=0;i<N_MovVert[0];i++)

    ReParam_axial3D_U(N_MovVert[0], Free_Cells,  IsoCellEdgeNos[1], IsoCellEdgeNos[0],
                      Velocity, NULL, NULL, FALSE);  
 
//    exit(0);
//     cout << "Surfact reparam has to be included in ReParam_axial3D_U" << endl;
    
    reparam = TRUE;   
    RefGridPos->GridToData();    
    Daxpy(2*GridLength, -1, ValuesX, RefValueX); // now reparamdisp in RefValueX

    //back to orig mesh (no reparam movement in calculation of free surf w)
    AuxGridPos->DataToGrid(); 
    
   //restore iso points
    IIso = 0;
    for(i=0;i<N_MovVert[0];i++)
    {
     cell = Free_Cells[i];
     j = IsoCellEdgeNos[1][i];      
     joint = cell->GetJoint(j);
     isojoint = (TIsoBoundEdge *)joint;
     k = isojoint->GetN_Vertices();
     Vertices = isojoint->GetVertices();
     
     for(l=0;l<k;l++)
      {
        Vertices[l]->SetCoords(IsoX[IIso], IsoY[IIso]);
        IIso++;
      }
    }// for(i=0;i<N_Cells;i++)    
    
    delete [] IsoX; delete [] IsoY; 
   } // if ( ((hmin < hlimit) || (hmax > 3.*h_tot/2.)) ||  reparam )
   //======================================================================       

//    MovBoundVert[2][0]->GetCoords(x, y); // right wetting point
//    y=0.;
//    h_tot = x;
//    h_tot /= (double)N_MovVert[0];
   for(i=0;i<N_MovVert[1];i++)
    {
     MovBoundVert[1][i]->GetCoords(x, y);
     MovBoundVert[1][i]->SetCoords(0., y);  
    }
 

// // axial boundary top
//    MovBoundVert[0][0]->GetCoords(x, Ay);
//    MovBoundVert[1][0]->GetCoords(x, y);
//  
//     cout<< " x " << x << " Ay" << Ay<<endl;  
//    N=N_MovVert[1];      
//    h_tot = (y-Ay)/(double)N;   
//    N--;
//    
//     for(i=1;i<N-1;i++)
//     {
//      MovBoundVert[1][i]->GetCoords(x, y);
//      cout<< " x " << x << " y " << y<<" new y " << y +((double)(N-i))*h_tot<<endl;      
//      y += ((double)(N-i))*h_tot;
//      MovBoundVert[1][i]->SetCoords(x, y);   
//     }      
//    
//     
//     
//     exit(0);
    
    
   AuxGridPos->GridToData();   
   GridPos->DataToGrid();
   
   
//       MovBoundVert[1][0]->GetCoords(x, y);
//     cout << " x " << x <<" y " << y <<endl;  
//     exit(0); 
    
   memset(Rhs, 0, 2*GridLength*SizeOfDouble);

   memcpy(d, NewValuesX, 2*GridLength*SizeOfDouble);
   Daxpy(2*GridLength, -1, ValuesX, d);
   memcpy(Rhs + (GridLength-N_BoundaryNodes), d+(GridLength-N_BoundaryNodes),
          N_BoundaryNodes*SizeOfDouble);
   memcpy(Rhs + (2*GridLength-N_BoundaryNodes), d+(2*GridLength-N_BoundaryNodes),
          N_BoundaryNodes*SizeOfDouble);

  memset(Sol, 0 , 2*GridLength*SizeOfDouble);
  memcpy(Sol + (GridLength-N_BoundaryNodes), d+(GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);
  memcpy(Sol + (2*GridLength-N_BoundaryNodes), d+(2*GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);
 
//     for(i=0;i<GridLength;i++)
//      cout<< i <<"  ---  "<< Rhs[i] << "  ---  " << Rhs[i+GridLength] << endl;
//     
  SolveGridEquation(Entries, Sol, Rhs, KCol, RowPtr, GridLength);

  gridvelo = GridVelocity->GetValues();
  memcpy(gridvelo, Sol, 2*GridLength*SizeOfDouble);
  Dscal(2*GridLength, 1./dt, gridvelo);
//     cout << " x " << x <<" y " << y <<endl;  
//     exit(0);  
  //======================================================================  
  //  Reparametrization of free surface 
  //======================================================================      
  if(reparam)
  {
   memset(Rhs, 0, 2*GridLength*SizeOfDouble);    
   memcpy(Rhs + (GridLength-N_BoundaryNodes), RefValueX+(GridLength-N_BoundaryNodes),
          N_BoundaryNodes*SizeOfDouble);
   memcpy(Rhs + (2*GridLength-N_BoundaryNodes), RefValueX+(2*GridLength-N_BoundaryNodes),
          N_BoundaryNodes*SizeOfDouble);
 
  memset(Sol, 0 , 2*GridLength*SizeOfDouble);
  memcpy(Sol + (GridLength-N_BoundaryNodes), RefValueX+(GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);
  memcpy(Sol + (2*GridLength-N_BoundaryNodes), RefValueX+(2*GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble); 
 
  SolveGridEquation(Entries, Sol, Rhs, KCol, RowPtr, GridLength);
  
  Dscal(2*GridLength, 1./dt, Sol);
   
  //only inner mesh velo, since U will be interpolated in free surf reparm     
  Daxpy(GridLength-N_BoundaryNodes, 1., Sol, gridvelo);
  Daxpy(GridLength-N_BoundaryNodes, 1., Sol+GridLength, gridvelo+GridLength);  
  
//     cout<< "reparam : "<<endl;
//     exit(0);
  }  
  
//     for(i=0;i<GridLength;i++)
//      cout<< i <<"  ---  "<< gridvelo[i] << "  ---  " << gridvelo[i+GridLength] << endl;

//   delete [] d;
   
  if(TDatabase::ParamDB->P5 > 0)
  { delete [] Nx;  delete [] Ny; }
  
  
} // GridVelo_imping

  

void MoveGrid_imping(double **Entries, double *Sol, double *d, double *Rhs,
                  int *KCol, int *RowPtr,
                  TFEVectFunct2D *GridPos,
                  TFEVectFunct2D *Velocity, double dt,
                  TFEVectFunct2D *GridVelocity, 
                  TFEVectFunct2D *NewGridPos, 
                  TVertex ***MovBoundVert, int *N_MovVert,
                  TBaseCell **Free_Cells, int **IsoCellEdgeNos,
                  TFEFunction2D *Surfactant, TFEFunction2D *SurfSurfactant,
                  bool &reparam, int &N_ReParam)
{
  int i,j,k,l,m, N;
  int *VeloBeginIndex, *VeloGlobalNumbers;
  int *GridBeginIndex, *GridGlobalNumbers;
  int N_Cells, N_Vertices, N_Edges, N_LocalDOFs;
  int N_Levels, *DOF, *JointDOF, GridLength, polydegree;
  int N_BoundaryNodes, N_LinePoints, IIso, N_Inner, N_;

  double xi[4], eta[4], X[4], Y[4], VX[4], VY[4];
  double FunctValues[4][MaxN_BaseFunctions2D];
  double FEValuesX[MaxN_BaseFunctions2D];
  double FEValuesY[MaxN_BaseFunctions2D];
  double *ValuesX, *ValuesY;
  double *ValuesVX, *ValuesVY;
  double *NewValuesX, *NewValuesY;
  double s, t, x, y, IsoX, IsoY, Ay;
  double x0, x1, y0, y1, h_tot, res, oldres;
  double *gridvelo, *Nx, *Ny, *LineWeights, *zeta;
  double normalx, normaly, tangenx, tangeny, nx, ny, tx, ty;
  double un, hE, t0,t1, temp2, x_max, x_min, temp, eps=1e-6;  
  
  TIsoBoundEdge *isojoint;
  TMGLevel2D *Level;  
  TFESpace2D *VelocitySpace, *GridSpace;
  TCollection *Coll;
  TBaseCell *cell;
  BoundTypes bdtype;
  TBoundEdge *BoundEdge;
  TBoundComp2D *BoundComp;
  FE2D FEId;
  TFE2D *Element;
  TFEDesc2D *FEDesc;
  BaseFunct2D BF;
  TBaseFunct2D *bf;
  boolean OnBoundary;
  TJoint *joint;
  TVertex **Vertices;
  TQuadFormula1D *qf1;
  QuadFormula1D LineQuadFormula;
  BF2DRefElements RefElement;
  TRefTrans2D *F_K;
  RefTrans2D RefTrans;  
  QuadFormula2D QuadFormula;
  BoundCond Cond0, Cond1;  
 
  VelocitySpace = Velocity->GetFESpace2D();
  VeloBeginIndex = VelocitySpace->GetBeginIndex();
  VeloGlobalNumbers = VelocitySpace->GetGlobalNumbers();
  ValuesVX = Velocity->GetValues();
  ValuesVY = ValuesVX + Velocity->GetLength();

  GridPos->GridToData();
  GridSpace = GridPos->GetFESpace2D();
  GridBeginIndex = GridSpace->GetBeginIndex();
  GridGlobalNumbers = GridSpace->GetGlobalNumbers();
  GridLength = GridPos->GetLength();
  ValuesX = GridPos->GetValues();
  ValuesY = ValuesX + GridLength;

  N_Inner = GridSpace->GetN_Inner();
  N_BoundaryNodes = GridLength - GridSpace->GetN_Inner();
  // cout << "N_BoundaryNodes: " << N_BoundaryNodes << endl;

//   d = new double[ 2*GridLength];
  
  if(TDatabase::ParamDB->P5 > 0)
  {  
   Nx = new double[2*N_BoundaryNodes]; // additional values for edge midpoints
   Ny = new double[2*N_BoundaryNodes]; // additional values for edge midpoints
   memset(Nx, 0, 2*N_BoundaryNodes*SizeOfDouble);
   memset(Ny, 0, 2*N_BoundaryNodes*SizeOfDouble);
  }
  //cout << "N_BoundaryNodes: " << N_BoundaryNodes << endl;
 // cout << GridLength << " " << N_DOF << endl;

  NewValuesX = NewGridPos->GetValues();
  NewValuesY = NewValuesX + GridLength;

  memcpy(NewValuesX, ValuesX, 2*GridLength*SizeOfDouble);

  Coll = VelocitySpace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  // determine outer normal vectors
  IIso = N_BoundaryNodes;

    // Outward normal no need if we move boundary with velocity
  if(TDatabase::ParamDB->P5 > 0)
  {
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_Edges = cell->GetN_Edges();

    for(j=0;j<N_Edges;j++)
    {
      if( !(cell->GetJoint(j)->InnerJoint()) )
      {
        FEId = VelocitySpace->GetFE2D(i, cell);
        l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
        LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
        qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
        qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);

        RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
        switch(RefElement)
        {
          case BFUnitTriangle:
            RefTrans = TriaIsoparametric;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
            QuadFormula = TFEDatabase2D::GetQFTriaFromDegree(2*polydegree-1);
//             ((TTriaIsoparametric *)F_K)->SetApproximationOrder(polydegree);
            ((TTriaIsoparametric *)F_K)->SetQuadFormula(QuadFormula);
            ((TTriaIsoparametric *)F_K)->SetCell(cell);
          break;

          case BFUnitSquare:
            RefTrans = QuadIsoparametric;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
            QuadFormula = TFEDatabase2D::GetQFQuadFromDegree(3*polydegree);
//             ((TQuadIsoparametric *)F_K)->SetApproximationOrder(polydegree);
            ((TQuadIsoparametric *)F_K)->SetQuadFormula(QuadFormula);
            ((TQuadIsoparametric *)F_K)->SetCell(cell);
          break;

          default:
            Error("only triangles and quadrilaterals are allowes" << endl);
            Error("file: " << __FILE__ << " line " << __LINE__ << endl);
            exit(-1);
        } // endswitch

        normalx = 0;
        normaly = 0;
        hE = 0;
        for(k=0;k<N_LinePoints;k++)
        {
          F_K->GetOuterNormal(j, zeta[k], nx, ny);
          F_K->GetTangent(j, zeta[k], tx, ty);
          t = sqrt(tx*tx+ty*ty);
          normalx += t * LineWeights[k] * nx;
          normaly += t * LineWeights[k] * ny;
          hE += t * LineWeights[k];
          // cout << "k= " << k << " " << nx << " " << ny << endl;
        } // endfor k

        DOF = GridGlobalNumbers + GridBeginIndex[i];
/*
        switch(N_Edges)
        {
          case 3:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[1] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[2] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;

              case 2:
                l = DOF[2] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[0] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;
            } // endswitch j
          break;

          case 4:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[1] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[3] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;

              case 2:
                l = DOF[3] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[2] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;

              case 3:
                l = DOF[2] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[0] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;
            } // endswitch j
          break;
        } // endswitch N_Edges
*/

// /*
        switch(N_Edges)
        {
          case 3:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 2:
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;
            } // endswitch j
          break;

          case 4:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[3] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 2:
                l = DOF[3] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 3:
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;
            } // endswitch j
          break;
        } // endswitch N_Edges
// */

        if(cell->GetJoint(j)->GetType() == IsoBoundEdge)
        {
          FEId = VelocitySpace->GetFE2D(i, cell);
          FEDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
          N_LocalDOFs = FEDesc->GetN_JointDOF();
          t = 2.0/(N_LocalDOFs-1);
          for(k=1;k<N_LocalDOFs-1;k++)
          {
            /*
            Nx[IIso] += normalx;
            Ny[IIso] += normaly;
            */
            // /*
            s = -1.0 + k*t;
            F_K->GetOuterNormal(j, s, nx, ny);
            Nx[IIso] += nx;
            Ny[IIso] += ny;
            // */
            IIso++;
          } // endfor
        } // endif
      } // !InnerJoint
    } // endfor j
  } // endfor i

  N_ = IIso;
  // normalize normal vector
  for(i=0;i<N_;i++)
  {
    x = Nx[i];
    y = Ny[i];
    t = sqrt(x*x+y*y);
    Nx[i] /= t;
    Ny[i] /= t;

    // cout << setw(5) << i << "n = (" << Nx[i] << ", " << Ny[i] << ")" << endl;
  }
}

  // determine new position of boundary vertices
  for(i=0;i<N_Cells;i++)
  {
    cell  = Coll->GetCell(i);
    N_Edges = cell->GetN_Edges();
    OnBoundary = FALSE;
    for(j=0;j<N_Edges;j++)
    {
      if( !(cell->GetJoint(j)->InnerJoint()) )
       {
        OnBoundary = TRUE;
        joint = cell->GetJoint(j);
       }
    } // endfor j

    if(OnBoundary)
    {
      FEId = VelocitySpace->GetFE2D(i, cell);
      Element = TFEDatabase2D::GetFE2D(FEId);
      bf = Element->GetBaseFunct2D();
      N_LocalDOFs = Element->GetN_DOF();

      switch(N_Edges)
      {
        case 3:
          xi[0]  = 0; xi[1]  = 1; xi[2]  = 0;
          eta[0] = 0; eta[1] = 0; eta[2] = 1;
          VX[0] = VX[1] = VX[2] = 0;
          VY[0] = VY[1] = VY[2] = 0;
          N_Vertices = N_Edges;
        break;

        case 4:
          // NOTE: sorted due to number of DOFs !!!
          xi[0]  = -1; xi[1]  =  1; xi[2]  = -1; xi[3]  =  1;
          eta[0] = -1; eta[1] = -1; eta[2] =  1; eta[3] =  1;
          VX[0] = VX[1] = VX[2] = VX[3] = 0;
          VY[0] = VY[1] = VY[2] = VY[3] = 0;
          N_Vertices = N_Edges;
        break;

        default:
          Error("only triangles and quadrilaterals are allowed!" << endl);
          exit(-1);
      } // endswitch

      for(j=0;j<N_Vertices;j++)
        bf->GetDerivatives(D00, xi[j], eta[j], FunctValues[j]);

      DOF = VeloGlobalNumbers + VeloBeginIndex[i];

      for(j=0;j<N_LocalDOFs;j++)
      {
        k = DOF[j];
        s = ValuesVX[k];
        t = ValuesVY[k];

        for(l=0;l<N_Vertices;l++)
        {
          VX[l] += FunctValues[l][j]*s;
          VY[l] += FunctValues[l][j]*t;
        } // endfor l
      } // endfor j

      FEId = GridSpace->GetFE2D(i, cell);
      Element = TFEDatabase2D::GetFE2D(FEId);
      BF = Element->GetBaseFunct2D_ID();
      if( (BF != BF_C_T_P1_2D) && (BF != BF_C_Q_Q1_2D) )
      {
        Error("Grid Space must be conforming and of first order!" << endl);
        exit(-1);
      }  // endif
      bf = Element->GetBaseFunct2D();
      N_LocalDOFs = Element->GetN_DOF();


      DOF = GridGlobalNumbers + GridBeginIndex[i];

      for(j=0;j<N_Vertices;j++)
      {
        l = DOF[j];
        k = l - N_Inner;
        if(k>=0)
        {
          if((TDatabase::ParamDB->P5 > 0) && (ValuesY[l] != 0) )
          {
            un = VX[j]*Nx[k] + VY[j]*Ny[k];
            NewValuesX[l] = ValuesX[l] + dt*un*Nx[k];
            NewValuesY[l] = ValuesY[l] + dt*un*Ny[k];
          }
          else
          {
 
// 	    if(ValuesX[l] == 0 || (ValuesX[l] == TDatabase::ParamDB->P7 && ValuesY[l]>=0) )
// 	    { NewValuesX[l] = ValuesX[l]; }
// 	    else
	    { NewValuesX[l] = ValuesX[l] + dt*VX[j]; }

/*            if( (ValuesX[l]==0 && ValuesY[l] >(0.5*TDatabase::ParamDB->P8)) || ValuesY[l] == TDatabase::ParamDB->P8 
	         || (ValuesX[l] == TDatabase::ParamDB->P7 && ValuesY[l]>=0) ) 
	     { NewValuesY[l] = ValuesY[l];  }
            else   */ 
	     { NewValuesY[l] = ValuesY[l] + dt*VY[j];  }                
          }
        }
 //    Due to spline approximation solid boundary end vertices may take negative y value
//          if(NewValuesY[l]<0.0 ) NewValuesY[l] = 0.0;
        if( fabs(NewValuesX[l]) < 1e-12 ) NewValuesX[l] = 0.0;
      } // endfor j
     } // endif
   } // endfor i

/*

   for(i=GridLength-N_BoundaryNodes;i<GridLength;i++)
 cout << i <<"  ---  " <<NewValuesX[i] << "  ---  " << NewValuesY[i]<<endl;

exit(0);
*/

   MovBoundVert[1][0]->GetCoords(x, Ay);   
    
   NewGridPos->DataToGrid();
    
   //axial boundary
   for(i=0;i<N_MovVert[1];i++)
    {
     MovBoundVert[1][i]->GetCoords(x, y);
     MovBoundVert[1][i]->SetCoords(0., y);  
    }
//    // axial boundary top
//    MovBoundVert[1][0]->GetCoords(x, y);
//    
//    N=N_MovVert[1];      
//    h_tot = (y-Ay)/(double)N;   
//    N--;
//    
//    for(i=1;i<N;i++)
//     {
//      MovBoundVert[1][i]->GetCoords(x, y);
// //      cout<< " y " << y <<" new y " << y +((double)(N-i))*h_tot<<endl;      
//      y += ((double)(N-i))*h_tot;
//      MovBoundVert[1][i]->SetCoords(x, y);   
//     }       
//  

  //======================================================================  
  //  Reparametrization of free surface - Begin
  //======================================================================     
   if(reparam)  
   {   
    // update the iso points and then reparmetrize the free surf
    IIso = 0;
    for(i=0;i<N_MovVert[0];i++)
    {
     cell = Free_Cells[i];
     j = IsoCellEdgeNos[1][i];
     joint = cell->GetJoint(j);
     isojoint = (TIsoBoundEdge *)joint;
     k = isojoint->GetN_Vertices();
     Vertices = isojoint->GetVertices();
     FEId = VelocitySpace->GetFE2D(IsoCellEdgeNos[0][i], cell);
     FEDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
     m = FEDesc->GetN_JointDOF();
     if(m == k+2)
      {
       JointDOF = FEDesc->GetJointDOF(j);
       DOF =  VeloGlobalNumbers+VeloBeginIndex[IsoCellEdgeNos[0][i]];
       for(l=0;l<k;l++)
        {
         m = DOF[JointDOF[l+1]];
         Vertices[l]->GetCoords(IsoX, IsoY);
         if(TDatabase::ParamDB->P5 > 0)
          {
           un = ValuesVX[m]*Nx[IIso+N_BoundaryNodes]
                + ValuesVY[m]*Ny[IIso+N_BoundaryNodes];
              IsoX += dt*un*Nx[IIso+N_BoundaryNodes];
              IsoY += dt*un*Ny[IIso+N_BoundaryNodes];
            }
            else
            {
             IsoX += dt*ValuesVX[m];
             IsoY += dt*ValuesVY[m];
            }
            
//            if(IsoY<=0) IsoY = 1e-5;
//            if(fabs(x)<1e-12) x = 0.;
           
           Vertices[l]->SetCoords(IsoX, IsoY);
           IIso++;
          } // endfor l
        }
       else
        {
          // approximation order of isoparametric boundary and velocity
          // element must be the same
          Error("No match in isoparametric case" << endl);
          exit(-1);
        } //  if(m == k+2)     
    }// for(i=0;i<N_MovVert[0];i++)     
       
    ReParam_axial3D_U(N_MovVert[0], Free_Cells,  IsoCellEdgeNos[1], IsoCellEdgeNos[0],
                      Velocity, Surfactant, SurfSurfactant, TRUE);     

    OutPut("ReParam CURRENT TIME: ");
    OutPut(TDatabase::TimeDB->CURRENTTIME << endl);   
   }  //if (reparam)  
   
   NewGridPos->GridToData();
   GridPos->DataToGrid();

   memset(Rhs, 0, 2*GridLength*SizeOfDouble);
   memcpy(d, NewValuesX, 2*GridLength*SizeOfDouble);
   Daxpy(2*GridLength, -1, ValuesX, d);
   memcpy(Rhs + (GridLength-N_BoundaryNodes), d+(GridLength-N_BoundaryNodes), N_BoundaryNodes*SizeOfDouble);
   memcpy(Rhs + (2*GridLength-N_BoundaryNodes), d+(2*GridLength-N_BoundaryNodes), N_BoundaryNodes*SizeOfDouble);

  memset(Sol, 0 , 2*GridLength*SizeOfDouble);
  memcpy(Sol + (GridLength-N_BoundaryNodes),
         d+(GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);
  memcpy(Sol + (2*GridLength-N_BoundaryNodes),
         d+(2*GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);


  SolveGridEquation(Entries, Sol, Rhs, KCol, RowPtr, GridLength);

  memcpy(d, ValuesX, 2*GridLength*SizeOfDouble);
  Daxpy(2*GridLength, 1, Sol, d);
  memcpy(NewValuesX, d, (GridLength-N_BoundaryNodes)*SizeOfDouble);
  memcpy(NewValuesY, d+GridLength, (GridLength-N_BoundaryNodes)*SizeOfDouble);

  // needed for scalar equations
  gridvelo = GridVelocity->GetValues();
  memcpy(gridvelo, Sol, 2*GridLength*SizeOfDouble);
  Dscal(2*GridLength, 1./dt, gridvelo);  
  
// for(i=GridLength-N_BoundaryNodes;i<GridLength;i++)
//  cout << i << "  ---  "<<NewValuesX[i] << "  ---  " << NewValuesX[i+GridLength] << endl;
//cout << "test " << endl;
  // put solution into grid position
  IIso = 0;
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
          cell->GetVertex(j)->SetCoords(NewValuesX[k], NewValuesY[k]); 
        }
      break;

      case 4:
        k = DOF[0];
        cell->GetVertex(0)->SetCoords(NewValuesX[k], NewValuesY[k]);
        k = DOF[1];
        cell->GetVertex(1)->SetCoords(NewValuesX[k], NewValuesY[k]);
        k = DOF[3];
        cell->GetVertex(2)->SetCoords(NewValuesX[k], NewValuesY[k]);
        k = DOF[2];
        cell->GetVertex(3)->SetCoords(NewValuesX[k], NewValuesY[k]);
      break;
    } // endswitch
 
  if (!reparam)  
   {  
    N_Edges = cell->GetN_Edges();

    for(j=0;j<N_Edges;j++)
    {
      joint = cell->GetJoint(j);
      if(joint->GetType() == IsoBoundEdge)
      {
        isojoint = (TIsoBoundEdge *)joint;
        k = isojoint->GetN_Vertices();
        Vertices = isojoint->GetVertices();
        FEId = VelocitySpace->GetFE2D(i, cell);
        FEDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
        m = FEDesc->GetN_JointDOF();
        if(m == k+2)
        {
          JointDOF = FEDesc->GetJointDOF(j);
          DOF =  VeloGlobalNumbers+VeloBeginIndex[i];
          for(l=0;l<k;l++)
          {
            m = DOF[JointDOF[l+1]];
            Vertices[l]->GetCoords(IsoX, IsoY);
            if(TDatabase::ParamDB->P5 > 0)
            {
              un = ValuesVX[m]*Nx[IIso+N_BoundaryNodes]
                  + ValuesVY[m]*Ny[IIso+N_BoundaryNodes];
              IsoX += dt*un*Nx[IIso+N_BoundaryNodes];
              IsoY += dt*un*Ny[IIso+N_BoundaryNodes];
              // cout << "U:   " << ValuesVX[m] << " " << ValuesVY[m] << endl;
              // cout << "N:   " << Nx[IIso+N_BoundaryNodes] << " "
              //                 << Ny[IIso+N_BoundaryNodes] << endl;
              // cout << "UNN: " << un*Nx[IIso+N_BoundaryNodes] << " "
              //                 << un*Ny[IIso+N_BoundaryNodes] << endl;
            }
            else
            {  
              IsoX += dt*ValuesVX[m];
              IsoY += dt*ValuesVY[m];
            }
//            if(IsoY<=0) IsoY = 1e-5;
//            if(fabs(IsoX)<1e-12) IsoX = 0.;
           Vertices[l]->SetCoords(IsoX, IsoY);
            IIso++;
          } // endfor l
        }
       else
        {
          // approximation order of isoparametric boundary and velocity
          // element must be the same
          Error("No match in isoparametric case" << endl);
          exit(-1);
        }
      } // endif
    } // endfor j  
   } // if(reparam)
  } // endfor i

   if(reparam)
    {
     N_ReParam++; 
     reparam = FALSE;
    }    
    
//   delete [] d;
  
  if(TDatabase::ParamDB->P5 > 0)
  { delete [] Nx;  delete [] Ny;}
  
} // MoveGrid_imping


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
//     cout<<i<<  " x : "<< x << " y : "<< y<<  " Angle of free Vertices "<<(180/Pi)*atan2(y,x)<<endl;
//    }
//   exit(0); 
}

 
void  GetMovingBoundData(TCollection *coll, int *N_MovVert, TBoundEdge *** &Bound_Joint, TVertex *** &MovBoundVert,
                         TIsoBoundEdge ** &Free_Joint, TBaseCell ** &Free_Cells, int ** &IsoCellEdgeNos, double x, double y)
{
 int i, j, k, l, m0, m1, m2, N_Cells, comp;
 int ORDER, VSP;
 
 double  TX[2], TY[2], x0, y0, x1, y1;
  
 TBaseCell *Me;
 TJoint *Joint;
 TBoundEdge *Solid_Joint;
 TBoundComp *BoundComp;  
 TVertex *temp_Mov;
 TBoundEdge *tempSlip_Joint;

  ORDER = 0;
  VSP = TDatabase::ParamDB->VELOCITY_SPACE;

  if (abs(VSP) > 20)
   {ORDER = abs(VSP) - 20;}
  else if ( abs(VSP) > 10)
    {ORDER = abs(VSP) - 10;}
  else ORDER = abs(VSP);
  
  
  N_Cells = coll->GetN_Cells();
  N_MovVert[0] = 0;
  N_MovVert[1] = 0;
  N_MovVert[2] = 0;
  
    for(j=0;j<N_Cells;j++)
      {
        Me = coll->GetCell(j);
        k = Me->GetN_Edges();
        for(l=0;l<k;l++)
        {
         Joint = Me->GetJoint(l);
          if(Joint->GetType() == BoundaryEdge)
            {
             Solid_Joint = (TBoundEdge *)Joint;
             BoundComp = Solid_Joint->GetBoundComp();
             comp=BoundComp->GetID();  
#ifdef  __WITHCYLINDERTICKNESS__             
             if(comp==4)
              { N_MovVert[1]++; }  
#else
             if(comp==3)
              { N_MovVert[1]++; }   
#endif
              
            }
           else if(Joint->GetType() == IsoBoundEdge)
            { N_MovVert[0]++;  }
            
          }// endfor l
        }// endfor j
  
//      for(i=0; i<2; i++)
//       cout<<"BDComp " << i << " N_Vert " << N_MovVert[i] << endl; 
// exit(0);


     // axial bound
     Bound_Joint[0] = new TBoundEdge* [N_MovVert[1]];
     MovBoundVert[1] = new TVertex* [N_MovVert[1]];

     // free bound
     Free_Joint = new TIsoBoundEdge*[N_MovVert[0]];
     MovBoundVert[0] = new TVertex*[N_MovVert[0]];
     Free_Cells = new TBaseCell*[N_MovVert[0]];
     IsoCellEdgeNos[0]  = new int [N_MovVert[0]];
     IsoCellEdgeNos[1]  = new int [N_MovVert[0]];
     
      m0 = 0;
      m1 = 0;
  
     for(j=0;j<N_Cells;j++)
      {
        Me = coll->GetCell(j);
        k = Me->GetN_Edges();
        for(l=0;l<k;l++)
        {
          Joint = Me->GetJoint(l);
          if(Joint->GetType() == BoundaryEdge)
          {
           Solid_Joint = (TBoundEdge *)Joint;
           BoundComp = Solid_Joint->GetBoundComp();
           comp=BoundComp->GetID();
#ifdef  __WITHCYLINDERTICKNESS__	   
           if(comp==4)
           {
             Bound_Joint[0][m1] = (TBoundEdge *)Joint;
             MovBoundVert[1][m1] = Me->GetVertex(l);
             m1++;
            }
#else
           if(comp==3)
           {
             Bound_Joint[0][m1] = (TBoundEdge *)Joint;
             MovBoundVert[1][m1] = Me->GetVertex(l);
             m1++;
            }
#endif
           }
         else if(Joint->GetType() == IsoBoundEdge)
           {
            Free_Joint[m0] = (TIsoBoundEdge *)Joint;
            MovBoundVert[0][m0] = Me->GetVertex(l);
            Free_Cells[m0] = Me;
            IsoCellEdgeNos[0][m0] = j;
            IsoCellEdgeNos[1][m0] = l;
            Me->GetVertex(l)->GetCoords(TX[0], TY[0]);
            Me->GetVertex((l+1) % k)->GetCoords(TX[1], TY[1]);
            //cout<< " x0: " << TX[0]<<" y0: " << TY[0]<< " x1: " << TX[1]<<" y1: " << TY[1]<< endl;
            Free_Joint[m0]->GeneratemidVert(ORDER-1, TX, TY);
            m0++;
           }
          }// endfor l
         }// endfor j    
 
//    exit(0);
         
  // sort         
    for(k=0;k<N_MovVert[1]-1;k++)
      {
      for(l=k+1;l<N_MovVert[1];l++)
       {MovBoundVert[1][k]->GetCoords(x0, y0);
	MovBoundVert[1][l]->GetCoords(x1, y1);
	if(y0 < y1)
	{
	  temp_Mov = MovBoundVert[1][k];
          MovBoundVert[1][k] = MovBoundVert[1][l];
          MovBoundVert[1][l] = temp_Mov;

	  tempSlip_Joint = Bound_Joint[0][k];
	  Bound_Joint[0][k] = Bound_Joint[0][l];
	  Bound_Joint[0][l] = tempSlip_Joint;
	 }
        }
       }         
//  for (k=0;k<N_MovVert[1];k++)
//       {
//        MovBoundVert[1][k]->GetCoords(x0, y0);
//        cout<< " x0 " << x0<<" y0 " << y0<<endl;
//        }
// exit(0);

//     for(k=0;k<N_MovVert[1]-1;k++)
//      {
//       for(l=k+1;l<N_MovVert[1];l++)
//        {MovBoundVert[1][k]->GetCoords(x0, y0);
// 	MovBoundVert[1][l]->GetCoords(x1, y1);
// 	if(y1 > y0)
// 	{
// 	  temp_Mov = MovBoundVert[1][k];
//           MovBoundVert[1][k] = MovBoundVert[1][l];
//           MovBoundVert[1][l] = temp_Mov;
// 
// 	  tempSlip_Joint = Bound_Joint[1][k];
// 	  Bound_Joint[1][k] = Bound_Joint[1][l];
// 	  Bound_Joint[1][l] = tempSlip_Joint;
// 	 }
//         }
//        }
//  for (k=0;k<N_MovVert[1];k++)
//       {
//        MovBoundVert[1][k]->GetCoords(x0, y0);
//        cout<< " x " << x0<<" y " << y0<<endl;
//        }
// exit(0);
   Sort_Imping(Free_Cells, MovBoundVert[0], IsoCellEdgeNos[0], IsoCellEdgeNos[1], N_MovVert[0], x, y);

}// GetMovingBoundData




void  Domain2DSurf_2Phase(TCollection *Coll, TDomain *SurfDomain, int **N_List)
  {
  int i, j, k, l, m, n,  m1, N_SurfCells, N_Cells, N_SurfVert, ID, N;
  int maxEpV, a, b, len1, len2, Neighb_tmp;
  int *Bd_Part, *Lines, N_G, *PointNeighb, CurrNeib,  Neib[2];
  int *Cell_No, *Joint_No;

  const int *TmpEV;
  TBaseCell *Me;
  TVertex **SurfVetrex;
  TJoint *Joint;
  TBaseCell  **SurfCellTree, *NullCell = NULL;
  boolean StartBD=FALSE, EndBD=FALSE;


  Cell_No = N_List[0];
  Joint_No = N_List[1];

  N_SurfCells = 0;
  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {
   Me = Coll->GetCell(i);
//    PhaseID = Me->GetPhase_ID();

//    if(PhaseID == 1) // outer phase cell is better for soluble surfactant coupling
   {
    k = Me->GetN_Edges();
    for(l=0;l<k;l++)
     {
       Joint = Me->GetJoint(l);
       if(Joint->GetType() == IsoBoundEdge   || Joint->GetType() == IsoInterfaceJoint
           || Joint->GetType() == InterfaceJoint  )
         N_SurfCells++;
     } // endfor l
   } //  if(PhaseID == 1)
  }// endfor i

//   OutPut("N_SurfCells: " << N_SurfCells << endl);

//   exit(0);
  
  Bd_Part = new int[N_SurfCells];
  Lines = new int[2*N_SurfCells];
  SurfVetrex = new TVertex*[2*N_SurfCells]; // maximum possible vertex

//   m = 0;
  m1 = 0;
  N    = 0; 
  N_SurfVert = 0;
  for(i=0;i<N_Cells;i++)
  {
   Me = Coll->GetCell(i);
//    PhaseID = Me->GetPhase_ID();

//    if(PhaseID == 1) // outer phase cell is better for soluble surfactant coupling
   {
    k = Me->GetN_Edges();
    for(l=0;l<k;l++)
     {
       Joint = Me->GetJoint(l);
       if(Joint->GetType() == IsoBoundEdge   || Joint->GetType() == IsoInterfaceJoint
           || Joint->GetType() == InterfaceJoint  )
        {
          Cell_No[N] = i;
          Joint_No[N] = l;
          Bd_Part[N++] =(((TBoundEdge *)Joint)->GetBoundComp())->GetID();
          Me->GetShapeDesc()->GetEdgeVertex(TmpEV);

           for(n=0;n<2;n++) // In 2D each edge contains only two vertex
            {
              ID = 0;
              for(m=0;m<N_SurfVert;m++)
              if( SurfVetrex[m]==Me->GetVertex(TmpEV[2*l + n ]) )
               {
//             vertex is already in the collection
                ID = 1;
                Lines[m1] = m;
                m1++;
                break; 
              }

              if(ID==0)
               {
                 Lines[m1] = N_SurfVert;
                 m1++;
                 SurfVetrex[N_SurfVert] = Me->GetVertex(TmpEV[2*l + n ]);
                 N_SurfVert++;
               }

           } //  for(n=0
        } // if(Joint->GetType()
     } // endfor l
   } // if(PhaseID == 1)
  }// endfor i

   OutPut("N_SurfVert: " << N_SurfVert << endl);

   SurfCellTree = new TBaseCell*[N_SurfCells];

   for (i=0;i<N_SurfCells;i++)
   {
    SurfCellTree[i] = new TMacroCell(TDatabase::RefDescDB[S_Line], 0);
    SurfCellTree[i]->SetVertex(0, SurfVetrex[ Lines[ 2*i       ]]);
    SurfCellTree[i]->SetVertex(1, SurfVetrex[ Lines[ 2*i + 1]]);

    ((TMacroCell *) SurfCellTree[i])->SetSubGridID(0);
    ((TBaseCell *) SurfCellTree[i])->SetBd_Part(Bd_Part[i] );
   }

   SurfDomain->SetTreeInfo(SurfCellTree, N_SurfCells);

   TDatabase::IteratorDB[It_EQ]->SetParam(SurfDomain);
   TDatabase::IteratorDB[It_LE]->SetParam(SurfDomain);
   TDatabase::IteratorDB[It_Finest]->SetParam(SurfDomain);
   TDatabase::IteratorDB[It_Between]->SetParam(SurfDomain);
   TDatabase::IteratorDB[It_OCAF]->SetParam(SurfDomain);

  // search neighbours
   maxEpV = 2; // max possible degree of each vertex
   N_G = N_SurfVert;
   maxEpV++;
   PointNeighb = new int[maxEpV * N_G];
   memset(PointNeighb, 0, maxEpV * N_G *SizeOfInt);

    // every row contains "maxEpV" columns
    // for every row at first colomn contains the number of cells containing this vertex
    // at further columns we set the index of the corresponding cells containing this vertex
   for(i=0;i<2*N_SurfCells;i++)
   {
     j = Lines[i]*maxEpV;
     PointNeighb[j]++;
     PointNeighb[j + PointNeighb[j]] = i / 2;
   }

  N_G = N_SurfCells;

  for(i=0;i<N_G;i++)
   {
    StartBD=FALSE;
    EndBD=FALSE;
    a =  Lines[2*i ];
    b =  Lines[2*i + 1];
    Neib[0] = -1;
    Neib[1] = -1;
    CurrNeib = 0;
    len1 = PointNeighb[a*maxEpV];
    len2 = PointNeighb[b*maxEpV];


    if(len2==1)
     {
      Neib[0] = PointNeighb[b*maxEpV + 1];
      Joint = new TJointEqN(SurfCellTree[Neib[0]]);
      SurfCellTree[Neib[0]]->SetJoint(1, Joint);
      cout<<" End cell " << Neib[0] <<" Neib[1] " << Neib[1] <<endl;
     }

  // find indexes of cells containing the current vertex (joint)
    if(len1==1)
      Neib[0] = PointNeighb[a*maxEpV + 1];
    else
      for(j=1;j<=len1;j++)
       {
        Neighb_tmp = PointNeighb[a*maxEpV + j];
        for (k=1;k<=len2;k++)
         {
          if(Neighb_tmp == PointNeighb[b*maxEpV + k])
           {
            Neib[0] = Neighb_tmp;

            if(j==1)
             Neib[1] = PointNeighb[a*maxEpV + j + 1];
            else if(j==2)
             Neib[1] = PointNeighb[a*maxEpV + j - 1];

            break;
           }
         }
       } //  for (j=1;j<=len1


    if(len1==2)
     {
      Joint = new TJointEqN(SurfCellTree[Neib[0]], SurfCellTree[Neib[1]]);

      SurfCellTree[Neib[0]]->SetJoint(0, Joint);
      SurfCellTree[Neib[1]]->SetJoint(1, Joint);
     }
    else if(len1==1)
     {
      Joint = new TJointEqN(SurfCellTree[Neib[0]]);
      SurfCellTree[Neib[0]]->SetJoint(0, Joint);
     cout<< " Start cell " << Neib[0] <<" Neib[1] " << Neib[1] <<endl;
     }

   } //  for(i=0;i<N_G;i++)

  delete [] SurfVetrex;
  delete [] PointNeighb;
  delete []  Bd_Part;
  delete []  Lines;
 }



void RemeshAxial3D(TDomain * &Domain, TFESpace2D ** &FESpaces_All,
                           TFEVectFunct2D ** &FEVectFuncts_All, TFEFunction2D ** &FEFunctions_All, int *N_MovVert, 
                           TBoundEdge *** &Bound_Joint, TVertex *** &MovBoundVert, TIsoBoundEdge ** &Free_Joint,
                           TBaseCell ** &Free_Cells, int ** &IsoCellEdgeNos, double ** &Sol_All,  double ** &Rhs_All, 
                           TSquareStructure2D ** &SquareStructure_All, TStructure2D ** &Structure_All,
                           TSquareMatrix2D ** &SqMat_All, TMatrix2D ** &Mat_All
#ifdef __SURFACT__   
                           , TDomain *SurfDomain, TFESpace1D **IFaceFeSpaces, int **N_List,
                           TFEFunction1D **IFaceFeFunct, TSquareMatrix1D **SqMat_IFace,
                           TSquareStructure1D **IFaceStruct, FE1D *FE1D_List
#endif                           
                     )
{
  int i, j, k, l, N_G, N_Cells, ORDER, VSP, N_DOF, N_ThermalDOF;
  int In_Index, CurrComp, Old_N_Cells, Old_N_RootCells, CurrVertex, N_Joints, N_Vertices, ID;
  int N_RootCells, *Triangles, *PointNeighb, maxEpV = 0, a, b, Neighb_tmp, Neib[2];
  int CurrNeib, len1, len2, pressure_space_code, N_U, N_P, N_Unknowns, comp;
  int *JointDOF, *DOF, *GlobalNumbers, *BeginIndex, N_refX, N;
  int N_TopBound_Vert, N_SlipBound_Vert, N_WallBound_Vert, N_Interf_Vertices, N_FreeBound_Vert;
  int N_IsurfactDOF, N_S, N_HorizBound_Vert;
  
  double area, *Coordinates, left, right, bottom, top, T_a, T_b, *sol, *ValuesU2;
  double x, y, y_begin, y_end, dx, y_AxialEnd, *S_BX, *S_BY;
  double CHAR_L=TDatabase::ParamDB->CHAR_L0, Wall_R;
   
  TBoundPart *BoundPart;
  TBoundComp *BoundComp;
  TBdLine *UpdateAxialBound;
  TCollection *coll, *Old_Coll, *mortarcoll = NULL;
  TBaseCell *cell, **CellTree, **Old_CellTree;
  TVertex **VertexDel, **NewVertices;
  TJoint *Joint;
  TFESpace2D *velocity_space, *pressure_space, *surfact_space;
  TFEVectFunct2D *u; 
  TFEFunction2D *p; 
  TBoundEdge *Solid_Joint;
  FE2D FeId;
  TFEDesc2D *FeDesc; 
    
  // strings
  char ReadinDat[] = "readin.dat";
  char CString[] = "C";
  char NameString[]  = "name";
  char UString[] = "u";
  char PString[] = "p";
  char WString[] = "w";
  char IFaceSString[] = "I_C";

  
#ifdef __SURFACT__   
  TCollection *SOld_Coll, *IFace_Coll; 
  TBaseCell **Surf_OldCellTree;
  TFEFunction2D *New_c, *SurfNew_c; 
  TFESpace2D *SurfSurfactSpace, *SurfactSpace;
     
  int N_IFaceCells, Old_S_N_RootCells;
  
  double *csol, *Surfcsol;
  
  TDatabase::IteratorDB[It_EQ]->SetParam(SurfDomain);
  TDatabase::IteratorDB[It_LE]->SetParam(SurfDomain);
  TDatabase::IteratorDB[It_Finest]->SetParam(SurfDomain);
  TDatabase::IteratorDB[It_Between]->SetParam(SurfDomain);
  TDatabase::IteratorDB[It_OCAF]->SetParam(SurfDomain);

  SOld_Coll = SurfDomain->GetCollection(It_Finest, 0);
  SurfDomain->GetTreeInfo(Surf_OldCellTree, Old_S_N_RootCells);

  // initialize iterators
  TDatabase::IteratorDB[It_EQ]->SetParam(Domain);
  TDatabase::IteratorDB[It_LE]->SetParam(Domain);
  TDatabase::IteratorDB[It_Finest]->SetParam(Domain);
  TDatabase::IteratorDB[It_Between]->SetParam(Domain);
  TDatabase::IteratorDB[It_OCAF]->SetParam(Domain);  
#endif     
  
  
  std::ostringstream opts;
  std::ostringstream os;
  os << " ";
  opts << " ";
   
// ======================================================================
// Triangular for grid generation begin
// ====================================================================== 
  boolean AllowEdgeRef = (boolean) TDatabase::ParamDB->MESHGEN_ALLOW_EDGE_REF;
  
  // free surface vertices
  S_BX = new double[N_MovVert[0]];
  S_BY = new double[N_MovVert[0]];
  
  
  for(i=0;i<N_MovVert[0];i++)
    MovBoundVert[0][i]->GetCoords(S_BX[i], S_BY[i]);

  
  
  struct triangulateio In, Out;
  
  BoundPart = Domain->GetBdPart(0);
//   UpdateFreeBound = (TBdCircle*)BoundPart->GetBdComp(0);
#ifdef  __WITHCYLINDERTICKNESS__  
  UpdateAxialBound = (TBdLine*)BoundPart->GetBdComp(4);
#else
  UpdateAxialBound = (TBdLine*)BoundPart->GetBdComp(3);
#endif
  
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

//   opts<<"a"<< area; // Imposes a maximum triangle area.
  opts<<'e'; // Outputs a list of edges of the triangulation
  opts<<'z'; // Numbers if items starting from 0
  //opts<<"VVVV"; // Gives detailed information about what Triangle is doing
  opts<<'Q'; // Supress all explanation of what Triangle is doing, unless an error occurs
//   opts<<'Y'; // Supress adding vertices on boundary edges
  opts<<ends;

//   N_FreeBound_Vert = int (TDatabase::ParamDB->P6);    //Freesurf except end point
  N_FreeBound_Vert = N_MovVert[0];
  N_TopBound_Vert = 10;
#ifdef  __WITHCYLINDERTICKNESS__    
  N_HorizBound_Vert = 5;
#else
  N_HorizBound_Vert = 0;// needed only if we consider the thickness of the pipe in meshing (needed for experimental comparison)
#endif  
  N_SlipBound_Vert = 25;     // Initially only three points on solid bound (except end point)
  N_WallBound_Vert =10;
  N_Interf_Vertices = N_FreeBound_Vert+N_WallBound_Vert+N_SlipBound_Vert + N_TopBound_Vert + N_HorizBound_Vert;
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
  
  for(i=0;i<N_FreeBound_Vert;i++) // without last point
   { 
      In.pointlist[2*In_Index] = S_BX[i];
      In.pointlist[2*In_Index+1] = S_BY[i];
//       cout<<" x : "<< S_BX[i] << " y : "<< S_BY[i]<<endl;
      //In.pointmarkerlist[In_Index] = CurrComp;
      In.segmentlist[2*In_Index] = In_Index;
      In.segmentlist[2*In_Index+1] = In_Index+1;
      if (AllowEdgeRef)
       { In.segmentmarkerlist[In_Index] = CurrComp; }
      else
       { In.segmentmarkerlist[In_Index] = 100000 + CurrComp; }

      In_Index++;  
   }
  CurrComp++;   
//      cout<<endl;
  
   Wall_R = 1./CHAR_L;
 #ifdef  __WITHCYLINDERTICKNESS__    
    /** horizontal BD (marker=1)*/
   x = TDatabase::ParamDB->P7; // last point from data file
   dx = (Wall_R - x)/(double)N_HorizBound_Vert;
   y=0.;
   for(i=0;i<N_HorizBound_Vert;i++) // without last point
    {
     In.pointlist[2*In_Index] = x;
     In.pointlist[2*In_Index+1] = y;
//      cout<<" x : "<< x << " y : "<< y<<endl;
     In.pointmarkerlist[In_Index] = CurrComp;
     In.segmentlist[2*In_Index] = In_Index;
     In.segmentlist[2*In_Index+1] = In_Index+1;
     In.segmentmarkerlist[In_Index] = CurrComp;
     In_Index++;
     x += dx;
    }
   CurrComp++;   
#endif  
//   cout << endl;
  
    x= Wall_R;
    y_begin = y = 0;
    y_end = TDatabase::ParamDB->P8/CHAR_L;  
    dx = y_end/(double)N_WallBound_Vert; 

   // points and segments on the horizontal boundary (marker=1)
   for(i=0;i<N_WallBound_Vert;i++) // without last point
    {
     In.pointlist[2*In_Index] = x;
     In.pointlist[2*In_Index+1] = y;
//      cout<<" x : "<< x << " y : "<< y<<endl;
     In.pointmarkerlist[In_Index] = CurrComp;
     In.segmentlist[2*In_Index] = In_Index;
     In.segmentlist[2*In_Index+1] = In_Index+1;
     In.segmentmarkerlist[In_Index] = CurrComp;
     In_Index++;
     y = y_begin + (double)(i+1)*dx;
    }
   CurrComp++;


//   cout << endl;
  
    x= Wall_R;
    y= TDatabase::ParamDB->P8/CHAR_L;
    dx = x/(double)N_TopBound_Vert; 

   // points and segments on the horizontal boundary (marker=1)
   for(i=0;i<N_TopBound_Vert;i++) // without last point
    {
     In.pointlist[2*In_Index] = x;
     In.pointlist[2*In_Index+1] = y;
//      cout<<" x : "<< x << " y : "<< y<<endl;
     In.pointmarkerlist[In_Index] = CurrComp;
     In.segmentlist[2*In_Index] = In_Index;
     In.segmentlist[2*In_Index+1] = In_Index+1;
     In.segmentmarkerlist[In_Index] = CurrComp;
     In_Index++;
     x -= dx;
    }
   CurrComp++;
//   cout << endl;
  

    x= 0.0;
    y_begin = y = TDatabase::ParamDB->P8/CHAR_L;
    y_AxialEnd = S_BY[0];  
    dx = (y_AxialEnd-y_begin)/(double)N_SlipBound_Vert; 

//     len1 = N_SlipBound_Vert/2;
//     
//    // points and segments on the horizontal boundary (marker=1)
//    for(i=0;i<len1;i++) // without last point
//     {
//      In.pointlist[2*In_Index] = x;
//      In.pointlist[2*In_Index+1] = y;
// //      cout<<" x : "<< x << " y : "<< y<<endl;
//      In.pointmarkerlist[In_Index] = CurrComp;
//      In.segmentlist[2*In_Index] = In_Index;
//      In.segmentlist[2*In_Index+1] = In_Index+1;
//      In.segmentmarkerlist[In_Index] = CurrComp;
//      In_Index++;
//      y = y_begin + double(i+1)*dx;
//     }
//    CurrComp++; 
//    cout<<endl;
    double BD_Y = y;
   // points and segments on the horizontal boundary (marker=1)
   for(i=0;i<N_SlipBound_Vert;i++) // without last point
    {
     In.pointlist[2*In_Index] = x;
     In.pointlist[2*In_Index+1] = y;
//      cout<<" x : "<< x << " y : "<< y<<endl;
     In.pointmarkerlist[In_Index] = CurrComp;
     In.segmentlist[2*In_Index] = In_Index;
     In.segmentlist[2*In_Index+1] = In_Index+1;
     In.segmentmarkerlist[In_Index] = CurrComp;
     In_Index++;
     y = y_begin + double(i+1)*dx;
    }
   CurrComp++;     
   
       
  In.segmentlist[2*(In.numberofsegments-1)+1] = 0;
  
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
    
  triangulate((char*)opts.str().c_str(), &In, &Out, (struct triangulateio *)NULL);    
 
  
  Old_Coll = Domain->GetCollection(It_Finest, 0);
  Old_N_Cells = Old_Coll->GetN_Cells();
  Domain->GetTreeInfo(Old_CellTree, Old_N_RootCells);  
  if(Old_N_Cells != Old_N_RootCells) 
   exit(-1);
  
  VertexDel = new TVertex*[3*Old_N_RootCells]; 
  CurrVertex = 0;
    
  for(i=0;i<Old_N_Cells;i++)
    {
      cell = Old_Coll->GetCell(i);
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
    
  N_RootCells = Out.numberoftriangles;   
  Coordinates = Out.pointlist;
  Triangles = Out.trianglelist;
  
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
  
  
  Domain->SetBoundBox(right-left,top-bottom);
  Domain->SetBoundBoxstart(left,bottom); 
  
  // Axial Bound startx, starty, x length and y length 
  UpdateAxialBound->SetParams(0., BD_Y, 0, y_AxialEnd - BD_Y);
//   cout <<  y_AxialEnd - BD_Y << " UpdateAxialBoundTop " << BD_Y << endl;  
  
//   UpdateAxialBoundTop->SetParams(0., y_begin, 0, BD_Y - y_begin);
//   cout << BD_Y - y_begin << " UpdateAxialBoundTop " << y_begin << endl;
 
//   cout <<  "N_RootCells " << N_RootCells << endl;
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

  for (i=0;i<N_G;i++)
    if (PointNeighb[i] > maxEpV) maxEpV = PointNeighb[i];

  delete [] PointNeighb;

  PointNeighb = new int[++maxEpV * N_G];

  memset(PointNeighb, 0, maxEpV * N_G *SizeOfInt);
  
  for(i=0;i<3*N_RootCells;i++)
  {
    j = Triangles[i]*maxEpV;
    PointNeighb[j]++;
    PointNeighb[j + PointNeighb[j]] = i/3;
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

   if(Out.edgemarkerlist[i]) // 0 for inner edges and Boundcomp+1 for Boundedge respect
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
          cout << " CurrComp " << CurrComp <<endl;
        //  exit(0);
       }

      if (CurrNeib == 2)    // 2 cells contain the current edge
        if(Domain->GetBdPart(0)->GetBdComp(CurrComp)->IsFreeBoundary())
          Joint = new TIsoInterfaceJoint(Domain->GetBdPart(0)->GetBdComp(CurrComp), T_a, T_b, 
                                         CellTree[Neib[0]], CellTree[Neib[1]]);
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
  
  } //  for (i=0;i<N_G;i++)

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

  
//   TriaReMeshGen(Domain, N_MovVert[0], S_BX, S_BY);
 
// ======================================================================
// Triangular for grid generation end
// ======================================================================
 
//    // write grid into an Postscript file
//    os.seekp(std::ios::beg);
//    os << "Domain" << ".ps" << ends;
//    Domain->PS(os.str().c_str(),It_Finest,0);
//    exit(0);
  
//======================================================================
// construct all finite element spaces
//======================================================================
  coll=Domain->GetCollection(It_Finest, 0);
  N_Cells = coll->GetN_Cells();
 
#ifdef __SURFACT__  
  TDatabase::IteratorDB[It_EQ]->SetParam(SurfDomain);
  TDatabase::IteratorDB[It_LE]->SetParam(SurfDomain);
  TDatabase::IteratorDB[It_Finest]->SetParam(SurfDomain);
  TDatabase::IteratorDB[It_Between]->SetParam(SurfDomain);
  TDatabase::IteratorDB[It_OCAF]->SetParam(SurfDomain);
    
  Domain2DSurf_2Phase(coll, SurfDomain, N_List);

  IFace_Coll = SurfDomain->GetCollection(It_Finest, 0);
  N_IFaceCells= IFace_Coll->GetN_Cells();
  cout<< " N_IFaceCells " << N_IFaceCells <<endl; 
  
//   delete [] FE1D_List; deleted in fespace1D
  FE1D_List = new FE1D[N_IFaceCells];
  for(j=0;j<N_IFaceCells;j++)
   FE1D_List[j] = FE1D(TDatabase::ParamDB->ANSATZ_ORDER);

  TDatabase::IteratorDB[It_EQ]->SetParam(Domain);
  TDatabase::IteratorDB[It_LE]->SetParam(Domain);
  TDatabase::IteratorDB[It_Finest]->SetParam(Domain);
  TDatabase::IteratorDB[It_Between]->SetParam(Domain);
  TDatabase::IteratorDB[It_OCAF]->SetParam(Domain);  
#endif     
  
  ORDER = 0;
  VSP = TDatabase::ParamDB->VELOCITY_SPACE;

  if (abs(VSP) > 20)
   {ORDER = abs(VSP) - 20;}
  else if ( abs(VSP) > 10)
    {ORDER = abs(VSP) - 10;}
  else ORDER = abs(VSP);

//======================================================================
// construct all finite element spaces
//======================================================================
  // get velocity and pressure spacess
  GetVelocityAndPressureSpace(coll,BoundCondition,
                              mortarcoll, velocity_space,
                              pressure_space, &pressure_space_code,
                              TDatabase::ParamDB->VELOCITY_SPACE,
                              TDatabase::ParamDB->PRESSURE_SPACE);
  N_U = velocity_space->GetN_DegreesOfFreedom();
  N_P = pressure_space->GetN_DegreesOfFreedom();
  GlobalNumbers = velocity_space->GetGlobalNumbers();
  BeginIndex = velocity_space->GetBeginIndex();

  // mesh velocity space 
  delete FESpaces_All[2];
  FESpaces_All[2] = new TFESpace2D(coll, NameString, WString, GridBoundCondition, 1, NULL);
  N_G = FESpaces_All[2]->GetN_DegreesOfFreedom();
  
  // surfact space
#ifdef __SURFACT__    
  SurfactSpace = new TFESpace2D(coll, NameString, CString, SurfactBoundCondition,
                                  TDatabase::ParamDB->ANSATZ_ORDER, NULL); 
   N_ThermalDOF = SurfactSpace->GetN_DegreesOfFreedom();
   OutPut("N_SurfactDOF    : "<< setw(10) << N_ThermalDOF  << endl); 
  
   delete IFaceFeSpaces[0];
   IFaceFeSpaces[0] = new TFESpace1D(IFace_Coll, IFaceSString, IFaceSString, FE1D_List);
 
   N_IsurfactDOF = IFaceFeSpaces[0]->GetN_DegreesOfFreedom();
   OutPut("N_SurfSurfactDOF    : "<< setw(10) << N_IsurfactDOF  << endl);   

   SurfSurfactSpace =  new TFESpace2D(coll, NameString, IFaceSString, SurfactBoundCondition,
                                     TDatabase::ParamDB->ANSATZ_ORDER, NULL);

   N_S =  SurfSurfactSpace->GetN_DegreesOfFreedom();      
#endif  

//===========================================================================
  delete [] Bound_Joint[0];
  delete [] Bound_Joint[1];                      

  delete [] MovBoundVert[0];
  delete [] MovBoundVert[1];   
  delete [] MovBoundVert[2];  

  delete [] Free_Joint;    
  delete [] Free_Cells;  
  delete [] IsoCellEdgeNos[0];  
  delete [] IsoCellEdgeNos[1];  
    
 
  GetMovingBoundData(coll, N_MovVert, Bound_Joint, MovBoundVert, Free_Joint,
                     Free_Cells, IsoCellEdgeNos, S_BX[0], S_BY[0]);
    
//======================================================================
// construct all finite element functions
//======================================================================
  N_Unknowns = 2*N_U + N_P;
  
  OutPut("dof velocity : "<< setw(10) << 2* N_U << endl);
  OutPut("dof pressure : "<< setw(10) << N_P << endl);
  OutPut("dof all      : "<<  setw(10) << N_Unknowns  << endl);  
  
  delete [] Rhs_All[0];
  Rhs_All[0] = new double[N_Unknowns];
  sol = new double[N_Unknowns];
  memset(sol, 0, N_Unknowns*SizeOfDouble);
 
  u =  new TFEVectFunct2D(velocity_space, UString, UString, sol, N_U, 2);
  p = new TFEFunction2D(pressure_space, PString,  PString,  sol+2*N_U, N_P);

  u->Interpolate(FEVectFuncts_All[0]);

//====================================================================== 
// impose no-penetration condition 
//====================================================================== 
  ValuesU2 = sol + N_U;
  
  for(i=0;i<N_Cells;i++)
   {
    cell = coll->GetCell(i);
    k = cell->GetN_Edges();
    
    for(l=0;l<k;l++)
     {
      Joint = cell->GetJoint(l);    
      
      if(Joint->GetType() == BoundaryEdge)
        {
         FeId = velocity_space->GetFE2D(i, cell);  
         FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FeId);
         JointDOF = FeDesc->GetJointDOF(l); 
         N_DOF = FeDesc->GetN_JointDOF();
         DOF = GlobalNumbers+BeginIndex[i];
    
         Solid_Joint = (TBoundEdge *)Joint;
         BoundComp = Solid_Joint->GetBoundComp();
         comp=BoundComp->GetID();

         if(comp==1 || comp==3 || comp==4) // side walls
          {
           for(j=0;j<N_DOF;j++)
            sol[DOF[JointDOF[j]]] = 0.;  
          }
 
        }// if(Joint->GetType() == Bo
     } // for(l=0;l<k;l++)
   } //  for(i=0;i<N_Cells;i++)   
//======================================================================   
// no need to interpolate pressure for DirectSolvers !!!
//======================================================================    
  delete FESpaces_All[0];
  delete FESpaces_All[1]; 
 
  FESpaces_All[0] = velocity_space;
  FESpaces_All[1] = pressure_space;    

  delete [] Sol_All[0];
  Sol_All[0] = sol;
 
  delete FEVectFuncts_All[0]; 
  FEVectFuncts_All[0] = u;
 
  FEFunctions_All[0] = FEVectFuncts_All[0]->GetComponent(0);
  FEFunctions_All[1] = FEVectFuncts_All[0]->GetComponent(1); 

  delete FEFunctions_All[2];
  FEFunctions_All[2] = p;
  
//======================================================================
// grid space finite element functions
//====================================================================== 
  delete [] Sol_All[1];
  Sol_All[1] = new double[2*N_G];
  delete [] Rhs_All[1]; 
  Rhs_All[1] = new double[2*N_G];  
  
  memset(Sol_All[1], 0, 2*N_G*SizeOfDouble); 
  delete FEVectFuncts_All[1];
  FEVectFuncts_All[1]  = new TFEVectFunct2D(FESpaces_All[2], WString, WString, Sol_All[1], N_G, 2);
  
  FEFunctions_All[3] = FEVectFuncts_All[1]->GetComponent(0);
  FEFunctions_All[4] = FEVectFuncts_All[1]->GetComponent(1); 

    
//======================================================================
// surfact space finite element functions
//======================================================================    
#ifdef __SURFACT__  
  csol  = new double[N_ThermalDOF];
  memset(csol, 0, N_ThermalDOF*SizeOfDouble);
  New_c = new TFEFunction2D(SurfactSpace, CString,  CString, csol, N_ThermalDOF);
    
// #ifdef  __SOLUBLE__  
  New_c->Interpolate(FEFunctions_All[5]);  
// #endif  
  
  delete [] Rhs_All[2];
  Rhs_All[2] = new double[N_ThermalDOF]; 
  memset(Rhs_All[2], 0, N_ThermalDOF*SizeOfDouble);
    
  delete [] Sol_All[2];
  Sol_All[2] = csol;
  
  // surfact fefunction
  delete FEFunctions_All[5];
  FEFunctions_All[5] =  New_c;
  
  delete FESpaces_All[3];
  FESpaces_All[3] = SurfactSpace;
  
   // surf surfact
  delete IFaceFeFunct[0];
  delete [] Sol_All[3];
  delete [] Rhs_All[3]; 
  
  Sol_All[3] = new double[N_IsurfactDOF]; 
  Rhs_All[3] = new double[N_IsurfactDOF]; 
  memset(Rhs_All[3], 0, N_IsurfactDOF*SizeOfDouble);

  
  IFaceFeFunct[0] = new TFEFunction1D(IFaceFeSpaces[0], IFaceSString, IFaceSString, Sol_All[3], N_IsurfactDOF);
  
  Surfcsol  =  new double[N_S];
  SurfNew_c  = new TFEFunction2D(SurfSurfactSpace, IFaceSString,  IFaceSString, Surfcsol, N_S); 
  SurfNew_c->Interpolate(FEFunctions_All[6]);   
 
  delete [] Sol_All[4];
  delete FEFunctions_All[6];
  delete FESpaces_All[4];     
  
  Sol_All[4] = Surfcsol;  
  FEFunctions_All[6] = SurfNew_c;  
  FESpaces_All[4] = SurfSurfactSpace;  
#endif

// ======================================================================
// allocate memory for all matrices
// ======================================================================  
  delete Structure_All[0];  delete Structure_All[1];
  
  Structure_All[0] = new TStructure2D(FESpaces_All[1], FESpaces_All[0]); // B
  Structure_All[1] = new TStructure2D(FESpaces_All[0], FESpaces_All[1]); // BT

  delete SquareStructure_All[0]; delete SquareStructure_All[1];
#ifdef __SURFACT__    
  delete SquareStructure_All[2];
  
  /* interface surfactant matrices */
  delete IFaceStruct[0];
#endif 

  //velo 
  SquareStructure_All[0] = new TSquareStructure2D(FESpaces_All[0]);  
  SquareStructure_All[0]->Sort();  

  // grid 
  SquareStructure_All[1] = new TSquareStructure2D(FESpaces_All[2]); 
  SquareStructure_All[1]->Sort();  

  //surfact
#ifdef __SURFACT__    
  SquareStructure_All[2] = new TSquareStructure2D(FESpaces_All[3]);
  SquareStructure_All[2]->Sort(); 
  
  IFaceStruct[0] = new TSquareStructure1D(IFaceFeSpaces[0]);
  IFaceStruct[0]->Sort();    
#endif
  // u
  for(i=0; i<10; i++)
   {
    delete SqMat_All[i];
    SqMat_All[i] = new TSquareMatrix2D(SquareStructure_All[0]);
   }

  // B
  for(i=0; i<2; i++)
   {
    delete Mat_All[i];
    Mat_All[i] = new TMatrix2D(Structure_All[0]);   
   }  

  // BT
  for(i=2; i<4; i++)
   {
    delete Mat_All[i];
    Mat_All[i] = new TMatrix2D(Structure_All[1]);   
   }  

  // for mesh
  for(i=10; i<14; i++)
   {
    delete SqMat_All[i];
    SqMat_All[i] = new TSquareMatrix2D(SquareStructure_All[1]);  
   }  
  
  // for heat
#ifdef __SURFACT__    
  for(i=14; i<16; i++)
   {
    delete SqMat_All[i];
    SqMat_All[i] = new TSquareMatrix2D(SquareStructure_All[2]);  
   }  
   
 /* interface surfactant matrices */
  for(i=0; i<2; i++)
   {
    delete SqMat_IFace[i];
    SqMat_IFace[i] = new TSquareMatrix1D(IFaceStruct[0]); // Gamma_A,Gamma_M   
   }    
#endif   

// ======================================================================
// delete old mesh
// ======================================================================  
    
  for(i=0;i<CurrVertex;i++)
    delete VertexDel[i];

  delete [] VertexDel;
  OutPut(CurrVertex<<" vertices were deleted"<<endl);
  
   // remove all existing cells and joints
  for(i=0;i<Old_N_RootCells;i++)
     delete (TGridCell*)Old_CellTree[i];
  OutPut(Old_N_RootCells<<" cells were deleted"<<endl);
   delete [] Old_CellTree;

  delete Old_Coll;

  delete [] S_BX;
  delete [] S_BY;
}// RemeshAxial3D 







