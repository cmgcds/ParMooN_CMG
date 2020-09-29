#include <TimeConvDiff2D.h>
#include <MacroCell.h>

void ExampleFile()
{
  OutPut("Example: 2D_1D_ADI_UnitSrq.h" << endl) ;

  TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD=1;
  #define __SIMPATURS__
}

void InitialU1(double x, double y, double *values)
{

  values[0] = 0.;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}


void InitialU2(double x, double y, double *values)
{
  values[0] = 0.;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}


// exact solution
void Exact(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}


void InitialCondition_Heat(double x, double y, double *values)
{
  double val;
  double T_D = TDatabase::ParamDB->REACTOR_P23;
  double T_W = TDatabase::ParamDB->REACTOR_P24;

  val = T_W/T_D;

  values[0] = val;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

// void Exact_Psd_Intl(double x, double y, double L, double *values)
void Exact_Psd_Intl(int N_Coord, double *X, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}
// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition_Heat(int i, double t, BoundCond &cond)
{

  switch(i)
  {
   case 0:
   case 2:
   case 4:
     cond = DIRICHLET;
   break;

   case 1:
   case 3: 
   case 5:
     cond = NEUMANN;
   break;

   default: 
      cout << "wrong boundary part number "  << i << endl;
      exit(0);
   break;
  }

}

void BoundValue_Heat(int BdComp, double Param, double &value)
{
  double val, dt, t,  y;
  double T_D = TDatabase::ParamDB->REACTOR_P23;
  double T_W = TDatabase::ParamDB->REACTOR_P24;

  val = T_W/T_D;

  switch(BdComp)
  {
   case 1:
   case 3:
   case 5:
     value=0.;
   break;

  case 0:
  case 2:
        value=val;
  break;
  

  case 4:
     y = 1 -  Param;
     dt = 1. - val;
//  cout << y << " t " << Param << endl;
     if( y>0 && y<(1./8.))
      {
	t = -6. + 12.*y*8.*(1. + 1./8. - y);
	t = val + dt*(1./(1. + exp(-t)));
//              cout << y << " t " << t << endl;
	value = t;  
      }
     else if(y<1. &&   y>(1. - 1./8.))
      {
	t = 6. + 12.*(1. - 1./8. - y)*8.*(1. + 1. -y);
	t = 1./(1. + exp(-t));
	t = val + dt*(1./(1. + exp(-t)));
//          cout << y << " t " << t << endl;
	value = t;      
      }
      else if(y==0. || y==1.)
      {
       value = val;  
      }   
      else
      {
	value = 1.;  
      }
     
  break;

  default:
     cout << "wrong boundary part number "  << BdComp << endl;
    break;
  }
}

// ========================================================================
// for scalar equations
// ========================================================================
void BilinearCoeffs_Heat(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  double eps;
  int i;
  double *coeff, *param;
  double x, y;

  if(TDatabase::ParamDB->REACTOR_P0)
    eps = 1.0/TDatabase::ParamDB->REACTOR_P0;
  else
    eps = 0.;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    x = X[i];
    y = Y[i];

    coeff[0] = eps;
    if(TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD)
     {
      coeff[1] = param[0];  // u1
      coeff[2] = param[1];  // u2
//       cout<< "coeff[0] eps " << eps << " u1 " << param[0]  << "u2 " << param[1] << endl;
     }
    else
     {
    coeff[1] = 0;  // u1
    coeff[2] = 0;  // u2
     }
    coeff[3] = 0.;

    coeff[4] = 0.; // f
  }
}



void InitialCondition_Conc(double x, double y, double *values)
{

  values[0] = 1.;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}


// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition_Conc(int i, double t, BoundCond &cond)
{
  switch(i)
  {
   case 0:
   case 1:
   case 2:
   case 3:
   case 5:
     cond = NEUMANN;
   break;

   case 4:
     cond = DIRICHLET;
   break;

   default: 
      cout << "wrong boundary part number "  << i << endl;
      exit(0);
   break;
  }
}

void BoundValue_Conc(int BdComp, double Param, double &value)
{
   double t,  y;
   
 switch(BdComp)
  {
  case 0:
  case 1:
  case 2:
  case 3:
  case 5:
     value=0.;
  break;

  case 4:
//      y = 1 -  Param;
// 
//      if( y>0 && y<(1./8.))
//       {
// 	t = -6. + 12.*y*8.*(1. + 1./8. - y);
// 	t = 1./(1. + exp(-t));
// //              cout << y << " t " << t << endl;
// 	value = t;  
//       }
//       else if(y<1. &&   y>(1. - 1./8.))
//       {
// 	t = 6. + 12.*(1. - 1./8. - y)*8.*(1. + 1. -y);
// 	t = 1./(1. + exp(-t));
// //          cout << y << " t " << t << endl;
// 	value = t;      
//       }
//       else if(y==0. || y==1.)
//       {
// 	value = 0.;  
//       }       
//       else
      {
	value = 1.;  
      }
   
      break;

      default: cout << "wrong boundary part number " << BdComp << endl;
    break;
  }
}


// ========================================================================
// for scalar equations
// ========================================================================
void BilinearCoeffs_Conc(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  double eps;
//   double a=0.2, b=0., c=0.;
  int i;
  double *coeff, *param;
  double x, y;

  if(TDatabase::ParamDB->REACTOR_P1)
    eps = 1.0/TDatabase::ParamDB->REACTOR_P1;
  else
    eps = 0.;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    x = X[i];
    y = Y[i];

    coeff[0] = eps;
    if(TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD)
     {
      coeff[1] = param[0];  // u1
      coeff[2] = param[1];  // u2
     }
    else
     {
      coeff[1] = 0;  // u1
      coeff[2] = 0;  // u2
     }
    coeff[3] = 0.;

    coeff[4] = 0.; // f
  }
}

void InitialCondition_Psd(double x, double y, double *values)
{
  values[0] = 0.;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition_Psd(int i, double t, BoundCond &cond)
{
  switch(i)
  {
   case 0:
   case 1:
   case 2:
   case 3:
   case 5:
     cond = NEUMANN;
   break;

   case 4:
     cond = DIRICHLET;
   break;

   default: 
      cout << "wrong boundary part number " << i << endl;
      exit(0);
   break;
  }
}

// void BoundValue_Psd(int BdComp, double Param, double &value)
// {
//   
//  double t, L = TDatabase::ParamDB->REACTOR_P29;
//   
//   
// //   if(L < 0.025) // condition at L_Min
//   if(L == 0 ) // condition at L_Min   
//    {
//     switch(BdComp)
//      {
//       case 0:
//       case 1:
//       case 2:
//       case 3:
//       case 5:
// 	value = 0.;
//       break;
//       case 4:
// 	value = 1.;
// //        value = exp(-L/0.0025);
//        
//       break;
//       
//       default: cout << "wrong boundary part number "  << BdComp  << endl;
// 	break;
//       }     
//    }
//   else if(L == TDatabase::ParamDB->REACTOR_P13) // condition at L_Max
//    {
//     switch(BdComp)
//      {
//       case 0:
//       case 1:
//       case 2:
//       case 3:
//       case 5:
// 	value = 0.;
//       break;
//       case 4:
// 	value = 0.;
//       break;
//       
//       default: cout << "wrong boundary part number "  << BdComp  << endl;
// 	break;
//       }     
//      }
//   else
//    {     
//     switch(BdComp)
//      {
//       case 0:
//       case 1:
//       case 2:
//       case 3:
//       case 5:
// 	value = 0.;
//       break;
//       case 4:
// 	value = 0.;
//       break;
//       
//       default: cout << "wrong boundary part number "  << BdComp  << endl;
// 	break;
//       }     
//      }
// 
//  }
//     

void BoundValue_Psd(int BdComp, double Param, double &value)
{  
 double t, L = TDatabase::ParamDB->REACTOR_P29;
 
    switch(BdComp)
     {
      case 0:
      case 1:
      case 2:
      case 3:
      case 5:
	value = 0.;
      break;
      case 4:
	if(L <0.025 && L>0.005)
	 {
	  t = (L - 0.005)/0.02;
	  value = 4.*t*(1.-t);
	 }
	else
	 value = 0.;
      break;
      
      default: cout << "wrong boundary part number "  << BdComp  << endl;
	break;
      }    

}

void BoundCondition_LminLMax(BoundCond &cond_Lmin, BoundCond &cond_Lmax)
{
//   cond_Lmin = DIRICHLET;
  cond_Lmin = NEUMANN;
  cond_Lmax = NEUMANN;
}


void BoundValue_LMin(double x, double y,  double *values)
 {
   double t;
  
   // no nucleation, so f = 0 at L=0
   
//  if((x==0) && y>(1./3) &&  y<(2./3))
//   {
//   if((x==0) && y>(1./3) &&  y<(1./3. + 1./24.))
//    {
//     t = -6. + 12.*(y - 1./3.)*24.*(1. + 1./3. + 1./24. - y);
//     t = 1./(1. + exp(-t));
// //              cout << y << " t " << t << endl;
//     values[0] = t;
//     values[1] = 0.;
//     values[2] = 0.;    
//    }
//   else if((x==0) && y>(2./3 - 1./24.)  &&  y<(2./3))
//    {
//      t = 6. + 12.*(2./3. - 1./24. - y)*24.*(1. + 2./3. -y);
//      t = 1./(1. + exp(-t));
// //          cout << y << " t " << t << endl;
//     values[0] = t;
//     values[1] = 0.;
//     values[2] = 0.;    
//    }  
//   else
//    {
//     values[0] = 1.;
//     values[1] = 0.;
//     values[2] = 0.;
//    }
//   }
//  else
  {
   values[0] = 0.;
   values[1] = 0.;
   values[2] = 0.;
  } 

 }


void BoundValue_LMax(double x, double y,  double *values)
 {
  values[0] = 0.;
  values[1] = 0.;
  values[2] = 0.;  
 }


// ========================================================================
// for scalar equations
// ========================================================================
void BilinearCoeffs_Psd(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  double eps;
  int i;
  double *coeff, *param;
  double x, y;

  if(TDatabase::ParamDB->REACTOR_P2)
    eps = 1.0/TDatabase::ParamDB->REACTOR_P2;
  else
    eps = 0.;


  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    x = X[i];
    y = Y[i];

    coeff[0] = eps;
    if(TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD)
     {
      coeff[1] = param[0];  // u1
      coeff[2] = param[1];  // u2
      //cout<< "coeff[0] eps " << eps << " u1 " << param[0]  << "u2 " << param[1] << endl;

     }
    else
     {
    coeff[1] = 0;  // u1
    coeff[2] = 0;  // u2
     }
    coeff[3] = 0.;

    coeff[4] = 0.; // f
  }
}


// // initial conditon
// void PsdAtLmin(double x, double y, double z, double &value)
// {
// //   //testing the implementation
// //   double t = TDatabase::TimeDB->CURRENTTIME;
// //   value = x*t/5.;
// 
// //   value = value:
// }


// initial conditon
// void InitialCondition_Psd_Intl(double x, double y, double z, double *values)
void InitialCondition_Psd_Intl(int N_Coord, double *X, double *values)
{

  //testing the implementation
//   double t = TDatabase::TimeDB->CURRENTTIME;
//    if(z < 0.025) 
//     values[0] = exp(-z/0.0025);  
//    else
     values[0] = 0;

//   values[0] = 0;
//    cout<< " test example "<<  values[0] <<endl;
}

// void BilinearCoeffs_Psd_Intl(int n_points, double *X, double *Y, double *Z,
//         double **parameters, double **coeffs)
void BilinearCoeffs_Psd_Intl(int n_points, int N_Dim, double **Coords, double **parameters, double **coeffs)
{
  int i;
  double eps, *coeff;                                  // *param;
  double x, y, z, c, a[3], b[3], s[3], h;
  double t = TDatabase::TimeDB->CURRENTTIME;

  double *X, *Y, *Z;
  
   X = coeffs[0];
   Y = coeffs[1];
   Z = coeffs[2];   
  
  b[0] = -1e8;// negative, so that C will be taken from the PBS growth term
  c = 0; 

  if(TDatabase::ParamDB->REACTOR_P3)
    eps = 1.0/TDatabase::ParamDB->REACTOR_P3;
  else
    eps = 0.;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    // param = parameters[i];

    x = X[i];
    y = Y[i];
    z = Z[i];

    // diffusion
    coeff[0] = eps;
    // convection in z direction
    coeff[1] = b[0];
    // reaction term
    coeff[2] = c;
    coeff[3] = 0.;
  }
}

//input and output : dimless
double GetCsat(double C, double T)
{
 double val;
 double T_ref = TDatabase::ParamDB->REACTOR_P23;
 double C_ref = TDatabase::ParamDB->REACTOR_P25;
  
//   val = (1.3045*(T_ref*T - 273.15) + 35.3642)/(C_ref*TDatabase::ParamDB->UREA_m_mol); 
   val = (1.3045*(T_ref*T - 273.15) + 35.3642)/(C_ref);  
  return(val);
}


double GetGrowth(double C, double T)
{
  double PbeGC = TDatabase::ParamDB->REACTOR_P27;
  double G, C_Sat, Mode = TDatabase::ParamDB->PBE_P0;

   C_Sat = GetCsat(C, T);
//             cout <<   "C: " <<C <<" C_Sat: " << C_Sat <<endl; 
    // growth rate based on concentration
    if(C>C_Sat && (Mode==1 || Mode==4|| Mode==6) )
     { 
      G = PbeGC*sqrt((C - C_Sat )/C_Sat);       
  
      if(G<0.) G = 0.;
     }
    else
     { G = 0.; }

//         if(fabs(G)>1e-4)
//          cout << "G: " << G << "C: " <<C <<"T: " << T <<endl;  
     
 return G;
 }
 
 

 // input and output in dimless form
void Get_GrowthAndB_Nuc(int N_Inputs, double *Inn, double *Out) 
{
//   double PbeGC = TDatabase::ParamDB->REACTOR_P27;
//   double C_Sat;
//   double T_ref = TDatabase::ParamDB->REACTOR_P23;
//   double C_ref = TDatabase::ParamDB->REACTOR_P25;
//   double alpha =  TDatabase::ParamDB->UREA_alfa_nuc;
//   double beta =  TDatabase::ParamDB->UREA_beta_nuc;
//   double f_ref = TDatabase::ParamDB->UREA_f_infty;
  double C, T, G, B_Nuc;
 
  if(N_Inputs!=2)
  {
    printf("N_Inputs != 2, needed Temp and conc to Get_GrowthAndB_Nuc: %d \n",N_Inputs); 
#ifdef _MPI
     MPI_Finalize();
#endif       
    exit(0);        
  }
  
  C = Inn[0];
  T = Inn[1];
  
    // growth rate based on concentration
    G = GetGrowth(C, T); 
    
    // set constant
    B_Nuc = TDatabase::ParamDB->UREA_alfa_nuc;  
    if( G>0.)
     { B_Nuc /= (G*TDatabase::ParamDB->UREA_f_infty); }
    else
     { B_Nuc = 0.; } // no growth      
    
 
//     B_Nuc = C/C_Sat;
//   
//     if (B_Nuc != 0)
//      {
//       B_Nuc=log(B_Nuc);
//       B_Nuc *= B_Nuc;
//       B_Nuc = -beta/B_Nuc;
//       B_Nuc = alpha*exp(B_Nuc);  
//       
//       // set constant
// //       B_Nuc = TDatabase::ParamDB->UREA_alfa_nuc;    
// 
//       if( G>0.)
//        { B_Nuc /= (f_ref*G); }
//       else
//        { B_Nuc /= f_ref; } // no growth
// 
//        if(B_Nuc>1.)
//          B_Nuc = 1.;
// 
//      }
//     else
//      { B_Nuc = 0.0; } 
     
    Out[0] = G;
    Out[1] = B_Nuc;       
    
//          cout << "G: " << G << " B_Nuc: " <<B_Nuc<< "  UREA_alfa_nuc: " <<TDatabase::ParamDB->UREA_alfa_nuc <<endl;      
}


void GetExampleFileData(BoundCondFunct2D **BoundaryConditions, BoundValueFunct2D **BoundValues, 
                        DoubleFunct2D **InitiaValues, CoeffFct2D **BilinearCoeffs, 
                        int &N_PBEqns, int &N_IndepntScalarEqns, int *Disctypes)
{

//   TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD=1;


   N_IndepntScalarEqns = 2;
   N_PBEqns = 1;
   #define __PBS__

   BilinearCoeffs[0] = BilinearCoeffs_Heat;
   BilinearCoeffs[1] = BilinearCoeffs_Conc;
   BilinearCoeffs[2] = BilinearCoeffs_Psd;

   BoundaryConditions[0] = BoundCondition_Heat;
   BoundaryConditions[1] = BoundCondition_Conc;
   BoundaryConditions[2] = BoundCondition_Psd;

   BoundValues[0] = BoundValue_Heat;
   BoundValues[1] = BoundValue_Conc;
   BoundValues[2] = BoundValue_Psd;

   InitiaValues[0] = InitialCondition_Heat;
   InitiaValues[1] = InitialCondition_Conc;
   InitiaValues[2] = InitialCondition_Psd;

   Disctypes[0] = GALERKIN;
   Disctypes[1] = GALERKIN;
//    Disctypes[2] = SDFEM;
   Disctypes[2] = GALERKIN;
}


void Generate1DMesh(TDomain *Domain, double Start, double End, int &N_Cells)
{
  int i, j, N_V;
  int *Lines;
  double len, h, x, y, *X;
  TVertex **Vetrex;
  TJoint *Joint;
  TBaseCell  **CellTree;

  N_V = N_Cells+1;
  X = new double[N_V];

  for(i=0; i<N_V; i++)
   X[i] = 1. + (1. - Start)*(tanh(2.75*(double(i)/double(N_Cells) - 1.)))/tanh(2.75);

//   h = (End - Start)/N_Cells;
//   for(i=1; i<N_V; i++)
//    X[i] =  h*(double)i;

//   X[0] = Start;
//   X[N_V-1] = End;
// 
// N_V = 88;
// N_Cells = N_V-1;
// 
// double X[88] = {0.001250000000000 ,
// 0.002110731579413 ,
// 0.002563914004088 ,
// 0.003166693847592 ,
// 0.003948458933296 ,
// 0.004948299443979 ,
// 0.006217671606001 ,
// 0.007823152131089 ,
// 0.009849848192817 ,
// 0.011272700427129 ,
// 0.012405802335839 ,
// 0.013362834027725 ,
// 0.014199484324047 ,
// 0.014947689230070 ,
// 0.015627666201516 ,
// 0.016833818754898 ,
// 0.017888194844297 ,
// 0.018831072546427 ,
// 0.019687946143196 ,
// 0.021207835206421 ,
// 0.022536431448077 ,
// 0.023724508066731 ,
// 0.024804199643075 ,
// 0.026719286105474 ,
// 0.028393316830700 ,
// 0.029890278429282 ,
// 0.031250666632101 ,
// 0.033663616515945 ,
// 0.035772828811153 ,
// 0.037658931922342 ,
// 0.039372952752786 ,
// 0.042413137150944 ,
// 0.045070619540390 ,
// 0.047446991851855 ,
// 0.049606547407142 ,
// 0.051592857816434 ,
// 0.053436976294678 ,
// 0.055161906663233 ,
// 0.056785220390605 ,
// 0.058320678197501 ,
// 0.059779281606086 ,
// 0.061169980847062 ,
// 0.062500166625532 ,
// 0.067326227648337 ,
// 0.071544767303250 ,
// 0.075317060474946 ,
// 0.078745170560444 ,
// 0.084825640943857 ,
// 0.090140678210433 ,
// 0.099212631825087 ,
// 0.113570087453620 ,
// 0.125000041585273 ,
// 0.134652203946564 ,
// 0.143089312023620 ,
// 0.150633920105047 ,
// 0.157490157382677 ,
// 0.169651123546876 ,
// 0.180281216202297 ,
// 0.189786828693602 ,
// 0.198425147902268 ,
// 0.206370468208258 ,
// 0.213747007445228 ,
// 0.220646784186759 ,
// 0.227140086575032 ,
// 0.239116409077369 ,
// 0.250000010253906 ,
// 0.269304345055468 ,
// 0.286178568401423 ,
// 0.301267789998808 ,
// 0.314980268830742 ,
// 0.339302207508500 ,
// 0.360562397349922 ,
// 0.379573625756073 ,
// 0.396850266867541 ,
// 0.427493989953308 ,
// 0.454280151067003 ,
// 0.478232798228372 ,
// 0.500000002278646 ,
// 0.538608674401517 ,
// 0.572357122891389 ,
// 0.602535567444792 ,
// 0.629960526177828 ,
// 0.678604405120686 ,
// 0.721124785936175 ,
// 0.759147243604363 ,
// 0.793700526500832 ,
// 0.908560296613240 ,
// 1.000000000000000 };

//   for(i=0; i<N_V; i++)
//    cout<<  i  << " X " << X[i] <<endl;
// exit(0);


  Lines = new int[2*N_Cells];
  Vetrex = new TVertex*[N_V]; 

  y=0.;

  for(i=0; i<N_Cells; i++)
   {
    Lines[2*i]=i;
    Lines[2*i+1]=i+1;
    Vetrex[i] = new TVertex(X[i], y);
   }

  Vetrex[N_Cells] = new TVertex(X[N_V-1], y);

  CellTree = new TBaseCell*[N_Cells];

   for (i=0;i<N_Cells;i++)
   {
//     Vetrex[ i ]->GetCoords(x, y);
//     cout<< " x " << x<< " y " << y<<endl;
    CellTree[i] = new TMacroCell(TDatabase::RefDescDB[S_Line], 0);
    CellTree[i]->SetVertex(0, Vetrex[ Lines[ 2*i       ]]);
    CellTree[i]->SetVertex(1, Vetrex[ Lines[ 2*i + 1]]);
    ((TMacroCell *) CellTree[i])->SetSubGridID(0);
   }
//     Vetrex[ i ]->GetCoords(x, y);
//     cout<< " x " << x<< " y " << y<<endl;
//     exit(0);

   Domain->SetTreeInfo(CellTree, N_Cells);

   TDatabase::IteratorDB[It_EQ]->SetParam(Domain);
   TDatabase::IteratorDB[It_LE]->SetParam(Domain);
   TDatabase::IteratorDB[It_Finest]->SetParam(Domain);
   TDatabase::IteratorDB[It_Between]->SetParam(Domain);
   TDatabase::IteratorDB[It_OCAF]->SetParam(Domain);

   // start joint(vertex)
   Joint = new TJointEqN(CellTree[0]);
   CellTree[0]->SetJoint(0, Joint);


   for(i=1;i<N_Cells;i++)
    {
     Joint = new TJointEqN(CellTree[i-1], CellTree[i]);

     CellTree[i-1]->SetJoint(1, Joint);
     CellTree[i]->SetJoint(0, Joint);
   } // for(i=0;i<N_Cells;i++)

   // end joint(vertex)
   Joint = new TJointEqN(CellTree[N_Cells-1]);
   CellTree[N_Cells-1]->SetJoint(1, Joint);

  delete []  Lines;
}


