#include <TimeConvDiff2D.h>
#include <MacroCell.h>



#define  __AXIAL3D__

void ExampleFile()
{
  OutPut("Example: 2D_1D_ADI_UnitSrq.h" << endl) ;

  TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD=1;
  #define __SIMPATURS__
}

//input and output : dimless
double GetCsat(double C, double T)
{
 double val;
//  double T_ref = TDatabase::ParamDB->REACTOR_P23;
//  double C_ref = TDatabase::ParamDB->REACTOR_P25;
  
//   val = (1.3045*(T_ref*T - 273.15) + 35.3642)/(C_ref*TDatabase::ParamDB->UREA_m_mol); 
   val = ( 35.3642 + 1.3045*(T - 273.15));

  if(val<0.) val =0.;

  return(val);
}

void InitialU1(double x, double y, double *values)
{

  values[0] = TDatabase::ParamDB->UREA_INFLOW_SCALE*(.25-y*y)/TDatabase::ParamDB->UREA_u_infty;
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
//   double T_W = TDatabase::ParamDB->REACTOR_P24;

//   val = T_W/T_D;

  values[0] = T_D;
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
   case 1:
     cond = NEUMANN;
   break;    
    
   case 2:
   case 3:
     cond = DIRICHLET;
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

//   val = T_W/T_D;

  switch(BdComp)
  {
   case 0:  
   case 1:
     value=0.;
   break;

   case 2:
        value=T_W;
   break;    

   
   case 3:
//      y = 0.5 -  Param;
//      dt = 1. - val;
// //  cout << y << " t " << Param << endl;
//      if( y>0 && y<(1./8.))
//       {
// 	t = -6. + 12.*y*8.*(1. + 1./8. - y);
// 	t = val + dt*(1./(1. + exp(-t)));
// //              cout << y << " t " << t << endl;
// 	value = t;  
//       }
//      else if(y<1. &&   y>(1. - 1./8.))
//       {
// 	t = 6. + 12.*(1. - 1./8. - y)*8.*(1. + 1. -y);
// 	t = 1./(1. + exp(-t));
// 	t = val + dt*(1./(1. + exp(-t)));
// //          cout << y << " t " << t << endl;
// 	value = t;      
//       }
//       else if(y==0. || y==0.5)
//       {
//        value = val;  
//       }   
//       else
      {
	value = T_D;  
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

    
    // x- axis as axis of symmetry
//     coeff[20] = y; 
    
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
  values[0] = 1.1972;
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
     cond = NEUMANN;
   break;

   case 3:
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
     value=0.;
  break;

  case 3:
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
	value = GetCsat(0, 301.5)/(TDatabase::ParamDB->UREA_c_infty * TDatabase::ParamDB->UREA_m_mol); 
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

    // x- axis as axis of symmetry
//     coeff[20] = y;     
    
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
     cond = NEUMANN;
   break;

   case 3:
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


double GetInletPsd()
{
 double  value, a = TDatabase::ParamDB->REACTOR_P29; 
 int i, found=0;
 double m, eps = 1e-8, L_max = TDatabase::ParamDB->UREA_D_P_MAX;
 double scale_f_max = 1./TDatabase::ParamDB->UREA_INFLOW_SCALE;
 
  double inlet_coord[100] ={2.500000e-06,1.356050e-05,3.068150e-05,4.780300e-05,6.492450e-05,8.204550e-05,9.916650e-05,1.162875e-04,1.334090e-04,1.505305e-04,1.676515e-04,1.847725e-04,2.018940e-04,2.190155e-04,2.361365e-04,2.532575e-04,2.703785e-04,2.875000e-04,3.046215e-04,3.217425e-04,3.388635e-04,3.559845e-04,3.731060e-04,3.902275e-04,4.073485e-04,4.244695e-04,4.415910e-04,4.587125e-04,4.758335e-04,4.929545e-04,5.100755e-04,5.271970e-04,5.443185e-04,5.614395e-04,5.785605e-04,5.956815e-04,6.128030e-04,6.299245e-04,6.470455e-04,6.641665e-04,6.812875e-04,6.984090e-04,7.155305e-04,7.326515e-04,7.497725e-04,7.668940e-04,7.840155e-04,8.011365e-04,8.182575e-04,8.353785e-04,8.525000e-04,8.696215e-04,8.867425e-04,9.038635e-04,9.209845e-04,9.381060e-04,9.552275e-04,9.723485e-04,9.894695e-04,1.006591e-03,1.023712e-03,1.040833e-03,1.057955e-03,1.075075e-03,1.092197e-03,1.109318e-03,1.126439e-03,1.143561e-03,1.160682e-03,1.177803e-03,1.194925e-03,1.212045e-03,1.229167e-03,1.246287e-03,1.263409e-03,1.280530e-03,1.297651e-03,1.314773e-03,1.331894e-03,1.349016e-03,1.366137e-03,1.383257e-03,1.400379e-03,1.417500e-03,1.434622e-03,1.451742e-03,1.468863e-03,1.485985e-03,1.503106e-03,1.520228e-03,1.537348e-03,1.554470e-03,1.571591e-03,1.588713e-03,1.605834e-03,1.622954e-03,1.640075e-03,1.657197e-03,1.674318e-03,1.691440e-03};

  double inlet_f_L_seed[100] = {0.000000e+00,1.539144e+09,2.082543e+09,1.056886e+09,1.032820e+09,9.174217e+08,7.339395e+08,6.383108e+08,5.630089e+08,4.776040e+08,3.952118e+08,3.244646e+08,2.668657e+08,2.186269e+08,1.774119e+08,1.430176e+08,1.142823e+08,9.072941e+07,7.094795e+07,5.432124e+07,4.189605e+07,3.280686e+07,2.562972e+07,2.013453e+07,1.620424e+07,1.315460e+07,1.040263e+07,7.846039e+06,5.728553e+06,4.495501e+06,3.800058e+06,3.005911e+06,2.229857e+06,1.675695e+06,1.307405e+06,1.003372e+06,7.486895e+05,6.805768e+05,7.435860e+05,7.812756e+05,7.396855e+05,6.082531e+05,3.731120e+05,1.582579e+05,8.976509e+04,1.072953e+05,9.867716e+04,5.302092e+04,3.384872e+04,7.221184e+04,1.380859e+05,1.713877e+05,1.450725e+05,9.890296e+04,6.954686e+04,3.619885e+04,1.297526e+04,1.931888e+04,5.546763e+04,1.066047e+05,1.070156e+05,5.457797e+04,1.540227e+04,3.078500e+03,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00};  
 
 
    scale_f_max /=TDatabase::ParamDB->UREA_inflow_time;   
    scale_f_max *=6.e7;
  
    
    Dscal(100, scale_f_max, inlet_f_L_seed);   
    
     for (i=100;i>0;i--)
     { 
      if((a*L_max<=inlet_coord[i])&&(a*L_max>inlet_coord[i-1]))
      {
        m=(inlet_f_L_seed[i]-inlet_f_L_seed[i-1])/(inlet_coord[i] - inlet_coord[i-1]);
        value = m *(a*L_max- inlet_coord[i-1]);
        value += inlet_f_L_seed[i-1];
	value/=TDatabase::ParamDB->UREA_f_infty;
//         OutPut(a <<" a - L_max is !!!"<<  L_max <<" i ist !!!" <<i << " i is !!!"<<inlet_coord[i]<<endl);
        found=1;
        break;
       }
      }      
 
    
    if(found==0)
      value= 0.;      
  
  return value;
}


void BoundValue_Psd(int BdComp, double Param, double &value)
{  
 double t = TDatabase::TimeDB->CURRENTTIME;
 int i;

    switch(BdComp)
     {
      case 0:
      case 1:
      case 2:
	value = 0.;
      break;
      case 3:
       
// 	if(t<=TDatabase::ParamDB->UREA_inflow_time)
	 {value = GetInletPsd();}
//         else
//          {value =0.;}
      break;
      
      default: cout << "wrong boundary part number "  << BdComp  << endl;
      break;
      }    

}

void BoundCondition_LminLMax(BoundCond &cond_Lmin, BoundCond &cond_Lmax)
{
  cond_Lmin = DIRICHLET;
//   cond_Lmin = NEUMANN;
  cond_Lmax = NEUMANN;
}



void BoundValue_LMin(double x, double y,  double *values)
 {
   double t;
  
   
   cout<<"Do not call here, set the nucleation value during assembling " <<endl;
   exit(0);
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
  { eps = 1.0/TDatabase::ParamDB->REACTOR_P2; }
  else
  { eps = 0.; }


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

     // Inlet 
     TDatabase::ParamDB->REACTOR_P29 = X[N_Coord -1];  // L value
       if(fabs(X[0])<1e-8 )
         {values[0] = GetInletPsd();}
        else
         {values[0] = 0;}

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


double GetGrowth(double C, double T)
{
  double kg = TDatabase::ParamDB->UREA_k_g;
  double G, C_Sat, Mode = TDatabase::ParamDB->PBE_P0;
  double C_orig;

   C_Sat = GetCsat(C, T);
   C_orig = C*TDatabase::ParamDB->UREA_c_infty * TDatabase::ParamDB->UREA_m_mol;

//             cout <<   "C_orig: " <<C_orig<<" C_Sat: " << C_Sat <<endl; 
    // growth rate based on concentration
    if(C_orig>C_Sat && (Mode==1 || Mode==4|| Mode==6) )
     { 
      G = kg*sqrt((C_orig - C_Sat)/C_Sat); 

      if(G<0.) G = 0.;
     }
    else
     { G = 0.; }

//         if(fabs(G)>1e-4)
//          cout << "G: " << G << "CS: " <<C <<"T: " << T <<endl;  
 return G;
 }

 // input and output in dimless form
void Get_GrowthAndB_Nuc(int N_Inputs, double *Inn, double *Out) 
{
//   double PbeGC = TDatabase::ParamDB->REACTOR_P27;
//   double C_Sat;
//   double T_ref = TDatabase::ParamDB->REACTOR_P23;
//   double C_ref = TDatabase::ParamDB->REACTOR_P25;
  double alpha =  TDatabase::ParamDB->UREA_alfa_nuc;
  double beta =  TDatabase::ParamDB->UREA_beta_nuc;
//   double f_ref = TDatabase::ParamDB->UREA_f_infty;
  double C, T, G, B_Nuc, C_Sat, val;
 
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
     { 
       C_Sat = GetCsat(C, T);
       val =  C*TDatabase::ParamDB->UREA_c_infty * TDatabase::ParamDB->UREA_m_mol/ C_Sat;   
       val=log(val);
       val = -beta/(val*val);
       val =  alpha*exp(val);     
       B_Nuc = val/(G*TDatabase::ParamDB->UREA_f_infty); 
    }
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
    
//          cout << "G: " << G << " B_Nuc: " <<B_Nuc<< "  C: " <<C <<" T: " <<T <<endl;
}


void GetExampleFileData(BoundCondFunct2D **BoundaryConditions, BoundValueFunct2D **BoundValues, 
                        DoubleFunct2D **InitiaValues, CoeffFct2D **BilinearCoeffs, 
                        int &N_PBEqns, int &N_IndepntScalarEqns, int *Disctypes)
{

  TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD=1;


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
  double len, h, x, y;
  TVertex **Vetrex;
  TJoint *Joint;
  TBaseCell  **CellTree;

//   N_V = N_Cells+1;
//   X = new double[N_V];

//   for(i=0; i<N_V; i++)
//    X[i] = 1. + (1. - Start)*(tanh(2.75*(double(i)/double(N_Cells) - 1.)))/tanh(2.75);

//   h = (End - Start)/N_Cells;
//   for(i=1; i<N_V; i++)
//    X[i] =  h*(double)i;

//   X[0] = Start;
//   X[N_V-1] = End;
// 
N_V = 94;
N_Cells = N_V-1;

  
  
 double X[94] =
  {
    0,    0.000488281,  0.000615194,    0.000775097,    0.000976561,    0.00123039,   
    0.00155019, 0.00195312,    0.00246078,    0.00310039,    0.00390625,    0.00447154,
    0.00492156,    0.0053016,    0.00563379,    0.00593084,    0.00620079,    0.00667959,
    0.00709813,    0.00747238,    0.0078125,    0.00841576,    0.00894308,
    0.00941462,    0.00984313,    0.0106032,    0.0112676,    0.0118617,    0.0124016,
    0.0133592,    0.0141962,    0.0149448,    0.015625,    0.0168315,    0.0178862,
    0.0188292,    0.0196863,    0.0204745,    0.0212064,    0.0218909,    0.0225352,
    0.0231445,    0.0237233,    0.0242752,    0.0248031,    0.0267184,    0.0283925,
    0.0298896,    0.03125,    0.0336631,    0.0357723,    0.0393725,    0.0450703,
    0.0496062,    0.0534368,    0.056785,    0.0597791,    0.0625,    0.0673261,
    0.0715446,    0.0753169,    0.0787451,    0.0818982,    0.0848255,    0.0875637,
    0.0901406,    0.0948934,    0.0992125,    0.106873,    0.11357,    0.119558,
    0.125,    0.134652,    0.143089,    0.150634,    0.15749,    0.169651,    0.180281,
    0.189787,    0.198425,    0.213747,    0.22714,    0.239117,    0.25,    0.269304,
    0.286179,    0.301268,    0.31498,    0.360562,    0.39685,    0.5,    0.629961,
    0.793701,   1.};
  
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


