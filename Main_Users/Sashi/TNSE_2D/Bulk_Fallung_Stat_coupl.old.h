// Navier-Stokes problem, SFB cavity
// 
// u(x,y) = unknown
// p(x,y) = unknown

void ExampleFile()
{
  int range;

  OutPut("Example: Bulk_Fallung_Stat.h ") ;
  OutPut("inflow (u_infty)" << TDatabase::ParamDB->P6);
  OutPut(" upper lid " << TDatabase::ParamDB->P5);
  OutPut(" left " << (int)TDatabase::ParamDB->P7);
  OutPut(" right " << (int)TDatabase::ParamDB->P8);
  OutPut(" lower " << (int)TDatabase::ParamDB->P9);
  OutPut(" d_p_max " << TDatabase::ParamDB->P2  << endl);

  range = (int)TDatabase::ParamDB->P7;
  if ((range<1)||(range>30))
  {
      OutPut("left boundary out of range !!!"<< endl);
      exit(4711);
  }

  range = (int)TDatabase::ParamDB->P8;
  if ((range<1)||(range>30))
  {
      OutPut("right boundary out of range !!!"<< endl);
      exit(4711);
  }

  range = (int)TDatabase::ParamDB->P9;
  if ((range<1)||(range>30))
  {
      OutPut("lower boundary out of range !!!"<< endl);
      exit(4711);
  }
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
  int lower = (int)TDatabase::ParamDB->P9;

  if ((i==0)&&((t>lower/32.0)&&(t<(lower+2)/32.0)))
      cond = NEUMANN;
  else
      cond = DIRICHLET;
  // cond = DIRICHLET;
  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0; 
}

void U1BoundValue(int BdComp, double Param, double &value)
{
    int range;
    double y, fac = 1024; // 32^2
    

  switch(BdComp)
  {
     case 0: 
        value = 0;
        break;
     case 1:
	range = (int)TDatabase::ParamDB->P8;
	if ((Param>range/32.0)&&(Param<(range+1)/32.0))
        {
           y = Param;
           value =  6*(y-range/32.0)*(y-(range+1)/32.0)*fac;
           //  OutPut(value << " ");
        }
	else
	    value = 0;
	break;
    case 2: if(Param<0.00001 || Param>0.99999) 
              value = 0;
            else
               value = TDatabase::ParamDB->P5/TDatabase::ParamDB->P6;
               //value = 0;
            break;
    case 3: 
      //range = (int)TDatabase::ParamDB->P7;
        range = (int)TDatabase::ParamDB->P8;
        y = 1-Param;
	if ((y>range/32.0)&&(y<(range+1)/32.0))
        {
           value = -6*(y-range/32.0)*(y-(range+1)/32.0)*fac;
        }
	else
	    value = 0;
	break;
    default: cout << "wrong boundary part number" << endl;
  }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  value = 0;
/*  int lower = (int)TDatabase::ParamDB->P9;

  if ((BdComp==0)&&((Param>lower/32.0)&&(Param<(lower+1)/32.0)))
      value = -2;
  else
      value = 0;
*/
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *x, double *y,
               double **parameters, double **coeffs)
{
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = 1e-6/TDatabase::ParamDB->P6;
    //coeff[0] = 1;
    coeff[1] = 0; // f1
    coeff[2] = 0; // f2
  }
}


// ========================================================================
// definitions for the substance A
// ========================================================================

// initial conditon
void InitialCondition_c_A(double x, double y, double *values)
{
   int range;
   double eps=1e-8;

   range = (int)TDatabase::ParamDB->P7;
   if ((fabs(x)<1e-7)&&(y>=range/32.0-eps)&&(y<=(range+1)/32.0+eps))
       values[0] = 1;
   else
       values[0] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition_c_A(int BdComp, double Param, BoundCond &cond)
{
   int range;
   double y,eps=1e-8;

   if (BdComp==3)
   {
      range = (int)TDatabase::ParamDB->P7;
      y = 1-Param;
      if ((y>=range/32.0-eps)&&(y<=(range+1)/32.0+eps))
      {
         cond = DIRICHLET;
      }
      else
	  cond = NEUMANN;
   }
   else
     cond = NEUMANN;
}

// value of boundary condition
void BoundValue_c_A(int BdComp, double Param, double &value)
{
   int range;
   double y,eps=1e-8;

   if (BdComp==3)
   {
      range = (int)TDatabase::ParamDB->P7;
      y = 1-Param;
      if ((y>=range/32.0-eps)&&(y<=(range+1)/32.0+eps))
      {
	  value = 1;
      }
      else
	  value = 0;
   }
   else
     value = 0;

}

void NoCoeffs(int n_points, double *X, double *Y,
		    double **parameters, double **coeffs)
{
    return;
}

void BilinearCoeffs(int n_points, double *X, double *Y,
		    double **parameters, double **coeffs)
{
  int i;
  double *coeff, *param;
  double x, y;
  double t = TDatabase::TimeDB->CURRENTTIME;
  double L_infty = 1;
  double U_infty = TDatabase::ParamDB->P6;
  double T_infty;
  double C_infty = 1;
  double D_A = 1.5e-9;
  double k_r = 1e-2;
  // double k_r = 10;
  double eps,c;
  
  T_infty = L_infty/U_infty;
  eps = D_A/(L_infty*U_infty);
  c = k_r*C_infty * L_infty /U_infty;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    x = X[i];
    y = Y[i];

    coeff[0] = eps;
    coeff[1] = param[1]; // u1
    coeff[2] = param[2]; // u2
    coeff[3] = c * param[0];
    coeff[4] = 0; 
    //OutPut(param[0] << " " << param[1] << " " << param[2] << endl); 
  }
}

void BilinearCoeffs_Cc(int n_points, double *X, double *Y,
		    double **parameters, double **coeffs)
{
  int i;
  double *coeff, *param;
  double x, y;
  double t = TDatabase::TimeDB->CURRENTTIME;
  double L_infty = 1;
  double U_infty = TDatabase::ParamDB->P6;
  double T_infty;
  double C_infty = 1;
  double D_A = 1.5e-9;
  double k_r = 1e-2;
  double C_nuc = 15.33;
  double d_p_0 = 1e-9;
  double C_sat = 1.37e-4;
  double k_nuc = 1e+24;
  double C_2 = 7.2e-9;
  double eps,c,c1,B_C_c;
  double C_c_infty;
  double f_infty = 1;
  double C_g = 45.98;
  double k_g = 1e-7;
  double d_p_max = 1e-3;
  double lambda_gamma;

  T_infty = L_infty/U_infty;
  eps = D_A/(L_infty*U_infty);
  C_c_infty = C_sat;
  c = k_r*C_infty*C_infty*L_infty /(U_infty*C_c_infty);
  c1 = L_infty*C_c_infty*C_c_infty*C_c_infty*C_c_infty/U_infty;

  //lambda_gamma = C_g * k_g * d_p_max * d_p_max * d_p_max * l_infty * f_infty * u_infty;
  // this is for an appropriate choice of f_infty, see 2d3d_FDM.C
  lambda_gamma = 1;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    x = X[i];
    y = Y[i];

    if (param[4] < exp(C_2/d_p_0))
      B_C_c = 0;
    else
      B_C_c = k_nuc*(param[4] - exp(C_2/d_p_0))*(param[4] - exp(C_2/d_p_0))*(param[4] - exp(C_2/d_p_0))*(param[4] - exp(C_2/d_p_0))*(param[4] - exp(C_2/d_p_0));

    coeff[0] = eps;
    coeff[1] = param[0]; // u1
    coeff[2] = param[1]; // u2
    coeff[3] = 0;
    coeff[4] = c*param[2]*param[3] - C_nuc*d_p_0*d_p_0*d_p_0*c1*B_C_c 
	- lambda_gamma*(param[4] - 1)*param[5];//r_chem - r_nuc - r_g
  }
}

// ========================================================================
// definitions for the substance B
// ========================================================================

// initial conditon
void InitialCondition_c_B(double x, double y, double *values)
{
    int range;
    double eps = 1e-8;

    range = (int)TDatabase::ParamDB->P8;
    if ((fabs(1-x) <1e-7)&&(y>=range/32.0-eps)&&(y<=(range+1)/32.0+eps))
	values[0] = 1;
    else
	values[0] = 0;
	
}

// kind of boundary condition (for FE space needed)
void BoundCondition_c_B(int BdComp, double Param, BoundCond &cond)
{
   int range;
   double eps = 1e-8;

   if (BdComp==1)
   {
      range = (int)TDatabase::ParamDB->P8;
      if ((Param>=range/32.0-eps)&&(Param<=(range+1)/32.0+eps))
      {
         cond = DIRICHLET;
      }
      else
	  cond = NEUMANN; 
   }
   else
     cond = NEUMANN;
}

// value of boundary condition
void BoundValue_c_B(int BdComp, double Param, double &value)
{
   int range;
   
   if (BdComp==1)
   {
      range = (int)TDatabase::ParamDB->P8;
      if ((Param>=range/32.0)&&(Param<=(range+1)/32.0))
      {
	  value = 1;
      }
      else
	  value = 0;
   }
   else
     value = 0;
}

// ========================================================================
// definitions for the substance C
// ========================================================================

// initial conditon
void InitialCondition_c_C(double x, double y, double *values)
{
  values[0] = 0;	
}

// kind of boundary condition (for FE space needed)
void BoundCondition_c_C(int BdComp, double Param, BoundCond &cond)
{
   cond = NEUMANN;
}

// value of boundary condition
void BoundValue_c_C(int BdComp, double Param, double &value)
{
   value = 0;
}
