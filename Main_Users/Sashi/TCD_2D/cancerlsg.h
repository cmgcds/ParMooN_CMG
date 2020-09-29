// ======================================================================
// instationary problem
// ======================================================================

/// ========================================================================
// example file
// ========================================================================

#define __CANCER__

void ExampleFile()
{
  OutPut("Example: cancerlsg.h" << endl); 
}

// exact solution
void Exact(double x, double y, double *values)
{
  values[0] =  0;
}

// =========================================================================
// kind of boundary condition (for FE space needed)
// =========================================================================
void BoundCondition(int BdComp, double t, BoundCond &cond)
{
 cond = NEUMANN;
  // cond = DIRICHLET;
}

void V_BoundCondition(int BdComp, double t, BoundCond &cond)
{
 cond = NEUMANN;
  // cond = DIRICHLET;
}

void W_BoundCondition(int BdComp, double t, BoundCond &cond)
{
  //cond = DIRICHLET;
 cond = NEUMANN;
}

// =========================================================================
// value of boundary condition
// =========================================================================
void BoundValue(int BdComp, double Param, double &value)
{
//    double t;
// 
//   t = TDatabase::TimeDB->CURRENTTIME;
//   switch(BdComp)
//   {
//     case 0: value = Param+t;
//             break;
//     case 1: value = 1+Param+t;
//             break;
//     case 2: value = 1+Param+t;
//             break;
//     case 3: value = Param+t;
//             break;
//   }
  value = 0;
}

void V_BoundValue(int BdComp, double Param, double &value)
{
//    double t;
// 
//   t = TDatabase::TimeDB->CURRENTTIME;
//   switch(BdComp)
//   {
//     case 0: value = -Param+t;
//             break;
//     case 1: value = Param-1+t;
//             break;
//     case 2: value = 1-Param+t;
//             break;
//     case 3: value = Param+t;
//             break;
//   }
  value = 0;
}

void W_BoundValue(int BdComp, double Param, double &value)
{
//    double t;
// 
//   t = TDatabase::TimeDB->CURRENTTIME;
//   switch(BdComp)
//   {
//     case 0: value = 0;
//             break;
//     case 1: value = Param*t;
//             break;
//     case 2: value = Param*t;
//             break;
//     case 3: value = 0;
//             break;
//   }
  value = 0;
}

// =========================================================================
// initial conditon
// =========================================================================
void InitialCondition(double x, double y, double *values)
{
   //double t;
//  
   //t = TDatabase::TimeDB->CURRENTTIME;
 // values[0] = exp(t)*(sin(2*Pi*x)*sin(2*Pi*y));
  
//   values[0]= 0.5*(1-tanh(x-3.0-0.2*cos(3.0*y)));
  
    double z;
//     z= (x-0.5)*(x-0.5)+(y-0.5)*(y-0.5);
//    if(z>0.1)
//     { values[0] = 0.; }
//    else
//     { 
//      values[0] =  exp(-z/0.01);
//     }
//     
  
  //circle Chaplain Paper
   z= x*x+y*y;
   if(z>0.25)
    { values[0] = 0.; }
   else
    { 
     values[0] =  exp(-z/0.01);
    }
//   \\cone
//   values[0] = 1-(sqrt((x-0.625)*(x-0.625)+(y-0.375)*(y-0.375)))/0.25;
    
    
 //circle 
//   values[0] =  exp(-((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5))/0.01);
  
 //asteriod
//   values[0] =  exp(-((x-0.5)*(x-0.5)*(y-0.5)*(y-0.5) )/0.01);
  
 
 //ellptic 
 //values[0] =  exp(-(((x-0.5)*(x-0.5))/0.10+((y-0.5)*(y-0.5)))/0.05);
 
}

void V_InitialCondition(double x, double y, double *values)
{/*
   double t = TDatabase::TimeDB->CURRENTTIME;
   values[0]=x-y+t;*/
  // values[0]=1-0,5*(x*x+y*y);
  //values[0] = exp(t)*(sin(2*Pi*x)*sin(2*Pi*y));
  
//     values[0]= 1- 0.5*(1-tanh(x-3.0-0.2*cos(3.0*y)));
  
 double z;
//     z= (x-0.5)*(x-0.5)+(y-0.5)*(y-0.5);
//     values[0] = 1-0.5*exp(- ((z*z-2*0.01*x)-(4*0.01*0.01*z))/0.01);


 //circle Chaplain Paper
   z= x*x+y*y;
   if(z>0.25)
    { values[0] = 1.; }
   else
    { 
     values[0] = 1. - 0.5*exp(-z/0.01);
    }
 
  //circle
  //values[0] = 1-( exp(-((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5))/0.01));
// values[0] = 1-(0.5* exp(-((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5))/0.01));
  
  //ellptic 
 // values[0] =  1-0.5*exp(-(((x-0.5)*(x-0.5))/0.10+((y-0.5)*(y-0.5)))/0.05);
  
   //asteriod
//   values[0] =  1-(0.5*exp(-((x-0.5)*(x-0.5)*(y-0.5)*(y-0.5) )/0.01));
}
void W_InitialCondition(double x, double y, double *values)
{
//    double t = TDatabase::TimeDB->CURRENTTIME;
  //values[0]= 0,5*(x*x+y*y);
//    values[0]=x*y*t;
  //values[0] = exp(t)*(sin(2*Pi*x)*sin(2*Pi*y));
  
//     values[0]= 0.75*(1-tanh(x-3.0-0.2*cos(3.0*y)));
    
    
  double z;
//     z= (x-0.5)*(x-0.5)+(y-0.5)*(y-0.5);
//     values[0] = 0.5*exp(- ((z*z-2*0.01*x)-(4*0.01*0.01*z))/0.01);
  
  
   //circle Chaplain Paper
   z= x*x+y*y;
   if(z>0.25)
    { values[0] = 0.; }
   else
    { 
     values[0] =  0.5*exp(-z/0.01);
    }
  
  //circle
// values[0] =  0.5*exp(-((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5))/0.01);
  
 //ellptic 
  // values[0] =  0.5*exp(-(((x-0.5)*(x-0.5))/0.10+((y-0.5)*(y-0.5)))/0.05);
  
   //asteriod
//   values[0] = (0.5*exp(-((x-0.5)*(x-0.5)*(y-0.5)*(y-0.5) )/0.01));
}
// =========================================================================
// linear coeffs
// =========================================================================
void BilinearCoeffs(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  double eps=1./TDatabase::ParamDB->PE_NR;
  double b1=0, b2=0, c=0;
  int i;
  double *coeff, *param;
  double x, y;
  double t = TDatabase::TimeDB->CURRENTTIME;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

//     x = X[i];
//     y = Y[i];

    coeff[0] = eps;
    coeff[1] = b1;
    coeff[2] = b2;
    coeff[3] = c;
    coeff[4] = 0.;
  }
}

void V_BilinearCoeffs(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  double eps=0./TDatabase::ParamDB->PE_NR;
  double b1=0, b2=0, c=0;
  int i;
  double *coeff, *param;
  double x, y;
  double t = TDatabase::TimeDB->CURRENTTIME;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

//     x = X[i];
//     y = Y[i];

    coeff[0] = eps;
    coeff[1] = b1;
    coeff[2] = b2;
    coeff[3] = c;
    coeff[4] = 0.;
  }
}

void W_BilinearCoeffs(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  double eps=5/TDatabase::ParamDB->PE_NR;
  double b1=0, b2=0, c=0;
  int i;
  double *coeff, *param;
  double x, y;
  double t = TDatabase::TimeDB->CURRENTTIME;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

//     x = X[i];
//     y = Y[i];

    coeff[0] = eps;
    coeff[1] = b1;
    coeff[2] = b2;
    coeff[3] = c;
    coeff[4] = 0.;
  }
}


// =========================================================================

