

void ExampleFile()
{
  OutPut("MassTransferTest.h" << endl) ;
  #define __MASSTRANSTEST__  
  #define __WITHSURFACTANT__    
}

// ========================================================================
// initial solution
// ========================================================================

void InitialU1(double x, double y, double *values)
{

  switch(int(values[0]))
   {
     case 0: // Phase 1
       values[0] = 0;
       values[1] = 0;
       values[2] = 0;
       values[3] = 0;
     break;
     
     case 1: // Phase 2 
       values[0] = 0;
       values[1] = 0;
       values[2] = 0;
       values[3] = 0;
     break;
     
  default:
    OutPut("only two phases are implemented check example and fefunction files " << endl);
    exit(1);

  }
}

void InitialU2(double x, double y, double *values)
{
  switch(int(values[0]))
   {
     case 0: // Phase 1
       values[0] = 0;
       values[1] = 0;
       values[2] = 0;
       values[3] = 0;
     break;
     case 1: // Phase 2
       values[0] = 0;
       values[1] = 0;
       values[2] = 0;
       values[3] = 0;
     break;
  default:
    OutPut("only two phases are implemented check example and fefunction files " << endl);
    exit(1);

  }

}

void InitialP(double x, double y, double *values)
{

  double We, U0;

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



void InitialS(double x, double y, double *values)
{
 double t = TDatabase::TimeDB->CURRENTTIME;
 double Re = TDatabase::ParamDB->RE_NR;
 double Pr = TDatabase::ParamDB->PR_NR;

  static double eps = 1.0/(Re*Pr);

  double phi = atan2(y-2., x);
  phi = Pi/2. - phi;

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

  values[0] = 0.;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void InitialS_Outer(double x, double y, double *values)
{
  values[0] = TDatabase::ParamDB->REACTOR_P19;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void ExactS_Outer(double x, double y, double *values)
{
 double Pe_c = TDatabase::ParamDB->REACTOR_P18;
 static double eps = 1.0/(Pe_c);

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

  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;

  switch(i)
      {
        case 0:
        case 1:
        case 3:
        case 5:
             cond = SLIP_FRICTION_PENETRATION_RESISTANCE;
             TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION = 1;
        break;

        case 2:
        case 4:
             cond = DIRICHLET;
        break;

        case 6:
             cond = FREESURF;
//        free surface boundary is inside the domain
             TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
        break;
  default:
            cout<<"BoundCond " << cond << endl;
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
  switch(i)
      {

        case 4:
             cond = NEUMANN;
        break;

        case 0:
        case 1:
        case 2:
        case 3:
        case 5:
        case 6:
             cond = NEUMANN;
        break;
  default:
            cout<<"BoundCond " << cond << endl;
            Error("Unknown Boundary component ref example file" << endl);
         exit(0);
       }
     }

void SurfactBoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
      {

        case 4:
             value = 0.0;
        break;

        case 0:
        case 1:
        case 2:
        case 3:
        case 5:
        case 6:
              value = 0;
        break;
  default:
            cout<<"BdComp " << BdComp << endl;
            Error("Unknown Boundary BdComp ref example file" << endl);
         exit(0);
       }
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
  int i, j, Cell_Phase_No;
  double *coeff;

 // for more details check DiscreteForm2D.C
  Cell_Phase_No = int(coeffs[0][0]);

// cout<< "Cell_Phase_No " << Cell_Phase_No << "eps1 " << eps1<< "eps2 " << eps2 <<endl; 
  switch(Cell_Phase_No)
   {
     case 0:  //   regular domain (inner)
      for(i=0;i<n_points;i++)
       {
        coeff = coeffs[i];
        coeff[0] = (1./TDatabase::ParamDB->P11)*eps;
        coeff[1] = 0; // f1
        if(TDatabase::ParamDB->FR_NR == 0)
         coeff[2] = 0;
       else
       coeff[2] = -1.;

// density ratio and viscosity ratio - no difference in outer phase
// dimensionless form w.r.t outer domain parameters
// so ratios will be in inner domain only
        coeff[3] = 1./TDatabase::ParamDB->P10;  // density ratio (inner/outer)
//         coeff[4] = TDatabase::ParamDB->P11;  // viscosity ratio (outer/inner)
//         coeff[5] = TDatabase::ParamDB->P13;  // density
//         coeff[6] = TDatabase::ParamDB->P15;  // viscosity
       }
     break;
     case 1: // Phase 1 Outer
      for(i=0;i<n_points;i++)
       { 
        coeff = coeffs[i];
        coeff[0] = eps;
        coeff[1] = 0; // f1

        if(TDatabase::ParamDB->FR_NR == 0)
         coeff[2] = 0;
        else
        coeff[2] = -1.;

// density ratio and viscosity ratio - no difference in outer phase
// scaled according to outer phase
        coeff[3] = 1.; // density ratio
//         coeff[4] = 1.; // viscosity ratio 
//         coeff[5] = 1.;  // density
//         coeff[6] = 1.;  // viscosity
       }
      break; 
 default:
    OutPut("only two phases are implemented check example and discreform files " << endl);
    exit(1);
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





