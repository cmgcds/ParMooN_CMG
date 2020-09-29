// Navier-Stokes problem, SFB cavity
// 
// u(x,y) = unknown
// p(x,y) = unknown

void ExampleFile()
{
  int range;

  OutPut("Example: SFBCavity_2.h ") ;
  OutPut("inflow " << TDatabase::ParamDB->P6);
  OutPut(" upper lid " << TDatabase::ParamDB->P5);
  OutPut(" left " << (int)TDatabase::ParamDB->P7);
  OutPut(" right " << (int)TDatabase::ParamDB->P8);
  OutPut(" lower " << (int)TDatabase::ParamDB->P9 << endl);

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
}

void U1BoundValue(int BdComp, double Param, double &value)
{
    int range;

  switch(BdComp)
  {
     case 0: value = 0;
        break;
    case 1:
	range = (int)TDatabase::ParamDB->P8;
	if ((Param>range/32.0)&&(Param<(range+1)/32.0))
	    value = -1;
	else
	    value = 0;
	break;
    case 2: if(Param<0.00001 || Param>0.99999) 
              value = 0;
            else
              value = TDatabase::ParamDB->P5/TDatabase::ParamDB->P6;
            break;
    case 3: 
	range = (int)TDatabase::ParamDB->P7;
	if (((1-Param)>range/32.0)&&((1-Param)<(range+1)/32.0))
	    value = 1;
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
    coeff[1] = 0; // f1
    coeff[2] = 0; // f2
  }
}

