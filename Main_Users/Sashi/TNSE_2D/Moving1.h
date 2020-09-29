// Navier-Stokes problem, Benchmark channel
// circle moves up and down
// 
// u(x,y) = unknown
// p(x,y) = unknown

// ======================================================================
// common declaration for time dependent Navier-Stokes problems
// ======================================================================

// Galerkin assembling with mass matrix
int MassN_Terms = 4;
MultiIndex2D MassDerivatives[4] = { D10, D01, D00, D00 };
int MassSpacesNumbers[4] = { 0, 0, 0, 1 };
int MassN_Matrices = 6;
int MassRowSpace[6]    = { 0, 0, 1, 1, 0, 0 };
int MassColumnSpace[6] = { 0, 0, 0, 0, 1, 1 };
int MassN_Rhs = 2;
int MassRhsSpace[2] = { 0, 0 };

// Laplacian assembling with mass matrix
int AllLapN_Terms = 4;
MultiIndex2D AllLapDerivatives[4] = { D10, D01, D00, D00 };
int AllLapSpacesNumbers[4] = { 0, 0, 0, 1 };
int AllLapN_Matrices = 6;
int AllLapRowSpace[6]    = { 0, 0, 1, 1, 0, 0 };
int AllLapColumnSpace[6] = { 0, 0, 0, 0, 1, 1 };
int AllLapN_Rhs = 2;
int AllLapRhsSpace[2] = { 0, 0 };

// assemble only right hand sides
int RhsN_Terms = 1;
MultiIndex2D RhsDerivatives[1] = { D00 };
int RhsSpacesNumbers[1] = { 0 };
int RhsN_Matrices = 0;
int *RhsRowSpace = NULL;
int *RhsColumnSpace = NULL;
int RhsN_Rhs = 2;
int RhsRhsSpace[2] = { 0, 0 };

// Galerkin assembling only
int GalerkinN_Terms = 3;
MultiIndex2D GalerkinDerivatives[3] = { D10, D01, D00 };
int GalerkinSpacesNumbers[3] = { 0, 0, 0 };
int GalerkinN_Matrices = 1;
int GalerkinRowSpace[1]    = { 0 };
int GalerkinColumnSpace[1] = { 0 };
int GalerkinN_Rhs = 0;
int *GalerkinRhsSpace = NULL;

// Laplacian assembling only
int LaplacianN_Terms = 2;
MultiIndex2D LaplacianDerivatives[2] = { D10, D01 };
int LaplacianSpacesNumbers[2] = { 0, 0 };
int LaplacianN_Matrices = 1;
int LaplacianRowSpace[1]    = { 0 };
int LaplacianColumnSpace[1] = { 0 };
int LaplacianN_Rhs = 0;
int *LaplacianRhsSpace = NULL;

// ========================================================================
// parameter routine
// ========================================================================
void NavierStokesParams(double *in, double *out)
{
  out[0] = in[2];
  out[1] = in[3];
}

// ========================================================================
// general settings
// ========================================================================
int fevalue_fctindex[2] = { 0, 1 };
MultiIndex2D fevalue_multiindex[2] = { D00, D00 };
ParamFct *paramfct[1] = { NavierStokesParams };
int beginpara[1] = { 0 };

MultiIndex2D AllDerivatives[3] = { D00, D10, D01 };

// ========================================================================
// Galerkin assembling with mass matrix
// ========================================================================
void TimeNavStoAllAssemble(double Mult, double *coeff, 
                    double *parameters, double hK, 
                    double **OrigValues, int *N_BaseFuncts,
                    double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, *MatrixRowA;
  double **MatrixM, *MatrixRowM;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *MatrixRow1, *MatrixRow2;
  double *Rhs1, *Rhs2, val;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  double c0, c1, c2, c3, c4, c5, c6;
  int i,j,k,l, N_U, N_P;

  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];

  MatrixB1 = LocMatrices[2];
  MatrixB2 = LocMatrices[3];
  MatrixB1T = LocMatrices[4];
  MatrixB2T = LocMatrices[5];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u
  Orig3 = OrigValues[3]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = parameters[0]; // u1
  c4 = parameters[1]; // u2
  c5 = parameters[2]; // grid move x
  c6 = parameters[3]; // grid move y

  for(i=0;i<N_U;i++)
  {
    MatrixRowA = MatrixA[i];
    MatrixRowM = MatrixM[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      
      val = c0*(test10*ansatz10+test01*ansatz01)
           +((c3-c5)*ansatz10+(c4-c6)*ansatz01)*test00;

      MatrixRowA[j] += val*Mult;

      val = Mult*ansatz00*test00;

      MatrixRowM[j] += val;
    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;
      
      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    } // endfor j

  } // endfor i

} // end TimeNavStoAllAssemble

// ========================================================================
// Laplacian assembling with mass matrix
// ========================================================================
void TimeNavStoAllLapAssemble(double Mult, double *coeff, 
                    double *parameters, double hK, 
                    double **OrigValues, int *N_BaseFuncts,
                    double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, *MatrixRowA;
  double **MatrixM, *MatrixRowM;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *MatrixRow1, *MatrixRow2;
  double *Rhs1, *Rhs2, val;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  double c0, c1, c2, c3, c4;
  int i,j,k,l, N_U, N_P;

  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];

  MatrixB1 = LocMatrices[2];
  MatrixB2 = LocMatrices[3];
  MatrixB1T = LocMatrices[4];
  MatrixB2T = LocMatrices[5];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u
  Orig3 = OrigValues[3]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = parameters[0]; // u1
  c4 = parameters[1]; // u2

  for(i=0;i<N_U;i++)
  {
    MatrixRowA = MatrixA[i];
    MatrixRowM = MatrixM[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      
      val = c0*(test10*ansatz10+test01*ansatz01);

      MatrixRowA[j] += val*Mult;

      val = Mult*ansatz00*test00;

      MatrixRowM[j] += val;
    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;
      
      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    } // endfor j

  } // endfor i

} // end TimeNavStoAllLapAssemble

// ========================================================================
// assemble only the right hand sides
// ========================================================================
void TimeNavStoRhsAssemble(double Mult, double *coeff, 
                    double *parameters, double hK, 
                    double **OrigValues, int *N_BaseFuncts,
                    double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2, val;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  double c0, c1, c2, c3, c4;
  int i,j,k,l, N_U, N_P;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u

  for(i=0;i<N_U;i++)
  {
    test00 = Orig0[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

  } // endfor i

} // end TimeNavStoAllAssemble

// ========================================================================
// Galerkin assembling only
// ========================================================================
void TimeNavStoGalerkinAssemble(double Mult, double *coeff, 
                    double *parameters, double hK, 
                    double **OrigValues, int *N_BaseFuncts,
                    double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, *MatrixRowA;
  double val;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  double c0, c3, c4, c5, c6;
  int i,j,k,l, N_U;

  MatrixA = LocMatrices[0];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u

  c0 = coeff[0]; // nu
  c3 = parameters[0]; // u1
  c4 = parameters[1]; // u2
  c5 = parameters[2]; // grid x movement
  c6 = parameters[3]; // grid y movement

  for(i=0;i<N_U;i++)
  {
    MatrixRowA = MatrixA[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      
      val = c0*(test10*ansatz10+test01*ansatz01)
           +((c3-c5)*ansatz10+(c4-c6)*ansatz01)*test00;

      MatrixRowA[j] += val*Mult;

    } // endfor j

  } // endfor i

} // end TimeNavStogalerkinAssemble

// ========================================================================
// Laplacian assembling only
// ========================================================================
void TimeNavStoLaplacianAssemble(double Mult, double *coeff, 
                    double *parameters, double hK, 
                    double **OrigValues, int *N_BaseFuncts,
                    double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, *MatrixRowA;
  double val;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  double c0, c3, c4;
  int i,j,k,l, N_U;

  MatrixA = LocMatrices[0];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y

  c0 = coeff[0]; // nu

  for(i=0;i<N_U;i++)
  {
    MatrixRowA = MatrixA[i];
    test10 = Orig0[i];
    test01 = Orig1[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      
      val = c0*(test10*ansatz10+test01*ansatz01);

      MatrixRowA[j] += val*Mult;

    } // endfor j

  } // endfor i

} // end TimeNavStoLaplacianAssemble

// ========================================================================
// put l infinity norm of u in coeff5
// ========================================================================
void linfu(int N_Points, double **Coeffs, double **Params, TBaseCell *cell)
{
  int i;
  double max, *coeff, u1, u2, u, *param;

  max = -1;

  for(i=0;i<N_Points;i++)
  {
    coeff = Coeffs[i];
    param = Params[i];
    u1 = param[0];
    u2 = param[1];

    u = MAX(fabs(u1), fabs(u2));

    if(u>max) max = u;
  }

  for(i=0;i<N_Points;i++)
    Coeffs[i][5] = max;

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
  cond = DIRICHLET;
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  switch(BdComp)
  {
    case 0: value = 0;
            break;
    case 1: value = 0;
            break;
    case 2: value = 0;
            break;
    case 3: value = 0;
            break;
    case 4: value = 0;
            break;
    default: cout << "wrong boundary part number: " << BdComp << endl;
  }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  double v=TDatabase::ParamDB->P7;

  switch(BdComp)
  {
    case 0: value = 0;
            break;
    case 1: value = 0;
            break;
    case 2: value = 0;
            break;
    case 3: value = 0;
            break;
    case 4: value = v*2*Pi*0.05*cos(v*2*Pi*t);
            break;
    default: cout << "wrong boundary part number: " << BdComp << endl;
  }
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *x, double *y,
               double **parameters, double **coeffs)
{
  static double eps = 1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = 0; // f1
    coeff[2] = 0; // f2
  }
}

// ======================================================================
// Laplacian problem for grid moving
// ======================================================================
void GridParams(double *in, double *out)
{
  out[0] = in[2];
  out[1] = in[3];
  out[2] = in[4];
  out[3] = in[5];
}

int GridN_Terms = 2;
MultiIndex2D GridDerivatives[2] = { D10, D01 };
int GridSpaceNumbers[2] = { 0, 0 };
int GridN_Matrices = 1;
int GridRowSpace[1] = { 0 };
int GridColumnSpace[1] = { 0 };
int GridN_Rhs = 0;
int *GridRhsSpace = NULL;

int Gridfevalue_fctindex[4] = { 0, 1, 2, 3 };
MultiIndex2D Gridfevalue_multiindex[4] = { D00, D00, D00, D00 };
ParamFct *Gridparamfct[1] = { GridParams };
int Gridbeginpara[] = { 0 };


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

    if( (r2 = (x[i]-0.2)*(x[i]-0.2) + (y[i]-0.2)*(y[i]-0.2)) < 0.01)
      coeff[0] = 10*sqrt(r2);
    else
      coeff[0] = 10*0.1;
    coeff[1] = 0;
    coeff[2] = 0;
    coeff[3] = 0;

    coeff[4] = 0;
  }
}

void GridAssemble(double Mult, double *coeff, double hK, 
                    double **OrigValues, int *N_BaseFuncts,
                    double ***LocMatrices, double **LocRhs)
{
  double **Matrix, *Rhs, val, *MatrixRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,k,l, N_;
  double c0, c1, c2, c3, c4; 

  Matrix = LocMatrices[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];

  c0 = coeff[0];
  c1 = coeff[1];
  c2 = coeff[2];
  c3 = coeff[3];
  c4 = coeff[4];

  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    test10 = Orig0[i];
    test01 = Orig1[i];

    for(j=0;j<N_;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      
      val = c0*(test10*ansatz10+test01*ansatz01);

      val *=Mult;

      MatrixRow[j] += val;
    } // endfor j
  } // endfor i
}

