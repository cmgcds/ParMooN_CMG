// instationary Rosensweig instability

void ExampleFile()
{
  OutPut("InstFFLayer" << endl);
}

// ========================================================================
// initial state is set to zero
// ========================================================================
void InitialU1(double x, double y, double z, double *values)
{
  values[0] = 0; // x;
}

void InitialU2(double x, double y, double z, double *values)
{
  values[0] = 0; // y;
}

void InitialU3(double x, double y, double z, double *values)
{
  values[0] = 0; // z;
}

void InitialP(double x, double y,  double z, double *values)
{
  values[0] = 0;
}

// ========================================================================
// exact solution is unknown, using zero
// ========================================================================
void ExactU1(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactU2(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactU3(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactP(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
  double height = TDatabase::ParamDB->DRIFT_Z;
  double Middle = 0.5*height;
  double LH = TDatabase::ParamDB->FS_LH;
  double bheight = (Middle-LH > 0) ? (Middle-LH) : 0;

  cond = FREESURF;

  if(fabs(z-bheight) < 1e-8)
    cond = DIRICHLET;
  else
  {
    if(
        // lower edge of hexagon
        (fabs(y+sqrt(3.0)) < 1e-8)
        // lower right edge of hexagon
     || (fabs(sqrt(3.0)*(x-2)-y) < 1e-8)
        // upper right edge of hexagon
     || (fabs(sqrt(3.0)*(x-2)+y) < 1e-8)
        // upper edge of hexagon
     || (fabs(y-sqrt(3.0)) < 1e-8)
        // upper left edge of hexagon
     || (fabs(sqrt(3.0)*(x+2)-y) < 1e-8)
        // lower left edge of hexagon
     || (fabs(sqrt(3.0)*(x+2)+y) < 1e-8) 
    )
    cond =  SLIP;
  }
}

// value of boundary condition
void U1BoundValue(double x, double y, double z, double &value)
{
  value = 0;     
}

// value of boundary condition
void U2BoundValue(double x, double y, double z, double &value)
{
  value = 0;
}

// value of boundary condition
void U3BoundValue(double x, double y, double z, double &value)
{
  value = 0;
}

// ========================================================================
// coefficients: RE_NR from database
// ========================================================================
void LinCoeffs(int n_points, double *x, double *y, double *z,
               double **parameters, double **coeffs)
{
  double eps = 1/TDatabase::ParamDB->RE_NR;
  double g = TDatabase::ParamDB->FS_G;
  double rho = TDatabase::ParamDB->FS_RHO;
  double L = TDatabase::ParamDB->FS_L; // length scale
  double U = TDatabase::ParamDB->FS_U; // length scale
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = 0; // f1
    coeff[2] = 0; // f2
    coeff[3] = -L*g/(U*U); // f3
  }
}

// ========================================================================
// part for magnetoproblem
// ========================================================================
// exact solution
void ExactPotential(double x, double y, double z, double *values)
{
  double height = TDatabase::ParamDB->DRIFT_Z;
  double Middle = 0.5*TDatabase::ParamDB->DRIFT_Z;
  double H1 = TDatabase::ParamDB->FS_H1;
  double H2 = TDatabase::ParamDB->FS_H2;
  double LH = TDatabase::ParamDB->FS_LH;
  double b = (0.5 - LH/height)>0 ? (0.5 - LH/height) : 0;

  if(z>Middle)
  {
    // above the layer
    values[0] = (z-Middle)*H2;
    values[1] = 0;
    values[2] = 0;
    values[3] = H2;
    values[4] = 0;
  }
  else
  {
    if(z>b*height)
    {
      // in the layer
      values[0] = (z-Middle)*H1;
      values[1] = 0;
      values[2] = 0;
      values[3] = H1;
      values[4] = 0;
    }
    else
    {
      // below the layer
      values[0] = height*(b-0.5)*H1 - (b*height + z)*H2;
      values[1] = 0;
      values[2] = 0;
      values[3] = H2;
      values[4] = 0;
    }
  }
}

// kind of boundary condition (for FE space needed)
void BoundConditionPotential(double x, double y, double z, BoundCond &cond)
{
  double height = TDatabase::ParamDB->DRIFT_Z;

  if(fabs(z)<1e-8 || fabs(z-height)<1e-8)
    cond = DIRICHLET;
  else
    cond = NEUMANN;
}

// value of boundary condition
void BoundValuePotential(double x, double y, double z, double &value)
{
  double height = TDatabase::ParamDB->DRIFT_Z;
  double Middle = 0.5*height;
  double H1 = TDatabase::ParamDB->FS_H1;
  double H2 = TDatabase::ParamDB->FS_H2;
  double LH = TDatabase::ParamDB->FS_LH;
  double b = (0.5 - LH/height)>0 ? (0.5 - LH/height) : 0;

  value = 0;

  if(fabs(z)<1e-8)
    value = height*( H1*(b-0.5) - H2*b );

  if(fabs(z-height)<1e-8)
    value = Middle*H2;
}

void BilinearCoeffsPotentialLangevin(int n_points,
        double *x, double *y, double *z,
        double **parameters, double **coeffs)
{
  double chi0 = TDatabase::ParamDB->FS_CHI0;
  double HM = TDatabase::ParamDB->FS_HM;
  double Ms = TDatabase::ParamDB->FS_MS;
  double gamma = TDatabase::ParamDB->FS_GAMMA;
  int i;
  double *coeff, *param, H, arg;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    H = param[0];

    coeff[0] = 1;
    coeff[1] = 0;
    coeff[2] = 0;
    coeff[3] = 0;
    coeff[4] = 0;
    coeff[5] = 0;

    arg = gamma*H;
    if(arg > 1e-2)
      coeff[6] = Ms/(H*HM) * ( 1/tanh(arg) - 1/arg ) + 1;
    else
      coeff[6] = chi0 + 1;

    coeff[7] = 1; // mu_vacuum
  }
}

void BilinearCoeffsPotentialVislovich(int n_points,
        double *x, double *y, double *z,
        double **parameters, double **coeffs)
{
  double chi0 = TDatabase::ParamDB->FS_CHI0;
  double HM = TDatabase::ParamDB->FS_HM;
  double Ms = TDatabase::ParamDB->FS_MS;
  double gamma = TDatabase::ParamDB->FS_GAMMA;
  double HT = TDatabase::ParamDB->FS_HT;
  int i;
  double *coeff, *param, H, arg;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    H = param[0];

    coeff[0] = 1;
    coeff[1] = 0;
    coeff[2] = 0;
    coeff[3] = 0;
    coeff[4] = 0;
    coeff[5] = 0;

    arg = gamma*H;
    if(arg > 1e-2)
      coeff[6] = 1 + Ms/(HM*H + HT);
    else
      coeff[6] = chi0 + 1;

    coeff[7] = 1; // mu_vacuum
  }
}

void SetMuR(int N_Points, double **Coeffs, double **Params, TBaseCell *cell)
{
  int i, ID;
  double *Coeff;

  ID = cell->GetSubGridID();

  for(i=0;i<N_Points;i++)
  {
    Coeff = Coeffs[i];
    if(ID == 1) // Fluid
      Coeff[0] = Coeff[6];
    else // air/vacuum
      Coeff[0] = Coeff[7];
  } // endfor i
}

void ParametersMF(double *in, double *out)
{
  out[0] = sqrt(in[3]*in[3] + in[4]*in[4] + in[5]*in[5]);
}

int N_FESpacesMF = 1;
int N_FEFunctionsMF = 1;
int N_ParamFctMF = 1;
int N_FEValuesMF = 3;
int FEValuesFctIndexMF[3] = { 0, 0, 0 };
MultiIndex3D FEValuesMultiIndexMF[3] = { D100, D010, D001 };
int N_ParametersMF = 1;
int BeginParametersMF[1] = { 0 };
ParamFct *ParameterFctMF[1] = { ParametersMF };

// !!! LANGEVIN only !!!
void GetMuAndM(TFEFunction3D *Potential, TFEFunction3D *Mu,
               TFEFunction3D *M, TFEFunction3D *HFct)
{
  int i,j,k,l;
  TCollection *Coll, *MuColl;
  int N_Cells, N_MuCells;
  double *potential, *mu, *m, *h;
  int N_PotentialDOF, N_MuDOF;
  TFESpace3D *PotentialSpace, *MuSpace;
  int *counter;
  TBaseCell *cell;
  int CellNumber;
  int *GlobalNumbersPot, *BeginIndexPot;
  int *GlobalNumbersMu, *BeginIndexMu;
  int *DOFPot, *DOFMu;
  FE3D Fe3DPot, Fe3DMu;
  TNodalFunctional3D *nf;
  TBaseFunct3D *bf;
  int N_Points;
  double *xi, *eta, *zeta;
  int N_BaseFuncts;
  double uorig[MaxN_BaseFunctions3D], uxorig[MaxN_BaseFunctions3D];
  double uyorig[MaxN_BaseFunctions3D], uzorig[MaxN_BaseFunctions3D];
  double uref[MaxN_BaseFunctions3D], uxiref[MaxN_BaseFunctions3D];
  double uetaref[MaxN_BaseFunctions3D], uzetaref[MaxN_BaseFunctions3D];
  RefTrans3D RefTrans;
  double u, ux, uy, uz;
  double vals[MaxN_BaseFunctions3D];
  double PointValuesMu[MaxN_PointsForNodal3D];
  double FctalValuesMu[MaxN_PointsForNodal3D];
  double PointValuesM[MaxN_PointsForNodal3D];
  double FctalValuesM[MaxN_PointsForNodal3D];
  double PointValuesH[MaxN_PointsForNodal3D];
  double FctalValuesH[MaxN_PointsForNodal3D];
  double muval, mval, H;
  int N_Fctals;

  double chi0 = TDatabase::ParamDB->FS_CHI0;
  double HM = TDatabase::ParamDB->FS_HM;
  double HT = TDatabase::ParamDB->FS_HT;
  double Ms = TDatabase::ParamDB->FS_MS;
  double gamma = TDatabase::ParamDB->FS_GAMMA;

  double arg;

  PotentialSpace = Potential->GetFESpace3D();
  Coll = PotentialSpace->GetCollection();
  N_Cells = Coll->GetN_Cells();
  potential = Potential->GetValues();
  N_PotentialDOF = Potential->GetLength();

  GlobalNumbersPot = PotentialSpace->GetGlobalNumbers();
  BeginIndexPot = PotentialSpace->GetBeginIndex();

  MuSpace = Mu->GetFESpace3D();
  MuColl = MuSpace->GetCollection();
  N_MuCells = MuColl->GetN_Cells();
  N_MuDOF = Mu->GetLength();

  mu = Mu->GetValues();
  memset(mu, 0, N_MuDOF*SizeOfDouble);
  m = M->GetValues();
  memset(m, 0, N_MuDOF*SizeOfDouble);
  h = HFct->GetValues();
  memset(h, 0, N_MuDOF*SizeOfDouble);

  counter = new int[N_MuDOF];
  memset(counter, 0, N_MuDOF*SizeOfInt);

  GlobalNumbersMu = MuSpace->GetGlobalNumbers();
  BeginIndexMu = MuSpace->GetBeginIndex();

  // every cell in whole collection gets its number
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  // loop over all cells in MuColl
  for(i=0;i<N_MuCells;i++)
  {
    cell = MuColl->GetCell(i);
    CellNumber = cell->GetClipBoard();

    DOFPot = GlobalNumbersPot + BeginIndexPot[CellNumber];
    DOFMu = GlobalNumbersMu + BeginIndexMu[i];

    Fe3DMu = MuSpace->GetFE3D(i, cell);
    nf = TFEDatabase3D::GetNodalFunctional3DFromFE3D(Fe3DMu);
    nf->GetPointsForAll(N_Points, xi, eta, zeta);

    Fe3DPot = PotentialSpace->GetFE3D(CellNumber, cell);
    bf = TFEDatabase3D::GetBaseFunct3DFromFE3D(Fe3DPot);
    N_BaseFuncts = bf->GetDimension();

    RefTrans = TFEDatabase3D::GetRefTrans3D_IDFromFE3D(Fe3DPot);
    TFEDatabase3D::SetCellForRefTrans(cell, RefTrans);

    for(j=0;j<N_BaseFuncts;j++)
      vals[j] = potential[DOFPot[j]];

    for(j=0;j<N_Points;j++)
    {
      bf->GetDerivatives(D000, xi[j], eta[j], zeta[j], uref);
      bf->GetDerivatives(D100, xi[j], eta[j], zeta[j], uxiref);
      bf->GetDerivatives(D010, xi[j], eta[j], zeta[j], uetaref);
      bf->GetDerivatives(D001, xi[j], eta[j], zeta[j], uzetaref);

      TFEDatabase3D::GetOrigValues(RefTrans,
                     xi[j], eta[j], zeta[j], N_BaseFuncts,
                     uref, uxiref, uetaref, uzetaref,
                     uorig, uxorig, uyorig, uzorig);

      u = 0; ux = 0; uy = 0; uz = 0;
      for(k=0;k<N_BaseFuncts;k++)
      {
        u  +=  uorig[k]*vals[k];
        ux += uxorig[k]*vals[k];
        uy += uyorig[k]*vals[k];
        uz += uzorig[k]*vals[k];
      } // endfor k

      H = sqrt(ux*ux + uy*uy + uz*uz);

      if(cell->GetSubGridID() == 1)
      {
        switch(TDatabase::ParamDB->FS_MAGNETLAW)
        {
          case 0: // Langevin
            arg = gamma*H;
            if(arg > 1e-2)
              muval = Ms/(H*HM) * ( 1/tanh(arg) - 1/arg ) + 1;
            else
              muval = chi0 + 1;
          break;

          case 1: // Vislovich
            muval = Ms/ (H*HM+HT)+1;
          break;
        } // end switch MAGNETLAW
      }
      else
      {
        muval = 1;
      }

      mval = (muval-1)*H*HM;

      PointValuesMu[j] = muval;
      PointValuesM[j] = mval;
      PointValuesH[j] = H*HM;
    } // endfor j

    nf->GetAllFunctionals(PointValuesMu, FctalValuesMu);
    nf->GetAllFunctionals(PointValuesM, FctalValuesM);
    nf->GetAllFunctionals(PointValuesH, FctalValuesH);

    N_Fctals = TFEDatabase3D::GetN_BaseFunctFromFE3D(Fe3DMu);
    for(j=0;j<N_Fctals;j++)
    {
      k = DOFMu[j];
      counter[k]++;
      mu[k] += FctalValuesMu[j];
      m[k] += FctalValuesM[j];
      h[k] += FctalValuesH[j];
    } // endfor j
  } // endfor i

  for(i=0;i<N_MuDOF;i++)
    if( (k = counter[i]) > 0)
    {
      mu[i] /= k;
      m[i] /= k;
      h[i] /= k;
    }

  delete counter;
} // end GetMuAndM

