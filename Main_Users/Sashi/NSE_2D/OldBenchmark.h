// DFG benchmark problem, stationary

// ========================================================================
// multi indices used for various things
// ========================================================================
MultiIndex2D StokesATestDerivatives[] = { D10, D01 };
MultiIndex2D StokesAAnsatzDerivatives[] = { D10, D01 };

MultiIndex2D NavierStokesAStabTestDerivatives[] = 
                { D10, D01, D00, D00, D10, D01 };
MultiIndex2D NavierStokesAStabAnsatzDerivatives[] = 
                { D10, D01, D10, D01, D00, D00 };

MultiIndex2D NavierStokesATestDerivatives[] = { D10, D01, D00, D00 };
MultiIndex2D NavierStokesAAnsatzDerivatives[] = { D10, D01, D10, D01 };

MultiIndex2D L2H1ErrorDerivatives[] = { D00, D10, D01 };

MultiIndex2D NavierStokesB1TestDerivatives[] = { D00 };
MultiIndex2D NavierStokesB1AnsatzDerivatives[] = { D10 };
MultiIndex2D NavierStokesB2TestDerivatives[] = { D00 };
MultiIndex2D NavierStokesB2AnsatzDerivatives[] = { D01 };

MultiIndex2D NavierStokesRhsDerivatives[] = { D00 };

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
  switch(BdComp)
  {
    case 0: value = 0;
            break;
    case 1: value=1.2*Param*(1-Param); // 4*0.3
            break;
    case 2: value = 0;
            break;
    case 3: value=1.2*Param*(1-Param); // 4*0.3
            break;
    case 4: value=0;
            break;
    default: cout << "wrong boundary part number" << endl;
  }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  value = 0;
  if(BdComp>4) cout << "wrong boundary part number" << endl;
}

// ========================================================================
// routine for evaluating errors
// ========================================================================
void L2H1Errors(int N_Points, double *X, double *Y, double *AbsDetjk, 
                double *Weights, double **Der, double **Exact,
                double *LocError)
{
  int i;
  double *deriv, *exactval, w, t;

  LocError[0] = 0.0;
  LocError[1] = 0.0;

  for(i=0;i<N_Points;i++)
  {
    deriv = Der[i];
    exactval = Exact[i];
    w = Weights[i]*AbsDetjk[i];

    t = deriv[0]-exactval[0];
    LocError[0] += w*t*t;

    t = deriv[1]-exactval[1];
    LocError[1] += w*t*t;
    t = deriv[2]-exactval[2];
    LocError[1] += w*t*t;
  } // endfor i

  // cout << "LocError[0]: " << LocError[0] << endl;
  // cout << "LocError[1]: " << LocError[1] << endl;
}

// ========================================================================
// bilinear forms for matrix blocks
// ========================================================================
void StokesACoeffs(int n_points, double *x, double *y,
                   double **parameters, double **coeffs)
{
  static double nu=1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = nu;
    coeff[1] = nu;
  }
}

void NavierStokesACoeffs(int n_points, double *x, double *y,
                   double **parameters, double **coeffs)
{
  static double nu=1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, *param;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    coeff[0] = nu;
    coeff[1] = nu;
    coeff[2] = param[0];
    coeff[3] = param[1];
    // coeff[2] = y[i];
    // coeff[3] = x[i];
    
/*
    cout << "param[0]: " << param[0] << endl;
    cout << "y[i]: " << y[i] << endl;
    cout << "param[1]: " << param[1] << endl;
    cout << "x[i] " << x[i] << endl;
    cout << endl;
*/
  }
}

void NavierStokesAStabCoeffs(int n_points, double *x, double *y,
                   double **parameters, double **coeffs)
{
  static double nu=1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, *param;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    coeff[0] = nu;
    coeff[1] = nu;
    coeff[2] = 0.5*param[0];
    coeff[3] = 0.5*param[1];
    coeff[4] = -0.5*param[0];
    coeff[5] = -0.5*param[1];
    // coeff[2] = y[i];
    // coeff[3] = x[i];
    
/*
    cout << "param[0]: " << param[0] << endl;
    cout << "y[i]: " << y[i] << endl;
    cout << "param[1]: " << param[1] << endl;
    cout << "x[i] " << x[i] << endl;
    cout << endl;
*/
  }
}

void NavierStokesB1Coeffs(int n_points, double *x, double *y,
            double **parameters, double **coeffs)
{
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = -1;
  }
}

void NavierStokesB2Coeffs(int n_points, double *x, double *y,
            double **parameters, double **coeffs)
{
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = -1;
  }
}

// ========================================================================
// linear forms for right-hand sides
// ========================================================================
void NavierStokesRhs1Coeffs(int n_points, double *X, double *Y,
                double **parameters, double **coeffs)
{
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    
    coeff[0] = 0;
  }
}

void NavierStokesRhs2Coeffs(int n_points, double *X, double *Y,
                double **parameters, double **coeffs)
{
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    
    coeff[0] = 0;
  }
}

// ========================================================================
// parameter routine
// ========================================================================
void NavierStokesParams(double *in, double *out)
{
  // in[0] = x, in[1] = y
  out[0] = in[2];
  out[1] = in[3];
}

void GetCdCl(double **rhs,
             TBilinearForm **bilinear, TLinearForm **linear,
             TAuxParam2D *Parameters,
             TFEFunction2D *u1, TFEFunction2D *u2,
             TFEFunction2D *p,
             double &cd, double &cl)
{
  int i,j,k,l,l1,l2,l3,n,m, N_UsedElements, N_LocalUsedElements;
  int N_Cells, N_Points, N_Parameters, N_Rows, N_, N_Hanging;
  int N_Test, N_Ansatz;
  int Used[N_FEs2D];
  int N_BaseFunct[N_FEs2D];
  BaseFunct2D BaseFuncts[N_FEs2D];
  TFESpace2D *fespace, *Uspace, *Pspace;
  FE2D *UsedElements, LocalUsedElements[N_FEs2D], CurrentElement;
  FE2D TestElement, AnsatzElement;
  QuadFormula_2D QuadFormula;
  TQuadFormula_2D *qf;
  BaseFunct2D BaseFunct;
  TBaseFunct2D *bf;
  TCollection *Coll;
  TBaseCell *cell;
  TJoint *joint;
  TBoundEdge *boundedge;
  int *UGlobalNumbers, *UBeginIndex;
  int *PGlobalNumbers, *PBeginIndex;
  TFE2D *ele;
  TFEDesc2D *FEDesc_Obj;
  double *weights, *xi, *eta;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  RefTrans2D RefTrans;
  double *Param[MaxN_QuadPoints_2D];
  double local_rhs[MaxN_BaseFunctions2D];
  double *Matrix[MaxN_QuadPoints_2D], *aux;
  double *AuxArray[MaxN_QuadPoints_2D];
  int *DOF, ActiveBound, DirichletBound, end, last;
  int *TestDOF, *AnsatzDOF;
  double *Entries;
  int *ColInd, *RowPtr;
  double *RHS, *MatrixRow;
  double **HangingEntries, **HangingRhs;
  double *CurrentHangingEntries, *CurrentHangingRhs;
  int *HangingRowPtr, *HangingColInd;
  THangingNode *hn, **HangingNodes;
  HNDesc2D HNDesc;
  THNDesc *HNDesc_Obj;
  TBoundComp *BoundComp;
  double t0, t1, t, s;
  int comp;
  BoundCond Cond0, Cond1;
  BoundCondFunct2D *BoundaryCondition;
  BoundValueFunct2D *BoundaryValue;
  TNodalFunctional2D *nf;
  int N_EdgePoints;
  double *EdgePoints;
  double PointValues[MaxN_PointsForNodal2D];
  double FunctionalValues[MaxN_BaseFunctions2D];
  int *EdgeDOF, N_EdgeDOF;
  int N_Vertices, N_Bound, v1, v2;
  TVertex *vert;
  double x, y;
  double v[4];
  double val;
  double *u1vect, *u2vect, *pvect;

  cd = 0.0;
  cl = 0.0;

// ########################################################################
// determine which finite elements are used
// ########################################################################
  memset(Used, 0, N_FEs2D*SizeOfInt);

  Uspace = u1->GetFESpace2D();
  UGlobalNumbers = Uspace->GetGlobalNumbers();
  UBeginIndex = Uspace->GetBeginIndex();
  u1vect = u1->GetValues();
  u2vect = u2->GetValues();
  n = Uspace->GetN_UsedElements(); 
  UsedElements = Uspace->GetUsedElements();
  for(j=0;j<n;j++)
  {
    CurrentElement = UsedElements[j];
    Used[CurrentElement] = 1;

    ele = TFEDatabase::GetFE2D(CurrentElement);
    BaseFuncts[CurrentElement] = ele->GetBaseFunct2D_ID();
    N_BaseFunct[CurrentElement] = ele->GetN_DOF();

    switch(ele->GetRefTransID())
    {
      case TriaAffin   : CurrentElement = C_P1_2D_T_A;
                         break;
      case QuadAffin   : CurrentElement = C_Q1_2D_Q_A;
                         break;
      case QuadBilinear: CurrentElement = C_Q1_2D_Q_M;
                         break;
    }
    Used[CurrentElement] = 1;

    ele = TFEDatabase::GetFE2D(CurrentElement);
    BaseFuncts[CurrentElement] = ele->GetBaseFunct2D_ID();
    N_BaseFunct[CurrentElement] = ele->GetN_DOF();

  }

  Pspace = p->GetFESpace2D();
  PGlobalNumbers = Pspace->GetGlobalNumbers();
  PBeginIndex = Pspace->GetBeginIndex();
  pvect = p->GetValues();
  n = Uspace->GetN_UsedElements(); 
  UsedElements = Pspace->GetUsedElements();
  for(j=0;j<n;j++)
  {
    CurrentElement = UsedElements[j];
    Used[CurrentElement] = 1;

    ele = TFEDatabase::GetFE2D(CurrentElement);
    BaseFuncts[CurrentElement] = ele->GetBaseFunct2D_ID();
    N_BaseFunct[CurrentElement] = ele->GetN_DOF();

    switch(ele->GetRefTransID())
    {
      case TriaAffin   : CurrentElement = C_P1_2D_T_A;
                         break;
      case QuadAffin   : CurrentElement = C_Q1_2D_Q_A;
                         break;
      case QuadBilinear: CurrentElement = C_Q1_2D_Q_M;
                         break;
    }
    Used[CurrentElement] = 1;

    ele = TFEDatabase::GetFE2D(CurrentElement);
    BaseFuncts[CurrentElement] = ele->GetBaseFunct2D_ID();
    N_BaseFunct[CurrentElement] = ele->GetN_DOF();

  }

  N_UsedElements = 0;
  for(i=0;i<N_FEs2D;i++)
    if(Used[i]) N_UsedElements++;

  UsedElements = new FE2D[N_UsedElements];
  j=0;
  for(i=0;i<N_FEs2D;i++)
    if(Used[i])
    {
      UsedElements[j] = (FE2D)i;
      j++;
    } // endif

  cout << "number of used elements: " << N_UsedElements << endl;
  for(i=0;i<N_UsedElements;i++)
    cout << "UsedElements[" << i << "]: " << UsedElements[i] << endl;

// ########################################################################
// calculate values of base functions and derivatives on ref element
// ########################################################################
  // QuadFormula = TQuadFormula_2D::FindQuadFormula_2D(UsedElements);
  for(i=0;i<N_UsedElements;i++)
  {
    CurrentElement = UsedElements[i];
    BaseFunct = BaseFuncts[CurrentElement];
    QuadFormula = TQuadFormula_2D::FindQF_2D(CurrentElement);

    bf = TFEDatabase::GetBaseFunct2D(BaseFunct);
    bf->MakeRefElementData(QuadFormula);
  } // endfor i

  N_Parameters = Parameters->GetN_Parameters();
  aux = new double [MaxN_QuadPoints_2D*N_Parameters];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    Param[j] = aux + j*N_Parameters;

  // 20 <= number of term in bilinear form
  aux = new double [MaxN_QuadPoints_2D*20]; 
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    AuxArray[j] = aux + j*20;

  aux = new double[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  for(j=0;j<MaxN_BaseFunctions2D;j++)
    Matrix[j] = aux+j*MaxN_BaseFunctions2D;

// ########################################################################
// loop over all cells
// ########################################################################
  // all spaces use same Coll
  Coll = Uspace->GetCollection(); 
  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);

    // check whether this cell has at least one vertex on circle
    N_Vertices = cell->GetN_Vertices();
    N_Bound = 0;
    for(j=0;j<N_Vertices;j++)
    {
      vert = cell->GetVertex(j);
      vert->GetCoords(x,y);
      s = x - 0.20;
      t = y - 0.20;
      if(ABS( t*t+s*s-0.0025 )<1e-6)
      {
        if(N_Bound==0) v1=j;
        if(N_Bound==1) v2=j;
        N_Bound++;
      }
    } // endfor j

    if(N_Bound>0)
    {
      if(cell->GetType() == Triangle)
      {
        v[0]=0; v[1]=0; v[2]=0;
        v[v1]=1;
        if(N_Bound==2) v[v2]=1;
      }
      else
      {
        v[0]=0; v[1]=0; v[2]=0; v[3]=0;
        switch(v1)
        {
          case 2: v1=3; break;
          case 3: v1=2; break;
        }
        v[v1]=1;
        
        if(N_Bound==2)
        {
          switch(v2)
          {
            case 2: v2=3; break;
            case 3: v2=2; break;
          }
          v[v2]=1;
        } // endif
      } // endif

      // ####################################################################
      // find local used elements on this cell
      // ####################################################################
      memset(Used, 0, N_FEs2D*SizeOfInt);
      Used[Uspace->GetFE2D(i, cell)] = 1;
      Used[Pspace->GetFE2D(i, cell)] = 1;
  
      switch(cell->GetType())
      {
        case Triangle:   Used[C_P1_2D_T_A] = 1;
                         break;
        case Quadrangle: Used[C_Q1_2D_Q_M] = 1;
                         break;
        case Rectangle:
        case Parallelogram: Used[C_Q1_2D_Q_A] = 1;
                            break;
      }
  
      N_LocalUsedElements = 0;
      memset(LocalUsedElements, 0, SizeOfInt*N_FEs2D);
      j = 0;
      for(k=0;k<N_FEs2D;k++)
        if(Used[k])
        {
          LocalUsedElements[j] = (FE2D)k;
          j++;
        }
      N_LocalUsedElements = j;

      QuadFormula = TQuadFormula_2D::FindLocalQuadFormula_2D
                  (N_LocalUsedElements, LocalUsedElements);
  
      // ####################################################################
      // calculate values on original element
      // ####################################################################
      qf = TFEDatabase::GetQuadFormula_2D(QuadFormula);
      qf->GetFormulaData(N_Points, weights, xi, eta);
      for(j=0;j<N_LocalUsedElements;j++)
      {
        CurrentElement = LocalUsedElements[j];
        ele = TFEDatabase::GetFE2D(CurrentElement);

        RefTrans=TFEDatabase::GetOrigValues(cell, ele, N_Points,
                          xi, eta, 
                          N_BaseFunct[CurrentElement], 
                          BaseFuncts[CurrentElement], QuadFormula);
      } // endfor j
  
      TFEDatabase::GetOrigFromRef(RefTrans,N_Points, xi, eta, 
                                  X, Y, AbsDetjk);
      Parameters->GetParameters(N_Points, cell, i, xi, eta, X, Y, Param); 
  
      // matrix block A
      AnsatzElement = Uspace->GetFE2D(i, cell);
      switch(cell->GetType())
      {
        case Triangle:   TestElement = C_P1_2D_T_A;
                         break;
        case Quadrangle: TestElement = C_Q1_2D_Q_M;
                         break;
        case Rectangle:
        case Parallelogram: TestElement = C_Q1_2D_Q_A;
                            break;
      }
  
      N_Test = N_BaseFunct[TestElement];
      N_Ansatz= N_BaseFunct[AnsatzElement];

      bilinear[0]->GetLocalMatrix(N_Points, weights, AbsDetjk, X, Y,
                                  N_Test, BaseFuncts[TestElement], 
                                  N_Ansatz, BaseFuncts[AnsatzElement], 
                                  Param, AuxArray, Matrix);
  
      DOF = UGlobalNumbers + UBeginIndex[i];

      val = 0.0;
      for(k=0;k<N_Ansatz;k++)
      {
        for(l=0;l<N_Test;l++)
        {
          val += v[l]*u1vect[DOF[k]]*Matrix[l][k]; 
        }
      }
      cd += val;
  
      val = 0.0;
      for(k=0;k<N_Ansatz;k++)
      {
        for(l=0;l<N_Test;l++)
        {
          val += v[l]*u2vect[DOF[k]]*Matrix[l][k]; 
        }
      }
      cl += val;
  
      // matrix blocks B1, B2
      for(j=0;j<2;j++)
      {
        TestElement = Pspace->GetFE2D(i, cell);
        switch(cell->GetType())
        {
          case Triangle:   AnsatzElement = C_P1_2D_T_A;
                           break;
          case Quadrangle: AnsatzElement = C_Q1_2D_Q_M;
                           break;
          case Rectangle:
          case Parallelogram: AnsatzElement = C_Q1_2D_Q_A;
                              break;
        }
  
        N_Test = N_BaseFunct[TestElement];
        N_Ansatz = N_BaseFunct[AnsatzElement];
  
        bilinear[j+1]
                  ->GetLocalMatrix(N_Points, weights, AbsDetjk, X, Y,
                    N_Test, BaseFuncts[TestElement], 
                    N_Ansatz, BaseFuncts[AnsatzElement], 
                    Param, AuxArray, Matrix);
  
        DOF = PGlobalNumbers + PBeginIndex[i];

        val = 0.0;
        for(k=0;k<N_Test;k++)
        {
          for(l=0;l<N_Ansatz;l++)
          {
            val += v[l]*pvect[DOF[k]]*Matrix[k][l]; 
          }
        }
        if(j==0)
        {
          cd += val;
        }
        else
        {
          cl += val;
        }
  
      } // endfor j
  
/*
      for(j=0;j<2;j++)
      {
        switch(cell->GetType())
        {
          case Triangle:   CurrentElement = C_P1_2D_T_A;
                           break;
          case Quadrangle: CurrentElement = C_Q1_2D_Q_M;
                           break;
          case Rectangle:
          case Parallelogram: CurrentElement = C_Q1_2D_Q_A;
                              break;
        }
  
        N_ = N_BaseFunct[CurrentElement];
        linear[j]->GetLocalRhs(N_Points, weights, AbsDetjk, X, Y,
                               N_, BaseFuncts[CurrentElement], 
                               Param, AuxArray, local_rhs);
  
      } // endfor j
*/
    } // N_Bound>0
  } // endfor i

  delete UsedElements;
  delete Param[0];
  delete AuxArray[0];
  delete Matrix[0];

} // end of GetCdCl
