// Navier-Stokes problem, Driven Cavity
// 
// u(x,y) = ?
// p(x,y) = ?

#define __MIXINGLAYERSLIP__

#define U_INFTY 1   

void ExampleFile()
{
  OutPut("Example: MixingLayerSlipOld.h" << endl) ;
}

// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double *values)
{
  double z;

  z = y*56;
  if (z>=0)
    values[0] = U_INFTY * (1-exp(-2*z))/(1+exp(-2*z));
  else
    values[0] = U_INFTY * (exp(2*z)-1)/(exp(2*z)+1);    
  values[0] -= 0.001* U_INFTY *exp(-z*z)*2*y*cos(8*Pi*x)*56*56;
  values[0] -= 0.001* U_INFTY *exp(-z*z)*2*y*cos(20*Pi*x)*56*56;
}

void InitialU2(double x, double y, double *values)
{
  double z;

  z = y*56;
  values[0] = 0.001*U_INFTY*exp(-z*z)*sin(8*Pi*x)*8*Pi;
  values[0] += 0.001*U_INFTY*exp(-z*z)*sin(20*Pi*x)*20*Pi;
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
  cond = SLIP_FRICTION_PENETRATION_RESISTANCE;
  TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION = 1;
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  double t = TDatabase::TimeDB->CURRENTTIME;

  switch(BdComp)
  {
    case 0: value=0;
            break;
    case 1: value=0;
            break;
    case 2: value=0;
            break;
    case 3: value=0;
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  value = 0;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  int i;
  double *coeff, x, y; 
  static double a=1;
  static double nu=1/TDatabase::ParamDB->RE_NR;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    x = X[i];
    y = Y[i];

    coeff[0] = nu;

    coeff[1] = 0;
    coeff[2] = 0;
  }
}

void ComputeVorticiyThickness(TFEFunction2D *Vorticity, double *thickness)
{
  int i, j, k, N_Cells, index, *DOF, N_loc_dofVort, N_Vort, found;
  int *GlobalNumbersVort, *BeginIndexVort;
  const int  MAXLINES=1025;
  TCollection *Coll;
  TBaseCell *cell;
  FE2D CurrentElementVort;
  TFE2D *Element;
  TNodalFunctional2D *nf;
  RefTrans2D RefTrans;
  TFESpace2D *vorticity_space;
  double X_orig[MaxN_PointsForNodal2D], Y_orig[MaxN_PointsForNodal2D];
  double *xi_ref, *eta_ref, AbsDetjkVort[MaxN_QuadPoints_2D];
  double *vort, *global_y_coord, max;
  double y_coord[MAXLINES], aver_vort[MAXLINES];
  int number[MAXLINES], max_lines;
    

  vorticity_space=Vorticity->GetFESpace2D();
  GlobalNumbersVort = vorticity_space->GetGlobalNumbers();
  BeginIndexVort = vorticity_space->GetBeginIndex();
  N_Vort = vorticity_space->GetN_DegreesOfFreedom();
  vort = Vorticity->GetValues(); 
  global_y_coord = new double[N_Vort];

  // get pointer to set of mesh cells which define the fe space
  Coll = vorticity_space->GetCollection();
  // get number of mesh cells
  N_Cells = Coll->GetN_Cells();
   // loop over all mesh cells
  for(i=0;i<N_Cells;i++)
  {
    // get current mesh cell
    cell = Coll->GetCell(i);        
   // compute geometric positions of the fe nodes
    // get id of finite element in current mesh cell
    CurrentElementVort = vorticity_space->GetFE2D(i, cell);
    // get fe from its id
    Element = TFEDatabase2D::GetFE2D(CurrentElementVort);
    // get reference transformation 
    RefTrans = Element->GetRefTransID();
    TFEDatabase2D::SetCellForRefTrans(cell,RefTrans);
    // get pointer to the nodal functionals (fe nodes) of the fe 
    // (in ref. cell)
    nf = Element->GetNodalFunctional2D();
    // get number and coordinates of local dof in ref cell
    // xi_ref, eta_ref are pointers
    nf->GetPointsForAll(N_loc_dofVort, xi_ref, eta_ref);
      
    // get coordinates of fe nodes in original cell
    // input: RefTrans, N_loc_dof, xi_ref, eta_ref
    // output : X_orig, Y_orig - pointers to arrays with coordinates
    //          AbsDetjk - same as above    
    TFEDatabase2D::GetOrigFromRef(RefTrans,N_loc_dofVort, xi_ref, 
                                  eta_ref, X_orig, Y_orig, AbsDetjkVort);
    
    DOF = GlobalNumbersVort + BeginIndexVort[i];
    for (j=0;j<N_loc_dofVort;j++)
    {
      // values computed with FindGradient have to averaged on a periodic boundary !!! 
      index = DOF[j];
      global_y_coord[index] =  Y_orig[j];
    }
  }


  max_lines = 0;
  // loop over all dofs
  for (i=0;i<N_Vort;i++)
  {
    found = 0;
    for (k=0;k<max_lines;k++)
    {
      // coordinate already found
      if (fabs(global_y_coord[i]-y_coord[k])<1e-7)
      {
        found++;
        aver_vort[k]+=vort[index];
        number[k]++;
      }
    }
    if (found) 
      continue;
    // new coordinate
    y_coord[max_lines] = global_y_coord[i];
    aver_vort[max_lines] = vort[i];
    number[max_lines] = 1;
    max_lines++;
    if (max_lines> MAXLINES)
    {
      OutPut("Increase MAXLINES in ComputeVorticiyThickness !!!");
      exit(4711);
    }
  }

  // compute averages and vorticity thickness
  max = -1;
  for (i=0;i<max_lines;i++)
  {
    aver_vort[i]/=number[i];
    if (fabs(aver_vort[i]) > max)
      max = fabs(aver_vort[i]);
  }

  thickness[0] =  2/max;
  delete global_y_coord;
  return;


}

