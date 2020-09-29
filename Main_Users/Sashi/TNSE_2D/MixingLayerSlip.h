// Navier-Stokes problem, Driven Cavity
// 
// u(x,y) = ?
// p(x,y) = ?

#define __MIXINGLAYERSLIP__
#define U_INFTY 1   

// ========================================================================
// example file
// ========================================================================

void ExampleFile()
{
  OutPut("Example: MixingLayerSlip.h (correct vorticity thickness,") 
  OutPut(" scale old one by 4), U_INFTY "<< U_INFTY<< endl);
}

// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double *values)
{
  double z, sigma;

  sigma = 1/TDatabase::ParamDB->P8;
  z = 2*y/sigma;
  if (z>=0)
    values[0] = U_INFTY * (1-exp(-2*z))/(1+exp(-2*z));
  else
    values[0] = U_INFTY * (exp(2*z)-1)/(exp(2*z)+1);    
  values[0] -= 0.001* U_INFTY *exp(-z*z)*cos(8*Pi*x)*8*y/(sigma*sigma);
  values[0] -= 0.001* U_INFTY *exp(-z*z)*cos(20*Pi*x)*8*y/(sigma*sigma);
}

void InitialU2(double x, double y, double *values)
{
  double z, sigma;

  sigma = 1/TDatabase::ParamDB->P8;
  z = 2*y/sigma;
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
  int i,j,k,l,found,max_lines,N_Cells,N_Edges,N_Vort;
  TCollection *Coll;
  TBaseCell *cell;
  double x0,x1,x,y0,y1,val[3],max;
  double *global_y_coord, *aver_vort;
  TFESpace2D *vorticity_space;
    

  vorticity_space=Vorticity->GetFESpace2D();
  N_Vort = vorticity_space->GetN_DegreesOfFreedom();
  global_y_coord = new double[N_Vort];
  aver_vort =  new double[N_Vort];

  // get pointer to set of mesh cells which define the fe space
  Coll = vorticity_space->GetCollection();
  // get number of mesh cells
  N_Cells = Coll->GetN_Cells();

  max_lines = 0;
  // loop over all mesh cells
  for(i=0;i<N_Cells;i++)
  {
    // get current mesh cell
    cell = Coll->GetCell(i);        
    // get number of edges
    N_Edges = cell->GetN_Edges();
    // for all vertices
    for (j=0;j<N_Edges;j++)
    {
      // compute coordinate of vertex
      cell->GetVertex(j)->GetCoords(x0, y0);
      // no vertices on lower and upper boundary
      if ((fabs(y0-1)< 1e-6)||(fabs(y0+1)< 1e-6))
        continue;
      // for all edges which has a larger number
      for (k=j+1; k < N_Edges;k++)
      {
        // compute coordinate of vertex
        cell->GetVertex(k)->GetCoords(x1, y1);
        // vertices do not posses the same y coordinate, continue
        if (fabs(y0-y1)> 1e-6)
          continue;
        // compute integral by edge midpoint rule
        x = (x0+x1)/2;
        // compute value of fe function in midpoint of edge
        Vorticity->FindGradientLocal(cell,i,x,y0,val);
        //OutPut(x << " " << y0 <<  " " << val[0] << endl);
        // add value to entry in array
        found = 0;
        for (l=0;l<max_lines;l++)
        {
          // not the correct entry
          if (fabs(y0-global_y_coord[l])> 1e-6)
            continue;
          // correct entry found
          found = 1;
          aver_vort[l] += fabs(x0-x1)*val[0];
        }
        // new entry has to be created
        if (!found)
        {
          global_y_coord[max_lines] = y0;
          aver_vort[max_lines] = fabs(x0-x1)*val[0];
          max_lines++;
          if (max_lines > N_Vort)
          {
            OutPut("Increase N_Vort in ComputeVorticiyThickness !!!"<<endl);
            exit(4711);
          }
        }
      }
    }
  }

  // compute averages and vorticity thickness
  max = -1;
  for (i=0;i<max_lines;i++)
  {
    // each integral is computed twice
    aver_vort[i]/=2.0;
    // OutPut(i << " " << aver_vort[i] << " " << global_y_coord[i] << endl);
    if (fabs(aver_vort[i]) > max)
      max = fabs(aver_vort[i]);
  }

  thickness[0] =  4*U_INFTY/max; // (2*U_INFTY)/(max/2) (divide max by length of
                               // the integration domain) 
  delete global_y_coord;
  delete aver_vort;
  return;


}

