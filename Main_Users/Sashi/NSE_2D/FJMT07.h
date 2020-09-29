// Navier-Stokes problem from Friedhelms Habil
// 

void ExampleFile()
{
  OutPut("Example: FJMT07.h" << endl) ;
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{

    values[0] = 2*x*x*(1-x)*(1-x)*y*(1-y)*(1-2*y);
    values[1] = 4*x*(1-x)*y*(1-y)*(1-2*x)*(1-2*y); 
    values[2] = 2*x*x*(1-x)*(1-x)*(1-6*y+6*y*y); 
    values[3] = 4*y*(1-y)*(1-2*y)*(1-6*x+6*x*x)+2*x*x*(1-x)*(1-x)*(12*y-6);
    return;
  values[0] = sin(x)*sin(y);
  values[1] = cos(x)*sin(y);
  values[2] = sin(x)*cos(y);
  values[3] = 0;
}

void ExactU2(double x, double y, double *values)
{

    values[0] = -2*y*y*(1-y)*(1-y)*x*(1-x)*(1-2*x);
    values[1] = -2*y*y*(1-y)*(1-y)*(1-6*x+6*x*x);
    values[2] = -4*x*(1-x)*(1-2*x)*y*(1-y)*(1-2*y);
  values[3] = -2*y*y*(1-y)*(1-y)*(12*x-6)-4*x*(1-x)*(1-2*x)*(1-6*y+6*y*y);
   return;
  values[0] = cos(x)*cos(y);
  values[1] = -sin(x)*cos(y);
  values[2] = -cos(x)*sin(y);
  values[3] =  0;
}

void ExactP(double x, double y, double *values)
{
  values[0] = (x*x*x+y*y*y-0.5)/1000.0;
  values[1] = 3*x*x/1000.0;
  values[2] = 3*y*y/1000.0;
 values[3] = 0;
  return;
  values[0] = 2*(cos(x)*sin(y)-(double)sin(1.0)+(double)sin(1.0)*(double)cos(1.0));
  values[1] = -2*sin(x)*sin(y);
  values[2] = 2*cos(x)*cos(y);
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
  return;
  switch(BdComp)
  {
  case 0: 
    value=0;
    break;
  case 1: 
    value=(double)sin(1.0)*sin(Param);
    break;
  case 2: 
    value=sin(1-Param)*(double)sin(1.0);
    break;
  case 3: 
    value=0;
    break;
  default: cout << "wrong boundary part number" << endl;
    break;
 }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
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
  return;
  switch(BdComp)
  {
  case 0: 
    value=cos(Param);
    break;
  case 1: 
    value=(double)cos(1.0)*cos(Param);
    break;
  case 2: 
    value=cos(1.0-Param)*(double)cos(1.0);
    break;
  case 3: 
    value=cos(1.0-Param);
    break;
  default: cout << "wrong boundary part number" << endl;
    break;
  }
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  double nu=1/TDatabase::ParamDB->RE_NR;
  double sigma = 1/TDatabase::TimeDB->TIMESTEPLENGTH;
  int i;
  double *coeff, x, y;
  double u1, u1x, u1y, u1lap, u2, u2x, u2y, u2lap, px, py;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    
    x = X[i];
    y = Y[i];

    coeff[0] = nu;
    // prescribed solution
    u1 =  2*x*x*(1-x)*(1-x)*y*(1-y)*(1-2*y);
    u1x = 4*x*(1-x)*y*(1-y)*(1-2*x)*(1-2*y);
    u1y = 2*x*x*(1-x)*(1-x)*(1-6*y+6*y*y);
    u1lap = 4*y*(1-y)*(1-2*y)*(1-6*x+6*x*x)+2*x*x*(1-x)*(1-x)*(12*y-6);
    u2 =  -2*y*y*(1-y)*(1-y)*x*(1-x)*(1-2*x);
    u2x = -2*y*y*(1-y)*(1-y)*(1-6*x+6*x*x);
    u2y = -4*x*(1-x)*(1-2*x)*y*(1-y)*(1-2*y);
    u2lap = -2*y*y*(1-y)*(1-y)*(12*x-6)-4*x*(1-x)*(1-2*x)*(1-6*y+6*y*y);
    px = 3*x*x/1000.0;	
    py = 3*y*y/1000.0;
    // convection field 
    coeff[3] = sin(x)*sin(y);
    coeff[4] = cos(x)*cos(y);
    //coeff[3] = x;
    //coeff[4] = -y;
    //coeff[3] = u1;
    //coeff[4] = u2;
   // rhs
  /*u1 = sin(x)*sin(y);
  u1x = cos(x)*sin(y);
  u1y = sin(x)*cos(y);
  u1lap = -2*sin(x)*sin(y);
  u2 = cos(x)*cos(y);
  u2x = -sin(x)*cos(y);
  u2y = -cos(x)*sin(y);
  u2lap = -2*cos(x)*cos(y);
  px = -2*sin(x)*sin(y);
  py = 2*cos(x)*cos(y);
*/
    coeff[1] = -nu*u1lap+coeff[3]*u1x+coeff[4]*u1y+sigma*u1+px;
    coeff[2] = -nu*u2lap+coeff[3]*u2x+coeff[4]*u2y+sigma*u2+py;
  }
}

void ComputeJumpTerm(TFESpace2D *fespace, 
		     TFEVectFunct2D *u,
		     double *rhs, 
		     double *sol, 
		     double *jump_error, 
		     CoeffFct2D *Coeffs)
{
  int i, j, k, ii, N_Cells, *ColInd, *RowPtr, *GlobalNumbers, *BeginIndex;
  int ActiveBound, *DOF, N_Edges, boundedge, locdof, N_U;
  int sold_parameter_type = TDatabase::ParamDB->SOLD_PARAMETER_TYPE;
  double *Entries, val_neigh[3], h, norm_t, x[3], y[3], area;
  double x0, x1, y0, y1, xs, ys, t1, t2, *coeff, *coeff1, jump, fac;
  double phi0_x, phi0_y, phi1_x, phi1_y, phi2_x, phi2_y, n1, n2, maxjump; 
  double sx, sy, tmp, val[3], error, sd_error, nx, ny;
  double valu1_v1[3], valu1_v2[3], valu2_v1[3], valu2_v2[3], hE;
  double valu1_v1_n[3], valu1_v2_n[3], valu2_v1_n[3], valu2_v2_n[3];
  double gamma, jump_u1, jump_u2, eps = 1e-10;
  TBaseCell *cell, *neigh;
  TCollection *coll;
  FE2D CurrentElement;
  TJoint *joint;
  TRefDesc *refdesc;
  TVertex *ver0,*ver1;
  TFEFunction2D *u1, *u2;
  const int *TmpEdVer;
  OutPut("Compute jump term"<<endl);
  // get arrays with the numbering of the dof
  GlobalNumbers = fespace->GetGlobalNumbers();
  BeginIndex = fespace->GetBeginIndex();
  N_U = fespace->GetN_DegreesOfFreedom();
  // get start of dirichlet nodes in dof array
  ActiveBound = fespace->GetActiveBound();
  // get functions
  u1 = u->GetComponent(0);
  u2 = u->GetComponent(1);
  // get collection and number of cells
  coll = fespace->GetCollection();
  N_Cells = coll->GetN_Cells();

  // assign a numbering to the cells
  for(i=0;i<N_Cells;i++)                        // do for all mesh cells
  {                                           // on the finest level 
    cell=coll->GetCell(i);
    cell->SetClipBoard(i);          
  }                                             // endfor i
  error = sd_error = 0.0;
  coeff = new double[5];
  coeff1 = new double[5];
  // loop over all cells for computing the edge stabilization
  for(i=0;i<N_Cells;i++)
  {
    // next cell
    cell = coll->GetCell(i);
    h = cell->GetDiameter();
    // pointer to global indices of dof connected with this cell
    DOF = GlobalNumbers + BeginIndex[i];
    
    // local dofs are arranged as follows
    // local dof 0 on edge0
    // local dof 1 on edge1
    // local dof 2 on edge2

    CurrentElement = fespace->GetFE2D(i, cell);
    if (CurrentElement!=N_P1_2D_T_A) 
    {
	OutPut("Jump terms for element " << CurrentElement <<
	       " not implemented !!!"<< endl);
	exit(4711);
    }
   // get refinement descriptor
    refdesc=cell->GetRefDesc();           
    refdesc->GetShapeDesc()->GetEdgeVertex(TmpEdVer); 
    // # of edges
    N_Edges = cell->GetN_Edges();
    sx = sy = 0;
    // compute center of the mesh cell
   for (j=0;j<N_Edges; j++)
    {
	sx += cell->GetVertex(j)->GetX();
	sy += cell->GetVertex(j)->GetY();
    }
    sx /= N_Edges;
    sy /= N_Edges;

    // compute values in the vertices of the edges
    for (j=0;j<N_Edges; j++)
    {
        // get vertices and length of edge
        ver0=cell->GetVertex(TmpEdVer[2*j]);  // get vertices of face j
        ver1=cell->GetVertex(TmpEdVer[2*j+1]);
	x[j] = ver0->GetX();
	y[j] = ver0->GetY();
        x[(j+1)%3] = ver1->GetX();
        y[(j+1)%3] = ver1->GetY();
        hE = sqrt((x[j]-x[(j+1)%3])*(x[j]-x[(j+1)%3])+(y[j]-y[(j+1)%3])*(y[j]-y[(j+1)%3]));
	nx = (y[(j+1)%3]-y[j])/hE;
	ny = -(x[(j+1)%3]-x[j])/hE;
        // get values of fe function in the vertices of the edge
 	u1->FindGradientLocal(cell, i, x[j], y[j], valu1_v1);
	u1->FindGradientLocal(cell, i, x[(j+1)%3], y[(j+1)%3], valu1_v2);
	u2->FindGradientLocal(cell, i, x[j], y[j], valu2_v1);
	u2->FindGradientLocal(cell, i, x[(j+1)%3], y[(j+1)%3], valu2_v2);
        //u1->FindGradientLocal(cell, i, (x[j]+x[(j+1)%3])/2, (y[j]+y[(j+1)%3])/2, val);
        //OutPut("a " << DOF[j] << " " << val[0] << endl); 
         // find neighbor cell
	joint=cell->GetJoint(j);
	// compute solution (including derivative) in the vertices
	// from point of view of neighbour mesh cell 
	// NO ADAPTIVE MESHES ALLOWED
	neigh=joint->GetNeighbour(cell); // neighbour cell
	if (neigh!=NULL)
	{
           ii =  neigh->GetClipBoard(); 
           u1->FindGradientLocal(neigh, ii, x[j], y[j], valu1_v1_n);
           u1->FindGradientLocal(neigh, ii, x[(j+1)%3], y[(j+1)%3], valu1_v2_n);
           u2->FindGradientLocal(neigh, ii, x[j], y[j], valu2_v1_n);
           u2->FindGradientLocal(neigh, ii, x[(j+1)%3], y[(j+1)%3], valu2_v2_n);
           fac = 0.5;
	 }
	 else
	 {
           // boundary edge
           // continue;
           valu1_v1_n[0] = valu1_v1_n[1] = valu1_v1_n[2] = 0; 
           valu1_v2_n[0] = valu1_v2_n[1] = valu1_v2_n[2] = 0; 
           valu2_v1_n[0] = valu2_v1_n[1] = valu2_v1_n[2] = 0; 
           valu2_v2_n[0] = valu2_v2_n[1] = valu2_v2_n[2] = 0;
           fac = 1;
           //continue; 
	 }
         //OutPut(" | " << valu1_v1[0] - valu1_v1_n[0] <<  " " << valu1_v2[0] - valu1_v2_n[0]  );
         //OutPut(" | " << valu2_v1[0] - valu2_v1_n[0] <<  " " << valu2_v2[0] - valu2_v2_n[0]  );
          // compute jumps of test functions
         // the test function which is 1 on the edge is continuous
	gamma = TDatabase::ParamDB->DELTA1/hE;

         // take test function which is 1 on the next edge (counterclockwise) 
         // this has value 1 on V2 and -1 on V1
         // compute integral on hE
         locdof = DOF[(j+1)%3];
         // check if correct dof, opposit of v1
         u1->FindGradientLocal(cell, i, sx+(sx-x[j])/2, sy+(sy-y[j])/2, val);
         if (fabs(val[0]-sol[locdof])>eps)
         {
            OutPut("Error in assigning dof (1) !!!"<<endl);
            exit(4711);
         }
         u2->FindGradientLocal(cell, i, sx+(sx-x[j])/2, sy+(sy-y[j])/2, val);
         if (fabs(val[0]-sol[locdof+N_U])>eps)
         {
            OutPut("Error in assigning dof (2) !!!"<<endl);
            exit(4711);
         }        
         jump_u1 = valu1_v2[0] - valu1_v2_n[0] - (valu1_v1[0] - valu1_v1_n[0]);
         //  OutPut(" | " << valu1_v2[0] - valu1_v2_n[0]  << " " << jump_u1 );
         jump_u1 *= gamma*hE/6.0;
         //OutPut(jump_u1 << " " );
         jump_u2 = valu2_v2[0] - valu2_v2_n[0] - (valu2_v1[0] - valu2_v1_n[0]);
         //  OutPut(" | " << valu2_v2[0] - valu2_v2_n[0]  << " " << jump_u2 );
         jump_u2 *= gamma*hE/6.0;
         //OutPut(jump_u2 << " " );
         rhs[locdof] += jump_u1;
         rhs[locdof+N_U] += jump_u2;
         //OutPut(" | " << jump_u1 << " " << jump_u2); 

          // take test function which is 1 on the former edge (counterclockwise) 
         // this has value 1 on V1 and -1 on V2
         // compute integral on hE
         locdof = DOF[(j+2)%3];
         u1->FindGradientLocal(cell, i, sx+(sx-x[(j+1)%3])/2, sy+(sy-y[(j+1)%3])/2, val);
         if (fabs(val[0]-sol[locdof])>eps)
         {
            OutPut("Error in assigning dof (1) !!!"<<endl);
            exit(4711);
         }
         u2->FindGradientLocal(cell, i, sx+(sx-x[(j+1)%3])/2, sy+(sy-y[(j+1)%3])/2, val);
         if (fabs(val[0]-sol[locdof+N_U])>eps)
         {
            OutPut("Error in assigning dof (2) !!!"<<endl);
            exit(4711);
         }        
         /*jump_u1 = valu1_v1[0] - valu1_v1_n[0] - (valu1_v2[0] - valu1_v2_n[0]);
         jump_u1 *= gamma*hE/6.0;
         jump_u2 = valu2_v1[0] - valu2_v1_n[0] - (valu2_v2[0] - valu2_v2_n[0]);
         jump_u2 *= gamma*hE/6.0;*/
         rhs[locdof] -= jump_u1;
         rhs[locdof+N_U] -= jump_u2;
         //OutPut(" | " << jump_u1 << " " << jump_u2);
         error += fac*gamma*hE*( (valu1_v2[0] - valu1_v2_n[0])* (valu1_v2[0] - valu1_v2_n[0])+
                  (valu1_v1[0] - valu1_v1_n[0])*(valu1_v1[0] - valu1_v1_n[0])+
                  (valu2_v2[0] - valu2_v2_n[0])* (valu2_v2[0] - valu2_v2_n[0])+ 
                  (valu2_v1[0] - valu2_v1_n[0])*(valu2_v1[0] - valu2_v1_n[0]))/6.0;

	 // jump term from convective term 
	 // get convection field
	 Coeffs(1, &x[j], &y[j] , NULL, &coeff);
	 Coeffs(1, &x[(j+1)%3], &y[(j+1)%3] , NULL, &coeff1);
	 // present edge
         locdof = DOF[j];
	 jump_u1 = (coeff[3]*nx+coeff[4]*ny)*(valu1_v1[0] - valu1_v1_n[0]);
	 jump_u1 += (coeff1[3]*nx+coeff1[4]*ny)*(valu1_v2[0] - valu1_v2_n[0]);
	 rhs[locdof] -= jump_u1*hE/12.0;
	 jump_u2 = (coeff[3]*nx+coeff[4]*ny)*(valu2_v1[0] - valu2_v1_n[0]);
	 jump_u2 += (coeff1[3]*nx+coeff1[4]*ny)*(valu2_v2[0] - valu2_v2_n[0]);
	 rhs[locdof+N_U] -= jump_u2*hE/12.0;
         // take test function which is 1 on the next edge (counterclockwise) 
         // this has value 1 on V2 and -1 on V1
         // compute integral on hE
         locdof = DOF[(j+1)%3];
	 jump_u1 = -(coeff[3]*nx+coeff[4]*ny)*(valu1_v1[0] - valu1_v1_n[0]);
	 jump_u1 += (coeff1[3]*nx+coeff1[4]*ny)*(valu1_v2[0] - valu1_v2_n[0]);
	 rhs[locdof] -= jump_u1*hE/12.0;
	 jump_u2 = -(coeff[3]*nx+coeff[4]*ny)*(valu2_v1[0] - valu2_v1_n[0]);
	 jump_u2 += (coeff1[3]*nx+coeff1[4]*ny)*(valu2_v2[0] - valu2_v2_n[0]);
	 rhs[locdof+N_U] -= jump_u2*hE/12.0;
          // take test function which is 1 on the former edge (counterclockwise) 
         // this has value 1 on V1 and -1 on V2
         // compute integral on hE
         locdof = DOF[(j+2)%3];
	 rhs[locdof] += jump_u1*hE/12.0;
	 rhs[locdof+N_U] += jump_u2*hE/12.0;

	 }// end loop over edges
        // compute sd error with edge mid point rule
     fac = 0;
   for (j=0;j<N_Edges; j++)
    {
	sx = cell->GetVertex(j)->GetX()+cell->GetVertex((j+1)%3)->GetX();
        sx /= 2;
	sy = cell->GetVertex(j)->GetY()+cell->GetVertex((j+1)%3)->GetY();
        sy /= 2;
        u1->FindGradientLocal(cell, i, sx, sy, valu1_v1);
        u2->FindGradientLocal(cell, i, sx, sy, valu2_v1);
        fac += (sin(sx)*sin(sy)*valu1_v1[1]+cos(sx)*cos(sy)*valu1_v1[2])*(sin(sx)*sin(sy)*valu1_v1[1]+cos(sx)*cos(sy)*valu1_v1[2]);
        fac += (sin(sx)*sin(sy)*valu2_v1[1]+cos(sx)*cos(sy)*valu2_v1[2])*(sin(sx)*sin(sy)*valu2_v1[1]+cos(sx)*cos(sy)*valu2_v1[2]);
    }
    sd_error += TDatabase::ParamDB->DELTA0*h*h*fac* cell->GetMeasure()/3.0;
        
    }// end loop over cells
   jump_error[0] = error;   
   jump_error[1] = sd_error;   
   delete coeff;
   delete coeff1;
}
