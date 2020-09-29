// Navier-Stokes problem 
// for example with analytical solution
// comparison of strong no-slip and weak no-slip
// strong no-slip case


#include <math.h>
#include <Constants.h>

void ExampleFile()
{
  OutPut("Example: DueseAxialSymm3D.00.h" << endl) ;
}
// ========================================================================
// exact solution
// ========================================================================
// first component of the velocity u1
// values[0] = u1
// values[1] = \partial u1/ \partial x
// values[2] = \partial u1/ \partial y
// values[3] = Laplacian of u1
// if the functions are not known or if they are not needed, set them to 0

void ExactU1(double x, double y, double *values)
{
  values[0] =1-y*y/0.25;
  values[1] =0;
  values[2] =-8*y;
  values[3] =0;
}

// second component of the velocity u2
// values[0] = u2
// values[1] = \partial u2/ \partial x
// values[2] = \partial u2/ \partial y
// values[3] = Laplacian of u2
// if the functions are not known or if they are not needed, set them to 0

void ExactU2(double x, double y, double *values)
{
  values[0] =0;
  values[1] =0;
  values[2] =0;
  values[3] =0;
}

// pressure p
// values[0] = p
// values[1] = \partial p/ \partial x
// values[2] = \partial p/ \partial y
// values[3] = Laplacian of p
// if the functions are not known or if they are not needed, set them to 0

void ExactP(double x, double y, double *values)
{
  double nu=1/TDatabase::ParamDB->RE_NR;
  double out = 10.0;

  values[0] =-16*nu*(x-out);
  values[1] =-16*nu;
  values[2] =0;
  values[3] =0;
}

// ========================================================================
// boundary conditions
// ========================================================================

// type of the boundary condition
// possible types:
// DIRICHLET, NEUMANN, SLIP_FRICTION_PENETRATION_RESISTANCE
void BoundConditionLower(int BdComp, double t, BoundCond &cond)
{
  // no-slip on the whole boundary
  switch(BdComp)
  {
     case 0:  cond =  SLIP_FRICTION_PENETRATION_RESISTANCE;
	 //  TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION = 1;
      break;
     case 1: cond = NEUMANN;
	 //TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;      
        break;
     case 2: cond = DIRICHLET;
//        OutPut("2" << endl);
        break;
     case 3: cond = DIRICHLET;
        //       OutPut("3" << endl);
        break;
	// b.c. on the interface
     case 4: cond = NEUMANN;
        //      OutPut("4" << endl);
        break;
  }
}

// boundary values of u1
// counting the boundary parts starts at 0
void U1BoundValue(int BdComp, double Param, double &value)
{
   double y;

  switch(BdComp)
  {
    case 0: value=0;
            break;
    case 1: value=0;
            break;
    case 2: value=0;
            break;
      case 3: //value=1-(1-Param)*(1-Param);
	  //break;
	y = (1-Param)*0.5;
	// fluid phase
	if (0.174515<=y)
	    value = 0.4398+6.10325*1e-9*y-1.75921*y*y;
	else
	    // gas phase
	    value = 2.74645 + 4.41189 *1e-10*y - 30.2233 * y*y;
	//OutPut("inflow " << y << " " << value << endl);
	break;
    case 4: value=0;
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
}

// boundary values of u2
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
    case 4: value=0;
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
}

// ========================================================================
// coefficients for the Navier--Stokes equations
// viscosity \nu and right hand side f 
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  // set nu
  double nu_lower=1.0/TDatabase::ParamDB->P5;
  double nu_upper=1.0/TDatabase::ParamDB->P6;
  double eps = 1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, x, y;
  double u1,u2, u1x, u1y, u2x, u2y, u1xx, u1yy, u2xx, u2yy, px, py;
 

  // the coefficients are needed for a set of quadrature points
  // loop over all points
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    x = X[i];
    y = Y[i];

    // coeff[0] is the viscosity
    if (7.1130e-04*x - y + 0.17455 > 0)
       coeff[0] = nu_lower;
    else
       coeff[0] = nu_upper;

    //coeff[0] = eps;
    coeff[1] = 0;
    coeff[2]= 0;
   
  }
}


void ComputeCollections(TCollection *coll_all, TCollection **coll_lower, TCollection **coll_upper)
{
   TBaseCell *cell, **Cells_lower, **Cells_upper, *parent;
   int i, j, n, N_Cells, N_Edges, lower, upper;
   int N_Cells_lower, N_Cells_upper;
   double x_bary, y_bary, x, y;

   lower = 0; 
   upper = 0;
   // compute first number of cells
   N_Cells = coll_all->GetN_Cells();
   for(i=0;i<N_Cells;i++)
   {
    cell = coll_all->GetCell(i);
    // search parent cell on coarsest grid
    parent = cell;
    while (parent->GetN_Parents())
	parent = parent->GetParent();
 
    N_Edges=parent->GetN_Edges();
    x_bary = 0;
    y_bary = 0;
    n = 0;
    for (j=0;j<N_Edges;j++)
    {
      x = parent->GetVertex(j)->GetX();
      y = parent->GetVertex(j)->GetY();
      x_bary += x;
      y_bary += y;
      n++;
    }
    x_bary /= n;
    y_bary /= n;
    if (7.1130e-04*x_bary- y_bary + 0.17455 > 0)
       lower++;
    else
       upper++;
   }
   OutPut("lower " << lower << " upper " << upper << endl);

   // create new collections
   Cells_lower = new TBaseCell*[lower];
   N_Cells_lower = 0;
   Cells_upper = new TBaseCell*[upper];
   N_Cells_upper = 0;
   for(i=0;i<N_Cells;i++)
   {
     cell = coll_all->GetCell(i);
     // search parent cell on coarsest grid
    parent = cell;
    while (parent->GetN_Parents())
	parent = parent->GetParent();
    
    N_Edges=parent->GetN_Edges();
     x_bary = 0;
     y_bary = 0;
     n = 0;
     for (j=0;j<N_Edges;j++)
     {
        x = parent->GetVertex(j)->GetX();
        y = parent->GetVertex(j)->GetY();
        x_bary += x;
        y_bary += y;
        n++;
     }
     x_bary /= n;
     y_bary /= n;
     if (7.1130e-04*x_bary- y_bary + 0.17455 > 0)
     {
        Cells_lower[N_Cells_lower] = cell;
        N_Cells_lower++; 
     }
     else
     {
        Cells_upper[N_Cells_upper] = cell;
        N_Cells_upper++;        
     }
   }
   *coll_lower = new TCollection(N_Cells_lower, Cells_lower);
   *coll_upper = new TCollection(N_Cells_upper, Cells_upper);
   return;
}


void MatchInterfaceDOFs(TFESpace2D *velocity_space_lower,TFESpace2D *velocity_space_upper,
   int *&interface_dof_lower, int *&interface_dof_upper, int *interface_dof)
{
   int i, j, k, l, N_Cells_lower, N_Cells_upper, N_local_dof, N_cells;
   int number, *lower_dof, *upper_dof, counter;
   int *GlobalNumbers_lower, *BeginIndex_lower, *GlobalNumbers_upper;
   int *BeginIndex_upper;
   int i1, j1, k1, N_Active_lower, N_Active_upper, count_cell = 0; 
   TCollection *lower, *upper;
   int *EdgeDOF;
   TBaseCell *cell, *neigh;
   TJoint *joint;
   FE2D FEType0, FEType1;
   TFE2D *FE0, *FE1;
   TFEDesc2D *FEDesc0_Obj, *FEDesc1_Obj;
   
   lower = velocity_space_lower->GetCollection();
   upper = velocity_space_upper->GetCollection();

   N_Active_lower = velocity_space_lower->GetActiveBound();
   N_Active_upper = velocity_space_upper->GetActiveBound();

   N_Cells_lower  = lower->GetN_Cells();
   N_Cells_upper  = upper->GetN_Cells();
   
   // set clip boards
   for(i=0;i<N_Cells_lower;i++)
   {
    cell = lower->GetCell(i);
    cell->SetClipBoard(i);
   } // endfor i
   for(i=0;i<N_Cells_upper;i++)
   {
    cell = upper->GetCell(i);
    cell->SetClipBoard(i);
   } // endfor i
   
   N_local_dof = 0;
   N_cells = 0;
   // loop over all cells for the lower collection
   for(i=0;i<N_Cells_lower;i++)
   {
      cell = lower->GetCell(i);
      // compute number of boundary nodes on an edge
      FEType0 = velocity_space_lower->GetFE2D(i, cell);
      FE0 = TFEDatabase2D::GetFE2D(FEType0);
      FEDesc0_Obj = FE0->GetFEDesc2D();   
      if (FEDesc0_Obj->GetN_JointDOF()> N_local_dof )
         N_local_dof = FEDesc0_Obj->GetN_JointDOF();
      k = cell->GetN_Edges();
      for(j=0;j<k;j++) // loop over all edges of cell
      {
         joint = cell->GetJoint(j);
         // found interface joint
         if(joint->GetType() == InterfaceJoint ||
            joint->GetType() == IsoInterfaceJoint)
            N_cells++;
      }
   }
  
   // allocate arrays for storing the matching information
   number = N_local_dof * N_cells;
   lower_dof = new int[number];
   upper_dof = new int[number];
   GlobalNumbers_lower = velocity_space_lower->GetGlobalNumbers();
   BeginIndex_lower =  velocity_space_lower->GetBeginIndex();
   GlobalNumbers_upper = velocity_space_upper->GetGlobalNumbers();
   BeginIndex_upper =  velocity_space_upper->GetBeginIndex();

   counter = 0;
   // loop over all cells for the lower collection
   for(i=0;i<N_Cells_lower;i++)
   {
      cell = lower->GetCell(i);
      FEType0 = velocity_space_lower->GetFE2D(i, cell);
      FE0 = TFEDatabase2D::GetFE2D(FEType0);
      FEDesc0_Obj = FE0->GetFEDesc2D();   
      N_local_dof = FEDesc0_Obj->GetN_JointDOF();
      k = cell->GetN_Edges();
      for(j=0;j<k;j++) // loop over all edges of cell
      {
         joint = cell->GetJoint(j);
         // found interface joint
         if(joint->GetType() == InterfaceJoint ||
            joint->GetType() == IsoInterfaceJoint)
         {
	    count_cell++;
            // find neighbour cell on upper collection 
            neigh = joint->GetNeighbour(cell);
            // OutPut(cell->GetClipBoard() << " " << neigh->GetClipBoard()<<endl);
            // find local dof on the interface
            EdgeDOF = FEDesc0_Obj->GetJointDOF(j);
            //   for (l=0;l<N_local_dof;l++)
            // OutPut(EdgeDOF[l] << " ");
            // OutPut(endl);
            // compute and store global dofs
            for (l=0;l<N_local_dof;l++)
            {
	      lower_dof[counter+l] = GlobalNumbers_lower[BeginIndex_lower[i]+EdgeDOF[l]];
	      //            OutPut(counter+l<< " " <<lower_dof[counter+l] << " ");
            }
            // OutPut(endl);

            // go to the neighbour above the interface
            // do the same
            i1 = neigh->GetClipBoard();
            FEType1 = velocity_space_upper->GetFE2D(i1,neigh);
            FE1 = TFEDatabase2D::GetFE2D(FEType1);
            FEDesc1_Obj = FE1->GetFEDesc2D(); 
            k1 = neigh->GetN_Edges();
            for(j1=0;j1<k1;j1++) // loop over all edges of cell
            {
               joint = neigh->GetJoint(j1);
               // found interface joint
               if(joint->GetType() == InterfaceJoint ||
                  joint->GetType() == IsoInterfaceJoint)
               {
                  EdgeDOF = FEDesc1_Obj->GetJointDOF(j1);
                  // for (l=0;l<N_local_dof;l++)
                  // OutPut(EdgeDOF[l] << " ");
                  // OutPut(endl);
                  for (l=0;l<N_local_dof;l++)
                  {
                     upper_dof[counter+N_local_dof-1-l] = GlobalNumbers_upper[BeginIndex_upper[i1]+EdgeDOF[l]];
		     //     OutPut(counter+N_local_dof-1-l << " " << upper_dof[counter+N_local_dof-1-l] << " ");
                  }
                  counter +=  N_local_dof;   
                  break;
               }
            }
            break; 
         }
      }
   }
   OutPut(count_cell << " cells on interface " << endl);

   // condens arrays
   counter = 0;
   for (i=0;i <number;i++)
   {
      i1 = 0;
      for (j=0;j<counter;j++)
      {
         if (lower_dof[i]==lower_dof[j])
         {
            i1 = 1;
            if (upper_dof[i]!=upper_dof[j])
            {
               OutPut(i <<" mismatched dof !!! lower " << lower_dof[i]
		      << " uppers " << upper_dof[i] << " " << upper_dof[j] << endl);
	       exit (4711);
            }
         }
      }
      // not yet found and active dof
      if ((!i1)&&
	  ( lower_dof[i] < N_Active_lower) && ( upper_dof[i] < N_Active_upper))
      {
         lower_dof[counter]= lower_dof[i];
         upper_dof[counter]= upper_dof[i];
         counter++;
      }
   }
   
   interface_dof[0] = counter;
   interface_dof_lower = lower_dof;
   interface_dof_upper = upper_dof;
   
   for (i=0;i <counter;i++)
   {
      OutPut(i << " " << lower_dof[i] << " " << upper_dof[i] << endl);
   }
   return;
}

void CheckStress(TFEFunction2D *u1_lower, TFEFunction2D *u2_lower,
  TFEFunction2D *p_lower, TFEFunction2D *u1_upper, TFEFunction2D *u2_upper,
  TFEFunction2D *p_upper)
{
  double t1=1, t2=7.113e-4,n,n1,n2,x,y;
  double h, values[18],sigma11,sigma12,sigma22,stress_lower,stress_upper;
  int number = 100, i;
  double nu_lower=1.0/TDatabase::ParamDB->P5;
  double nu_upper=1.0/TDatabase::ParamDB->P6;
  double eps=1e-4;

  // compute normal and tangential
  n=sqrt(t1*t1+t2*t2);
  t1/=n;
  t2/=n;
  n1=-t2;
  n2=t1;
  h = 10.0/number;

  for (i=0;i<=number;i++)
  {
      x = h*i;
      y = 0.17455 + 7.1130e-04*x;
      y-=eps;
     u1_lower->FindGradient(x,y,values);
      u2_lower->FindGradient(x,y,values+3);
      p_lower->FindGradient(x,y,values+6);
      y+=2*eps;
      u1_upper->FindGradient(x,y,values+9);
      u2_upper->FindGradient(x,y,values+12);
      p_upper->FindGradient(x,y,values+15);
      sigma11 = 2 * nu_lower * values[1] - values[6];
      sigma12 = nu_lower * ( values[1] + values[2]);
      sigma22 = 2*nu_lower * values[2] - values[6];
      //sigma11= -values[6];
      //sigma12 = 0;
      //sigma22 = -values[6];
      stress_lower = t1 * sigma11 * n1 + t1 * sigma12*n2
	  + t2 * sigma12 * n1 + t2 * sigma22 * n2;
      //OutPut(sigma11 << " " << sigma12 << " " << sigma22 << " ");
      sigma11 = 2 * nu_upper * values[10] - values[15];
      sigma12 = nu_upper * ( values[10] + values[11]);
      sigma22 = 2*nu_upper * values[11] - values[15];
      //sigma11= -values[15];
      //sigma12 = 0;
      //sigma22 = -values[15];
      stress_upper = t1 * sigma11 * n1 + t1 * sigma12*n2
	  + t2 * sigma12 * n1 + t2 * sigma22 * n2;
      //OutPut(sigma11 << " " << sigma12 << " " << sigma22 << " ");
      OutPut("stress lower " << stress_lower << " upper " 
      << stress_upper << " diff " << stress_lower - stress_upper << endl);

  }

}

