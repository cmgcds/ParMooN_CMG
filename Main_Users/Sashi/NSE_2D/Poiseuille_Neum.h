// Navier-Stokes problem, Poiseuille-Problem
// 
// u(x,y) = (4*y*(1-y), 0)
// p(x,y) = x-1/2
#define  __DOWNWIND__

void ExampleFile()
{
  OutPut("Example: Poiseuille_Neum.h" << endl) ;
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  values[0] = y*(1-y);
  values[1] = 0;
  values[2] = 1-2*y;
  values[3] = -2;
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
  static double eps = 1/TDatabase::ParamDB->RE_NR;

  values[0] = -2*eps*(x-1);
  values[1] = -2*eps;
  values[2] = 0;
  values[3] = 0;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int i, double t, BoundCond &cond)
{
  if (i==1)
  {
    cond = NEUMANN;
    TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
  }
  else
    cond = DIRICHLET;
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
      // case 1:
    case 3: value=Param*(1-Param);
            break;
    default: value = 0;
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


void DownwindNumberingCells(TCollection *Coll, int *downwind)
{
   int i,j,N_Cells, N_V, changes, change_now,tmp;
   double x, y, *sx, *sy;
   TBaseCell *cell;

   OutPut("downwind numbering started"<< endl);
      
   N_Cells = Coll->GetN_Cells();

   if (TDatabase::ParamDB->SC_DOWNWIND_TYPE == 0)
   {
      for(i=0;i<N_Cells;i++)
         downwind[i] = i;
      return;
   }

   sx = new double[N_Cells];
   sy = new double[N_Cells];
 
   // compute center of mesh cells
   for(i=0;i<N_Cells;i++)
   {
      cell = Coll->GetCell(i);
      cell->SetClipBoard(-1); 
      downwind[i] = i;
      N_V = cell->GetN_Vertices();
      x = 0;
      y = 0;
      for (j=0;j<N_V;j++)
      {
         x += cell->GetVertex(j)->GetX();
         y += cell->GetVertex(j)->GetY();
      }
      sx[i] = x/N_V;
      sy[i] = y/N_V;
      //      OutPut(i << " " << sx[i] << " " << sy[i] << " " << sz[i] << endl);
   }

   changes = 1;
   i =0;
   while (changes)
   {
      i++;
      changes = 0;
      for (j=0;j<N_Cells-1;j++)
      {
         change_now = 0;
         switch(TDatabase::ParamDB->SC_DOWNWIND_TYPE)
         {
            case 1:
               // sort for x
               if (sx[j+1] - sx[j] < -1e-8)
                  change_now = 1;
               break;
	       /* case 2:
                // sort for x and y (first variant), cells in the center of the channel come first
               if (sx[j+1] - sx[j] < -1e-8)
                  change_now = 1;
               if ((fabs(sx[j+1]-sx[j]) < 1e-8) && (fabs(sy[j+1]-0.5) < fabs(sy[j]-0.5)-1e-8))
                  change_now = 1;
               if ((fabs(sx[j+1]-sx[j]) < 1e-8) && (fabs(fabs(sy[j+1]-0.5) - fabs(sy[j]-0.5))<1e-8)
                   &&  (fabs(sz[j+1]-0.5) < fabs(sz[j]-0.5)-1e-8))
                  change_now = 1;
               break;
            case 3:
                // sort for x and y (first variant), cells in the center of the channel come last
               if (sx[j+1] - sx[j] < 0)
                  change_now = 1;
               if ((fabs(sx[j+1]-sx[j]) ==0) && (fabs(sy[j+1]-0.5) > fabs(sy[j]-0.5)))
                  change_now = 1;
               if ((fabs(sx[j+1]-sx[j]) ==0) && (fabs(fabs(sy[j+1]-0.5) - fabs(sy[j]-0.5))==0)
                   &&  (fabs(sz[j+1]-0.5) > fabs(sz[j]-0.5)))
                  change_now = 1;
               break;
           case 4:
               // sort for x, only at the beginning of the channel
                if ((sx[j+1] - sx[j] < -1e-8)&&(sx[j+1] <= 0.5 ))
                   change_now = 1;
               break;
                  
	       */
            default :
            {
               OutPut(" no routine for SC_DOWNWIND_TYPE = " << TDatabase::ParamDB->SC_DOWNWIND_TYPE 
                      << " defined" << endl);
               exit(4711);
            }
         }
         if (change_now)
         {
            tmp = downwind[j];
            downwind[j] = downwind[j+1];
            downwind[j+1] = tmp;
            x = sx[j];
            sx[j] = sx[j+1];
            sx[j+1] = x;
            y = sy[j];
            sy[j] = sy[j+1];
            sy[j+1] = y;
            changes++;
         }
      } // endfor j
      OutPut("changes " << i << " " << changes << endl);
   }
   //  OutPut("changes " << i << " " << changes << endl);
//   for (i=0;i<N_Cells;i++)
//   {
//      OutPut("sort " << downwind[i] << endl);
//      OutPut("sort " << i << " " << downwind[i] << " " << sx[i] << " " << sy[i] << " " << sz[i] << endl);
//   }
   

   delete sx;
   delete sy;
   // exit(1);
   OutPut("downwind numbering finished"<< endl);
   return;
}



