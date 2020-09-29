// Navier-Stokes problem
// channel flow in 3D
//

#define  __DOWNWIND__

void ExampleFile()
{
  OutPut("Example: Channel3D.h" << endl);
}

// ========================================================================
// exact solution
// ========================================================================
void InitialU1(double x, double y, double z, double *values)
{
  values[0] = y*(1-y)*z*(1-z)*16;
  values[0] = 0;
}

void InitialU2(double x, double y, double z, double *values)
{
  values[0] = 0;
}

void InitialU3(double x, double y, double z, double *values)
{
  values[0] = 0;
}

void InitialP(double x, double y,  double z, double *values)
{
  values[0] = 0;
}

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
  if (fabs(10-x)<1e-5)
  {
    cond = DIRICHLET;
  }
  else
  {
    cond = DIRICHLET;
  }
//  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
}

// value of boundary condition
void U1BoundValue(double x, double y, double z, double &value)
{
  if (fabs(10-x)<1e-5)
    value = 0;
  else
  {  
    if (fabs(x)<1e-5)
      value = y*(1-y)*z*(1-z)*16;
    else
      value = 0;
  }
  value = y*(1-y)*z*(1-z)*16;
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
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y, double *Z,
               double **parameters, double **coeffs)
{
  double eps = 1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, x, y, z, u1, u2, u3, ux, uy, uz;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    coeff[0] = eps;
    coeff[1] = 0;
    coeff[2] = 0;
    coeff[3] = 0;
  }
}



void DownwindNumberingCells(TCollection *Coll, int *downwind)
{
   int i,j,N_Cells, N_V, changes, change_now;
   double x, y, z, *sx, *sy, *sz, tmp;
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
   sz = new double[N_Cells];

   // compute center of mesh cells
   for(i=0;i<N_Cells;i++)
   {
      cell = Coll->GetCell(i);
      cell->SetClipBoard(-1); 
      downwind[i] = i;
      N_V = cell->GetN_Vertices();
      x = 0;
      y = 0;
      z = 0;
      for (j=0;j<N_V;j++)
      {
         x += cell->GetVertex(j)->GetX();
         y += cell->GetVertex(j)->GetY();
         z += cell->GetVertex(j)->GetZ();
      }
      sx[i] = x/N_V;
      sy[i] = y/N_V;
      sz[i] = z/N_V;
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
            case 2:
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
            z = sz[j];
            sz[j] = sz[j+1];
            sz[j+1] = z;
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
   delete sz;
   // exit(1);
   OutPut("downwind numbering finished"<< endl);
   return;
}



