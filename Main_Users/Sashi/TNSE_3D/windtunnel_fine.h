// time-dependent NSE
// simulation of wind tunnel experiments
// u(x,y) = unknown
// p(x,y) = unknown

#define __WINDTUNNEL__

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstddef>
using namespace std;
#include <Constants.h>
#include <Windtunnel_3d4d.h>
//#include <Examples/TNSE_3D/windtunnel_log_normal_parameters.h>

void ExampleFile()
{
  TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY = 1506;
  TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION_IDENTITY = 5;

  OutPut("Example: windtunnel_fine.h " << endl);

  TDatabase::ParamDB->DRIFT_Z = 0.18;
  OutPut("Setting DRIFT_Z to " << TDatabase::ParamDB->DRIFT_Z << endl);

  TDatabase::ParamDB->N_CELL_LAYERS = 18;
  OutPut("Setting N_CELL_LAYERS to " << TDatabase::ParamDB->N_CELL_LAYERS  << endl);

  OutFile << "WINDTUNNEL_DIM_Y: " <<  TDatabase::ParamDB->WINDTUNNEL_DIM_Y << endl;
  OutFile << "WINDTUNNEL_DIM_Z: " <<  TDatabase::ParamDB->WINDTUNNEL_DIM_Z << endl;
  OutFile << "WINDTUNNEL_DIM_R: " <<  TDatabase::ParamDB->WINDTUNNEL_DIM_R << endl;
  OutFile << "WINDTUNNEL_ENVIR_COND: " << TDatabase::ParamDB->WINDTUNNEL_ENVIR_COND << endl;
  OutFile << "WINDTUNNEL_SUPERSAT: " << TDatabase::ParamDB->WINDTUNNEL_SUPERSAT << endl;
  OutFile << "WINDTUNNEL_U_INFTY: " <<  TDatabase::ParamDB->WINDTUNNEL_U_INFTY<< endl;
  OutFile << "WINDTUNNEL_L_INFTY: " << TDatabase::ParamDB->WINDTUNNEL_L_INFTY << endl;
  OutFile << "WINDTUNNEL_R_MIN: " << TDatabase::ParamDB->WINDTUNNEL_R_MIN << endl;
  OutFile << "WINDTUNNEL_R_INFTY: " << TDatabase::ParamDB->WINDTUNNEL_R_INFTY << endl;
  OutFile << "WINDTUNNEL_F_INFTY: " << TDatabase::ParamDB->WINDTUNNEL_F_INFTY << endl;
  OutFile << "WINDTUNNEL_kinematic_viscosity: " << TDatabase::ParamDB->WINDTUNNEL_kinematic_viscosity << endl;
}


void InitialU1(double x, double y, double z, double *values)
{
  values[0] = 3;
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


// ========================================================================
// exact solution
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
  double eps = 1e-6;

  if (fabs(x) < eps)
    cond = DIRICHLET;
  else
  {
    if (fabs(x-0.5) < eps)
    {
      cond = NEUMANN;
      TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
    }
    else
    {
      cond  =  SLIP_FRICTION_PENETRATION_RESISTANCE;
      TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION = 1;
    }
  }
}


// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// value of boundary condition
//boundary values of u1
void U1BoundValue(double x, double y, double z, double &value)
{
  double eps = 1e-6;
  int coord_y, coord_z;
  double dist_y, dist_z, y1, z1;
  double mean_value, stand_dev;

  value = 0;
  stand_dev = 0;
  // on the inlet
  if(fabs(x)< eps)
  {
    // values y,z are in a range of -0.225 - 0.225 or -0.18 - 0
    // experimental y,z are in the range of -22.5 22.5 or 18.0 - 0 (centimeters)
    // the inlet values are order 
    // (y_0,z_0) (y_0,z_1) ...
    // (y_1,z_0) (y_1,z_1) ... 
    //for this reason they must be converted in cordinates between 0-45 or 0-18
    //since the experimental data are given in centimeters
    y1 = y + 0.225;
    z1 = z + 0.18;
   
    coord_y = (int) (y1*100 + 0.0001);             //left neighbour
    coord_z = (int) (z1*100 + 0.0001);
 
    dist_y = fabs((double) coord_y - (y1 * 100.)); //distance between coordinate and left neighbour
    dist_z = fabs((double) coord_z - (z1 * 100.));
  
    //bilinear interpolation
    mean_value = (1-dist_y) * (1-dist_z) * TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL[coord_y][coord_z][0]
      +  dist_y    * (1-dist_z) *  TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL[coord_y+1][coord_z][0]
      + (1-dist_y) *  dist_z    *  TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL[coord_y][coord_z+1][0]
      +  dist_y    *  dist_z    *  TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL[coord_y+1][coord_z+1][0];

   /*  variance = (1-dist_y) * (1-dist_z) * (1-dist_y) * (1-dist_z) * TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL[coord_y][coord_z][1]
      +  dist_y    * (1-dist_z) *  dist_y    * (1-dist_z) *  TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL[coord_y+1][coord_z][1]
      + (1-dist_y) *  dist_z    * (1-dist_y) *  dist_z    *  TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL[coord_y][coord_z+1][1]
      +  dist_y    *  dist_z    *  dist_y    *  dist_z    *  TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL[coord_y+1][coord_z+1][1];*/
stand_dev = (1-dist_y) * (1-dist_z)* TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL[coord_y][coord_z][1]
      +  dist_y    * (1-dist_z) *   TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL[coord_y+1][coord_z][1]
      + (1-dist_y) *  dist_z    *   TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL[coord_y][coord_z+1][1]
      +  dist_y    *  dist_z    *   TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL[coord_y+1][coord_z+1][1];
    //OutPut(y << " " << z << " " << variance << ":");

    //stand_dev =sqrt(variance);
    //normal distrubution in each point: the mean value as expected value
    //and variance
    //OutPut("rand " << mean_value << " " );
    value = mean_value + normal_rand()* stand_dev;
    //if (value > 3.2)
    //OutPut(" inlet " << y << " "  <<  y  << " "  << z   <<  " " << z <<  " : " << mean_value << " " << stand_dev << 
    //	   " " << value <<  endl);
    value /= TDatabase::ParamDB->WINDTUNNEL_U_INFTY;
    //OutPut(normal_rand() << " " << value << endl);
  }
}


// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

// value of boundary condition
//boundary value of u2
void U2BoundValue(double x, double y, double z, double &value)
{
  value = 0;
}


// value of boundary condition
//boundary value of u3
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
  double eps;
  int i;
  double *coeff, x, y, z;

  eps = TDatabase::ParamDB->WINDTUNNEL_kinematic_viscosity/
    (TDatabase::ParamDB->WINDTUNNEL_L_INFTY*TDatabase::ParamDB->WINDTUNNEL_U_INFTY);

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    x = X[i];
    y = Y[i];
    z = Z[i];

    coeff[0] = eps;
    // f1
    coeff[1] = 0;
    // f2
    coeff[2] = 0;
    // f3
    coeff[3] = 0;
  }
}


void ReadExperimentalBoundaryConditions()
{
  string line, entry;
  size_t found;
  int dim_y = TDatabase::ParamDB->WINDTUNNEL_DIM_Y;
  int dim_z = TDatabase::ParamDB->WINDTUNNEL_DIM_Z;
  int data_number = dim_y*dim_z;
  double A[data_number][12];
  int counter, i,n;
  int coord_y, coord_z; 
 
  string name_data_file = "../WINDTUNNEL_EXP_DATA/M2/M2x0u.txt";
 // string name_data_file = "../Examples/TNSE_3D/x=0_u.txt";
// cout <<  name_data_file << endl;
  ifstream data_file((name_data_file).c_str());
  if (!data_file)
  {
    cout << " cannot open data-file: " << name_data_file << endl;
    exit(-1);
  }
  char read_line[600];

  //loop counts the number of lines in the txt-file
  /* while( data_file.getline(read_line,200))
   {
     line = read_line;

     if ( line != "\0"  )
     {
       NumberOFLines++;
     }

     line = "";
   }
  data_file.close();
  //number of data sets: Number of Lines minus header (5lines)
  NumberOFLines= NumberOFLines-5;
  //open data file
  ifstream data_file_1((name_data_file).c_str());

  OutPut("NOL: "  <<  NumberOFLines << endl);*/

  n=0;
  //read datas from the file in the array A
  while( data_file.getline(read_line,600))
  {
    i=0;
    line=read_line;
     //ignore the header of the file
    if(n>4)
    { 
      // replace ": and , " by  "."
      found=line.find_first_of(":,");
      while (found!=string::npos)
      {//cout << read_line <<endl;
        line[found]='.';
        found=line.find_first_of(":,",found+1);
      }
      //devide the string in its entries
      found = line.find("\t");
      while (found!=string::npos)
      {
        entry = line;
        counter = entry.length();
        entry.erase(found,counter);
        line.erase(0,found+1);
        stringstream str;
        str << entry;
        str >> A[n-5][i];
        found=line.find("\t");
        i++;
      }
      //last entry of the line
      stringstream str;
      str << line;
      str >> A[n-5][i];
    }n++;
  }
  data_file.close();

  for(n=0;n<data_number;n++)
  {                                               //coordinates of the sorted array
    coord_y=  (int) ((A[n][1]+0.01+225)/10);
    coord_z = (int) (((180-A[n][2])+0.01)/10);

    TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL[coord_y][coord_z][0] = A[n][7];
                                                  //stand_deviation
    TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL[coord_y][coord_z][1] = A[n][8];
  }

  //boundary
  for (i=0;i<dim_y;i++)
  {
    TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL[i][dim_z][0] = TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL[i][dim_z-1][0];
    TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL[i][dim_z][1] = TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL[i][dim_z-1][1];
  }

  for (i=0;i<dim_z+1;i++)
  {
    TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL[dim_y][i][0] = TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL[dim_y-1][i][0];
    TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL[dim_y][i][1] = TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL[dim_y-1][i][1];
  }

   /* for(i=0;i<dim_y;i++)
     {for( int j=0;j< dim_z;j++)
      {
      cout << TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL[i][j][1] <<" ";
    } cout << endl;
   }*/

}
//read the velocity of the
void ReadExperimentalVelocityDrops(double ***diff_velo_air_drops)
{
  string line, entry;
  size_t found;
  int dim_x = TDatabase::ParamDB->WINDTUNNEL_LAYER_NUMBER_X;
  int dim_y = TDatabase::ParamDB->WINDTUNNEL_DIM_Y;
  int dim_z = TDatabase::ParamDB->WINDTUNNEL_DIM_Z;
  int data_number = dim_y*dim_z;
  double A[data_number][19];
  int counter, i,j,n, iter;
  int coord_y, coord_z;
 // name of the 3 data files 
 string name_data_file[3]={"../WINDTUNNEL_EXP_DATA/M2/M2x0Mom.txt","../WINDTUNNEL_EXP_DATA/M2/M2x200Mom.txt","../WINDTUNNEL_EXP_DATA/M2/M2x400Mom.txt"};
//name_data_file[3]={"../Examples/TNSE_3D/M1x0Mom.txt","../Examples/TNSE_3D/M1x200Mom.txt","../Examples/TNSE_3D/M1x400Mom.txt"};
  for(iter=0;iter<dim_x;iter++)
  {
 
  ifstream data_file((name_data_file[iter]).c_str());
  if (!data_file)
  {
    cout << " cannot open data-file: " << name_data_file[iter] << endl;
    exit(-1);
  }
  char read_line[500];

 
  n=0;
  //read datas from the file in the array A
  while( data_file.getline(read_line,500))
  {
    i=0;
    line=read_line;
   // cout << line <<endl;
    //ignore the header of the file
    if(n>4)
    {
      // replace ": and , " by  "."
      found=line.find_first_of(":,");
      while (found!=string::npos)
      {
        line[found]='.';
        found=line.find_first_of(":,",found+1);
      }
      //devide the string in its entries
      found = line.find("\t");
      while (found!=string::npos)
      {
        entry = line;
        counter = entry.length();
        entry.erase(found,counter);
        line.erase(0,found+1);
        stringstream str;
        str << entry;
        str >> A[n-5][i];
        found=line.find("\t");
        i++;
      }
      //last entry of the line
      stringstream str;
      str << line;
      str >> A[n-5][i];
    }n++;
  }

  data_file.close();

  for(n=0;n<data_number;n++)
  {                                               //coordinates of the sorted array
    coord_y=  (int) ((A[n][1]+0.01+225)/10);
    coord_z = (int) (((180-A[n][2])+0.01)/10);
//store the data in the array
 TDatabase::ParamDB->WINDTUNNEL_DROP_VELO[iter][coord_y][coord_z]= A[n][7];
   }

//mirroring of the boundaries to garantie a correct treatment during the interpolation process
  //"top"
 for (i=0;i<dim_y;i++)
  {
    for (j=0;j<dim_z;j++)
    {
     TDatabase::ParamDB->WINDTUNNEL_DROP_VELO[dim_x][i][j] = TDatabase::ParamDB->WINDTUNNEL_DROP_VELO[dim_x-1][i][j];
     
    }
  }

  for (i=0;i<dim_y;i++)
  {
    for (j=0;j<dim_x+1;j++)
    {
      TDatabase::ParamDB->WINDTUNNEL_DROP_VELO[j][i][dim_z] =
        TDatabase::ParamDB->WINDTUNNEL_DROP_VELO[j][i][dim_z-1];
        
    }
  }

  for (i=0;i<dim_z+1;i++)
  {
    for (j=0;j<dim_x+1;j++)
    {
      TDatabase::ParamDB->WINDTUNNEL_DROP_VELO[j][dim_y][i] =
        TDatabase::ParamDB->WINDTUNNEL_DROP_VELO[j][dim_y-1][i]; 
      }
  }
}//end iter
//difference between airflow and droplet velocity 
//difference is stored in diff_velo_air_drops
for(i=0;i<=dim_x;i++)
    for(j=0;j<=dim_y;j++)
        for(n=0;n<=dim_z;n++)
            {
             diff_velo_air_drops[i][j][n]=
               TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL[j][n][0]-
               TDatabase::ParamDB->WINDTUNNEL_DROP_VELO[i][j][n];
            }
              
/*double a=0.;
int b=0;
    for(i=0;i<dim_y;i++)
     {for(  j=0;j<dim_z;j++)
      {a=a+TDatabase::ParamDB->WINDTUNNEL_DROP_VELO[1][i][j]; b++;
      cout << TDatabase::ParamDB->WINDTUNNEL_DROP_VELO[1][i][j]<<" ";
    } cout << endl;
   }

exit(-1);*/
 }

void ReadExperimentalBoundaryConditionsDrops()
{
  char read_line[200];
  string line, entry, prevline;
  size_t found;
  int dim_y, dim_z, dim_r;
  dim_y = TDatabase::ParamDB->WINDTUNNEL_DIM_Y;
  dim_z = TDatabase::ParamDB->WINDTUNNEL_DIM_Z;
  dim_r = TDatabase::ParamDB->WINDTUNNEL_DIM_R;
  int i,j,r, counter;
  int coord_y, coord_z, coord_r;
  double help[3], A[dim_y][dim_z][dim_r][2], number_density[dim_y][dim_z];
  double eps = 1e-8, val;

  TDatabase::ParamDB->WINDTUNNEL_R_INFTY_EXPERIMENT = 0;

  for (i=0; i< dim_y; i++)
    for(j=0; j<dim_z; j++)
      for(r=0; r<dim_r; r++)
      {
        A[i][j][r][0] = 0;
        A[i][j][r][1] = 0;
      }

  // file for the time-averaged boundary data
  string name_data_file = "../Examples/TNSE_3D/concentration_d.txt";
  // open file
  ifstream data_file((name_data_file).c_str());

  if (!data_file)
  {
    OutPut("cannot open data-file for boundary data: " << name_data_file << endl);
    exit(-1);
  }

  prevline = "init";
  //loop counts the number of lines in the txt-file
  while( data_file.getline(read_line,200))
  {
    i=0;
    line = read_line;
    found = line.find("\t");
    while (found!=string::npos)
    {
      entry = line;
      counter = entry.length();
      entry.erase(found,counter);
      line.erase(0,found+1);
      //storage into help
      stringstream str;
      str << entry;
      str >> help[i];
      found=line.find("\t");
      i++;
    }
    //last entry of the line
    stringstream str;
    str << line;
    str >> help[i];

    //distinguish: coordinate and normal entry
    //coordinate

    if (prevline.compare(0,8,"--------") == 0)
    {                                             //cout << "hallo" <<endl;
      coord_y=  (int) ((help[1]+225)/10+eps);
      //coord_z = (int) ((help[2])/10+eps);
      coord_z = (int) ((180-help[2])/10+eps);
      if(coord_y > dim_y || coord_z > dim_z)
      {
        OutPut("error in dimension of experimental boundary condition ")
          exit(4711);
      }
    }
    else
    {
      // normal entry
      if(i==2 && prevline.compare(0,8,"--------") != 0)
      {
        //conversion r-coordinate
        // be careful in the file the values are given as diameters not as radius
	val = (double) (help[0]);
	if (val >  TDatabase::ParamDB->WINDTUNNEL_R_INFTY*1e6)
	  continue;
	if (val > TDatabase::ParamDB->WINDTUNNEL_R_INFTY_EXPERIMENT)
	    TDatabase::ParamDB->WINDTUNNEL_R_INFTY_EXPERIMENT = val;
        coord_r = (int)(help[0] + 0.1 -5)/10;
        if(coord_r > dim_r)
        {
          OutPut("error in dimension of experimental boundary condition ")
            exit(4711);
        }
        A[coord_y][coord_z][coord_r][0]= help[1];
        A[coord_y][coord_z][coord_r][1]= help[2];
      }
    }
    prevline=read_line;
  }
  data_file.close();

  //probably not useful
  /* for(i=0;i<dim_y;i++)
    {
      for( int j=0;j<dim_z;j++)
      {
        number_density[i][j]=0;
        for (int r=0;r<dim_r;r++)
        {
          number_density[i][j]+= A[i][j][r];
        }
      }
    }*/
 //mirroring of the boundaries to garantie a correct treatment during the interpolation process
  //"top"
 for (i=0;i<dim_y;i++)
  {
    for (j=0;j<dim_z;j++)
    {
     TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL_DROPS[i][j][dim_r][0] = TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL_DROPS[i][j][dim_r-1][0];
     TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL_DROPS[i][j][dim_r][1] = TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL_DROPS[i][j][dim_r-1][1];
    }
  }

  for (i=0;i<dim_y;i++)
  {
    for (j=0;j<dim_r+1;j++)
    {
      TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL_DROPS[i][dim_z][j][0] =
        TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL_DROPS[i][dim_z-1][j][0];
        TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL_DROPS[i][dim_z][j][1] =
        TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL_DROPS[i][dim_z-1][j][1];
    }
  }

  for (i=0;i<dim_z+1;i++)
  {
    for (j=0;j<dim_r+1;j++)
    {
      TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL_DROPS[dim_y][i][j][0] =
        TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL_DROPS[dim_y-1][i][j][0]; 
       TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL_DROPS[dim_y][i][j][1] =
        TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL_DROPS[dim_y-1][i][j][1];
    }
  }
 /* for (i=0; i< dim_y; i++)
    for(j=0; j<dim_z; j++)
    {
        TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL_DROPS[i][j][0][0]= 0;
	TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL_DROPS[i][j][0][1]= 0;
      for(r=0; r<dim_r; r++)
      {
        TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL_DROPS[i][j][r+1][0]= A[i][j][r][0];
	//variance instead of standard deviation
        TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL_DROPS[i][j][r+1][1]= A[i][j][r][1]*A[i][j][r][1];
      }
    }*/
  // convert micrometer to meters
  TDatabase::ParamDB->WINDTUNNEL_R_INFTY_EXPERIMENT/=1e+6;
  OutPut("largest experimental drop diameter " << TDatabase::ParamDB->WINDTUNNEL_R_INFTY_EXPERIMENT 
	 << " m"<<endl);
}


void accumulate_bound_condition_drops(int Ny, int Nz, int Nr)
{
  int dim_y = TDatabase::ParamDB->WINDTUNNEL_DIM_Y;
  int dim_z = TDatabase::ParamDB->WINDTUNNEL_DIM_Z;
  int dim_r = TDatabase::ParamDB->WINDTUNNEL_DIM_R;
  
  dim_r = 17;

  double acc[Ny+1][dim_z][dim_r][2],acc1[Ny+1][Nz+1][dim_r][2],acc2[Ny+1][Nz+1][Nr+1][2];
  int comp, firstcomp, i, j,r ,delta;
  comp=0;
  //accumulation of values in y-direction
  for(r=0;r<dim_r;r++)
    for(j=0;j<dim_z;j++)
  {
    comp=0;
    for(i=0;i<Ny+1;i++)
    {
      acc[i][j][r][0]=0.;
      acc[i][j][r][1]=0.;
      firstcomp=comp;
                                                  //computiation of the new partition, linear transform from dim_y intervals to N_y+1 intervals
      while(comp<(double)(i+1)*(double)(dim_y/(Ny+1.)))
      {
        acc[i][j][r][0]+=TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL_DROPS[comp][j][r][0];
        acc[i][j][r][1]+=TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL_DROPS[comp][j][r][1];
        comp++;
      }
      delta=comp-firstcomp;                       //delta:counts the summands in the sum
      acc[i][j][r][0]=acc[i][j][r][0]/delta;      // arithmetic mean
      acc[i][j][r][1]=acc[i][j][r][1]/delta;
      //cout <<"acc:"<< acc[i][j][r]<< "delta " <<delta<<" "<< "i "<<i<< "r " << "z" << j<< endl; ;
    }
  }
  //accumulation of values in z-direction
  for(r=0;r<dim_r;r++)
    for(i=0;i<Ny+1;i++)
  {
    comp=0;
    for(j=0;j<Nz+1;j++)
    {
      acc1[i][j][r][0]=0;
      acc1[i][j][r][1]=0;
      firstcomp=comp;
      while(comp<(double)(j+1)*(double)(dim_z/(Nz+1.)))
      {
        acc1[i][j][r][0]+=acc[i][comp][r][0];
        acc1[i][j][r][1]+=acc[i][comp][r][1];
        comp++;
      }
      delta=comp-firstcomp;
      acc1[i][j][r][0]=acc1[i][j][r][0]/delta;
      acc1[i][j][r][1]=acc1[i][j][r][1]/delta;
    }
  }
  //accumulation of values in r-direction
  for(i=0;i<Ny+1;i++)
    for(j=0;j<Nz+1;j++)
  {
    comp=0;
    for(r=0;r<Nr+1;r++)
    {
      acc2[i][j][r][0]=0;
      acc2[i][j][r][1]=0;
      firstcomp=comp;
      while(comp<(double)(r+1)*(double)(dim_r/(Nr+1.)))
      {
        acc2[i][j][r][0]+=acc1[i][j][comp][0];
        acc2[i][j][r][1]+=acc1[i][j][comp][1];
        comp++;
      }
      delta=comp-firstcomp;
      acc2[i][j][r][0] = acc2[i][j][r][0]/delta;
      acc2[i][j][r][1] = acc2[i][j][r][1]/delta;
    }
  }

  /* for(i=0;i<Ny+1;i++)
  for(j=0;j<Nz+1;j++)
     for(r=0;r<Nr+1;r++)
     {
   TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL_DROPS[i][j][r][0]=acc2[i][j][r][0];
   TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL_DROPS[i][j][r][1]=acc2[i][j][r][1];
   //if (i==Ny)
   //  OutPut(i << " " << acc2[i][j][r][0] << " " << acc2[i][j][r][1] << endl);
     } */
  /*  for(i=0;i<Ny+1;i++)
  for(j=0;j<Nz+1;j++)
           for(r=Nr+1; r<dim_r; r++)
           {
           TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL_DROPS[i][j][r][0]=0.;
           TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL_DROPS[i][j][r][1]=0.;
            } */

  /*for(i=0;i<Ny+1;i++)
       {OutPut(endl);
    for(j=0;j<Nz+1;j++)
             for(r=0; r<Nr+1; r++)
              {if (i==2 && j==1)
           OutPut( " " <<acc2[i][j][r][0] );
                }}*/

}


void accumulate_bound_condition_drops_neu(int N_y, int N_z, int N_r, double *a_layers_coord)
{
  int dim_y = TDatabase::ParamDB->WINDTUNNEL_DIM_Y;
  int dim_z = TDatabase::ParamDB->WINDTUNNEL_DIM_Z;
  int dim_r = TDatabase::ParamDB->WINDTUNNEL_DIM_R;
  double r_infty = TDatabase::ParamDB->WINDTUNNEL_R_INFTY;
  double r_min =  TDatabase::ParamDB->WINDTUNNEL_R_MIN;
  double r_max_exp = TDatabase::ParamDB->WINDTUNNEL_R_INFTY_EXPERIMENT;
  double interval_bound, interval_bound_old, interval, h;
  double acc[N_y+1][dim_z][dim_r][2],acc1[N_y+1][N_z+1][dim_r][2],acc2[N_y+1][N_z+1][N_r+1][2];
  int comp, firstcomp, i, j,r ,delta;
  comp=0;
  r_max_exp /= r_infty;
  OutPut(dim_y << endl);
  OutPut(N_y+1 << " " <<  N_z+1 << " " << N_r+1 << endl);
  if (N_r+1 >  WINDTUNNEL_DIM_R_CONST)
  {
      OutPut("dimension of WINDTUNNEL_BOUND_VAL_DROPS for internal coordinate too small !!!"<<endl);
      OutPut("increase WINDTUNNEL_DIM_R_CONST in Constants.h to " << N_r+1 << endl);
      exit(4711);
  }
  //accumulation of values in y-direction
  h = (double)(dim_y/(N_y+1.));
  OutPut(h<<endl);
  for(r=0;r<dim_r;r++)
  {
      for(j=0;j<dim_z;j++)
      {
	  comp=0;
	  for(i=0;i<N_y+1;i++)
	  {
	      acc[i][j][r][0]=0.;
	      acc[i][j][r][1]=0.;
	      firstcomp=comp;
	      //computation of the new partition, linear transform from dim_y intervals to N_y+1 intervals
	      while(comp< (i+1)*h)
	      {
		  acc[i][j][r][0]+=TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL_DROPS[comp][j][r][0];
		  acc[i][j][r][1]+=TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL_DROPS[comp][j][r][1];
		  comp++;
	      }
	      delta=comp-firstcomp;                       //delta:counts the summands in the sum
	      if (delta > 0)
	      {
		  acc[i][j][r][0]=acc[i][j][r][0]/delta;    // arithmetic mean
		  acc[i][j][r][1]=acc[i][j][r][1]/delta;
	      }
	  }
      }
  }

  //accumulation of values in z-direction
  h = (double)(dim_z/(N_z+1.));
  OutPut(h<<endl);
  for(r=0;r<dim_r;r++)
  {
      for(i=0;i<N_y+1;i++)
      {
	  comp=0;
	  for(j=0;j<N_z+1;j++)
	  {
	      acc1[i][j][r][0]=0;
	      acc1[i][j][r][1]=0;
	      firstcomp=comp;
	      while(comp<(double)(j+1)*h)
	      {
		  acc1[i][j][r][0]+=acc[i][comp][r][0];
		  acc1[i][j][r][1]+=acc[i][comp][r][1];
		  comp++;
	      }
	      delta=comp-firstcomp;
	      if (delta > 0)
	      {
		  acc1[i][j][r][0]=acc1[i][j][r][0]/delta;
		  acc1[i][j][r][1]=acc1[i][j][r][1]/delta;
		  //if (i==10)
		  //OutPut("acc1 " << i << " " << j << " " << r << " " << acc1[i][j][r][0] << endl);
	      }
	  }
      }
  }
  
  h = (double)(dim_r/(N_r+1.));
  OutPut(h<<endl);
  //accumulation of values in r-direction
  for(i=0;i<N_y+1;i++)
  {
      for(j=0;j<N_z+1;j++)
      {
	  comp=0;
	  // loop over layers of internal coordinate
	  for(r=0;r<N_r+1;r++)
	  {
	      acc2[i][j][r][0]=0;
	      acc2[i][j][r][1]=0;
	      firstcomp=comp;
	      while(comp<(double)(r+1)*h)
	      {
		  acc2[i][j][r][0]+=acc1[i][j][comp][0];
		  acc2[i][j][r][1]+=acc1[i][j][comp][1];
		  comp++;
	      }
	      delta=comp-firstcomp;
	      if (delta>0)
	      {
		  acc2[i][j][r][0] = acc2[i][j][r][0]/delta;
		  acc2[i][j][r][1] = acc2[i][j][r][1]/delta;
		  if (i==10)
		      OutPut(i << " "  << j << " " << r << " " << acc2[i][j][r][0] << endl);
	      }
	  }
      }
  }
  
  // reset array WINDTUNNEL_BOUND_VAL_DROPS
  for(i=0;i<WINDTUNNEL_DIM_Y_CONST;i++)
    for(j=0;j<WINDTUNNEL_DIM_Z_CONST;j++)
      for(r=0;r<WINDTUNNEL_DIM_R_CONST;r++)
      {
	  TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL_DROPS[i][j][r][0]=0;
	  TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL_DROPS[i][j][r][1]=0;
      }

  // fill a part of WINDTUNNEL_BOUND_VAL_DROPS with the accumulated data
  for(i=0;i<N_y+1;i++)
    for(j=0;j<N_z+1;j++)
      for(r=0;r<N_r+1;r++)
      {
	  //if ((i==11) && (j==N_z-1))
	  //  OutPut(r << " " << acc2[i][j][r][0] << endl);
        TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL_DROPS[i][j][r][0]=acc2[i][j][r][0];
        TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL_DROPS[i][j][r][1]=acc2[i][j][r][1];
      }

  /*for(i=0;i<N_y+1;i++)
  for(j=0;j<N_z+1;j++)
         for(r=N_r+1; r<dim_r; r++)
         {
         TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL_DROPS[i][j][r][0]=0.;
         TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL_DROPS[i][j][r][1]=0.;
          } */

  /* for(i=0;i<N_y+1;i++)
    {OutPut(endl);
  for(j=0;j<N_z+1;j++)
          for(r=0; r<N_r+1; r++)
           {if (i==2 &&j ==1)
        OutPut( " " <<TDatabase::ParamDB->WINDTUNNEL_BOUND_VAL_DROPS[i][j][r][0] );
             }}*/
}
