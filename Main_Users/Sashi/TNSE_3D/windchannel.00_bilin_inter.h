// time-dependent NSE
// simulation of wind channel experiments
// u(x,y) = unknown
// p(x,y) = unknown

#define __WINDCHANNEL__

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
using namespace std;
#include <Constants.h>

void ExampleFile()
{
    OutPut("Example: windchannel.00_bilin_inter.h " << TDatabase::ParamDB->P7/100.0
         << " % noise (only for U1 !!!)" << endl);
    TDatabase::ParamDB->DRIFT_Z = 0.18;
    OutPut("Setting DRIFT_Z to " << TDatabase::ParamDB->DRIFT_Z << endl);
    TDatabase::ParamDB->N_CELL_LAYERS = 2;
    OutPut("Setting N_CELL_LAYERS to " << TDatabase::ParamDB->N_CELL_LAYERS  << endl);    
}


void InitialU1(double x, double y, double z, double *values)
{
  double w, t1, t2, n,eps=1e-3;

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


// value of boundary condition
void U1BoundValue(double x, double y, double z, double &value)
{   
    double eps = 1e-6;
    int coord_y, coord_z;
    double dec_plac_y, dec_plac_z;
    double noise = TDatabase::ParamDB->P7/100.0;

    value = 0;
    //values y,z are in a range of 0-0.45 or 0-0.18
    //for this reason they must be converted in values between 0-45 or 0-18 
    coord_y= (int) (y*100 + 0.0001);  //the decimal places are cut off 
    coord_z= (int) (z*100 + 0.0001);

    dec_plac_y = fabs((double) coord_y - (y * 100.)); //decimal places
    dec_plac_z = fabs((double) coord_z - (z * 100.)); //are used for interpolation

    // on the inflow boundary
    if(fabs(x)< eps)
    {
	// not on the edges
	if(coord_y != 45 && coord_z != 18)
	{ 
	    value = (1-dec_plac_y) * (1-dec_plac_z) * TDatabase::ParamDB->WINDCHANNEL_BOUND_VAL[coord_y][coord_z]
		+  dec_plac_y    * (1-dec_plac_z) *  TDatabase::ParamDB->WINDCHANNEL_BOUND_VAL[coord_y+1][coord_z]
		+ (1-dec_plac_y) *  dec_plac_z    *  TDatabase::ParamDB->WINDCHANNEL_BOUND_VAL[coord_y][coord_z+1]
		+  dec_plac_y    *  dec_plac_z    *  TDatabase::ParamDB->WINDCHANNEL_BOUND_VAL[coord_y+1][coord_z+1]; 
	    //OutPut(y<< " : " <<coord_y<< " : "<<dec_plac_y<<" : "<< z<< " : " << coord_z << " :: " << value << endl);
	}
	else
	{
	    if (coord_y == 45 && coord_z != 18)
		value = (1-dec_plac_z) *  TDatabase::ParamDB->WINDCHANNEL_BOUND_VAL[coord_y][coord_z]
		    +  dec_plac_z    *  TDatabase::ParamDB->WINDCHANNEL_BOUND_VAL[coord_y][coord_z+1];
	    if(coord_y != 45 && coord_z == 18)
		value = (1-dec_plac_y) *  TDatabase::ParamDB->WINDCHANNEL_BOUND_VAL[coord_y][coord_z]
		    +  dec_plac_y    *  TDatabase::ParamDB->WINDCHANNEL_BOUND_VAL[coord_y+1][coord_z];
	    if(coord_y == 45 && coord_z == 18)   
		value = TDatabase::ParamDB->WINDCHANNEL_BOUND_VAL[45][18];

	    //OutPut(y<< " : " <<coord_y<< " : "<<dec_plac_y<<" : "<< z<< " : " << coord_z << " ::: " << value << endl);
	}
    }
        
    // MAN BEKOMMT HIER (x,y,z)
    // FALLS fabs(x-0.5) < eps MUSS ANHAND DER (y,z) KOORDINATEN 
    // BESTIMMT WREDEN, WELCHEN WERT MAN AUS DEM ARRAY WINDCHANNEL_BOUND_VAL
    // NIMMT
    // 1. VARIANTE: MAN NIMMT DEN WERT, BEI DEM DIE KOORDINATEN AM NAECHSTEN
    //              AN (y,z) LIEGEN
    // 2. VARIANTE: MAN NIMMT EINE GEEIGNETE MITTELUNG AUS DEN UMLIEGENDEN WERTEN

    // nur zum Test
    //value = TDatabase::ParamDB->WINDCHANNEL_BOUND_VAL[7][5];
   // if (value < 0.5 )
    // {
    //OutPut(value << " " << x <<" " << z<< " " << endl);
     //}

    // add the noise
    //OutPut(value << " ")
    value *= (1 + noise * ((double)rand()/RAND_MAX-0.5));
    //OutPut(value << " " << endl)
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
  static double eps = 1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, x, y, z;

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
 string name_data_file = "../Examples/TNSE_3D/Inlet1.2bar.txt";

  ifstream data_file((name_data_file).c_str());
   if (!data_file) 
  {                           
    cout << " cannot open data-file: " << name_data_file << endl;
    exit(-1);
  }
  char read_line[200];
  string line;
  int  NumberOFLines = 0;

  //loop counts the number of lines in the txt-file
  while( data_file.getline(read_line,200))
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

  OutPut("NOL: "  <<  NumberOFLines << endl);
  int counter, i=0;
  double A[NumberOFLines][12];
  string entry;
  size_t found;
  int n=0;
//read datas from the file in the array A
  while( data_file_1.getline(read_line,200))
  { 
    i=0;
    line=read_line;
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
  data_file_1.close();

/*  for(i=0;i<NumberOFLines;i++)
  {
    for(int j=0;j<12;j++)
    {
     cout << A[i][j] <<" ";
    } 
    cout << endl;
  }*/

  
  double max_y =-1000.;
  double max_z =-1000.;
  double min_y = 1000.;
  double min_z = 1000.; 

  //stepsize
  int ny=10;
  int nz=10;
//convert the array into a square matrix 
  for(i=0;i<NumberOFLines;i++)
  {
    if(A[i][1]>max_y) max_y=A[i][1];
    if(A[i][2]>max_z) max_z=A[i][2];
    if(A[i][1]<min_y) min_y=A[i][1];
    if(A[i][2]<min_z) min_z=A[i][2];
  }

  int dim_y = (int)((max_y - min_y)/ny +1 +1e-6);
  int dim_z = (int)((max_z - min_z)/nz +1 +1e-6);
  if (( TDatabase::ParamDB->WINDCHANNEL_DIM_Y!=dim_y)||
      ( TDatabase::ParamDB->WINDCHANNEL_DIM_Z!=dim_z))
  {
      OutPut("error in dimension of experimental boundary condition "
	     << dim_y << " " << dim_z << endl);
      exit(4711);
  }
 
  for(n=0;n<NumberOFLines;n++)
  { 
      // AN DIESER STELLE MUESSEN NOCH DIE ARRAY WINDCHANNEL_X UND
      // WINDCHANNEL_Y AUS DER DATABASE BELEGT WERDEN
  //x-coordinate of the sorted array 
   int coordy=(int) (((A[n][1]-min_y)/ny)+ 1e-6);
   //z-coordinate of the sorted array
   int coordz=(int) (((A[n][2]-min_z)/nz) + 1e-6);
    //write the datas in the array
    TDatabase::ParamDB->WINDCHANNEL_BOUND_VAL[coordy][coordz] = A[n][7]; 
  }
 for(i=0;i<dim_y;i++)
  {for( int j=0;j<dim_z;j++)
   {
   cout << TDatabase::ParamDB->WINDCHANNEL_BOUND_VAL[i][j] <<" ";
 } cout << endl;
}
}


