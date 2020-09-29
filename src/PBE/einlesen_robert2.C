#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
using namespace std;

int main () 
{
  //cout.precision(12);

  string name_data_file = "/home/standard/Magdeburg/Inlet1.2bar.txt";

  ifstream data_file((name_data_file).c_str());

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

  cout << "NOL: "  <<  NumberOFLines << endl;
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

  
  double max_x =-1000.;
  double max_y =-1000.;
  double min_x = 1000.;
  double min_y = 1000.; 

  //stepsize
  int nx=10;
  int ny=10;
//convert the array into a square matrix 
  for(i=0;i<NumberOFLines;i++)
  {
    if(A[i][1]>max_x) max_x=A[i][1];
    if(A[i][2]>max_y) max_y=A[i][2];
    if(A[i][1]<min_x) min_x=A[i][1];
    if(A[i][2]<min_y) min_y=A[i][2];
  }

  int dim_x = (int)((max_x - min_x)/nx +1 +1e-6);
  int dim_y = (int)((max_y - min_y)/ny +1 +1e-6);
  double C[dim_x][dim_y];

  for(n=0;n<NumberOFLines;n++)
  { //x-coordinate of the sorted array 
   int cordx=(int) (((A[n][1]-min_x)/nx)+ 1e-6);
   //x-coordinate of the sorted array
   int cordy=(int) (((A[n][2]-min_y)/ny) + 1e-6);
    //write the datas in the array
    C[cordx][cordy] = A[n][7]; 
  }
 for(i=0;i<dim_x;i++)
  {for( int j=0;j<dim_y;j++)
   {
   cout << C[i][j] <<" ";
 } cout << endl;
}
  return 0;
}


