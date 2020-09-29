#include <FESpace2D.h>
#include <Database.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <LinAlg.h>
// #include <utility>
// #include <tuple>
#include <iostream>
// #include <ofstream>
#include <fstream>
#define PI 3.14159265

using namespace std;

  struct bodydes  // body description - forCGx and forCGy are x and y coordinates of CMs of all cells. Update this bodydes at each time step.
  {
   double *forCGx, *forCGy, *bodyprm, *cellareas;   // bodyprm = [CGx, CGy, totalarea]
  };
// computes total force in 2 directions and total moment given forCGx and forCGy. totFx, totFy and totmoment are passed by reference.
void compforce(double *fforCGx, double *fforCGy, double fCGx, double fCGy, double *fforceX, double *fforceY, int fN_Cells, double &totFx, double &totFy, double &totmoment)
{
  int i;
  double ftotFx = 0, ftotFy = 0, ftotmoment = 0;
  for(i=0;i<fN_Cells;i++)
  {
  fforceX[i] = sin(fforCGx[i]) + cos(fforCGy[i]);// + sin(fforCGx[i])/(cos(fforCGy[i])+2);
  fforceY[i] = cos(fforCGx[i]) + sin(fforCGy[i]);// + cos(fforCGx[i])/(sin(fforCGy[i])+2);
//     fforceX[i] = sin(fforCGx[i]) + cos(fforCGy[i]);
//     fforceY[i] = pow(fforCGx[i],0.5) + pow(fforCGy[i],0.5);
//     fforceX[i] = cos(fforCGx[i]);// + cos(fforCGy[i]);// + sin(fforCGx[i])/(cos(fforCGy[i])+2);
//     fforceY[i] = cos(fforCGy[i]);// + sin(fforCGy[i]);// + cos(fforCGx[i])/(sin(fforCGy[i])+2);
//     fforceX[i] = fforCGx[i];
//     fforceY[i] = fforCGy[i];    
  ftotFx += fforceX[i];  
  ftotFy += fforceY[i];
  ftotmoment += fforceX[i]*(fCGy - fforCGy[i]) + fforceY[i]*(fCGx - fforCGx[i]);
  }
  totFx = ftotFx; totFy = ftotFy; totmoment = ftotmoment;
//   cout<<"total force in X dir is "<<totFx << "\n";
//   cout<<"total force in Y dir is "<<totFy << "\n";
//   cout<<"total moment is "<<totmoment << "\n";
}
// (overloaded)computes total force in 2 directions and total moment
void compforce(bodydes *frigidb, double *fforceX, double *fforceY, int N_Cells, double &totFx, double &totFy, double &totmoment)
{
   int i;
   double fCGx = frigidb->bodyprm[0], fCGy = frigidb->bodyprm[1];
   double ftotFx = 0, ftotFy = 0, ftotmoment = 0;
   for(i=0;i<N_Cells;i++)
  {
  fforceX[i] = sin(frigidb->forCGx[i]) + cos(frigidb->forCGy[i]);// + sin(frigidb->forCGx[i])/(cos(frigidb->forCGy[i])+2);
  fforceY[i] = cos(frigidb->forCGx[i]) + sin(frigidb->forCGy[i]);// + cos(frigidb->forCGx[i])/(sin(frigidb->forCGy[i])+2);
//   fforceX[i] = sin(frigidb->forCGx[i]) + cos(frigidb->forCGy[i]);
//   fforceY[i] = pow(frigidb->forCGx[i],0.5) + pow(frigidb->forCGy[i],0.5);
  ftotFx += fforceX[i];
  ftotFy += fforceY[i];
  ftotmoment += fforceX[i]*(fCGy - frigidb->forCGy[i]) + fforceY[i]*(fCGx - frigidb->forCGx[i]);
  }
  totFx = ftotFx; totFy = ftotFy; totmoment = ftotmoment;
  cout<<"total force in X dir is "<<totFx << "\n";
  cout<<"total force in Y dir is "<<totFy << "\n";
  cout<<"total moment is "<<totmoment << "\n";
}
//std::tuple<totalarea, CGx, CGy, forCGx, forCGy>
//std::tuple<double, double, double, double*, double*> centerofgravity(TFESpace2D *FeSpace)
void centerofgravity(TFESpace2D *FeSpace, bodydes *rigidb)  // called this function first in main. passed empty bodydes to update with CMs etc
{
 int i, j;
 int N_cells, nvertices;
 double *x, *y, *cellareas, *forCGx, *forCGy, totalarea=0, CGx=0,CGy=0;
 TCollection *Coll;
 TGridCell *Cell;
 TVertex  **Vertices;
  
  Coll = FeSpace->GetCollection();
  N_cells = Coll->GetN_Cells();
  cellareas = new double[N_cells];
  forCGx = new double[N_cells];
  forCGy = new double[N_cells];
  memset(cellareas, 0, N_cells*SizeOfDouble);
  memset(forCGx, 0, N_cells*SizeOfDouble);
  memset(forCGy, 0, N_cells*SizeOfDouble);
  for(i=0;i<N_cells;i++)
  {
    Cell = (TGridCell*)Coll->GetCell(i);
    nvertices = Cell->GetN_Vertices();
    Vertices = Cell->GetVertices();
    x = new double[nvertices];
    y = new double[nvertices];
   
    for(j=0;j<nvertices;j++)
    {
      x[j] = Vertices[j]->GetX();
      y[j] = Vertices[j]->GetY();
    } 
    for(j=0;j<nvertices-1;j++)
    {
      cellareas[i]+=(x[j]*y[j+1] - x[j+1]*y[j]);
      forCGx[i]+=(x[j]+x[j+1])*(x[j]*y[j+1] - x[j+1]*y[j]);  //Centre of mass of a polygon given its x and y coordinates of vertices
      forCGy[i]+=(y[j]+y[j+1])*(x[j]*y[j+1] - x[j+1]*y[j]);
    }
    cellareas[i]+=(x[nvertices-1]*y[0] - x[0]*y[nvertices-1]);
    forCGx[i]+=(x[nvertices-1]+x[0])*(x[nvertices-1]*y[0] - x[0]*y[nvertices-1]);
    forCGy[i]+=(y[nvertices-1]+y[0])*(x[nvertices-1]*y[0] - x[0]*y[nvertices-1]);
    delete[] x,y;
    cellareas[i] = cellareas[i]/2;
    forCGx[i] = forCGx[i]/(6*cellareas[i]);
    forCGy[i] = forCGy[i]/(6*cellareas[i]);
    totalarea+=cellareas[i];
    CGx+=forCGx[i]*cellareas[i];
    CGy+=forCGy[i]*cellareas[i];
  }
  CGx = CGx/totalarea;
  CGy = CGy/totalarea;
  rigidb->bodyprm[0] = CGx;   // update bodydes (rigidb is object of type bodydes struct)
  rigidb->bodyprm[1] = CGy;
  rigidb->bodyprm[2] = totalarea;
  rigidb->forCGx = forCGx;
  rigidb->forCGy = forCGy;
  rigidb->cellareas = cellareas;
  
  delete[] cellareas, forCGx, forCGy;
  cout<<"totalarea is "<<totalarea<<"\n";
  cout<<"CGx is "<<CGx<<"\n";
  cout<<"CGy is "<<CGy<<"\n";  
//   return std::make_tuple(totalarea, CGx, CGy, forCGx, forCGy);
}// end of centerofgravity
//compute moment of inertia (/Sigma(m.r^2)) (accurate when mesh is fine. Not accurate for coarse mesh)
double compmominertia(bodydes *frigidb, int fN_Cells)
{
 int i;
 double fmominertia = 0;
 for(i=0;i<fN_Cells;i++)
 {
 fmominertia += (frigidb->cellareas[i])* (pow((frigidb->bodyprm[0] - frigidb->forCGx[i]),2) + pow((frigidb->bodyprm[1] - frigidb->forCGy[i]),2)); 
 }
 return fmominertia;
}
/** y := alpha*x */
void vecscal(int n, double alpha, double *x, double *y)
{
  register int i;
  register double scal;
  register double *a;

  scal = alpha;
  a = x;
  for(i=0; i<n; i++)
  {
    *y = *a*scal;
    a++;
    y++;
  }
}
/** x := alpha*x */
void vecscal(int n, double alpha, double *x)
{
  register int i;
  register double scal;
  register double *a;

  scal = alpha;
  a = x;
  for(i=0; i<n; i++)
  {
    *a *=scal;
    a++;
  }
}
/** x := alpha+x */
void vecadd(int n, double alpha, double *x)
{
  register int i;
  register double scal;
  register double *a;

  scal = alpha;
  a = x;
  for(i=0; i<n; i++)
  {
    *a += scal;
    a++;
  }
}
/** y := alpha+x */
void vecadd(int n, double alpha, double *x, double *y)
{
  register int i;
  register double scal;
  register double *a, *b;

  scal = alpha;
  a = x;
  b = y;
  for(i=0; i<n; i++)
  {
    *b = *a + scal;
    a++;
    b++;
  }
}
/** y := alpha*x + y */
void vecsum(int n, double alpha, double *x, double *y)
{
  register int i;
  register double *a, *b;
  register double scal;

  a = x;
  b = y;
  scal = alpha;
  for(i=0;i<n;i++)
  {
    *b += scal * *a;
    a++;
    b++;
  }
}
/** y := x + y */
void vecsum(int n, double *x, double *y)
{
  register int i;
  register double *a, *b;

  a = x;
  b = y;
  for(i=0;i<n;i++)
  {
    *b += *a;
    a++;
    b++;
  }
}
// z = x + y
void vecsum(int n, double *x, double *y, double *z)
{
  register int i;
  register double *a, *b, *c;

  a = x;
  b = y;
  c = z;
  for(i=0;i<n;i++)
  {
    *c = *a + *b;
    a++;
    b++;
    c++;
  }
}
//generate F vector given the input U vector
void generateFvec(int fN_Cells, double *ffinputkn, double *fFvec, bodydes *frigidb, double mass, double mominertia, double *fUn0)
{
  double CGx, CGy, *fforceX, *fforceY, totFx, totFy, totmoment, *inpforCGx, *inpforCGy;
  fforceX = new double[fN_Cells];
  fforceY = new double[fN_Cells];
  inpforCGx = new double[fN_Cells];
  inpforCGy = new double[fN_Cells];
  fFvec[0] = ffinputkn[3]; fFvec[1] = ffinputkn[4]; fFvec[2] = ffinputkn[5];
  vecadd(fN_Cells, ffinputkn[0] - fUn0[0], frigidb->forCGx, inpforCGx);  // shift x coordinates according to k1/2 or k2/2 or k3
  vecadd(fN_Cells, ffinputkn[1] - fUn0[1], frigidb->forCGy, inpforCGy);  // shift y coordinates according to k1/2 or k2/2 or k3
  CGx = ffinputkn[0]; CGy = ffinputkn[1];
  compforce(inpforCGx, inpforCGy, CGx, CGy, fforceX, fforceY, fN_Cells, totFx, totFy, totmoment);  // compute forces from these shifted points
  fFvec[3] = totFx/mass; fFvec[4] = totFy/mass; fFvec[5] = totmoment/mominertia;
}
//======================================================================
int main(int argc, char* argv[])
{ 
  int i, j, k, l, N_Cells, ret;
  double totalarea, CGx, CGy; // forCGx, forCGy;
  double *forceX, *forceY;
  char *PRM, *GEO;
  char ReadinDat[] = "readin.dat";
  char UString[] = "U";
  char Description[] = "description";
  
  TDomain *Domain = new TDomain();
  TDatabase *Database = new TDatabase();
  TCollection *Coll;
  TBaseCell *Cell;
  TFESpace2D *FeSpace;

  std::ostringstream os;
  os << " ";
// read parameter file
  if(argc>=2)
    ret=Domain->ReadParam(argv[1]);
  else
    ret=Domain->ReadParam(ReadinDat);
  if(ret==-1)
  {
    exit(-1);
  }
  OpenFiles();
  OutFile.setf(std::ios::scientific);
  Database->WriteParamDB(argv[0]);
//============== get the constants from ParamDB=============
  PRM = TDatabase::ParamDB->BNDFILE;
  GEO = TDatabase::ParamDB->GEOFILE;
// *****************************************************************************
    Domain->Init(PRM, GEO);
// write grid into an Postscript file
    os.seekp(std::ios::beg);
    os << "Domain_Coarse" << ".ps" << ends;
    Domain->PS(os.str().c_str(),It_Finest,0);
// refine grid
  for(i=0;i<TDatabase::ParamDB->UNIFORM_STEPS;i++)
	      Domain->RegRefineAll();
// write grid into an Postscript file
    os.seekp(std::ios::beg);
    os << "Domain" << ".ps" << ends;
    Domain->PS(os.str().c_str(),It_Finest,0);
    
  Coll = Domain->GetCollection(It_Finest, 0); 
  N_Cells = Coll->GetN_Cells();
    
    bodydes rigidba;  
    bodydes *rigidb;
    rigidb = &rigidba;
    rigidb->bodyprm = new double[3];  // malloc bodydes
    rigidb->forCGx = new double[N_Cells];
    rigidb->forCGy = new double[N_Cells];
    rigidb->cellareas = new double[N_Cells];

  FeSpace = new TFESpace2D(Coll, UString, Description);
  centerofgravity(FeSpace, rigidb);  // updates rigidb->bodyprm[0] = CGx;rigidb->bodyprm[1] = CGy;rigidb->bodyprm[2] = totalarea; rigidb->forCGx = forCGx; rigidb->forCGy = forCGy; rigidb->cellareas = cellareas;

  forceX = new double[N_Cells]; // forces in individual cells
  forceY = new double[N_Cells];
  CGx = rigidb->bodyprm[0];
  CGy = rigidb->bodyprm[1];
  double arealden = 100; // 100 kg/m^-2  
  double totFx=0.0, totFy=0.0, totmoment =0.0, mominertia = 0;
  mominertia = arealden * compmominertia(rigidb, N_Cells);
  cout << "mominertia is "<< mominertia<<"\n";
  double mass = (rigidb->bodyprm[2])*arealden;
  cout << "mass is "<< mass<<"\n";
  double inx = rigidb->bodyprm[0], iny = rigidb->bodyprm[1], inphi = 0, xdot=0, ydot=0, phidot=0; //initial x, y, phi, xdot, ydot, phidot
  double Fvec[6], k1[6], k2[6], k3[6], k4[6], Un0[6], Un1[6];
  Un0[0] = inx; Un0[1] = iny; Un0[2] = inphi; Un0[3] = xdot; Un0[4] = ydot; Un0[5] = phidot;
  double ksum[6], finputk2[6], Fveck2[6], finputk3[6], Fveck3[6], finputk4[6], Fveck4[6];
  memset(ksum, 0, 6*SizeOfDouble);
  memset(finputk2, 0, 6*SizeOfDouble);
  memset(finputk3, 0, 6*SizeOfDouble);
  memset(finputk4, 0, 6*SizeOfDouble);
  double *rotx, *roty;
  rotx = new double[N_Cells];
  roty = new double[N_Cells];
  double t = 0, deltat, timeend;
  deltat = TDatabase::TimeDB->TIMESTEPLENGTH;
  timeend = TDatabase::TimeDB->ENDTIME;
  
  fstream myfile;
  myfile.open("rigidout.txt");
  myfile << Un0[0]<< " "<< Un0[1]<<" "<<Un0[2]<<"\n";
  exit(0);
  
  //time stepping and RK4
  do
  {
//  compforce(rigidb, forceX, forceY, N_Cells, totFx, totFy, totmoment);
    compforce(rigidb->forCGx, rigidb->forCGy, rigidb->bodyprm[0], rigidb->bodyprm[1], forceX, forceY, N_Cells, totFx, totFy, totmoment);
    Fvec[0] = xdot; Fvec[1] = ydot; Fvec[2] = phidot; Fvec[3] = totFx/mass; Fvec[4] = totFy/mass; Fvec[5] = totmoment/mominertia;
    vecscal(6, deltat, Fvec);  // Fvec is now k1
    memcpy(ksum, Fvec, 6*sizeof(double)); // ksum has k1 now
    vecscal(6, 0.5, Fvec); //Fvec is now k1/2
    memcpy(finputk2, Fvec, 6*sizeof(double)); //finputk2 has k1/2 (Fvec) now
    vecsum(6, Un0, finputk2); //finputk2 has Un0 + k1/2    
    generateFvec(N_Cells, finputk2, Fveck2, rigidb, mass, mominertia, Un0);
    vecscal(6, deltat, Fveck2); // Fveck2 is now k2
    memcpy(finputk3, Fveck2, 6*sizeof(double)); //finputk3 has k2 (Fveck2) now
    vecscal(6, 2, Fveck2); //Fveck2 is now 2*k2
    vecsum(6, Fveck2, ksum); //ksum has k1 + 2*k2    
    vecscal(6, 0.5, finputk3); //finputk3 has k2/2 now
    vecsum(6, Un0, finputk3); //finputk3 is now Un0 + k2/2
    generateFvec(N_Cells, finputk3, Fveck3, rigidb, mass, mominertia, Un0);
    vecscal(6, deltat, Fveck3); // Fveck3 is now k3
    memcpy(finputk4, Fveck3, 6*sizeof(double)); //finputk4 has k3 (Fveck3) now
    vecscal(6, 2, Fveck3); //Fveck3 is now 2*k3
    vecsum(6, Fveck3, ksum); //ksum has k1 + 2*k2 + 2*k3 
    vecsum(6, Un0, finputk4); //finputk4 is now Un0 + k3
    generateFvec(N_Cells, finputk4, Fveck4, rigidb, mass, mominertia, Un0);
    vecscal(6, deltat, Fveck4); // Fveck4 is now k4
    vecsum(6, Fveck4, ksum); //ksum has k1 + 2*k2 + 2*k3 + k4
    vecscal(6, 1.0/6.0, ksum); // ksum is now (k1 + 2*k2 + 2*k3 + k4)/6
    vecsum(6, Un0, ksum); // ksum is now Un0 + (k1 + 2*k2 + 2*k3 + k4)/6
    //update rigidb (description of the rigid body) with the solution and copy Un1 to Un0
    memcpy(Un1, ksum, 6*sizeof(double));
    myfile << Un1[0]<< " "<< Un1[1]<<" "<<Un1[2]<<"\n";
    cout <<"x coordinate of centre of mass at t = "<< t<<" is "<<Un1[0]<<"\n";
    cout <<"y coordinate of centre of mass at t = "<< t<<" is "<<Un1[1]<<"\n";
    cout <<"Rotation of the rigid body about its centre of mass from initial position is "<<Un1[2]<<"\n";
    vecadd(N_Cells, Un1[0]-Un0[0], rigidb->forCGx); // translation of x coordinates of centres of mass of cells
    vecadd(N_Cells, Un1[1]-Un0[1], rigidb->forCGy); // translation of y coordinates of centres of mass of cells
    rigidb->bodyprm[0] += Un1[0]-Un0[0]; // translation of x coordinate of centre of mass of 2D domain
    rigidb->bodyprm[1] += Un1[1]-Un0[1]; // translation of y coordinate of centre of mass of 2D domain
    for(j=0;j<N_Cells;j++)  // rotation of x and y coordinates of centres of mass of cells
    {
      rotx[j] = rigidb->forCGx[j]*cos(Un1[2]-Un0[2]) + rigidb->forCGy[j]*sin(Un1[2]-Un0[2]);
      roty[j] = -1*rigidb->forCGx[j]*sin(Un1[2]-Un0[2]) + rigidb->forCGy[j]*cos(Un1[2]-Un0[2]);
    }  
    memcpy(rigidb->forCGx, rotx, N_Cells*sizeof(double));
    memcpy(rigidb->forCGy, roty, N_Cells*sizeof(double));
//     cout << "CGx is "<< rigidb->bodyprm[0] <<"\n";
//     cout << "CGy is "<< rigidb->bodyprm[1] <<"\n";
    totFx = 0; totFy = 0; totmoment = 0;
    xdot = Un1[3]; ydot = Un1[4]; phidot = Un1[5];
    memcpy(Un0, Un1, 6*sizeof(double));
    t+=deltat;
    } while(t<timeend);
    delete[] forceX, forceY, rigidb->bodyprm, rotx, roty, rigidb;  // delete malloced arrays after end of all timesteps
    cout << "N cells is "<< N_Cells <<"\n";
    myfile.close();
}
