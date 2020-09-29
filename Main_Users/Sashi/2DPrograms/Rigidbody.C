// =======================================================================
// 
// Purpose:     Main program for rigid body motion
//
// Author   :     Bhanu Teja, Sashikumaar Ganesan
//  date    :     29.05.2014
// modified :   
// ======================================================================= 
#include <Domain.h>
#include <FEDatabase2D.h>
#include <SquareStructure2D.h>
#include <Structure2D.h>
#include <AuxParam2D.h>
#include <DirectSolver.h>
#include <Assemble2D.h>
#include <Output2D.h>
#include <DiscreteForm2D.h>
#include <QuadAffin.h>
#include <QuadBilinear.h>
#include <QuadIsoparametric.h>
#include <TriaAffin.h>
#include <TriaIsoparametric.h>
#include <IsoBoundEdge.h>
#include <IsoInterfaceJoint.h>

#include <FESpace2D.h>
#include <Database.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <LinAlg.h>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
// #define PI 3.14159265

// ============================================================
//  include example file
#include "../Examples/TNSE_2D/RotatingBody.h"
// ============================================================

double GetMOI(TFESpace2D *FESpace, double *CGxy, double &CGx_Body, double &CGy_Body, double &mass)
{
  int i, j, N_Cells, N_Edges, polydegree; 
  
  double MoI = 0.;
  double locvol, r2; 
  
  TCollection *Coll;
  TBaseCell *cell;
  FE2D FEId;
  RefTrans2D RefTrans;
  TRefTrans2D *rt;
  QuadFormula2D QuadFormula;
  TQuadFormula2D *qf2;
  boolean IsIsoparametric;
  TJoint *joint;
  JointType jointtype;
  BoundTypes bdtype;

  Coll = FESpace->GetCollection();
  N_Cells = Coll->GetN_Cells();
  
  CGx_Body = 0.;
  CGy_Body = 0.; 
  
  for(i=0;i<N_Cells;i++)
   {
    CGx_Body +=CGxy[i];
    CGy_Body +=CGxy[N_Cells+i];     
   }
   CGx_Body /=(double)N_Cells;
   CGy_Body /=(double)N_Cells;    
  
  mass=0.;
  for(i=0;i<N_Cells;i++)
   {
    cell = Coll->GetCell(i);
    FEId = FESpace->GetFE2D(i, cell);

    RefTrans = TFEDatabase2D::GetRefTrans2D_IDFromFE2D(FEId);
    N_Edges = cell->GetN_Edges(); 

    IsIsoparametric = FALSE;
    if (TDatabase::ParamDB->USE_ISOPARAMETRIC)
    {
      for(j=0;j<N_Edges;j++)
      {
        joint = cell->GetJoint(j);
        jointtype = joint->GetType();
        if(jointtype == BoundaryEdge)
        {
          bdtype = ((TBoundEdge *)(joint))->GetBoundComp()->GetType();
          if(bdtype != Line)
            IsIsoparametric = TRUE;
        }
        if(jointtype == InterfaceJoint)
        {
          bdtype = ((TInterfaceJoint *)(joint))->GetBoundComp()->GetType();
          if(bdtype != Line)
            IsIsoparametric = TRUE;
        }
        if(jointtype == IsoInterfaceJoint ||
           jointtype == IsoBoundEdge)
          IsIsoparametric = TRUE;
      }
    } // endif 

    if(IsIsoparametric)
    {
      switch(N_Edges)
      {
        case 4:
          RefTrans = QuadIsoparametric;
        break;

        case 3:
          RefTrans = TriaIsoparametric;
        break;
      }
    } // endif IsIsoparametric

    rt = TFEDatabase2D::GetRefTrans2D(RefTrans);
    switch(RefTrans)
    {
      case TriaAffin:
        ((TTriaAffin *)rt)->SetCell(cell);
        locvol = ((TTriaAffin *)rt)->GetVolume();
      break;

      case TriaIsoparametric:
        polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
        QuadFormula = TFEDatabase2D::GetQFTriaFromDegree(2*polydegree-1);
        ((TTriaIsoparametric *)rt)->SetApproximationOrder(polydegree);
        ((TTriaIsoparametric *)rt)->SetQuadFormula(QuadFormula);
        ((TTriaIsoparametric *)rt)->SetCell(cell);
        locvol = ((TTriaIsoparametric *)rt)->GetVolume();
      break;

      case QuadAffin:
        ((TQuadAffin *)rt)->SetCell(cell);
        locvol = ((TQuadAffin *)rt)->GetVolume();
      break;

      case QuadBilinear:
        ((TQuadBilinear *)rt)->SetCell(cell);
        locvol = ((TQuadBilinear *)rt)->GetVolume();
      break;

      case QuadIsoparametric:
        polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
        QuadFormula = TFEDatabase2D::GetQFQuadFromDegree(3*polydegree);
        ((TQuadIsoparametric *)rt)->SetApproximationOrder(polydegree);
        ((TQuadIsoparametric *)rt)->SetQuadFormula(QuadFormula);
        ((TQuadIsoparametric *)rt)->SetCell(cell);
        locvol = ((TQuadIsoparametric *)rt)->GetVolume();
      break;
    }

   r2 =  (CGx_Body- CGxy[i])*(CGx_Body- CGxy[i]) + (CGy_Body- CGxy[N_Cells+i])*(CGy_Body- CGxy[N_Cells+i]);
   MoI += locvol*r2;
   mass +=locvol;
  } // endfor i

  //multiply with the density
  MoI *=TDatabase::ParamDB->P0;
  mass*=TDatabase::ParamDB->P0;
  
  return  MoI;
} // GetMOI

int main(int argc, char* argv[])
{ 
  TDomain *Domain = new TDomain();
  TDatabase *Database = new TDatabase();
  TFEDatabase2D *FEDatabase = new TFEDatabase2D();  
  TCollection *Coll;
  TBaseCell *Cell;
  TFESpace2D *Grid_Space, *Cell_Space;
  TFEVectFunct2D *GridPos, *GridCGxy;
  TOutput2D *Output;
  
  int i, j, k, l, N_Cells, ret, img=0;
  int N_GridDOFs, N_C, m=0;
  
  double *gridpos, *gridpos_old, MoI, mass;
  double *CGxy, mominertia, CGx_Body, CGy_Body;
  double t = 0, tau, timeend;
  
  char *PRM, *GEO;
  char ReadinDat[] = "readin.dat";
  char GridString[] = "X";
  char Description[] = "description";
  const char vtkdir[] = "VTK";  
  char *VtkBaseName;
    
  std::ostringstream os;
  os << " ";
  
//======================================================================
// read parameter file
//======================================================================
   if(argc>=2)
    { ret=Domain->ReadParam(argv[1]); }
    else
      { ret=Domain->ReadParam(ReadinDat); }  
  if(ret==-1)
   {
    exit(-1);
   }


  OpenFiles();
  OutFile.setf(std::ios::scientific);
  Database->WriteParamDB(argv[0]);
  Database->WriteTimeDB();
  ExampleFile();
    
  
  VtkBaseName = TDatabase::ParamDB->VTKBASENAME;
  mkdir(vtkdir, 0777);  
//======================================================================
// read boundary parameterization and initialize coarse grid
//======================================================================
   PRM = TDatabase::ParamDB->BNDFILE;
   GEO = TDatabase::ParamDB->GEOFILE;
   
   Domain->Init(PRM, GEO);
   
// write grid into an Postscript file
//    os.seekp(std::ios::beg);
//    os << "Domain_Coarse" << ".ps" << ends;
//    Domain->PS(os.str().c_str(),It_Finest,0);
    
// refine grid
   for(i=0;i<TDatabase::ParamDB->UNIFORM_STEPS;i++)
     Domain->RegRefineAll();
    
 
   // write grid into an Postscript file
//    os.seekp(std::ios::beg);
//    os << "Domain" << ".ps" << ends;
//    Domain->PS(os.str().c_str(),It_Finest,0);
    
   Coll = Domain->GetCollection(It_Finest, 0); 
   N_Cells = Coll->GetN_Cells();
  
   //=========================================================================
   // construct all finite element spaces
   //========================================================================= 
   Grid_Space  =  new TFESpace2D(Coll, GridString, Description, GridBoundCondition, 1, NULL);
   
   N_GridDOFs = Grid_Space->GetN_DegreesOfFreedom();
   OutPut("Dof       : "<< setw(10) << N_GridDOFs  << endl); 
   
   // piecewise constant space
   Cell_Space  =  new TFESpace2D(Coll, GridString, Description, GridBoundCondition, 0, NULL);
   N_C = Cell_Space->GetN_DegreesOfFreedom();
   
   if(N_C!=N_Cells)
    {
     OutPut(" Cell_Space must be Q_0 or P_0 !!!" <<endl);
     exit(-1);
    }
   //=========================================================================
   // memory allocate all vectors and construction of all fefunction
   //=========================================================================      
   gridpos = new double[2*N_GridDOFs];
   gridpos_old = new double[2*N_GridDOFs];   

   memset(gridpos, 0, 2*N_GridDOFs*SizeOfDouble);
   GridPos = new TFEVectFunct2D(Grid_Space, GridString, GridString, gridpos, N_GridDOFs, 2);  
   GridPos->GridToData();      
   memcpy(gridpos_old, gridpos, 2*N_GridDOFs*SizeOfDouble); 
   
   // CGx and CGy of each cell in the array CGxy
   CGxy = new double[2*N_Cells];
   memset(CGxy, 0, 2*N_Cells*SizeOfDouble);
   GridCGxy = new TFEVectFunct2D(Cell_Space, GridString, GridString, CGxy, N_Cells, 2);
   //compute CGx and CGy each cell 
   GridCGxy->GridToData();        
   

   
   //=========================================================================
   //create output 
    Output = new TOutput2D(1, 1, 1, 1, Domain);   
          
    Output->AddFEFunction(GridPos); 
    os.seekp(std::ios::beg);
    Output->AddParameter(TDatabase::TimeDB->CURRENTTIME,os.str().c_str());      
  //=========================================================================

//   data for time stepping
    MoI = GetMOI(Grid_Space, CGxy, CGx_Body, CGy_Body, mass);
    OutPut("CGx_Body   "<<CGx_Body<< " CGy_Body "<<CGy_Body<<" MoI   "<<MoI<<" mass   "<<mass<<endl);    
    
    
    tau = TDatabase::TimeDB->TIMESTEPLENGTH;
    timeend = TDatabase::TimeDB->ENDTIME;

    
   if(TDatabase::ParamDB->WRITE_VTK)
    {
        os.seekp(std::ios::beg);
        if(img<10) os << "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
         else if(img<100) os << "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
         else if(img<1000) os <<"VTK/"<< VtkBaseName<<".00"<<img<<".vtk" << ends;
         else if(img<10000) os << "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
         else  os << "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
        Output->WriteVtk(os.str().c_str());
     img++;
    }    
    
 //=========================================================================   
  //time stepping and RK4
  do
  {
   m++; 
   TDatabase::TimeDB->CURRENTTIME += tau;
    

   
 // ===================================================================  
 // x[i] = gridpos[i] and y[i] = gridpos[N_GridDOFs+i]   
 //  add arbitrary displacement
    for(i=0;i<N_GridDOFs;i++)
      ModifyCoords(gridpos_old[N_GridDOFs+i], gridpos[N_GridDOFs+i]);
// ===================================================================    
    
     // update the mesh position
     GridPos->DataToGrid();      
     
     
  if(m==1 || m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0) 
   if(TDatabase::ParamDB->WRITE_VTK)
    {
        os.seekp(std::ios::beg);
        if(img<10) os << "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
         else if(img<100) os << "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
         else if(img<1000) os <<"VTK/"<< VtkBaseName<<".00"<<img<<".vtk" << ends;
         else if(img<10000) os << "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
         else  os << "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
        Output->WriteVtk(os.str().c_str());
     img++;
    }    
    
     t+=tau;
  } while(t<timeend);   
    
  
  
  
  if(TDatabase::ParamDB->WRITE_VTK)
   {
        os.seekp(std::ios::beg);
        if(img<10) os << "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
         else if(img<100) os << "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
         else if(img<1000) os <<"VTK/"<< VtkBaseName<<".00"<<img<<".vtk" << ends;
         else if(img<10000) os << "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
         else  os << "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
        Output->WriteVtk(os.str().c_str());
     img++;
    }      
      
  CloseFiles();
  return 0;  
  
} 

