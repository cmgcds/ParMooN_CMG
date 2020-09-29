// =======================================================================
// Purpose  :   Main program for rigid body motion
// Author   :   Bhanu Teja, Sashikumaar Ganesan
// Date     :   29.05.2014
// Modified :   
// ======================================================================= 
#include <Domain.h>
#include <FEDatabase3D.h>
#include <SquareStructure3D.h>
#include <Structure3D.h>
#include <AuxParam3D.h>
#include <DirectSolver.h>
#include <Assemble3D.h>
#include <Output3D.h>
#include <DiscreteForm3D.h>
#include <QuadAffin.h>
#include <QuadBilinear.h>
#include <QuadIsoparametric.h>
#include <TetraAffin.h>
#include <TetraIsoparametric.h>
#include <HexaAffin.h>
#include <HexaIsoparametric.h>
#include <HexaTrilinear.h>
#include <TriaAffin.h>
#include <TriaIsoparametric.h>
#include <IsoBoundEdge.h>
#include <IsoInterfaceJoint.h>
#include <FESpace3D.h>
#include <Database.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <LinAlg.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <BoundFace.h>
#include <InterfaceJoint3D.h>
#include <MainUtilities.h>
// #include <TimeUtilities.h>
#include <tetgen.h>
#include <GridCell.h>
#include <MacroCell.h>
#include <BdPlane.h>
#include <BdSphere.h>
#include <IsoBoundFace.h>
#include <InterfaceJoint3D.h>
#include <IsoInterfaceJoint3D.h>
#include "../../Examples/RB_3D/import3Dmesh4.h"
// ============================================================

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

int main(int argc, char* argv[])
{ 
  TDomain *CompleteDomain = new TDomain();
  TDatabase *Database = new TDatabase();
  TFEDatabase3D *TotalFEDatabase = new TFEDatabase3D();
  TCollection *TotalColl;
  TBaseCell *Cell;
  TFESpace3D *TotalSpace;
  TFEVectFunct3D *GridVertices;
  TOutput3D *TotalOutput;
  
//   BoundCondFunct3D *BoundaryConditions[1];
  
  int TotalCells=0, ret, img=0;
  int TotalDOFs=0,i;  
  double *gridvertices;
  char ReadinDat[] = "readin.dat";
  char GridString[] = "X";
  char SolidString[] = "S";
  char FluidString[] = "F";
  char Description[] = "description";
  const char completevtkdir[] = "completeVTK";
  char *VtkBaseName;  
  std::ostringstream os;
  os << " ";
// read parameter file
   if(argc>=2)
    { ret=CompleteDomain->ReadParam(argv[1]); }
    else
      { ret=CompleteDomain->ReadParam(ReadinDat); }  
  if(ret==-1)
   {
    exit(-1);
   }
  OpenFiles();
  OutFile.setf(std::ios::scientific);
  Database->WriteParamDB(argv[0]);
  Database->WriteTimeDB();
  ExampleFile(); 
   TetrameshCreate(CompleteDomain);   
//      for(i=0;i<1;i++)
//     CompleteDomain->RegRefineAll();    
   TotalColl = CompleteDomain->GetCollection(It_Finest, 0);
   TotalCells = TotalColl->GetN_Cells();
   cout<<"Number of FE Cells "<<TotalCells<<endl;
   TotalSpace  =  new TFESpace3D(TotalColl, GridString, Description, GridBoundCondition, 1);
    os.seekp(std::ios::beg);
    VtkBaseName = TDatabase::ParamDB->VTKBASENAME;
  //=====================================================================
   TotalDOFs = TotalSpace->GetN_DegreesOfFreedom();
   cout << "total dof "<<TotalDOFs<<endl;
   gridvertices = new double[3*TotalDOFs];
   GridVertices = new TFEVectFunct3D(TotalSpace, GridString, GridString, gridvertices, TotalDOFs, 3);
   GridVertices->GridToData();
    TotalOutput = new TOutput3D(1, 1, 1, 1, CompleteDomain, TotalColl, completevtkdir);
    TotalOutput->AddFEVectFunct(GridVertices);    
    TotalOutput->AddParameter(TDatabase::TimeDB->CURRENTTIME,os.str().c_str());
    VtkBaseName = TDatabase::ParamDB->VTKBASENAME;
    mkdir(completevtkdir, 0777);
  //===========================  
   if(TDatabase::ParamDB->WRITE_VTK)     
    {
	os << "completeVTK/"<<VtkBaseName<<"."<<std::fixed<<std::setprecision(4)<<img/10000.<<".vtk" << ends;
	TotalOutput->WriteVtk(os.str().c_str());
        os.str("");   
     img++;
    }
}
