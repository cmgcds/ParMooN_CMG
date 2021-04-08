/** =======================================================================
* @class     TSystem1D
* @brief     stores the information of system 1D
* @author    Sashikumaar Ganesan
* @date      12.12.2020
* @History 
* ======================================================================= */
#include <System1D.h>
#include <Database.h>
#include <FEDatabase2D.h>
// #include <DirectSolver.h>
// #include <SquareStructure1D.h>
// #include <SquareMatrix1D.h>
// #include <FEFunction1D.h>
#include <LineAffin.h>
#include <LinAlg.h>

#include <MacroCell.h>
#include <JointEqN.h>
// #include <NodalFunctional1D.h>
// #include <SquareStructure1D.h>
// #include <SquareMatrix1D.h>
// #include <SquareMatrix.h>
// #include <Matrix.h>

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdio.h>
#include <stdlib.h>


TSystem1D::TSystem1D(int N_L, double start, double end, BoundCond1D *boundConLminLMax, DoubleFunctND *BdValues, char *ParamFile)
{
 BundValues = BdValues;
 N_Coord = 0;

 //mesh
 Domain_Intl = new TDomain(ParamFile);
 this->Generate1DMesh(start, end, N_L);
 Coll_Intl = Domain_Intl->GetCollection(It_Finest, 0);

 cout<< "N_Cells_Internal " << Coll_Intl->GetN_Cells() <<endl;

 // Finite difference in internal 
 if(TDatabase::ParamDB->INTL_DISCTYPE==FD)
  { 
   TDatabase::ParamDB->ANSATZ_ORDER_INTL = 1;
  }

 char IString[] = "I";
 FESpace1D = new TFESpace1D(Coll_Intl, IString, IString, TDatabase::ParamDB->ANSATZ_ORDER);

 N_Dof = FESpace1D->GetN_DegreesOfFreedom();

 Sol = new double[N_Dof];
 Rhs = new double[N_Dof];
 B = new double[N_Dof];
 defect = new double[N_Dof];

 if(TDatabase::ParamDB->ANSATZ_ORDER_INTL<0)
  {
   FESpace1D->SetAsDGSpace(); 
   dGDisc = TRUE;
   TDatabase::ParamDB->INTL_DISCTYPE=DG;  // no supg method
  }
  else
  {
    dGDisc = FALSE;
  }

  //read boundary conditions 
  boundConLminLMax(cond_Lmin, cond_Lmax); 
}

TSystem1D::~TSystem1D()
{
 delete [] Sol;
 delete [] Rhs;
 delete [] B;
 delete [] defect;
}


void TSystem1D::Generate1DMesh(double Start, double End, int N_Cells)
{
  int i, j, N_Vert;
  int *Lines;
  double len, h, x, y, *X;
  double hmin, hmax;
  TVertex **Vetrex;
  TJoint *Joint;
  TBaseCell  **CellTree;

  N_Vert = N_Cells+1;
  X = new double[N_Vert];

  h = (End-Start)/(double)N_Cells;

  X[0] = Start;

  for(i=1; i<N_Vert; i++)
   X[i] = X[i-1] + (double)h;

  X[N_Vert-1] = End;

    hmin = 1.e8; hmax = -1.e8; 
    for(i=0; i<N_Vert-1; i++)
     {
      len = sqrt ((X[i+1] - X[i])*(X[i+1] - X[i]));
      if(len< hmin) hmin = len;
      if(len> hmax) hmax = len;        
     }
     OutPut("L h_min : " << hmin << " L h_max : " << hmax << endl);
    
//    for(i=0; i<N_Vert; i++)
//     cout<< i << " X[i] " << X[i] <<endl;
 
//  exit(0);

  Lines = new int[2*N_Cells];
  Vetrex = new TVertex*[N_Vert]; 

  y=0.;

  for(i=0; i<N_Cells; i++)
   {
    Lines[2*i]=i;
    Lines[2*i+1]=i+1;
#ifdef __2D__
    Vetrex[i] = new TVertex(X[i], y);
#else
  cout << "Not Yet Implemented " <<endl;
#endif
   }
   
#ifdef __2D__
  Vetrex[N_Cells] = new TVertex(X[N_Vert-1], y);
#endif
  CellTree = new TBaseCell*[N_Cells];

   for (i=0;i<N_Cells;i++)
   {
 //     Vetrex[ i ]->GetCoords(x, y);
 //     cout<< " x " << x<< " y " << y<<endl;
    CellTree[i] = new TMacroCell(TDatabase::RefDescDB[S_Line], 0);
    CellTree[i]->SetVertex(0, Vetrex[ Lines[ 2*i       ]]);
    CellTree[i]->SetVertex(1, Vetrex[ Lines[ 2*i + 1]]);
    ((TMacroCell *) CellTree[i])->SetSubGridID(0);
   }
 //     Vetrex[ i ]->GetCoords(x, y);
 //     cout<< " x " << x<< " y " << y<<endl;
 //     exit(0);

   Domain_Intl->SetTreeInfo(CellTree, N_Cells);

   TDatabase::IteratorDB[It_EQ]->SetParam(Domain_Intl);
   TDatabase::IteratorDB[It_LE]->SetParam(Domain_Intl);
   TDatabase::IteratorDB[It_Finest]->SetParam(Domain_Intl);
   TDatabase::IteratorDB[It_Between]->SetParam(Domain_Intl);
   TDatabase::IteratorDB[It_OCAF]->SetParam(Domain_Intl);

   // start joint(vertex)
   Joint = new TJointEqN(CellTree[0]);
   CellTree[0]->SetJoint(0, Joint);


   for(i=1;i<N_Cells;i++)
    {
     Joint = new TJointEqN(CellTree[i-1], CellTree[i]);

     CellTree[i-1]->SetJoint(1, Joint);
     CellTree[i]->SetJoint(0, Joint);
   } // for(i=0;i<N_Cells;i++)

   // end joint(vertex)
   Joint = new TJointEqN(CellTree[N_Cells-1]);
   CellTree[N_Cells-1]->SetJoint(1, Joint);

  delete []  Lines;
}

