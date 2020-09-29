// =======================================================================
// 
// Purpose:     Main program for impinging droplet
//
// Author:      Sashikumaar Ganesan
// modified    10.06.2010 
//             29.01.2014 (surfactant)
// ======================================================================= 
#include <Domain.h>
#include <Database.h>
#include <DiscreteForm2D.h>
#include <FEDatabase2D.h>
#include <FESpace2D.h>
#include <SquareStructure1D.h>
#include <SquareStructure2D.h>
#include <SquareMatrix1D.h>
#include <FEFunction1D.h>
#include <SquareStructure2D.h>
#include <Structure2D.h>
#include <AuxParam2D.h>
#include <DirectSolver.h>
#include <Assemble2D.h>
#include <Output2D.h>
#include <DiscreteForm2D.h>
#include <LinAlg.h>
#include <TNSE2D_ParamRout.h>
#include <AuxParam2D.h>
 
#include <Collection.h>
#include <NodalFunctional2D.h>
#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
// #include <malloc.h>

#include <Upwind.h>
#include <NSE_MultiGrid.h>
#include <NSE_MGLevel1.h>
#include <NSE_MGLevel2.h>
#include <NSE_MGLevel3.h>
#include <NSE_MGLevel4.h>
#include <Convolution.h>
#include <ItMethod.h>
#include <FixedPointIte.h>
#include <MultiGridIte.h>
#include <FgmresIte.h>

#include <MultiGrid2D.h>
#include <MGLevel2D.h>
#include <FreeSurface2D.h>

#include <MainUtilities.h>
#include <TimeDiscRout.h>


#include <MacroCell.h>
#include <BoundEdge.h>
#include <IsoBoundEdge.h>
#include <IsoInterfaceJoint.h>
#include <BdLine.h>
#include <BdCircle.h>
#include <GridCell.h>

#include <QuadAffin.h>
#include <QuadBilinear.h>
#include <QuadIsoparametric.h>
#include <TriaAffin.h>
#include <TriaIsoparametric.h>

#  include <sys/stat.h>
#  include <sys/types.h>

// #include "../Examples/TNSE_2D/Droponsolid.h"
#include "../TNSE_2D/Drop_Imping_Axial3D.h"
// #include "../Examples/TNSE_2D/DropHeat_imping_axial3D.h"
// #include "../TNSE_2D/Drop_Imping_Axial3D_DiffusionTest.h"
// #include "../TNSE_2D/Drop_Imping_Axial3D_ConvTest.h"

extern "C"
{
  #include <gridgen.h>    
  void triangulate(char*, struct triangulateio*,
                   struct triangulateio*, struct triangulateio*);
}
// sorting curved surfrace vertices - general case
void Sort_Imping(TBaseCell **cell, TVertex **Vertex, int *CellNo, int *EdgeNo, int N,
                     double X0, double Y0 )
{
 int i, j, k, temp, test=0;
 
 double x, y, x1, y1;
 
 TVertex *temp_vert;
 TBaseCell *temp_cell;

 // finding the right(starting) vertex
  for(i=0;i<N;i++)
   {
    Vertex[i]->GetCoords(x, y);
    if(  sqrt((x-X0)*(x-X0) +(y-Y0)* (y-Y0))<1e-5  )
      {
//        cout << " sorting " << x << ' ' << y<<endl;
       temp_vert = Vertex[0];
       Vertex[0] = Vertex[i];
       Vertex[i] = temp_vert;

       temp_cell = cell[0];
       cell[0] = cell[i];
       cell[i] = temp_cell;
       
       temp = EdgeNo[0];
       EdgeNo[0] = EdgeNo[i];
       EdgeNo[i] = temp;
       
       temp = CellNo[0];
       CellNo[0] = CellNo[i];
       CellNo[i] = temp;        

       test++;
      }
    if(test) break;
   }

  if(i==N)
  {
   cout<< "Error in finding start vert " << X0 << endl;   
   exit(0);
  }

  for(i=0;i<N-1;i++)
   {
    test = 0; 
    k = cell[i]->GetN_Edges();
    cell[i]->GetVertex((EdgeNo[i]+1) % k)->GetCoords(x, y);

     for(j=i+1;j<N;j++)
      {
       Vertex[j]->GetCoords(x1, y1);
       if((x==x1) && (y==y1))
        {
         temp_vert = Vertex[j];
         Vertex[j] = Vertex[i+1];
         Vertex[i+1] = temp_vert;
       
         temp_cell = cell[j];
         cell[j] = cell[i+1];
         cell[i+1] = temp_cell;
       
         temp = EdgeNo[j];
         EdgeNo[j] = EdgeNo[i+1];
         EdgeNo[i+1] = temp; 
         
         temp = CellNo[j];
         CellNo[j] = CellNo[i+1];
         CellNo[i+1] = temp;          
         
         test++;
        }
      if(test) break;  
     }
   }

//    print   
//   for(i=0;i<N;i++)
//    {
//     Vertex[i]->GetCoords(x, y);
//     cout<<i<<  " x : "<< x << " y : "<< y<<  " Angle of free Vertices "<<(180/Pi)*atan2(y,x)<<endl;
//    }
//   exit(0); 
}


void FreeSurf_axial3D_new(TSquareMatrix2D *A11, TSquareMatrix2D *A22,
                          double *rhs1, double *rhs2,
                          BoundCondFunct2D *BoundaryCondition,
                          double dt, double *Ucl, TFEFunction2D *Surfact, double *param)
{
  int i, j, k, l, DOF_R, SDOF_R, DOF_L, m, mm;
  int *KCol, *RowPtr, *JointDOF, N_DOF;
  int N_LinePoints;
  int *BeginIndex, *GlobalNumbers, *DOF, TestDOF, AnsatzDOF;
  int index1, index2;
  int count=0, count1=0, count2=0;
  int N_BaseFunct, *N_BaseFuncts;
  int JointNumbers[MAXN_JOINTS], IJoint, N_IsoJoints;  
  int comp, N_U, test_L=1, test_R=1;
  int N_Cells, N_Vertices, N_Edges, Semi_implicit=0;  
    
  double r2, r;  
  double  X_B[100], Y_B[100], r_axial, d1, d2, e1, e2, ngrad_test, ngrad_ansatz, tmp;
  double val, theta, factor1, factor2, angle;
  double *ValuesA11, *ValuesA12, *ValuesA21, *ValuesA22;  
  double *LineWeights, *zeta;
  double x0, y0, x1, y1,tx,ty,mod_t, x, y;
  double **uref, **uxiref, **uetaref;
  double uorig[MaxN_BaseFunctions2D], uxorig[MaxN_BaseFunctions2D];
  double uyorig[MaxN_BaseFunctions2D];
  double t0, t1, n0, n1, normn, line_wgt;
  double  Gamma, ngrad_Gamma;    
  double D = TDatabase::ParamDB->REACTOR_P21 / TDatabase::ParamDB->REACTOR_P16 ;

  
  TBaseCell *cell;
  TFEDesc2D *FeDesc;
  BaseFunct2D *BaseFuncts;
  TCollection *Coll;
  TJoint *joint;
  TIsoBoundEdge *isoboundedge;
  TBoundComp *BoundComp;
  BoundCond Cond0, Cond1;
  FE2D FEId;
  TFE2D *ele;
  RefTrans2D RefTrans;
  TRefTrans2D *F_K;
  TFESpace2D *fespace;
  BF2DRefElements RefElement;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1;  
  
#ifdef __SURFACT__ 
  TFESpace2D *surfactantspace;  
  FE2D TFEId;
  TFEDesc2D *TFeDesc;  
  
  double **Turef, **Tuxiref, **Tuetaref;  
  double Tuorig[MaxN_BaseFunctions2D], Tuxorig[MaxN_BaseFunctions2D];
  double Tuyorig[MaxN_BaseFunctions2D];
  double T_val[3], *S_Values, S;
  
  
  int  *TGlobalNumbers, *TBeginIndex, local_dof;
  int TN_BaseFunct, *TJointDOF, TN_DOF_Local, *TDOF;

// surfactant elasticity E
  double E = TDatabase::ParamDB->REACTOR_P10;
//Equation of state, 0 linear, 1 non-linear
  int EOS = int(TDatabase::ParamDB->REACTOR_P11);
//\Gamma_1/Gamma_\infty

  double st = TDatabase::ParamDB->SURF_TENSION;
  
  surfactantspace = Surfact->GetFESpace2D();
  S_Values=Surfact->GetValues();
  TGlobalNumbers = surfactantspace->GetGlobalNumbers();
  TBeginIndex = surfactantspace->GetBeginIndex();  
  
//   Gamma_Max = 0;
//   N_Surf = Surfact->GetLength();

//   for(i=0;i<N_Surf;i++)
//    if(Gamma_Max<S_Values[i]) Gamma_Max=S_Values[i];
// 
//    OutPut("Gamma_Max " << Gamma_Max<<endl);  
  
#endif   
  
 
  
  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFuncts = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  fespace = A11->GetFESpace();
  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  RowPtr = A11->GetRowPtr();
  KCol = A11->GetKCol();

  ValuesA11 = A11->GetEntries();
  ValuesA22 = A22->GetEntries();

  double Re = TDatabase::ParamDB->RE_NR;
  double We = TDatabase::ParamDB->WB_NR, U;
  double Ca = We/Re, D_Angle;
  double beta = TDatabase::ParamDB->FRICTION_CONSTANT;

  double EQ_Angle = TDatabase::ParamDB->EQ_CONTACT_ANGLE;
  
  EQ_Angle = (3.141592654/180)*EQ_Angle;

 for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_Edges = cell->GetN_Edges();
    IJoint = 0;
    for(j=0;j<N_Edges;j++)
     {
      joint = cell->GetJoint(j);
      if(joint->GetType() == IsoBoundEdge)
      {
        isoboundedge = (TIsoBoundEdge *)joint;
        BoundComp = isoboundedge->GetBoundComp();
        isoboundedge->GetParameters(t0, t1);
        comp=BoundComp->GetID();
        BoundaryCondition(comp, t0, Cond0);
        BoundaryCondition(comp, t1, Cond1);

        if(Cond0 == FREESURF)
        {
          JointNumbers[IJoint] = j;
          IJoint++;
        }
      } // endif
     } // endfor j

    N_IsoJoints = IJoint;
    if(N_IsoJoints > 0)
    {
      FEId = fespace->GetFE2D(i, cell);
      
#ifdef __SURFACT__ 
      TFEId = surfactantspace->GetFE2D(i, cell);
#endif
       
      for(j=0;j<N_IsoJoints;j++)
      {
//      cout << "Cell " << i << " has free surface." << endl;
        IJoint = JointNumbers[j];
        // cout << "joint number: " << IJoint << endl;
        cell->GetVertex(IJoint)->GetCoords(x0, y0);
        cell->GetVertex((IJoint+1) % N_Edges)->GetCoords(x1, y1);
        //   if(y0==0||y1==0)
        //   cout<< " y0= " <<y0<<" y1= "<<y1<<"  x0= "<<x0<<"  x1= "<<x1<<endl;
        //   cout<< " N_LinePoints= " <<N_LinePoints<<endl;

     // entries for wetting DOF
      if(y0==0) // right wett point edge (bottom)
       {
        FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
        JointDOF = FeDesc->GetJointDOF(IJoint);
        N_DOF = FeDesc->GetN_JointDOF();
#ifdef __SURFACT__ 
        TFeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(TFEId);
        TJointDOF = TFeDesc->GetJointDOF(IJoint);
        TN_DOF_Local = TFeDesc->GetN_JointDOF();
        TDOF = TGlobalNumbers + TBeginIndex[i];
       
        for(m=0;m<TN_DOF_Local;m++)	
         {
          SDOF_R = TDOF[TJointDOF[m]];
          surfactantspace->GetDOFPosition(SDOF_R, x, y);

          if(y==0) // right wett point
           {break;}
         }
         
         if(m==TN_DOF_Local)
           OutPut("Error in finding the Surfact Wetting Pt. " <<endl);     
#endif      
        for(m=0;m<N_DOF;m++)
         {
          DOF_R =  GlobalNumbers[BeginIndex[i]+JointDOF[m]];
          fespace->GetDOFPosition(DOF_R, x, y);
  
          if(y==0) // right wett point
          {
           U = Ucl[DOF_R];

#ifdef __SURFACT__ 
           S = S_Values[SDOF_R];

           if(EOS==0)
            {
             Gamma =(1. + E*(D - S) ); 
            }
           else
           {
	     if(S<0.9)
	     { Gamma =(1. + E*log(1. - S) );
	       
	      if(Gamma<0.05) Gamma=0.05;
	    }
	     else
	     {
	       Gamma = 0.05;
	    }
	    
           //see SolubleSurf JCP paper
           if(Gamma<0.1)
            {
             Gamma = 0.1;
            }
          }             
//        cout<< "Surfact at Wetting point :" << S << " Gamma " << Gamma <<endl;     
       if(fabs(cos(EQ_Angle)/Gamma)>1.)
        {

         if( (cos(EQ_Angle)/Gamma) < -1.)
          { EQ_Angle = acos(-1.); }
         else
          { EQ_Angle = acos(1.);  }

         OutPut("  x= "<< x <<"  y= "<< y << " Gamma " << S<<  " EQ_Angle: " <<  (180./Pi)*EQ_Angle<< endl); 
// 	 exit(0);
        }
      else
        {
         EQ_Angle = acos(cos(EQ_Angle)/Gamma); 
        }
       
//         exit(0);
#endif
         switch((int)TDatabase::ParamDB->CONTACT_ANGLE_TYPE)
          {
           case 0:
              D_Angle = EQ_Angle;
           break;

           case 1:
            if(U>0.)
	     { D_Angle = (3.141592654/180.)*TDatabase::ParamDB->AD_CONTACT_ANGLE; } // advanving angle
            else if(U<0.)
	     { D_Angle = (3.141592654/180.)*TDatabase::ParamDB->RE_CONTACT_ANGLE; }// receding angle
            else
             {
              D_Angle = EQ_Angle;
//                        + tanh(50.*U)*(TDatabase::ParamDB->AD_CONTACT_ANGLE
//                                        - TDatabase::ParamDB->RE_CONTACT_ANGLE);
             }  
           break;

           case 2:
//   Hocking's expression
              D_Angle = pow(EQ_Angle, (double)3.0) 
                           + 9.0 * Ca * (fabs(U))* log(beta);

              D_Angle = pow(fabs(D_Angle), 1./3.);
           break;

           case 3:
              Ca *= fabs(U); // capillary number w.r.t contact line velocity
//   Jiang et al (1979) expression
              D_Angle = acos( cos(EQ_Angle) - (cos(EQ_Angle) + 1.)*tanh( 4.96*pow(Ca,0.702) )   );
           break;
           case 4:
              Ca *= fabs(U); // capillary number w.r.t contact line velocity
//   Bracke et al (1989) expression
              D_Angle = acos( cos(EQ_Angle) - (cos(EQ_Angle) + 1.)*( 2.*pow(Ca,0.5) )   );
           break;
           case 5:
              Ca *= fabs(U); // capillary number w.r.t contact line velocity
//   Berg et al (1992) expression
              D_Angle = acos( cos(EQ_Angle) - (cos(EQ_Angle) + 1.)*( 2.24*pow(Ca,0.54) )   );
           break;

           case 6:
//  Berg et al (1992) expression
              Ca *= fabs(U); // capillary number w.r.t contact line velocity
              D_Angle = acos( cos(EQ_Angle) - (cos(EQ_Angle) + 1.)*( 4.47*pow(Ca,0.42) )   );
           break;

          }
// OutPut("  x= "<< x <<"  y= "<< y << " U " << U<<  " D_Angle: " << (180./Pi)*D_Angle<< endl);
// exit(0);
          param[0] = x;
          param[1] = y;
          param[2] = U;
          r_axial = x;       // r value in the axial symmetric integral
          rhs1[DOF_R] +=  Gamma*r_axial*((cos(D_Angle))/We);   break;
         }
        }
       }

      DOF = GlobalNumbers + BeginIndex[i];
      N_BaseFunct = N_BaseFuncts[FEId];
      ele = TFEDatabase2D::GetFE2D(FEId);
      RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);

      l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
      LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
      qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
      qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);
      TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(LineQuadFormula);
   
#ifdef __SURFACT__ 
       TFEDatabase2D::GetBaseFunct2DFromFE2D(TFEId)->MakeRefElementData(LineQuadFormula);
       TFeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(TFEId);
       TN_BaseFunct = N_BaseFuncts[TFEId];
       TJointDOF = TFeDesc->GetJointDOF(IJoint);
       TN_DOF_Local = TFeDesc->GetN_JointDOF();
       TDOF = TGlobalNumbers + TBeginIndex[i];
#endif

      switch(RefElement)
      {
        case BFUnitSquare:
          RefTrans = QuadIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TQuadIsoparametric *)F_K)->SetCell(cell);
          ((TQuadIsoparametric *)F_K)->GetOrigBoundFromRef(IJoint, N_LinePoints, zeta, X_B, Y_B);
        break;

        case BFUnitTriangle:
          RefTrans = TriaIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TTriaIsoparametric *)F_K)->SetCell(cell);
          ((TTriaIsoparametric *)F_K)->GetOrigBoundFromRef(IJoint, N_LinePoints, zeta, X_B, Y_B);

        break;
      } // endswitch


       uref = TFEDatabase2D::GetJointValues2D(BaseFuncts[FEId], LineQuadFormula, IJoint);
       uxiref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId], LineQuadFormula, IJoint, D10);
       uetaref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId], LineQuadFormula, IJoint, D01);

#ifdef __SURFACT__ 
       Turef = TFEDatabase2D::GetJointValues2D(BaseFuncts[TFEId], LineQuadFormula, IJoint);
//        Tuxiref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[TFEId], LineQuadFormula, IJoint, D10);
//        Tuetaref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[TFEId], LineQuadFormula, IJoint, D01);       
#endif
       for(k=0;k<N_LinePoints;k++)
        {
         switch(RefElement)
          {
            case BFUnitSquare:
              ((TQuadIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
// #ifdef __SURFACT__                         
//               ((TQuadIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
//                         TN_BaseFunct, Turef[k], Tuxiref[k], Tuetaref[k],
//                         Tuorig, Tuxorig, Tuyorig);
// #endif 
            break;

            case BFUnitTriangle:
              ((TTriaIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
                        
// #ifdef __SURFACT__                         
//               ((TTriaIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
//                         TN_BaseFunct, Turef[k], Tuxiref[k], Tuetaref[k],
//                         Tuorig, Tuxorig, Tuyorig);
// #endif                         
            break;
          } // endswitch

          // modify matrices
         F_K->GetTangent(IJoint, zeta[k], t0, t1);  // old line
         r_axial = fabs(X_B[k]);   // r value in the axial symmetric integral
         normn = sqrt(t0*t0+t1*t1);
         n0 =  t1/normn;
         n1 = -t0/normn;
         t1 /= normn;
         t0 /= normn;     

#ifdef __SURFACT__   
      for(mm=0;mm<TN_BaseFunct;mm++)
           Tuorig[mm] = Turef[k][mm];  
      
       T_val[0] = 0.;
//        T_val[1] = 0.; T_val[2] = 0.;
       
       for(l=0;l<TN_DOF_Local;l++)
        {
          // assumed that the velo space and 2D surfactant space are same fe space
          local_dof   = TJointDOF[l];
          m = TDOF[local_dof];

//           if(S_Values[m]<0)
// 	  {
// 	   OutPut(i<< "x : "<<X_B[k]<< " y: " << Y_B[k] <<"  Surfactant exceeds the reference value, T_val= " <<S_Values[m]<<endl);
// 	    S_Values[m] = 0.;
// 	  }
	  
	  
          val = S_Values[m];
          T_val[0] += val*Tuorig[local_dof];  // Surfactant C
        } // for(l=0;l<TN_

        if(T_val[0]<0. )
	 {
          OutPut(i<< "x : "<<X_B[k]<< " y: " << Y_B[k] <<"  Surfactant exceeds the reference value, T_val= " <<T_val[0]<<endl);
//        for(l=0;l<TN_DOF_Local;l++)
//         {
//           // assumed that the velo space and 2D surfactant space are same fe space
//           local_dof   = TJointDOF[l];
//           m = TDOF[local_dof];
//           val = S_Values[m];
// 	  
// 	   OutPut(l << " sval  "<<val<<endl);
// 	
// 	   
// 	}
          //numerical correction
          T_val[0]=0.; 
         }  

       // Marangoni effect in weak formulation
         if(EOS==0)
          {
           Gamma =(1. + E*(D - T_val[0]) ); 
           }
         else
          {
           Gamma =(1. + E*log(1. - T_val[0]) );
   
           //see SolubleSurf JCP paper
           if(Gamma<0.1)
            {
             Gamma = 0.1;
            }
          }              
#else
    Gamma = 1;
#endif       
          // Multiply with time step dt in the main program not here
          r = normn/We;
          for(l=0;l<N_BaseFunct;l++)
          {
           TestDOF = DOF[l];

           // updating rhs
            ngrad_test= n0*uxorig[l] + n1*uyorig[l];
            d1 = uxorig[l] - ngrad_test*n0;
            d2 = uyorig[l] - ngrad_test*n1;

// rhs1
//             val = r_axial*( (1.-n0*n0)*(d1 - GammaE1) - n0*n1*(d2 -GammaE2) );
	    
            val = r_axial*( (1.-n0*n0)*d1 - n0*n1*d2 );
            val += uorig[l]; // due to axialsymmetric
            val *= LineWeights[k]*r*Gamma;
            rhs1[TestDOF] -= val;

// rhs2
            val =  r_axial*( -n1*n0*d1 + (1.-n1*n1)*d2 );   
            val *= LineWeights[k]*r*Gamma;
            rhs2[TestDOF] -= val;
 
            index2 = RowPtr[TestDOF+1];
//               cout << TestDOF  << " RhsA  " << rhs1[TestDOF] << " RhsB  " << rhs2[TestDOF] << endl;
	    
            for(m=0;m<N_BaseFunct;m++)
            {
              AnsatzDOF = DOF[m];
              // cout << AnsatzDOF << " -- " << TestDOF << endl;
              index1 = RowPtr[TestDOF];
              if(index1+1 == index2) continue;
              while(KCol[index1] != AnsatzDOF) index1++;

              ngrad_ansatz= n0*uxorig[m] + n1*uyorig[m];
              e1 = uxorig[m] - ngrad_ansatz*n0;
              e2 = uyorig[m] - ngrad_ansatz*n1;

              val =d1*e1 + d2*e2 + (uorig[l]*uorig[m]/(r_axial*r_axial));
              val *= dt*LineWeights[k]*r*Gamma*r_axial;
              // cout << "A11: " << TestDOF << " ";
              // cout << AnsatzDOF << " " << val << endl;
              ValuesA11[index1] += val;

              val = d1*e1 + d2*e2;
              val *= dt*LineWeights[k]*r*Gamma*r_axial;

              // cout << "A22: " << TestDOF << " ";
              // cout << AnsatzDOF << " " << val << endl;
              ValuesA22[index1] += val;
            } // endfor m
          } // endfor l
        } // endfor k
      } // endfor j

    } // end (N_IsoJoints > 0)
  } // endfor i
 } //FreeSurf_axial3D_new(TSquareMatr

 

void GetSurfErrors(TFEFunction2D *fefunction, TFEFunction1D *fefunct_low, double *errors,
                   int *Cell_array, int *Joint_array)
{
  int i, j, k, l, m, n, N, N_Cells_low, N_LocalUsedElements, local_i, local_dof, ORDER;
  int N_BaseFunct,  N_Points, N_JointDOF, begin, end, *N_BaseFuncts, IJoint;
  TCollection *Coll, *Coll_low;
  TBaseCell *Me, *Me_low;
  FE2D FEId;
  TFE2D *ele;
  FE1D FEId_low;
  TFE1D *Element;
  TFESpace2D *fespace;
  TFESpace1D *fespace_low;
  BaseFunct2D LocBF[N_BaseFuncts2D];
  BaseFunct2D *BaseFuncts;
  int *DOF, *JointDOF, *BeginIndex, *GlobalNumbers, N_BaseFunct_low;
  TFEDesc2D *FeDesc;
  TFEDesc1D *FeDesc_low;
  BF2DRefElements RefElement;
  RefTrans2D RefTrans;
  TRefTrans2D *F_K;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1;
  int N_LinePoints;
  double *LineWeights, t0, t1, normn, n0, n1, X_P, Y_P;
  double *u, **uref, **uxiref, **uetaref, **uzetaref;
  double *Weights, *p1, *p2, *zeta, LocL2U, LocH1U,  Surf[100], Surf_DOF[10], Surf_val[10];
  double r_axial;
  double  X_B[MaxN_QuadPoints_2D], Y_B[MaxN_QuadPoints_2D], Exact_Surf[10], v, Mult;
  double uorig[MaxN_BaseFunctions2D], uxorig[MaxN_BaseFunctions2D];
  double uyorig[MaxN_BaseFunctions2D];
  double test00, test10, test01;
  double U, ux, uy, d1, d2, ngrad, h_K =0., h_K_min=1e8, h_K_max=-1e8;

  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFuncts = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  fespace=fefunction->GetFESpace2D();
  Coll = fespace->GetCollection();
  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();
  u = fefunction->GetValues();

  fespace_low=fefunct_low->GetFESpace1D();
  Coll_low = fespace_low->GetCollection();
  N_Cells_low = Coll_low->GetN_Cells();

  errors[0] =0.; // L2-norm
  errors[1] =0.; // H1-seminorm

// ########################################################################
// loop over all surf cells
// ########################################################################
  for(i=0;i<N_Cells_low;i++)
  {
    h_K =0.;
    Me_low = Coll_low->GetCell(i);
    FEId_low = fespace_low->GetFE1D(i, Me_low);
    Element = TFEDatabase2D::GetFE1D(FEId_low);
    N_BaseFunct_low = Element->GetN_DOF();

    N = Cell_array[i];
    Me = Coll->GetCell(N);
    IJoint = Joint_array[i];

    FEId = fespace->GetFE2D(N, Me);
    ele = TFEDatabase2D::GetFE2D(FEId);

    RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
    FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
    N_JointDOF =  FeDesc->GetN_JointDOF();
    JointDOF = FeDesc->GetJointDOF(IJoint);
    N_BaseFunct = FeDesc->GetN_DOF();
    DOF = GlobalNumbers + BeginIndex[N];

    l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
    ORDER = TFEDatabase2D::GetAccuracyFromFE2D(FEId);
    LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
    qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);

    switch(RefElement)
     {
      case BFUnitSquare:
          RefTrans = QuadIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TQuadIsoparametric *)F_K)->SetCell(Me);
          ((TQuadIsoparametric *)F_K)->GetOrigBoundFromRef(IJoint,
	                            N_LinePoints, zeta, X_B, Y_B);
          ((TTriaIsoparametric *)F_K)->SetApproximationOrder(ORDER);
      break;

      case BFUnitTriangle:
          RefTrans = TriaIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TTriaIsoparametric *)F_K)->SetCell(Me);
          ((TTriaIsoparametric *)F_K)->GetOrigBoundFromRef(IJoint,
	                            N_LinePoints, zeta, X_B, Y_B);
          ((TTriaIsoparametric *)F_K)->SetApproximationOrder(ORDER);

      break;
    } // endswitch

    TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)
                  ->MakeRefElementData(LineQuadFormula);


      uref = TFEDatabase2D::GetJointValues2D(BaseFuncts[FEId],
                        LineQuadFormula, IJoint);
      uxiref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId],
                        LineQuadFormula, IJoint, D10);
      uetaref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId],
                        LineQuadFormula, IJoint, D01);


    for(k=0;k<N_LinePoints;k++)
     {

         switch(RefElement)
          {
            case BFUnitSquare:
              ((TQuadIsoparametric *)F_K)->SetCell(Me);
              ((TQuadIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
              ((TQuadIsoparametric *) F_K)->GetTangent(IJoint, zeta[k], t0, t1);
            break;

            case BFUnitTriangle:
              ((TTriaIsoparametric *)F_K)->SetCell(Me); 
              ((TTriaIsoparametric *)F_K)->SetApproximationOrder(ORDER);
              ((TTriaIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
              ((TTriaIsoparametric *) F_K)->GetTangent(IJoint, zeta[k], t0, t1);
            break;
          } // endswitch


      normn = sqrt(t0*t0+t1*t1);
      n0 =  t1/normn;
      n1 = -t0/normn;

//  //  Projection of all points to the free boundary in normal direction  !! only for circle!!!
//       v =  sqrt(X_B[k]* X_B[k] + Y_B[k]*Y_B[k] );// /r;
      X_P =  X_B[k]  ;
      Y_P =  Y_B[k]  ;
      r_axial = fabs(X_P);

//   get solution and its gradients
     U=0.; ux=0.;  uy=0.; 
     for(l=0;l<N_BaseFunct_low;l++)
      {
       local_dof   = JointDOF[l];

       test00 = uorig[local_dof];
//        test10 = uxorig[local_dof];
//        test01 = uyorig[local_dof];

//        ngrad = test10*n0 + test01*n1;
//        d1 = test10 - ngrad*n0;
//        d2 = test01 - ngrad*n1;

       m = DOF[local_dof];

       U  += u[m]*test00;
//        ux += u[m]*d1;
//        uy += u[m]*d2;
      }

      ExactS(X_P, Y_P, Exact_Surf);

      Mult = LineWeights[k]*normn*r_axial; // since the line integral is -1 to +1
      h_K +=normn;
      errors[0] +=Mult*(Exact_Surf[0]-U)*(Exact_Surf[0]-U);
// //       errors[0] +=Mult;
     } //  for(k=0;k<N_LinePoints;k++)

     if(h_K_min>h_K)
         h_K_min = h_K;
     if(h_K_max<h_K)
         h_K_max = h_K;

  } // for(i=0;i<N
//     OutPut("h_K_min and h_K_max of free surface: "<< h_K_min << " " << h_K_max<<endl; );
   errors[0] = sqrt(errors[0]);
//    cout << "errors[0] " << errors[0]<< " "<<endl;

// exit(0);
}



// ====================================================================
// Get the inner angles of the cells in whole domain
// ====================================================================
void Getcellangle(TFESpace2D *Space, double *MinMaxAngle)
{
 int i,j,k,l, N_Cells, N_Edges;
 int found,  N_LinePoints;

 double TX[4], TY[4], hE[4], Theta, tx, ty, Test, MQI=0.;
 TBaseCell *cell;
 FE2D FEId;
 BF2DRefElements RefElement;
 TRefTrans2D *F_K;
 RefTrans2D RefTrans;
 TCollection *Cells;

  MinMaxAngle[0] = 180;  // Min_Angel = 180
  MinMaxAngle[1] = 0;  // Max_Angel = 0
  Cells = Space->GetCollection();
  N_Cells = Cells->GetN_Cells();
     
//      TX      = new double[4];  // Max no edges in 2d
//      TY      = new double[4];  // Max no edges in 2d
//      hE      = new double[4];  // Max no edges in 2d

  for(i=0;i<N_Cells;i++)
   {
     cell    = Cells->GetCell(i);
     N_Edges = cell->GetN_Edges();

     FEId = Space->GetFE2D(i, cell);
     RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);

     switch(RefElement)
        {
         case BFUnitTriangle:

            RefTrans = TriaAffin;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            ((TTriaAffin *)F_K)->SetCell(cell);

          break;

          case BFUnitSquare:

            RefTrans = QuadAffin;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            ((TQuadAffin *)F_K)->SetCell(cell);

          break;

          default:
            Error("only triangles and quadrilaterals are allowed" << endl);
            Error("file: " << __FILE__ << " line " << __LINE__ << endl);
            exit(-1);
          } // endswitch

     for(j=0;j<N_Edges;j++)
      {
        F_K->GetTangent(j, 0, tx, ty);
        TX[j] = tx;
        TY[j] = ty;
        hE[j] = sqrt(tx*tx+ty*ty);

   // cout <<"cell : " <<i << "  j= " << j << ": " <<TX[j]<< "------ " << TY[j] << endl;
       } // endfor j

//      Test = 0;
      k = N_Edges -1;
      for(j=0;j<N_Edges;j++)
      {
       if(hE[j]==0.0 || hE[k]== 0.0 )
        Theta = 0.0;
       else
        Theta = acos(-(TX[j]*TX[k]+TY[j]*TY[k])/(hE[j]*hE[k]))*(180/3.141592654);

       k = j;
//        Test +=Theta;
       if(MinMaxAngle[0]>Theta) MinMaxAngle[0] = Theta;
       if(MinMaxAngle[1]<Theta) MinMaxAngle[1] = Theta;
//        cout <<"cell : " <<i << "  j= " << j << ": " << " Theta : " << Theta << endl;
//  *****************************************************
//  Grid test

      MQI += (60. - Theta)*(60. - Theta);
//  *****************************************************

     }
//       cout <<"cell : " <<i <<  " sum of 3 angels : " << Test << endl;
     //  cout<<endl;

   } // endfor i

   MQI /=double(3*N_Cells);
   MQI = sqrt(MQI);

// OutPut("Mesh Quality Indicator: "<< MQI<< endl);
//    delete [] TX;
//    delete [] TY;
//    delete [] hE;
 //cout<< " Min_Angel: "<< MinMaxAngle[0]<< "  Max_Angel : "<<MinMaxAngle[1]<< endl;
// exit(0);
}


void SolveGridEquation(double **Entries, double *sol, double *rhs,
                       int *KCol, int *RowPtr, int N_DOF)
{
  int i,j,k, col, Diognal;
  double *Entries11, *Entries12, *Entries21, *Entries22;
  double sum1, sum2, max_sum1, max_sum2;
  int start, end;

  double max_error, error=1.e-12;
  int iter;

  Entries11 = Entries[0];
  Entries12 = Entries[1];
  Entries21 = Entries[2];
  Entries22 = Entries[3];
  
  max_error = 1.; iter = 0;
  while(max_error>error)
  {
    max_error = 0.0; iter++;
    for(i=0;i<N_DOF;i++)
    {
      start = RowPtr[i];
      end = RowPtr[i+1];
      sum1 = rhs[i];
      sum2 = rhs[i+N_DOF];
      for(k=start;k<end;k++)
      {
        col = KCol[k];
        if (col==i) Diognal = k;
        sum1 -= Entries11[k] * sol[col]
              +Entries12[k] * sol[col+N_DOF];
        sum2 -= Entries21[k] * sol[col]
              +Entries22[k] * sol[col+N_DOF];
      } // endfor k
     // sol[i] += sum1/Entries11[start];
     // sol[i+N_DOF] += sum2/Entries22[start];
        sol[i] += sum1/Entries11[Diognal];
        sol[i+N_DOF] += sum2/Entries22[Diognal];
      if(max_error<fabs(sum1/Entries11[Diognal])) max_error = fabs(sum1/Entries11[Diognal]);
      if(max_error<fabs(sum2/Entries22[Diognal])) max_error = fabs(sum2/Entries22[Diognal]);
    } // endfor i
    if(iter == 1000) break;
  } // end while
//OutPut("Grid Solver: Number iteration "<<iter<<endl);
}

void GridVelo_imping(double **Entries, double *Sol, double *d, double *Rhs,
                     int *KCol, int *RowPtr,
                     TFEVectFunct2D *GridPos,
                     TFEVectFunct2D *AuxGridPos,
                     TFEVectFunct2D *Velocity, double dt,
                     TFEVectFunct2D *GridVelocity, 
                     TVertex ***MovBoundVert, int *N_MovVert,
                     TBaseCell **Free_Cells, int **IsoCellEdgeNos,
                     bool &reparam, TFEVectFunct2D *RefGridPos)
{
  int i,j,k,l,m,comp, N;
  int *VeloBeginIndex, *VeloGlobalNumbers;
  int *GridBeginIndex, *GridGlobalNumbers;
  int N_Cells, N_Vertices, N_Edges, N_LocalDOFs;
  int N_Levels, *DOF, *JointDOF, GridLength;
  int N_BoundaryNodes, N_LinePoints, IIso, N_Inner, N_;
  
  double xi[4], eta[4], X[4], Y[4], VX[4], VY[4];
  double FunctValues[4][MaxN_BaseFunctions2D];
  double FEValuesX[MaxN_BaseFunctions2D];
  double FEValuesY[MaxN_BaseFunctions2D];
  double *ValuesX, *ValuesY, *RefValueX, *RefValueY;
  double *ValuesVX, *ValuesVY, *NewValuesX, *NewValuesY;
  double s, t, x, y, h_tot, x0, x1, y0, y1;
  double res, oldres, *gridvelo, *Nx, *Ny, Ay;
  double *LineWeights, *zeta;
  double normalx, normaly, tangenx, tangeny, nx, ny, tx, ty;
  double un, hE, t0,t1, temp2, eps=1e-6;
  double h, hmin, hmax, hlimit, *IsoX, *IsoY;
   
  BF2DRefElements RefElement;
  TRefTrans2D *F_K;
  RefTrans2D RefTrans;
  TIsoBoundEdge *isojoint;  
  TMGLevel2D *Level;
  TFESpace2D *VelocitySpace, *GridSpace;
  TCollection *Coll;
  TBaseCell *cell;
  BoundTypes bdtype;
  TBoundEdge *BoundEdge;
  TBoundComp2D *BoundComp;
  FE2D FEId;
  TFE2D *Element;
  TFEDesc2D *FEDesc;
  BaseFunct2D BF;
  TBaseFunct2D *bf;
  boolean OnBoundary;
  TJoint *joint;
  TVertex **Vertices;
  TQuadFormula1D *qf1;
  QuadFormula1D LineQuadFormula;
  BoundCond Cond0, Cond1;
    
  VelocitySpace = Velocity->GetFESpace2D();
  VeloBeginIndex = VelocitySpace->GetBeginIndex();
  VeloGlobalNumbers = VelocitySpace->GetGlobalNumbers();
  ValuesVX = Velocity->GetValues();
  ValuesVY = ValuesVX + Velocity->GetLength();

  GridPos->GridToData();
  GridSpace = GridPos->GetFESpace2D();
  GridBeginIndex = GridSpace->GetBeginIndex();
  GridGlobalNumbers = GridSpace->GetGlobalNumbers();
  GridLength = GridPos->GetLength();
  ValuesX = GridPos->GetValues();
  ValuesY = ValuesX + GridLength;

  N_Inner = GridSpace->GetN_Inner();
  N_BoundaryNodes = GridLength - GridSpace->GetN_Inner();
  
  if(TDatabase::ParamDB->P5 > 0)
  {
   Nx = new double[N_BoundaryNodes];
   Ny = new double[N_BoundaryNodes];
   memset(Nx, 0, N_BoundaryNodes*SizeOfDouble);
   memset(Ny, 0, N_BoundaryNodes*SizeOfDouble);
  }

 // cout << "N_BoundaryNodes: " << N_BoundaryNodes << endl;
//   cout << GridLength << " " << N_DOF << endl;

  NewValuesX = AuxGridPos->GetValues();
  NewValuesY = NewValuesX + GridLength;

  memcpy(NewValuesX, ValuesX, 2*GridLength*SizeOfDouble);

  RefValueX = RefGridPos->GetValues();
  RefValueY = RefValueX + GridLength;   
  
  Coll = VelocitySpace->GetCollection();
  N_Cells = Coll->GetN_Cells();
//  cout << "N_Cells: " <<N_Cells<< endl;
  // determine outer normal vectors

  IIso = N_BoundaryNodes;
  // Outward normal no need if we move boundary with velocity
 if(TDatabase::ParamDB->P5 > 0)
  {
  // determine outer normal vectors
   for(i=0;i<N_Cells;i++)
    {
    // cout << "cell: " << i << endl;
    cell = Coll->GetCell(i);
    cell->SetClipBoard(i);
    N_Edges = cell->GetN_Edges();

    for(j=0;j<N_Edges;j++)
    {
      if( !(cell->GetJoint(j)->InnerJoint()) )
      {
        // cout << "joint: " << j << endl;
        cell->GetVertex(j)->GetCoords(x0, y0);
        cell->GetVertex((j+1)%N_Edges)->GetCoords(x1, y1);
        t = x1-x0;
        s = y1-y0;
        FEId = VelocitySpace->GetFE2D(i, cell);
        l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
        LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
        qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
        qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);

        RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
        switch(RefElement)
        {
          case BFUnitTriangle:
            RefTrans = TriaIsoparametric;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            ((TTriaIsoparametric *)F_K)->SetCell(cell);
          break;

          case BFUnitSquare:
            RefTrans = QuadIsoparametric;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            ((TQuadIsoparametric *)F_K)->SetCell(cell);
          break;

          default:
            Error("only triangles and quadrilaterals are allowes" << endl);
            Error("file: " << __FILE__ << " line " << __LINE__ << endl);
            exit(-1);
        } // endswitch

        normalx = 0;
        normaly = 0;
        hE = 0;
        for(k=0;k<N_LinePoints;k++)
        {
          F_K->GetOuterNormal(j, zeta[k], nx, ny);
          F_K->GetTangent(j, zeta[k], tx, ty);
          t = sqrt(tx*tx+ty*ty);
          normalx += t*LineWeights[k]*nx;
          normaly += t*LineWeights[k]*ny;
          hE += t*LineWeights[k];
          // cout << "k= " << k << " " << nx << " " << ny << endl;
        } // endfor k

        DOF = GridGlobalNumbers + GridBeginIndex[i];

        switch(N_Edges)
        {
          case 3:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 2:
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;
            } // endswitch j
          break;

          case 4:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[3] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 2:
                l = DOF[3] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 3:
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;
            } // endswitch j
          break;
        } // endswitch N_Edges
      } // !InnerJoint
    } // endfor j
  } // endfor i

  N_ = IIso;
  // normalize normal vector
  for(i=0;i<N_;i++)
  {
    x = Nx[i];
    y = Ny[i];
    t = sqrt(x*x+y*y);
    Nx[i] /= t;
    Ny[i] /= t;
  }
 }
  // determine new position of boundary vertices
  for(i=0;i<N_Cells;i++)
  {
    cell  = Coll->GetCell(i);
    N_Edges = cell->GetN_Edges();
    OnBoundary = FALSE;
    for(j=0;j<N_Edges;j++)
    {
      if( !(cell->GetJoint(j)->InnerJoint()) || 
           (cell->GetJoint(j)->GetType() == InterfaceJoint)  || 
           (cell->GetJoint(j)->GetType() == IsoInterfaceJoint) )
       {
        OnBoundary = TRUE;
        joint = cell->GetJoint(j);
       }
    } // endfor j

    if(OnBoundary)
    {

      FEId = VelocitySpace->GetFE2D(i, cell);
      Element = TFEDatabase2D::GetFE2D(FEId);
      bf = Element->GetBaseFunct2D();
      N_LocalDOFs = Element->GetN_DOF();

      switch(N_Edges)
      {
        case 3:
          xi[0]  = 0; xi[1]  = 1; xi[2]  = 0;
          eta[0] = 0; eta[1] = 0; eta[2] = 1;
          VX[0] = VX[1] = VX[2] = 0;
          VY[0] = VY[1] = VY[2] = 0;
          N_Vertices = N_Edges;
        break;

        case 4:
          // NOTE: sorted due to number of DOFs !!!
          xi[0]  = -1; xi[1]  =  1; xi[2]  = -1; xi[3]  =  1;
          eta[0] = -1; eta[1] = -1; eta[2] =  1; eta[3] =  1;
          VX[0] = VX[1] = VX[2] = VX[3] = 0;
          VY[0] = VY[1] = VY[2] = VY[3] = 0;
          N_Vertices = N_Edges;
        break;

        default:
          Error("only triangles and quadrilaterals are allowed!" << endl);
          exit(-1);
      } // endswitch

      for(j=0;j<N_Vertices;j++)
        bf->GetDerivatives(D00, xi[j], eta[j], FunctValues[j]);

      DOF = VeloGlobalNumbers + VeloBeginIndex[i];

      for(j=0;j<N_LocalDOFs;j++)
      {
        k = DOF[j];
        s = ValuesVX[k];
        t = ValuesVY[k];

        for(l=0;l<N_Vertices;l++)
        {
          VX[l] += FunctValues[l][j]*s;
          VY[l] += FunctValues[l][j]*t;
        } // endfor l
      } // endfor j

      FEId = GridSpace->GetFE2D(i, cell);
      Element = TFEDatabase2D::GetFE2D(FEId);
      BF = Element->GetBaseFunct2D_ID();
      if( (BF != BF_C_T_P1_2D) && (BF != BF_C_Q_Q1_2D) )
      {
        Error("Grid Space must be conforming and of first order!" << endl);
        exit(-1);
      }  // endif
      bf = Element->GetBaseFunct2D();
      N_LocalDOFs = Element->GetN_DOF();

      switch(N_Edges)
      {
        case 3:
          xi[0]  = 0; xi[1]  = 1; xi[2]  = 0;
          eta[0] = 0; eta[1] = 0; eta[2] = 1;
          X[0] = X[1] = X[2] = 0;
          Y[0] = Y[1] = Y[2] = 0;
          N_Vertices = N_Edges;
        break;

        case 4:
          // NOTE: sorted due to number of DOF !!!
          xi[0]  = -1; xi[1]  =  1; xi[2]  = -1; xi[3]  =  1;
          eta[0] = -1; eta[1] = -1; eta[2] =  1; eta[3] =  1;
          X[0] = X[1] = X[2] = X[3] = 0;
          Y[0] = Y[1] = Y[2] = Y[3] = 0;
          N_Vertices = N_Edges;
        break;

        default:
          Error("only triangles and quadrilaterals are allowed!" << endl);
          exit(-1);
      } // endswitch

      DOF = GridGlobalNumbers + GridBeginIndex[i];

      for(j=0;j<N_Vertices;j++)
      {
        l = DOF[j];
        k = l - N_Inner;
        if(k>=0)
        {
          if((TDatabase::ParamDB->P5 > 0) )
          { 
           un = VX[j]*Nx[k] + VY[j]*Ny[k];
            
           if(ValuesX[l] == 0 )
            { NewValuesX[l] = ValuesX[l]; }
           else
            { NewValuesX[l] = ValuesX[l] + dt*un*Nx[k]; }

            if(ValuesY[l] == 0) 
             { NewValuesY[l] = ValuesY[l];  }
            else    
            { NewValuesY[l] = ValuesY[l] + dt*un*Ny[k]; }
          }
          else
          {
	    if(ValuesX[l] == 0 )
	    { NewValuesX[l] = ValuesX[l]; }
	    else
	    { NewValuesX[l] = ValuesX[l] + dt*VX[j]; }

            if(ValuesY[l] == 0) 
	     { NewValuesY[l] = ValuesY[l];  }
            else    
	     { NewValuesY[l] = ValuesY[l] + dt*VY[j];  }
	    
	  }
       } //  if(k>=0)
 //    Due to spline approximation solid boundary end vertices may take negative y value
        if(NewValuesY[l]<0.0 ) NewValuesY[l] = 0.0;
        if( fabs(NewValuesX[l]) < 1e-10 ) NewValuesX[l] = 0.0;
      } // endfor j
    } // endif
  } // endfor i

// cout << " dt " << dt <<endl;
/*
   for(i=0;i<GridLength;i++)
 cout << i <<"  ---  " <<ValuesX[i] << "  ---  " << NewValuesX[i] << endl;
*/
// exit(0);
  
   MovBoundVert[0][0]->GetCoords(x, Ay);   

   AuxGridPos->DataToGrid();
   
  //======================================================================  
  //  Reparametrization of free surface - Begin
  //======================================================================  
  if(!reparam)
  {
   h_tot = 0;
   hmin = 100;
   hmax = 0.0;

   for(k=0;k<N_MovVert[2];k++)
    {
     MovBoundVert[2][k]->GetCoords(x1, y1);

     if(k==N_MovVert[2]-1)
     { MovBoundVert[1][0]->GetCoords(x, y);}
     else
     { MovBoundVert[2][k+1]->GetCoords(x, y); }

     h = sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1));
     h_tot +=h;
     if (h < hmin) hmin = h;
     if (h > hmax) hmax = h;
    } // for(k=0;k<N_MovVert[2];k++) 

    h_tot /= (double)N_MovVert[2];
    hlimit =  0.8*h_tot;   
   }

   if ( ((hmin < hlimit) || (hmax > 3.*h_tot/2.)) ||  reparam )
   {
    //before reparam update iso points, since we use interpolated cubic spline
    //which pass through iso points also
    IsoX = new double[IIso];
    IsoY = new double[IIso];

    IIso = 0;
    for(i=0;i<N_MovVert[2];i++)
    {
     cell = Free_Cells[i];
      j = IsoCellEdgeNos[1][i];
      joint = cell->GetJoint(j);
      isojoint = (TIsoBoundEdge *)joint;
      k = isojoint->GetN_Vertices();
      Vertices = isojoint->GetVertices();
      FEId = VelocitySpace->GetFE2D(IsoCellEdgeNos[0][i], cell);
      FEDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
      m = FEDesc->GetN_JointDOF();
      if(m == k+2)
        {
          JointDOF = FEDesc->GetJointDOF(j);
          DOF =  VeloGlobalNumbers+VeloBeginIndex[IsoCellEdgeNos[0][i]];
          for(l=0;l<k;l++)
          {
            m = DOF[JointDOF[l+1]];
            Vertices[l]->GetCoords(IsoX[IIso], IsoY[IIso]);
            if(TDatabase::ParamDB->P5 > 0)
            {
              un = ValuesVX[m]*Nx[IIso+N_BoundaryNodes]
                  + ValuesVY[m]*Ny[IIso+N_BoundaryNodes];
              x  = IsoX[IIso] + dt*un*Nx[IIso+N_BoundaryNodes];
              y  = IsoY[IIso] + dt*un*Ny[IIso+N_BoundaryNodes];
            }
            else
            {
             x  = IsoX[IIso] + dt*ValuesVX[m];
             y  = IsoY[IIso] + dt*ValuesVY[m];
            }

           if(y<=0) y = 1e-5;
//            if(fabs(x)<1e-12) x = 0.;

           Vertices[l]->SetCoords(x, y);
           IIso++;
          } // endfor l
        }
       else
        {
          // approximation order of isoparametric boundary and velocity
          // element must be the same
          Error("No match in isoparametric case" << endl);
          exit(-1);
        } //  if(m == k+2)
    }// for(i=0;i<N_MovVert[2];i++)

    ReParam_axial3D_U(N_MovVert[2], Free_Cells,  IsoCellEdgeNos[1], IsoCellEdgeNos[0],
                      Velocity, NULL, NULL, FALSE);  

    reparam = TRUE;
    RefGridPos->GridToData();
    Daxpy(2*GridLength, -1, ValuesX, RefValueX); // now reparamdisp in RefValueX

    //back to orig mesh (no reparam movement in calculation of free surf w)
    AuxGridPos->DataToGrid(); 

   //restore iso points
    IIso = 0;
    for(i=0;i<N_MovVert[2];i++)
    {
     cell = Free_Cells[i];
     j = IsoCellEdgeNos[1][i];
     joint = cell->GetJoint(j);
     isojoint = (TIsoBoundEdge *)joint;
     k = isojoint->GetN_Vertices();
     Vertices = isojoint->GetVertices();
     FEId = VelocitySpace->GetFE2D(IsoCellEdgeNos[0][i], cell);
     FEDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
     m = FEDesc->GetN_JointDOF();
     for(l=0;l<k;l++)
      {
        Vertices[l]->SetCoords(IsoX[IIso], IsoY[IIso]);
        IIso++;
      }
    }// for(i=0;i<N_Cells;i++)

    delete [] IsoX; delete [] IsoY; 
   } // if ( ((hmin < hlimit) || (hmax > 3.*h_tot/2.)) ||  reparam )
   //======================================================================       

   MovBoundVert[2][0]->GetCoords(x, y); // right wetting point

   y=0.;
   h_tot = x;
   h_tot /= (double)N_MovVert[0];
   for(i=1;i<N_MovVert[0];i++)
     MovBoundVert[0][i]->SetCoords(h_tot*(double)i, y);

// axial boundary
   MovBoundVert[1][0]->GetCoords(x, y);

   N=N_MovVert[1];      
   h_tot = (y-Ay)/(double)N;   
//    N--;
//      cout<< " y " << y <<" nAy " <<Ay<< " h_tot " <<h_tot<<endl; 
    for(i=1;i<N;i++)
    {
//      MovBoundVert[1][i]->GetCoords(x, y);
     y = ((double)(N-i))*h_tot; 
//      cout<< " y " << y <<" new y " << y <<endl;

     MovBoundVert[1][i]->SetCoords(x, y);   
    }

//    exit(0);

//    x=0.;
//    h_tot = -y;
//    h_tot /= (double)N_MovVert[1];
//    for(i=1;i<N_MovVert[1];i++)
//     MovBoundVert[1][i]->SetCoords(x,  y +  h_tot*(double)i );
//      
   AuxGridPos->GridToData();   
   GridPos->DataToGrid();
   
   
//       MovBoundVert[1][0]->GetCoords(x, y);
//     cout << " x " << x <<" y " << y <<endl;  
//     exit(0); 
    
   memset(Rhs, 0, 2*GridLength*SizeOfDouble);

   memcpy(d, NewValuesX, 2*GridLength*SizeOfDouble);
   Daxpy(2*GridLength, -1, ValuesX, d);
   memcpy(Rhs + (GridLength-N_BoundaryNodes), d+(GridLength-N_BoundaryNodes),
          N_BoundaryNodes*SizeOfDouble);
   memcpy(Rhs + (2*GridLength-N_BoundaryNodes), d+(2*GridLength-N_BoundaryNodes),
          N_BoundaryNodes*SizeOfDouble);

  memset(Sol, 0 , 2*GridLength*SizeOfDouble);
  memcpy(Sol + (GridLength-N_BoundaryNodes), d+(GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);
  memcpy(Sol + (2*GridLength-N_BoundaryNodes), d+(2*GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);
	 
//     for(i=0;i<GridLength;i++)
//      cout<< i <<"  ---  "<< Rhs[i] << "  ---  " << Rhs[i+GridLength] << endl;
//     
  SolveGridEquation(Entries, Sol, Rhs, KCol, RowPtr, GridLength);

  gridvelo = GridVelocity->GetValues();
  memcpy(gridvelo, Sol, 2*GridLength*SizeOfDouble);
  Dscal(2*GridLength, 1./dt, gridvelo);

  //======================================================================  
  // move inner points accordingly with the anticipation of Reparam of free surface
  // in this timestep, however the free surface points move only with the fluid velocity 
  // but the velo values are interpolated during reparam to satisfy kinematic condition
  //======================================================================      
  if(reparam)
  {
   memset(Rhs, 0, 2*GridLength*SizeOfDouble);    
   memcpy(Rhs + (GridLength-N_BoundaryNodes), RefValueX+(GridLength-N_BoundaryNodes),
          N_BoundaryNodes*SizeOfDouble);
   memcpy(Rhs + (2*GridLength-N_BoundaryNodes), RefValueX+(2*GridLength-N_BoundaryNodes),
          N_BoundaryNodes*SizeOfDouble);
 
  memset(Sol, 0 , 2*GridLength*SizeOfDouble);
  memcpy(Sol + (GridLength-N_BoundaryNodes), RefValueX+(GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);
  memcpy(Sol + (2*GridLength-N_BoundaryNodes), RefValueX+(2*GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble); 
 
  SolveGridEquation(Entries, Sol, Rhs, KCol, RowPtr, GridLength);
  
  Dscal(2*GridLength, 1./dt, Sol);
   
  //only inner mesh velo, since U will be interpolated in free surf reparm     
  Daxpy(GridLength-N_BoundaryNodes, 1., Sol, gridvelo);
  Daxpy(GridLength-N_BoundaryNodes, 1., Sol+GridLength, gridvelo+GridLength);  
  
//     cout<< "reparam : "<<endl;
//     exit(0);
  }  
  
//     for(i=0;i<GridLength;i++)
//      cout<< i <<"  ---  "<< gridvelo[i] << "  ---  " << gridvelo[i+GridLength] << endl;

//   delete [] d;
   
  if(TDatabase::ParamDB->P5 > 0)
  { delete [] Nx;  delete [] Ny; }
  
  
} // GridVelo_imping


void Get_KE(TFEVectFunct2D *Velocity, double *parameters)
 {
  int i,j,k,l, polydegree, Phase_No;
  int N_Cells, N_Joints, N_Vertices;
  TBaseCell *cell;
  TCollection *coll;
  int *BeginIndex, *GlobalNumbers, *DOF;
  TJoint *joint;
  double KE=0., volume=0., volume_1=0., volume_2=0.;
  double KE_QP, KE_1=0., x_mass=0., y_mass=0.;
  double u1_rise, u2_rise, U1_Rise=0., U2_Rise=0.;
  TRefTrans2D *F_K;
  double U1, U2;
  FE2D FEid;
  TBaseFunct2D *bf;
  int N_QFPoints;
  double *weights, *xi, *eta;
  double values[MaxN_QuadPoints_2D][MaxN_BaseFunctions2D];
  double AbsDetjk[MaxN_QuadPoints_2D], X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double *ValuesVX, *ValuesVY;
  int N_BF;
  double Mult, r_axial;
  TFESpace2D *VelocitySpace;
  JointType jointtype;
  BoundTypes bdtype;
  RefTrans2D RefTrans;
  boolean IsIsoparametric;
  QuadFormula2D QuadFormula;
  TQuadFormula2D *qf2;

  VelocitySpace = Velocity->GetFESpace2D();
  BeginIndex = VelocitySpace->GetBeginIndex();
  GlobalNumbers = VelocitySpace->GetGlobalNumbers();
  ValuesVX = Velocity->GetValues();
  ValuesVY = ValuesVX + Velocity->GetLength();

  coll = VelocitySpace->GetCollection();
  N_Cells = coll->GetN_Cells();

  for(i=0;i<N_Cells;i++)
   {
    cell = coll->GetCell(i);
    FEid = VelocitySpace->GetFE2D(i, cell);

    RefTrans = TFEDatabase2D::GetRefTrans2D_IDFromFE2D(FEid);
    N_Joints = cell->GetN_Joints();

    IsIsoparametric = FALSE;

    if(TDatabase::ParamDB->USE_ISOPARAMETRIC)
     {
     for(j=0;j<N_Joints;j++)
      {
       joint = cell->GetJoint(j);
       jointtype = joint->GetType();
       if(jointtype == BoundaryEdge)
        {
         bdtype = ((TBoundEdge *)(joint))->GetBoundComp()->GetType();
         if(bdtype != Line)  IsIsoparametric = TRUE;
        }
       if(jointtype == InterfaceJoint)
        {
          bdtype = ((TInterfaceJoint *)(joint))->GetBoundComp()->GetType();
          if(bdtype != Line)
            IsIsoparametric = TRUE;
        }
        if(jointtype == IsoInterfaceJoint || jointtype == IsoBoundEdge)
         IsIsoparametric = TRUE;

      } // for(j=0;j<
     } // if(TDatabase::ParamDB->USE_ISOPARAMETRIC)

   if(IsIsoparametric)
    {
      switch(N_Joints)
      {
        case 4:
          RefTrans = QuadIsoparametric;
        break;

        case 3:
          RefTrans = TriaIsoparametric;
        break;
      }
    } // endif IsIsoparametric

    F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);

    switch(RefTrans)
    {
      case TriaAffin:
        polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEid);
        QuadFormula = TFEDatabase2D::GetQFTriaFromDegree(9);
        qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
        qf2->GetFormulaData(N_QFPoints, weights, xi, eta);
        ((TTriaAffin *)F_K)->SetCell(cell);
        ((TTriaAffin *)F_K)->GetOrigFromRef(N_QFPoints, xi, eta,
                                         X, Y, AbsDetjk);
      break;

      case TriaIsoparametric:
        polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEid);
        QuadFormula = TFEDatabase2D::GetQFTriaFromDegree(9);
        qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
        qf2->GetFormulaData(N_QFPoints, weights, xi, eta);
        ((TTriaIsoparametric *)F_K)->SetApproximationOrder(polydegree);
        ((TTriaIsoparametric *)F_K)->SetQuadFormula(QuadFormula);
        ((TTriaIsoparametric *)F_K)->SetCell(cell);
        ((TTriaIsoparametric *)F_K)->GetOrigFromRef(N_QFPoints, xi, eta,
                                         X, Y, AbsDetjk);
      break;

      case QuadAffin:
        polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEid);
        QuadFormula = TFEDatabase2D::GetQFQuadFromDegree(3*polydegree);
        qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
        qf2->GetFormulaData(N_QFPoints, weights, xi, eta);
        ((TQuadAffin *)F_K)->SetCell(cell);
        ((TQuadAffin *)F_K)->GetOrigFromRef(N_QFPoints, xi, eta,
                                         X, Y, AbsDetjk);
      break;

      case QuadBilinear:
        polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEid);
        QuadFormula = TFEDatabase2D::GetQFQuadFromDegree(3*polydegree);
        qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
        qf2->GetFormulaData(N_QFPoints, weights, xi, eta);
        ((TQuadBilinear *)F_K)->SetCell(cell);
        ((TQuadBilinear *)F_K)->GetOrigFromRef(N_QFPoints, xi, eta,
                                         X, Y, AbsDetjk);
      break;

      case QuadIsoparametric:
        polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEid);
        QuadFormula = TFEDatabase2D::GetQFQuadFromDegree(3*polydegree);
        qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
        qf2->GetFormulaData(N_QFPoints, weights, xi, eta);
        ((TQuadIsoparametric *)F_K)->SetApproximationOrder(polydegree);
        ((TQuadIsoparametric *)F_K)->SetQuadFormula(QuadFormula);
        ((TQuadIsoparametric *)F_K)->SetCell(cell);
        ((TQuadIsoparametric *)F_K)->GetOrigFromRef(N_QFPoints, xi, eta,
                                         X, Y, AbsDetjk);
      break;
    }

    // find basis functions on cell i
    bf = TFEDatabase2D::GetBaseFunct2DFromFE2D(FEid);
    N_BF = bf->GetDimension();
    DOF = GlobalNumbers + BeginIndex[i];

    for(k=0;k<N_QFPoints;k++)
     {
      bf->GetDerivatives(D00, xi[k], eta[k], values[k]);
      r_axial = fabs(X[k]);
      if(X[k]<=0)
       {
        cout <<"X[k] negative in Get_KE change Quad rule " <<  X[k] <<endl;
//         exit(0);
       }

      Mult = r_axial*weights[k]*AbsDetjk[k];

      KE_QP = 0.;
      u1_rise = 0.;
      u2_rise = 0.;

      for(l=0;l<N_BF;l++)
       {
        j = DOF[l];
        U1 = ValuesVX[j];
        U2 = ValuesVY[j];


 //       u.u
        KE_QP += (values[k][l]*U1 * values[k][l]*U1
	        + values[k][l]*U2 * values[k][l]*U2 );
        u1_rise += values[k][l]*U1;
        u2_rise += values[k][l]*U2;
       }

     KE += (KE_QP*Mult);
     volume += Mult;
     x_mass   += (X[k]*Mult);
     y_mass   += (Y[k]*Mult);
     U1_Rise  += (u1_rise * Mult);
     U2_Rise  += (u2_rise * Mult);

    }  // for(k=0;k<N_QF
   } // endfor i

    KE      /=volume;
    x_mass  /=volume;
    y_mass  /=volume;
    U1_Rise /=volume;
    U2_Rise /=volume;

 parameters[0] = volume;
 parameters[1] = KE;
 parameters[2] = x_mass;
 parameters[3] = y_mass;
 parameters[4] = U1_Rise;
 parameters[5] = U2_Rise;

//   cout<< " Volume_1: "<< volume_1<< " Volume_2: "<< volume_2
//       << " Volume: "<< volume<< endl;
//   cout<< " x_mass: "<< x_mass << " y_mass: "<< y_mass<< endl;
//   cout<< " KE_1: "<< KE_1<< " KE: "<< KE<< endl;
//   cout<< " U1_Rise: "<< U1_Rise << " U2_Rise: "<< U2_Rise<< endl;
//   exit(0);
 }



void MoveGrid_imping(double **Entries, double *Sol, double *d, double *Rhs,
                  int *KCol, int *RowPtr,
                  TFEVectFunct2D *GridPos,
                  TFEVectFunct2D *Velocity, double dt,
		  TFEVectFunct2D *GridVelocity,
                  TFEVectFunct2D *NewGridPos, 
                  TVertex ***MovBoundVert, int *N_MovVert,
                  TBaseCell **Free_Cells, int **IsoCellEdgeNos,
                  TFEFunction2D *Surfactant, TFEFunction2D *SurfSurfactant,
                  bool &reparam, int &N_ReParam)
{
  int i,j,k,l,m,comp, N;
  int *VeloBeginIndex, *VeloGlobalNumbers;
  int *GridBeginIndex, *GridGlobalNumbers;
  int N_Cells, N_Vertices, N_Edges, N_LocalDOFs;
  int N_Levels, *DOF, *JointDOF, GridLength, polydegree;
  int N_BoundaryNodes, N_LinePoints, IIso, N_Inner, N_;

  double xi[4], eta[4], X[4], Y[4], VX[4], VY[4];
  double FunctValues[4][MaxN_BaseFunctions2D];
  double FEValuesX[MaxN_BaseFunctions2D];
  double FEValuesY[MaxN_BaseFunctions2D];
  double *ValuesX, *ValuesY;
  double *ValuesVX, *ValuesVY;
  double *NewValuesX, *NewValuesY;
  double s, t, x, y, IsoX, IsoY, Ay;
  double x0, x1, y0, y1, h_tot, res, oldres;
  double *gridvelo, *Nx, *Ny, *LineWeights, *zeta;
  double normalx, normaly, tangenx, tangeny, nx, ny, tx, ty;
  double un, hE, t0,t1, temp2, x_max, x_min, temp, eps=1e-6;  
  
  TIsoBoundEdge *isojoint;
  TMGLevel2D *Level;  
  TFESpace2D *VelocitySpace, *GridSpace;
  TCollection *Coll;
  TBaseCell *cell;
  BoundTypes bdtype;
  TBoundEdge *BoundEdge;
  TBoundComp2D *BoundComp;
  FE2D FEId;
  TFE2D *Element;
  TFEDesc2D *FEDesc;
  BaseFunct2D BF;
  TBaseFunct2D *bf;
  boolean OnBoundary;
  TJoint *joint;
  TVertex **Vertices;
  TQuadFormula1D *qf1;
  QuadFormula1D LineQuadFormula;
  BF2DRefElements RefElement;
  TRefTrans2D *F_K;
  RefTrans2D RefTrans;  
  QuadFormula2D QuadFormula;
  BoundCond Cond0, Cond1;  
 
  VelocitySpace = Velocity->GetFESpace2D();
  VeloBeginIndex = VelocitySpace->GetBeginIndex();
  VeloGlobalNumbers = VelocitySpace->GetGlobalNumbers();
  ValuesVX = Velocity->GetValues();
  ValuesVY = ValuesVX + Velocity->GetLength();

  GridPos->GridToData();
  GridSpace = GridPos->GetFESpace2D();
  GridBeginIndex = GridSpace->GetBeginIndex();
  GridGlobalNumbers = GridSpace->GetGlobalNumbers();
  GridLength = GridPos->GetLength();
  ValuesX = GridPos->GetValues();
  ValuesY = ValuesX + GridLength;

  N_Inner = GridSpace->GetN_Inner();
  N_BoundaryNodes = GridLength - GridSpace->GetN_Inner();
  // cout << "N_BoundaryNodes: " << N_BoundaryNodes << endl;

//   d = new double[ 2*GridLength];
  
  if(TDatabase::ParamDB->P5 > 0)
  {  
   Nx = new double[2*N_BoundaryNodes]; // additional values for edge midpoints
   Ny = new double[2*N_BoundaryNodes]; // additional values for edge midpoints
   memset(Nx, 0, 2*N_BoundaryNodes*SizeOfDouble);
   memset(Ny, 0, 2*N_BoundaryNodes*SizeOfDouble);
  }
  //cout << "N_BoundaryNodes: " << N_BoundaryNodes << endl;
 // cout << GridLength << " " << N_DOF << endl;

  NewValuesX = NewGridPos->GetValues();
  NewValuesY = NewValuesX + GridLength;

  memcpy(NewValuesX, ValuesX, 2*GridLength*SizeOfDouble);

  Coll = VelocitySpace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  // determine outer normal vectors
  IIso = N_BoundaryNodes;

    // Outward normal no need if we move boundary with velocity
  if(TDatabase::ParamDB->P5 > 0)
  {
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_Edges = cell->GetN_Edges();

    for(j=0;j<N_Edges;j++)
    {
      if( !(cell->GetJoint(j)->InnerJoint()) )
      {
        FEId = VelocitySpace->GetFE2D(i, cell);
        l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
        LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
        qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
        qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);

        RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
        switch(RefElement)
        {
          case BFUnitTriangle:
            RefTrans = TriaIsoparametric;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
            QuadFormula = TFEDatabase2D::GetQFTriaFromDegree(2*polydegree-1);
//             ((TTriaIsoparametric *)F_K)->SetApproximationOrder(polydegree);
            ((TTriaIsoparametric *)F_K)->SetQuadFormula(QuadFormula);
            ((TTriaIsoparametric *)F_K)->SetCell(cell);
          break;

          case BFUnitSquare:
            RefTrans = QuadIsoparametric;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
            QuadFormula = TFEDatabase2D::GetQFQuadFromDegree(3*polydegree);
//             ((TQuadIsoparametric *)F_K)->SetApproximationOrder(polydegree);
            ((TQuadIsoparametric *)F_K)->SetQuadFormula(QuadFormula);
            ((TQuadIsoparametric *)F_K)->SetCell(cell);
          break;

          default:
            Error("only triangles and quadrilaterals are allowes" << endl);
            Error("file: " << __FILE__ << " line " << __LINE__ << endl);
            exit(-1);
        } // endswitch

        normalx = 0;
        normaly = 0;
        hE = 0;
        for(k=0;k<N_LinePoints;k++)
        {
          F_K->GetOuterNormal(j, zeta[k], nx, ny);
          F_K->GetTangent(j, zeta[k], tx, ty);
          t = sqrt(tx*tx+ty*ty);
          normalx += t * LineWeights[k] * nx;
          normaly += t * LineWeights[k] * ny;
          hE += t * LineWeights[k];
          // cout << "k= " << k << " " << nx << " " << ny << endl;
        } // endfor k

        DOF = GridGlobalNumbers + GridBeginIndex[i];
/*
        switch(N_Edges)
        {
          case 3:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[1] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[2] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;

              case 2:
                l = DOF[2] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[0] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;
            } // endswitch j
          break;

          case 4:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[1] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[3] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;

              case 2:
                l = DOF[3] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[2] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;

              case 3:
                l = DOF[2] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[0] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;
            } // endswitch j
          break;
        } // endswitch N_Edges
*/

// /*
        switch(N_Edges)
        {
          case 3:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 2:
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;
            } // endswitch j
          break;

          case 4:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[3] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 2:
                l = DOF[3] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 3:
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;
            } // endswitch j
          break;
        } // endswitch N_Edges
// */

        if(cell->GetJoint(j)->GetType() == IsoBoundEdge)
        {
          FEId = VelocitySpace->GetFE2D(i, cell);
          FEDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
          N_LocalDOFs = FEDesc->GetN_JointDOF();
          t = 2.0/(N_LocalDOFs-1);
          for(k=1;k<N_LocalDOFs-1;k++)
          {
            /*
            Nx[IIso] += normalx;
            Ny[IIso] += normaly;
            */
            // /*
            s = -1.0 + k*t;
            F_K->GetOuterNormal(j, s, nx, ny);
            Nx[IIso] += nx;
            Ny[IIso] += ny;
            // */
            IIso++;
          } // endfor
        } // endif
      } // !InnerJoint
    } // endfor j
  } // endfor i

  N_ = IIso;
  // normalize normal vector
  for(i=0;i<N_;i++)
  {
    x = Nx[i];
    y = Ny[i];
    t = sqrt(x*x+y*y);
    Nx[i] /= t;
    Ny[i] /= t;

    // cout << setw(5) << i << "n = (" << Nx[i] << ", " << Ny[i] << ")" << endl;
  }
}

  // determine new position of boundary vertices
  for(i=0;i<N_Cells;i++)
  {
    cell  = Coll->GetCell(i);
    N_Edges = cell->GetN_Edges();
    OnBoundary = FALSE;
    for(j=0;j<N_Edges;j++)
    {
      if( !(cell->GetJoint(j)->InnerJoint()) )
       {
        OnBoundary = TRUE;
        joint = cell->GetJoint(j);
       }
    } // endfor j

    if(OnBoundary)
    {
      FEId = VelocitySpace->GetFE2D(i, cell);
      Element = TFEDatabase2D::GetFE2D(FEId);
      bf = Element->GetBaseFunct2D();
      N_LocalDOFs = Element->GetN_DOF();

      switch(N_Edges)
      {
        case 3:
          xi[0]  = 0; xi[1]  = 1; xi[2]  = 0;
          eta[0] = 0; eta[1] = 0; eta[2] = 1;
          VX[0] = VX[1] = VX[2] = 0;
          VY[0] = VY[1] = VY[2] = 0;
          N_Vertices = N_Edges;
        break;

        case 4:
          // NOTE: sorted due to number of DOFs !!!
          xi[0]  = -1; xi[1]  =  1; xi[2]  = -1; xi[3]  =  1;
          eta[0] = -1; eta[1] = -1; eta[2] =  1; eta[3] =  1;
          VX[0] = VX[1] = VX[2] = VX[3] = 0;
          VY[0] = VY[1] = VY[2] = VY[3] = 0;
          N_Vertices = N_Edges;
        break;

        default:
          Error("only triangles and quadrilaterals are allowed!" << endl);
          exit(-1);
      } // endswitch

      for(j=0;j<N_Vertices;j++)
        bf->GetDerivatives(D00, xi[j], eta[j], FunctValues[j]);

      DOF = VeloGlobalNumbers + VeloBeginIndex[i];

      for(j=0;j<N_LocalDOFs;j++)
      {
        k = DOF[j];
        s = ValuesVX[k];
        t = ValuesVY[k];

        for(l=0;l<N_Vertices;l++)
        {
          VX[l] += FunctValues[l][j]*s;
          VY[l] += FunctValues[l][j]*t;
        } // endfor l
      } // endfor j

      FEId = GridSpace->GetFE2D(i, cell);
      Element = TFEDatabase2D::GetFE2D(FEId);
      BF = Element->GetBaseFunct2D_ID();
      if( (BF != BF_C_T_P1_2D) && (BF != BF_C_Q_Q1_2D) )
      {
        Error("Grid Space must be conforming and of first order!" << endl);
        exit(-1);
      }  // endif
      bf = Element->GetBaseFunct2D();
      N_LocalDOFs = Element->GetN_DOF();


      DOF = GridGlobalNumbers + GridBeginIndex[i];

      for(j=0;j<N_Vertices;j++)
      {
        l = DOF[j];
        k = l - N_Inner;
        if(k>=0)
        {
          if((TDatabase::ParamDB->P5 > 0) && (ValuesY[l] != 0) )
          {
           un = VX[j]*Nx[k] + VY[j]*Ny[k];

           if(ValuesX[l] == 0 )
            { NewValuesX[l] = ValuesX[l]; }
           else
            { NewValuesX[l] = ValuesX[l] + dt*un*Nx[k]; }

            if(ValuesY[l] == 0) 
             { NewValuesY[l] = ValuesY[l];  }
            else    
            { NewValuesY[l] = ValuesY[l] + dt*un*Ny[k]; } 
            
//             NewValuesX[l] = ValuesX[l] + dt*un*Nx[k];
//             NewValuesY[l] = ValuesY[l] + dt*un*Ny[k];
          }
          else
          {
 
	    if(ValuesX[l] == 0 )
	    { NewValuesX[l] = ValuesX[l]; }
	    else
	    { NewValuesX[l] = ValuesX[l] + dt*VX[j]; }

            if(ValuesY[l] == 0) 
	     { NewValuesY[l] = ValuesY[l];  }
            else    
	     { NewValuesY[l] = ValuesY[l] + dt*VY[j];  }                
          }
        }
 //    Due to spline approximation solid boundary end vertices may take negative y value
         if(NewValuesY[l]<0.0 ) NewValuesY[l] = 0.0;
        if( fabs(NewValuesX[l]) < 1e-12 ) NewValuesX[l] = 0.0;
      } // endfor j
     } // endif
   } // endfor i

/*

   for(i=GridLength-N_BoundaryNodes;i<GridLength;i++)
 cout << i <<"  ---  " <<NewValuesX[i] << "  ---  " << NewValuesY[i]<<endl;

exit(0);
*/

   MovBoundVert[0][0]->GetCoords(x, Ay);

   NewGridPos->DataToGrid();

   MovBoundVert[2][0]->GetCoords(x, y); // right wetting point
   y=0.;
   h_tot = x;
   h_tot /= (double)N_MovVert[0];
   for(i=1;i<N_MovVert[0];i++)
     MovBoundVert[0][i]->SetCoords(h_tot*(double)i, y);


   // axial boundary
   MovBoundVert[1][0]->GetCoords(x, y);

   N=N_MovVert[1];      
   h_tot = (y-Ay)/(double)N;   
//    N--;
   
   for(i=1;i<N;i++)
    {
//      MovBoundVert[1][i]->GetCoords(x, Ay);
     y  = ((double)(N-i))*h_tot;
//      cout<< " y " << Ay <<" new y " << y <<endl;      

     MovBoundVert[1][i]->SetCoords(x, y);   
    }       
 
// exit(0);
  //======================================================================  
  //  Reparametrization of free surface - Begin
  //======================================================================     
   if(reparam)  
   {   
    // update the iso points and then reparmetrize the free surf
    IIso = 0;
    for(i=0;i<N_MovVert[2];i++)
    {
     cell = Free_Cells[i];
     j = IsoCellEdgeNos[1][i];
     joint = cell->GetJoint(j);
     isojoint = (TIsoBoundEdge *)joint;
     k = isojoint->GetN_Vertices();
     Vertices = isojoint->GetVertices();
     FEId = VelocitySpace->GetFE2D(IsoCellEdgeNos[0][i], cell);
     FEDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
     m = FEDesc->GetN_JointDOF();
     if(m == k+2)
      {
       JointDOF = FEDesc->GetJointDOF(j);
       DOF =  VeloGlobalNumbers+VeloBeginIndex[IsoCellEdgeNos[0][i]];
       for(l=0;l<k;l++)
        {
         m = DOF[JointDOF[l+1]];
         Vertices[l]->GetCoords(IsoX, IsoY);
         if(TDatabase::ParamDB->P5 > 0)
          {
           un = ValuesVX[m]*Nx[IIso+N_BoundaryNodes]
                + ValuesVY[m]*Ny[IIso+N_BoundaryNodes];
              IsoX += dt*un*Nx[IIso+N_BoundaryNodes];
              IsoY += dt*un*Ny[IIso+N_BoundaryNodes];
            }
            else
            {
             IsoX += dt*ValuesVX[m];
             IsoY += dt*ValuesVY[m];
            }
            
           if(IsoY<=0) IsoY = 1e-5;
//            if(fabs(x)<1e-12) x = 0.;
           
           Vertices[l]->SetCoords(IsoX, IsoY);
           IIso++;
          } // endfor l
        }
       else
        {
          // approximation order of isoparametric boundary and velocity
          // element must be the same
          Error("No match in isoparametric case" << endl);
          exit(-1);
        } //  if(m == k+2)     
    }// for(i=0;i<N_MovVert[2];i++)     
     
    ReParam_axial3D_U(N_MovVert[2], Free_Cells,  IsoCellEdgeNos[1], IsoCellEdgeNos[0],
                      Velocity, Surfactant, SurfSurfactant, TRUE);   

    OutPut("ReParam CURRENT TIME: ");
    OutPut(TDatabase::TimeDB->CURRENTTIME << endl);   
   }  //if (reparam)  
   
   NewGridPos->GridToData();
   GridPos->DataToGrid();

   memset(Rhs, 0, 2*GridLength*SizeOfDouble);
   memcpy(d, NewValuesX, 2*GridLength*SizeOfDouble);
   Daxpy(2*GridLength, -1, ValuesX, d);
   memcpy(Rhs + (GridLength-N_BoundaryNodes),
          d+(GridLength-N_BoundaryNodes),
          N_BoundaryNodes*SizeOfDouble);
   memcpy(Rhs + (2*GridLength-N_BoundaryNodes),
          d+(2*GridLength-N_BoundaryNodes),
          N_BoundaryNodes*SizeOfDouble);

  memset(Sol, 0 , 2*GridLength*SizeOfDouble);
  memcpy(Sol + (GridLength-N_BoundaryNodes),
         d+(GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);
  memcpy(Sol + (2*GridLength-N_BoundaryNodes),
         d+(2*GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);


  SolveGridEquation(Entries, Sol, Rhs, KCol, RowPtr, GridLength);

  memcpy(d, ValuesX, 2*GridLength*SizeOfDouble);
  Daxpy(2*GridLength, 1, Sol, d);
  memcpy(NewValuesX, d, (GridLength-N_BoundaryNodes)*SizeOfDouble);
  memcpy(NewValuesY, d+GridLength, (GridLength-N_BoundaryNodes)*SizeOfDouble);

  gridvelo = GridVelocity->GetValues();
  memcpy(gridvelo, Sol, 2*GridLength*SizeOfDouble);
  Dscal(2*GridLength, 1./dt, gridvelo);  

  
  
// for(i=GridLength-N_BoundaryNodes;i<GridLength;i++)
//  cout << i << "  ---  "<<NewValuesX[i] << "  ---  " << NewValuesX[i+GridLength] << endl;
//cout << "test " << endl;
  // put solution into grid position
  IIso = 0;
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    DOF = GridGlobalNumbers + GridBeginIndex[i];
    N_Vertices = cell->GetN_Vertices();

    switch(N_Vertices)
    {
      case 3:
        for(j=0;j<N_Vertices;j++)
        {
          k = DOF[j];
          cell->GetVertex(j)->SetCoords(NewValuesX[k], NewValuesY[k]); 
        }
      break;

      case 4:
        k = DOF[0];
        cell->GetVertex(0)->SetCoords(NewValuesX[k], NewValuesY[k]);
        k = DOF[1];
        cell->GetVertex(1)->SetCoords(NewValuesX[k], NewValuesY[k]);
        k = DOF[3];
        cell->GetVertex(2)->SetCoords(NewValuesX[k], NewValuesY[k]);
        k = DOF[2];
        cell->GetVertex(3)->SetCoords(NewValuesX[k], NewValuesY[k]);
      break;
    } // endswitch
 
  if (!reparam)  
   {  
    N_Edges = cell->GetN_Edges();

    for(j=0;j<N_Edges;j++)
    {
      joint = cell->GetJoint(j);
      if(joint->GetType() == IsoBoundEdge)
      {
        isojoint = (TIsoBoundEdge *)joint;
        k = isojoint->GetN_Vertices();
        Vertices = isojoint->GetVertices();
        FEId = VelocitySpace->GetFE2D(i, cell);
        FEDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
        m = FEDesc->GetN_JointDOF();
        if(m == k+2)
        {
          JointDOF = FEDesc->GetJointDOF(j);
          DOF =  VeloGlobalNumbers+VeloBeginIndex[i];
          for(l=0;l<k;l++)
          {
            m = DOF[JointDOF[l+1]];
            Vertices[l]->GetCoords(IsoX, IsoY);
            if(TDatabase::ParamDB->P5 > 0)
            {
              un = ValuesVX[m]*Nx[IIso+N_BoundaryNodes]
                  + ValuesVY[m]*Ny[IIso+N_BoundaryNodes];
              IsoX += dt*un*Nx[IIso+N_BoundaryNodes];
              IsoY += dt*un*Ny[IIso+N_BoundaryNodes];
              // cout << "U:   " << ValuesVX[m] << " " << ValuesVY[m] << endl;
              // cout << "N:   " << Nx[IIso+N_BoundaryNodes] << " "
              //                 << Ny[IIso+N_BoundaryNodes] << endl;
              // cout << "UNN: " << un*Nx[IIso+N_BoundaryNodes] << " "
              //                 << un*Ny[IIso+N_BoundaryNodes] << endl;
            }
            else
            {  
              IsoX += dt*ValuesVX[m];
              IsoY += dt*ValuesVY[m];
            }
           if(IsoY<=0) IsoY = 1e-5;
//            if(fabs(IsoX)<1e-12) IsoX = 0.;
           Vertices[l]->SetCoords(IsoX, IsoY);
            IIso++;
          } // endfor l
        }
       else
        {
          // approximation order of isoparametric boundary and velocity
          // element must be the same
          Error("No match in isoparametric case" << endl);
          exit(-1);
        }
      } // endif
    } // endfor j  
   } // if(reparam)
  } // endfor i

   if(reparam)
    {
     N_ReParam++; 
     reparam = FALSE;
    }    
    
//   delete [] d;
  
  if(TDatabase::ParamDB->P5 > 0)
  { delete [] Nx;  delete [] Ny;}
  
} // MoveGrid_imping



void MapNSESol(TFESpace2D *velocity_space_nse, double *sol_nse, int N_U_nse, TFESpace2D *pressure_space_nse, int N_P_nse,
	       TFESpace2D *velocity_space,  double *sol, int N_U, TFESpace2D *pressure_space, int N_P, int *NSE_GlobalCllNo)
{
  int i, j, k, l, N, N_Cells_NSE, n1, n2;
  int N_DOF, *UGlobalNumbers_NSE, *UBeginIndex_NSE, *PGlobalNumbers_NSE, *PBeginIndex_NSE;
  int *UGlobalNumbers, *UBeginIndex, *PGlobalNumbers, *PBeginIndex;
  int *DOF_NSE, *DOF;
  
  TCollection *Coll_NSE, *Coll;  
  TBaseCell *Me;
  FE2D FeId;
  TFEDesc2D *FeDesc;  
  
  
  Coll_NSE = velocity_space_nse->GetCollection();
  N_Cells_NSE = Coll_NSE->GetN_Cells();
  UGlobalNumbers_NSE = velocity_space_nse->GetGlobalNumbers();
  UBeginIndex_NSE = velocity_space_nse->GetBeginIndex(); 
  PGlobalNumbers_NSE = pressure_space_nse->GetGlobalNumbers();
  PBeginIndex_NSE = pressure_space_nse->GetBeginIndex(); 
  
//   Coll = velocity_space->GetCollection();
  UGlobalNumbers = velocity_space->GetGlobalNumbers();
  UBeginIndex = velocity_space->GetBeginIndex(); 
  PGlobalNumbers = pressure_space->GetGlobalNumbers();
  PBeginIndex = pressure_space->GetBeginIndex(); 
  
  memset(sol, 0, 2*N_U + N_P*SizeOfDouble);    
  for(i=0;i<N_Cells_NSE;i++)
   {
    N = NSE_GlobalCllNo[i];
    Me = Coll_NSE->GetCell(i); 
    
    //velocity
    FeId = velocity_space_nse->GetFE2D(i, Me);
    FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FeId);
    N_DOF = FeDesc->GetN_DOF(); 
    
    DOF_NSE = UGlobalNumbers_NSE + UBeginIndex_NSE[i];
    DOF = UGlobalNumbers + UBeginIndex[N];  
    
     for(j=0;j<N_DOF;j++)
      {
       n1 = DOF_NSE[j];
       n2 = DOF[j];
       sol[n2] = sol_nse[n1];
       sol[N_U + n2] = sol_nse[N_U_nse + n1];      
      }
  
    //pressure
    FeId = pressure_space_nse->GetFE2D(i, Me);
    FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FeId);
    N_DOF = FeDesc->GetN_DOF(); 
    DOF_NSE = PGlobalNumbers_NSE + PBeginIndex_NSE[i];
    DOF = PGlobalNumbers + PBeginIndex[N]; 
     for(j=0;j<N_DOF;j++)
      {
       n1 = DOF_NSE[j];
       n2 = DOF[j];
       sol[2*N_U + n2] = sol_nse[2*N_U_nse + n1];   
      }  
      
//      exit(0);  
   } 
//    exit(0);
}

void  GetMovingBoundData(TCollection *coll, int *N_MovVert, TBoundEdge *** &Bound_Joint, TVertex *** &MovBoundVert,
                         TIsoBoundEdge ** &Free_Joint, TBaseCell ** &Free_Cells, int ** &IsoCellEdgeNos, double x, double y)
{
 int i, j, k, l, m0, m1, m2, N_Cells, comp;
 int ORDER, VSP;
 
 double  TX[2], TY[2], x0, y0, x1, y1;
  
 TBaseCell *Me;
 TJoint *Joint;
 TBoundEdge *Solid_Joint;
 TBoundComp *BoundComp;  
 TVertex *temp_Mov;
 TBoundEdge *tempSlip_Joint;

  ORDER = 0;
  VSP = TDatabase::ParamDB->VELOCITY_SPACE;

  if (abs(VSP) > 20)
   {ORDER = abs(VSP) - 20;}
  else if ( abs(VSP) > 10)
    {ORDER = abs(VSP) - 10;}
  else ORDER = abs(VSP);
  
  
  N_Cells = coll->GetN_Cells();
  N_MovVert[0] = 0;
  N_MovVert[1] = 0;
  N_MovVert[2] = 0;
  
    for(j=0;j<N_Cells;j++)
      {
        Me = coll->GetCell(j);
        k = Me->GetN_Edges();
        for(l=0;l<k;l++)
        {
         Joint = Me->GetJoint(l);
          if(Joint->GetType() == BoundaryEdge)
            {
             Solid_Joint = (TBoundEdge *)Joint;
             BoundComp = Solid_Joint->GetBoundComp();
             comp=BoundComp->GetID();
             
             if(comp==0)
              { N_MovVert[0]++; }  
             else if(comp==2)
              {N_MovVert[1]++; }
             else
              {
               cout<<"Error comp " <<endl;
               exit(0);      
              }                            
            }
           else if(Joint->GetType() == IsoBoundEdge)
            { N_MovVert[2]++;  }
            
          }// endfor l
        }// endfor j
  
//      for(i=0; i<3; i++)
//       cout<<"BDComp " << i << " N_Vert " << N_MovVert[i] << endl; 
     

     // solid bound
     Bound_Joint[0] = new TBoundEdge* [N_MovVert[0]];
     Bound_Joint[1] = new TBoundEdge* [N_MovVert[1]];

     MovBoundVert[0] = new TVertex* [N_MovVert[0]];
     MovBoundVert[1] = new TVertex*[N_MovVert[1]];

     // free bound
     Free_Joint = new TIsoBoundEdge*[N_MovVert[2]];
     MovBoundVert[2] = new TVertex*[N_MovVert[2]];
     Free_Cells = new TBaseCell*[N_MovVert[2]];
     IsoCellEdgeNos[0]  = new int [N_MovVert[2]];
     IsoCellEdgeNos[1]  = new int [N_MovVert[2]];
     
      m0 = 0;
      m1 = 0;
      m2 = 0;  
  
     for(j=0;j<N_Cells;j++)
      {
        Me = coll->GetCell(j);
        k = Me->GetN_Edges();
        for(l=0;l<k;l++)
        {
          Joint = Me->GetJoint(l);
          if(Joint->GetType() == BoundaryEdge)
          {
           Solid_Joint = (TBoundEdge *)Joint;
           BoundComp = Solid_Joint->GetBoundComp();
           comp=BoundComp->GetID();
           if(comp==0)
           {
             Bound_Joint[0][m0] = (TBoundEdge *)Joint;
             MovBoundVert[0][m0] = Me->GetVertex(l);
             m0++;
            }
         else if(comp==2)
            {
             Bound_Joint[1][m1] = (TBoundEdge *)Joint;
             MovBoundVert[1][m1] = Me->GetVertex(l);
             m1++;
            }
         
           }
         else if(Joint->GetType() == IsoBoundEdge)
           {
            Free_Joint[m2] = (TIsoBoundEdge *)Joint;
            MovBoundVert[2][m2] = Me->GetVertex(l);
            Free_Cells[m2] = Me;
            IsoCellEdgeNos[0][m2] = j;
            IsoCellEdgeNos[1][m2] = l;

            Me->GetVertex(l)->GetCoords(TX[0], TY[0]);
            Me->GetVertex((l+1) % k)->GetCoords(TX[1], TY[1]);
#ifdef __SURFACTDIFFTEST__  
            Free_Joint[m2]->GenerateVertices(ORDER-1);
#else    
            Free_Joint[m2]->GeneratemidVert(ORDER-1, TX, TY);
#endif    
            m2++;
           }
          }// endfor l
         }// endfor j    
         
         
  // sort         
    for(k=0;k<N_MovVert[0]-1;k++)
      {
      for(l=k+1;l<N_MovVert[0];l++)
       {MovBoundVert[0][k]->GetCoords(x0, y0);
	MovBoundVert[0][l]->GetCoords(x1, y1);
	if(x0 > x1)
	{
	  temp_Mov = MovBoundVert[0][k];
          MovBoundVert[0][k] = MovBoundVert[0][l];
          MovBoundVert[0][l] = temp_Mov;

	  tempSlip_Joint = Bound_Joint[0][k];
	  Bound_Joint[0][k] = Bound_Joint[0][l];
	  Bound_Joint[0][l] = tempSlip_Joint;
	 }
        }
       }         
//  for (k=0;k<N_MovVert[0];k++)
//       {
//        MovBoundVert[0][k]->GetCoords(x0, y0);
//        cout<< " x0 " << x0<<" y0 " << y0<<endl;
//        }

    for(k=0;k<N_MovVert[1]-1;k++)
     {
      for(l=k+1;l<N_MovVert[1];l++)
       {MovBoundVert[1][k]->GetCoords(x0, y0);
	MovBoundVert[1][l]->GetCoords(x1, y1);
	if(y1 > y0)
	{
	  temp_Mov = MovBoundVert[1][k];
          MovBoundVert[1][k] = MovBoundVert[1][l];
          MovBoundVert[1][l] = temp_Mov;

	  tempSlip_Joint = Bound_Joint[1][k];
	  Bound_Joint[1][k] = Bound_Joint[1][l];
	  Bound_Joint[1][l] = tempSlip_Joint;
	 }
        }
       }
//  for (k=0;k<N_MovVert[1];k++)
//       {
//        MovBoundVert[1][k]->GetCoords(x0, y0);
//        cout<< " x " << x0<<" y " << y0<<endl;
//        }

   Sort_Imping(Free_Cells, MovBoundVert[2], IsoCellEdgeNos[0], IsoCellEdgeNos[1], N_MovVert[2], x, y);

}// GetMovingBoundData

void  Domain2DSurf_2Phase(TCollection *Coll, TDomain *SurfDomain, int **N_List)
  {
  int i, j, k, l, m, n,  m1, N_SurfCells, N_Cells, N_SurfVert, ID, N;
  int maxEpV, a, b, len1, len2, Neighb_tmp;
  int *Bd_Part, *Lines, N_G, *PointNeighb, CurrNeib,  Neib[2];
  int *Cell_No, *Joint_No, comp, N_SolidEdges=0;

  double x1, y1;
  const int *TmpEV;
  TBaseCell *Me;
  TVertex **SurfVetrex;
  TJoint *Joint;
  TBaseCell  **SurfCellTree, *NullCell = NULL;
  boolean StartBD=FALSE, EndBD=FALSE;


  Cell_No = N_List[0];
  Joint_No = N_List[1];

  N_SurfCells = 0;
  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {
   Me = Coll->GetCell(i);
//    PhaseID = Me->GetPhase_ID();

//    if(PhaseID == 1) // outer phase cell is better for soluble surfactant coupling
   {
    k = Me->GetN_Edges();
    for(l=0;l<k;l++)
     {
       Joint = Me->GetJoint(l);
       
       //liquid-solid interface also for surfactants
       comp = -1;
#ifndef __SURFACTCONVTEST__         
#ifndef __SURFACTDIFFTEST__         
        if(Joint->GetType() == BoundaryEdge)
         comp=(((TBoundEdge *)Joint)->GetBoundComp())->GetID();
#endif 
#endif
       if(Joint->GetType() == IsoBoundEdge   || Joint->GetType() == IsoInterfaceJoint
           || Joint->GetType() == InterfaceJoint || comp==0  )
         N_SurfCells++;
     } // endfor l
   } //  if(PhaseID == 1)
  }// endfor i

//   OutPut("N_SurfCells: " << N_SurfCells << endl);
//   exit(0);
  
  Bd_Part = new int[N_SurfCells];
  Lines = new int[2*N_SurfCells];
  SurfVetrex = new TVertex*[2*N_SurfCells]; // maximum possible vertex

//   m = 0;
  m1 = 0;
  N    = 0; 
  N_SurfVert = 0;
  for(i=0;i<N_Cells;i++)
  {
   Me = Coll->GetCell(i);
//    PhaseID = Me->GetPhase_ID();

//    if(PhaseID == 1) // outer phase cell is better for soluble surfactant coupling
   {
    k = Me->GetN_Edges();
    for(l=0;l<k;l++)
     {
       Joint = Me->GetJoint(l);
       
       comp = -1;
#ifndef __SURFACTCONVTEST__    
#ifndef __SURFACTDIFFTEST__         
        if(Joint->GetType() == BoundaryEdge)
         comp=(((TBoundEdge *)Joint)->GetBoundComp())->GetID();       
#endif      
#endif
       
       if(Joint->GetType() == IsoBoundEdge   || Joint->GetType() == IsoInterfaceJoint
           || Joint->GetType() == InterfaceJoint || comp==0 )
        {

	  
          Cell_No[N] = i;
          Joint_No[N] = l;
          Bd_Part[N++] =(((TBoundEdge *)Joint)->GetBoundComp())->GetID();
          Me->GetShapeDesc()->GetEdgeVertex(TmpEV);
	  
	  if(comp==0)
	  {N_SolidEdges++;}
	  
	  
           for(n=0;n<2;n++) // In 2D each edge contains only two vertex
            {
              ID = 0;
              for(m=0;m<N_SurfVert;m++)
              if( SurfVetrex[m]==Me->GetVertex(TmpEV[2*l + n ]) )
               {
//             vertex is already in the collection
                ID = 1;
                Lines[m1] = m;
                m1++;
                break; 
              }

              if(ID==0)
               {
                 Lines[m1] = N_SurfVert;
                 m1++;
                 SurfVetrex[N_SurfVert] = Me->GetVertex(TmpEV[2*l + n ]);
                 N_SurfVert++;
		 
// 	        if(comp==0)
// 	        {		 
// 	         SurfVetrex[N_SurfVert-1]->GetCoords(x1, y1);	
// 		 cout << "x : " << x1 << " y1: " << y1 <<endl;
// 		}
               }

           } //  for(n=0
        } // if(Joint->GetType()
     } // endfor l
   } // if(PhaseID == 1)
  }// endfor i

   OutPut("N_SurfVert: " << N_SurfVert << endl);
   OutPut("N_SolidEdges: " << N_SolidEdges << endl);
//    exit(0);
   
   SurfCellTree = new TBaseCell*[N_SurfCells];

   for (i=0;i<N_SurfCells;i++)
   {
    SurfCellTree[i] = new TMacroCell(TDatabase::RefDescDB[S_Line], 0);
    SurfCellTree[i]->SetVertex(0, SurfVetrex[ Lines[ 2*i       ]]);
    SurfCellTree[i]->SetVertex(1, SurfVetrex[ Lines[ 2*i + 1]]);

    ((TMacroCell *) SurfCellTree[i])->SetSubGridID(0);
    ((TBaseCell *) SurfCellTree[i])->SetBd_Part(Bd_Part[i] );
   }

   SurfDomain->SetTreeInfo(SurfCellTree, N_SurfCells);

   TDatabase::IteratorDB[It_EQ]->SetParam(SurfDomain);
   TDatabase::IteratorDB[It_LE]->SetParam(SurfDomain);
   TDatabase::IteratorDB[It_Finest]->SetParam(SurfDomain);
   TDatabase::IteratorDB[It_Between]->SetParam(SurfDomain);
   TDatabase::IteratorDB[It_OCAF]->SetParam(SurfDomain);

  // search neighbours
   maxEpV = 2; // max possible degree of each vertex
   N_G = N_SurfVert;
   maxEpV++;
   PointNeighb = new int[maxEpV * N_G];
   memset(PointNeighb, 0, maxEpV * N_G *SizeOfInt);

    // every row contains "maxEpV" columns
    // for every row at first colomn contains the number of cells containing this vertex
    // at further columns we set the index of the corresponding cells containing this vertex
   for(i=0;i<2*N_SurfCells;i++)
   {
     j = Lines[i]*maxEpV;
     PointNeighb[j]++;
     PointNeighb[j + PointNeighb[j]] = i / 2;
   }

  N_G = N_SurfCells;

  for(i=0;i<N_G;i++)
   {
    StartBD=FALSE;
    EndBD=FALSE;
    a =  Lines[2*i ];
    b =  Lines[2*i + 1];
    Neib[0] = -1;
    Neib[1] = -1;
    CurrNeib = 0;
    len1 = PointNeighb[a*maxEpV];
    len2 = PointNeighb[b*maxEpV];


    if(len2==1)
     {
      Neib[0] = PointNeighb[b*maxEpV + 1];
      Joint = new TJointEqN(SurfCellTree[Neib[0]]);
      SurfCellTree[Neib[0]]->SetJoint(1, Joint);
      cout<<" End cell " << Neib[0] <<" Neib[1] " << Neib[1] <<endl;
     }

  // find indexes of cells containing the current vertex (joint)
    if(len1==1)
      Neib[0] = PointNeighb[a*maxEpV + 1];
    else
      for(j=1;j<=len1;j++)
       {
        Neighb_tmp = PointNeighb[a*maxEpV + j];
        for (k=1;k<=len2;k++)
         {
          if(Neighb_tmp == PointNeighb[b*maxEpV + k])
           {
            Neib[0] = Neighb_tmp;

            if(j==1)
             Neib[1] = PointNeighb[a*maxEpV + j + 1];
            else if(j==2)
             Neib[1] = PointNeighb[a*maxEpV + j - 1];

            break;
           }
         }
       } //  for (j=1;j<=len1


    if(len1==2)
     {
      Joint = new TJointEqN(SurfCellTree[Neib[0]], SurfCellTree[Neib[1]]);

      SurfCellTree[Neib[0]]->SetJoint(0, Joint);
      SurfCellTree[Neib[1]]->SetJoint(1, Joint);
     }
    else if(len1==1)
     {
      Joint = new TJointEqN(SurfCellTree[Neib[0]]);
      SurfCellTree[Neib[0]]->SetJoint(0, Joint);
     cout<< " Start cell " << Neib[0] <<" Neib[1] " << Neib[1] <<endl;
     }

   } //  for(i=0;i<N_G;i++)

  delete [] SurfVetrex;
  delete [] PointNeighb;
  delete []  Bd_Part;
  delete []  Lines;
  
//   cout << "Domain2DSurf_2Phase done !" << endl;
//   exit(0);
 }

void MapDomainToSurf(TFEFunction2D *Fe2D, TFEFunction1D *Fe1D,
                     int *Cell_array, int *Joint_array)
{
  int i,j,k,l,n1,n2, N_Cells1D, N_Cells2D, N_DOF, N_DOF1D, N;
  int *GlobalNumbers1D, *GlobalNumbers2D,  *BeginIndex1D, *BeginIndex2D, *JointDOF;
  double  *Values1D, *Values2D;
  TCollection *Coll1D, *Coll2D;
  TFESpace1D *SurfSpace;
  TFESpace2D *FeSpace;
  TBaseCell *Me;
  TJoint *joint;
  FE2D FeId;
  TFEDesc2D *FeDesc;
  FE1D FeId1D;
  TFE1D *Element;

  SurfSpace = Fe1D->GetFESpace1D();
  FeSpace = Fe2D->GetFESpace2D();
  Values1D= Fe1D->GetValues();
  Values2D= Fe2D->GetValues();

  GlobalNumbers1D = SurfSpace->GetGlobalNumbers();
  BeginIndex1D = SurfSpace->GetBeginIndex();

  GlobalNumbers2D = FeSpace->GetGlobalNumbers();
  BeginIndex2D = FeSpace->GetBeginIndex();


  Coll1D = SurfSpace->GetCollection();
  Coll2D = FeSpace->GetCollection();
  N_Cells1D = Coll1D->GetN_Cells();
  N_Cells2D = Coll2D->GetN_Cells();

//   cout <<N_Cells1D<< "  test MapSurfToDomain " << N_Cells2D << endl;
  for(i=0;i<N_Cells1D;i++)
   {
     N = Cell_array[i];
     Me = Coll2D->GetCell(N);
     l = Joint_array[i];
     joint = Me->GetJoint(l);
     FeId = FeSpace->GetFE2D(N, Me);
     FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FeId);
     JointDOF = FeDesc->GetJointDOF(l);
     N_DOF = FeDesc->GetN_JointDOF();

     FeId1D = SurfSpace->GetFE1D(i, Coll1D->GetCell(i));
     Element = TFEDatabase2D::GetFE1D(FeId1D);
     N_DOF1D = Element->GetN_DOF();

     if(N_DOF1D !=N_DOF )
      {
        cout << "Number of degrees of freedon on surface line should be matched with triangular edge" << endl;
        cout << " Check MapDomainToSurf !!!!" << endl;
        cout << "N_DOF1D " << N_DOF1D << " N_DOF2D " << N_DOF << endl;
        exit(0);
      }

     for(j=0;j<N_DOF1D;j++)
      {
       n1 = GlobalNumbers1D[BeginIndex1D[i] + j];
       n2 = GlobalNumbers2D[BeginIndex2D[N]+JointDOF[j]];
       Values1D[n1] = Values2D[n2];
      }// endfor j
   } //  for(i=0;i<N_Cel

//   cout <<N_Cells1D<< "  2 test MapSurfToDomain " << N_Cells2D << endl;

}


 void MapSurfToDomain(TFEFunction1D *Fe1D, TFEFunction2D *Fe2D,
                     int *Cell_array, int *Joint_array)
{
  int i,j,k,l,n1,n2, N_Cells1D, N_Cells2D, N_DOF, N_DOF1D, N;
  int *GlobalNumbers1D, *GlobalNumbers2D,  *BeginIndex1D, *BeginIndex2D, *JointDOF, N_G;
  double  *Values1D, *Values2D;
  TCollection *Coll1D, *Coll2D;
  TFESpace1D *SurfSpace;
  TFESpace2D *FeSpace;
  TBaseCell *Me;
  TJoint *joint;
  FE2D FeId;
  TFEDesc2D *FeDesc;
  FE1D FeId1D;
  TFE1D *Element;

  SurfSpace = Fe1D->GetFESpace1D();
  FeSpace = Fe2D->GetFESpace2D();
  Values1D= Fe1D->GetValues();
  Values2D= Fe2D->GetValues();
  N_G = Fe2D->GetLength();
  
  //set it to zero
  memset(Values2D, 0, N_G *SizeOfInt);  
  
  GlobalNumbers1D = SurfSpace->GetGlobalNumbers();
  BeginIndex1D = SurfSpace->GetBeginIndex();

  GlobalNumbers2D = FeSpace->GetGlobalNumbers();
  BeginIndex2D = FeSpace->GetBeginIndex();

  Coll1D = SurfSpace->GetCollection();
  Coll2D = FeSpace->GetCollection();
  N_Cells1D = Coll1D->GetN_Cells();
  N_Cells2D = Coll2D->GetN_Cells();
  
//   cout <<N_Cells1D<< "  test MapSurfToDomain " << N_Cells2D << endl;
//   exit(0);
  for(i=0;i<N_Cells1D;i++)
   {
     N = Cell_array[i];
     Me = Coll2D->GetCell(N);
     l = Joint_array[i];
     joint = Me->GetJoint(l);
     FeId = FeSpace->GetFE2D(N, Me);
     FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FeId);
     JointDOF = FeDesc->GetJointDOF(l);
     N_DOF = FeDesc->GetN_JointDOF();

     FeId1D = SurfSpace->GetFE1D(i, Coll1D->GetCell(i));
     Element = TFEDatabase2D::GetFE1D(FeId1D);
     N_DOF1D = Element->GetN_DOF();

     if(N_DOF1D !=N_DOF )
      {
        cout << "Number of degrees of freedon on surface line should be matched with triangular edge" << endl;
        cout << " Check MapDomainToSurf !!!!" << endl;
        cout << "N_DOF1D " << N_DOF1D << " N_DOF2D " << N_DOF << endl;
        exit(0);
      }

     for(j=0;j<N_DOF1D;j++)
      {
       n1 = GlobalNumbers1D[BeginIndex1D[i] + j];
       n2 = GlobalNumbers2D[BeginIndex2D[N]+JointDOF[j]];
       Values2D[n2] = Values1D[n1];
      }// endfor j


   } //  for(i=0;i<N_Cel

//   cout <<N_Cells1D<< "  test MapSurfToDomain " << N_Cells2D << endl;
//   exit(0);   
}



void PrintSurfSurfactant(int N, TVertex **Vertex, TFEFunction2D *Surfact, int &N_BData, char *Surf_Char)
{
 int i, j, k, Cell_No, IJoint, N_DOF_Local;
 double  T_val[3], x1, y1, x2, y2, ArcLength=0.;
 char *VtkBaseName;
 VtkBaseName = TDatabase::ParamDB->VTKBASENAME;

 std::ostringstream os;
 os << " ";

 os.seekp(std::ios::beg);
//       os << "surfact"<< i << ".dat" << ends;
  if(N_BData<10) os << "BDData/"<<VtkBaseName<<Surf_Char<<"_0000"<<N_BData<<".data" << ends;
  else if(N_BData<100) os <<"BDData/"<<VtkBaseName<<Surf_Char<<"_000"<<N_BData<<".data" << ends;
  else if(N_BData<1000) os <<"BDData/"<<VtkBaseName<<Surf_Char<<"_00"<<N_BData<<".data" << ends;
  else if(N_BData<10000) os <<"BDData/"<<VtkBaseName<<Surf_Char<<"_0"<<N_BData<<".data" << ends;
  else  os <<"BDData/"<<VtkBaseName<<Surf_Char<<"Gamma_"<<N_BData<<".data" << ends;

  std::ofstream dat(os.str().c_str());

  if (!dat)
   {
    cerr << "cannot open file for output" << endl;
    exit(0);
   }
  dat << "%% Surfactant data created by MooNMD" << endl;
  dat << "%% Current Reference Time :" << TDatabase::TimeDB->CURRENTTIME << endl;
  dat << "%% x, y, ArcLength,  Surfact" << endl;

  Vertex[0]->GetCoords(x1, y1);
  for(i=0;i<N;i++)
   {
    Vertex[i]->GetCoords(x2, y2);
//    cout<<i<< " Angle of free Vertices "<<(180./Pi)*atan2(y,x)<<endl;

    ArcLength += sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) );

    Surfact->FindGradient( x1,  y1, T_val);
    dat << x1<< " " << y1 << " "<< ArcLength << "  " << T_val[0] <<endl;

    x1=x2; y1=y2;
   }

      dat.close();
      cout << endl;
      cout << "Surfactant data wrote into file " << endl;
 N_BData++;

}


void RemeshAxial3D_ImpDrop(TDomain * &Domain, TFESpace2D ** &FESpaces_All,
                           TFEVectFunct2D ** &FEVectFuncts_All, TFEFunction2D ** &FEFunctions_All, int *N_MovVert, 
                           TBoundEdge *** &Bound_Joint, TVertex *** &MovBoundVert, TIsoBoundEdge ** &Free_Joint,
                           TBaseCell ** &Free_Cells, int ** &IsoCellEdgeNos, double ** &Sol_All,  double ** &Rhs_All, 
                           TSquareStructure2D ** &SquareStructure_All, TStructure2D ** &Structure_All,
                           TSquareMatrix2D ** &SqMat_All, TMatrix2D ** &Mat_All,
                           TDomain *SurfDomain, TFESpace1D **IFaceFeSpaces, int **N_List,
                           TFEFunction1D **IFaceFeFunct, TSquareMatrix1D **SqMat_IFace,
                           TSquareStructure1D **IFaceStruct, FE1D *FE1D_List
                           )
{
  int i, j, k, l, N_G, N_Cells, ORDER, VSP, N_DOF, N_ThermalDOF;
  int In_Index, CurrComp, Old_N_Cells, Old_N_RootCells, CurrVertex, N_Joints, N_Vertices, ID;
  int N_RootCells, *Triangles, *PointNeighb, maxEpV = 0, a, b, Neighb_tmp, Neib[2];
  int CurrNeib, len1, len2, pressure_space_code, N_U, N_P, N_Unknowns, comp;
  int *JointDOF, *DOF, *GlobalNumbers, *BeginIndex, N_refX, N, N_surfactDOF, N_IsurfactDOF, N_S;
  
  double d, h, t, tx, ty, x1, x2, y1, y2, Lx, Ly, Rx, Ry, *S_BX, *S_BY, *A_BX, *A_BY, refX;
  double area, *Coordinates, left, right, bottom, top, T_a, T_b, *sol, *ValuesU2, *Sx, *Sy;
  double TX[2], TY[2], temp0, temp1;
  
  
  TBoundPart *BoundPart;
  TBoundComp *BoundComp;
  TBdLine *UpdateSlipBound, *UpdateAxialBound;
  TCollection *coll, *Old_Coll, *mortarcoll = NULL;
  TBaseCell *cell, **CellTree, **Old_CellTree;
  TVertex **VertexDel, **NewVertices;
  TJoint *Joint;
  TFESpace2D *velocity_space, *pressure_space, *thermal_space;
  TFEVectFunct2D *u; 
  TFEFunction2D *p; 
  TBoundEdge *Solid_Joint;
  FE2D FeId;
  TFEDesc2D *FeDesc; 
  TBaseCell *Me;
  TVertex **IsoVertices;
  
  // strings
  char ReadinDat[] = "readin.dat";
  char TString[] = "T";
  char NameString[]  = "name";
  char UString[] = "u";
  char PString[] = "p";
  char WString[] = "w";
  char CString[] = "C";
  char IFaceSString[] = "C_I";

    
  std::ostringstream opts;
  std::ostringstream os;
  os << " ";
  opts << " ";

  struct triangulateio In, Out;

#ifdef __SURFACT__   
  TCollection *SOld_Coll, *IFace_Coll; 
  TBaseCell **Surf_OldCellTree;
  TFEFunction2D *New_c, *SurfNew_c; 
  TFESpace2D *SurfSurfactSpace, *SurfactSpace;
     
  int N_IFaceCells, Old_S_N_RootCells;
  
  double *csol, *Surfcsol;
  
  TDatabase::IteratorDB[It_EQ]->SetParam(SurfDomain);
  TDatabase::IteratorDB[It_LE]->SetParam(SurfDomain);
  TDatabase::IteratorDB[It_Finest]->SetParam(SurfDomain);
  TDatabase::IteratorDB[It_Between]->SetParam(SurfDomain);
  TDatabase::IteratorDB[It_OCAF]->SetParam(SurfDomain);

  SOld_Coll = SurfDomain->GetCollection(It_Finest, 0);
  SurfDomain->GetTreeInfo(Surf_OldCellTree, Old_S_N_RootCells);

  // initialize iterators
  TDatabase::IteratorDB[It_EQ]->SetParam(Domain);
  TDatabase::IteratorDB[It_LE]->SetParam(Domain);
  TDatabase::IteratorDB[It_Finest]->SetParam(Domain);
  TDatabase::IteratorDB[It_Between]->SetParam(Domain);
  TDatabase::IteratorDB[It_OCAF]->SetParam(Domain);  
  
#endif    
  
  // free surface vertices
  d = 0;
  for(i=0;i<N_MovVert[2]-1;i++) // without last point
   {
    MovBoundVert[2][i]->GetCoords(x1, y1);
    MovBoundVert[2][i+1]->GetCoords(x2, y2);
//     x1 = FreePts[0][i];
//     x2 = FreePts[0][i+1];
//     y1 = FreePts[1][i];
//     y2 = FreePts[1][i+1];
     d +=sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)  );
   }
   h = d/((double)i-1.);
  
  MovBoundVert[0][0]->GetCoords(Lx, Ly);
  MovBoundVert[2][0]->GetCoords(Rx, Ry);
 
  //liquid-solid BD
  d = Rx-Lx;
  k = (int)(d/h); // No of intervals with step length 0.01
  if(k<2) k=2;     // minimum two intervals
  t = d/(double)k;
  N_MovVert[0] = k;
  S_BX = new double[N_MovVert[0]];
  S_BY = new double[N_MovVert[0]];
     
  for(i=0;i<N_MovVert[0];i++)
   {
    S_BX[i]= Lx + (double)i*t;
    S_BY[i]=0.0;
//      cout<<i<< " x :" << S_BX[i] << " -----------------y: " <<S_BY[i]<< endl;
   }
 
  //axial bound
  MovBoundVert[1][0]->GetCoords(x1, y1);
  MovBoundVert[0][0]->GetCoords(x2, y2);
  
   d  = y2 - y1; 
   k = (int)(fabs(d)/h);
   if(k<2) k=2;     // minimum two intervals     
   t = d/(double)k;
   N_MovVert[1] = k;
   A_BY = new double[N_MovVert[1]]; 

  for(i=0;i<N_MovVert[1];i++)
  {
   A_BY[i] = y1 + t*(double)i;
//    cout<<i<< " x :" << 0 << " -----------------y: " <<A_BY[i]<< endl;
  }

  
  area = TDatabase::ParamDB->Area; 
//======================================================================
// Triangular for grid generation begin
//======================================================================
  BoundPart = Domain->GetBdPart(0);
  UpdateSlipBound = (TBdLine*)BoundPart->GetBdComp(0);
  UpdateAxialBound = (TBdLine*)BoundPart->GetBdComp(2);

  Out.pointlist = NULL;
  Out.pointattributelist = NULL;
  Out.pointmarkerlist = NULL;
  Out.trianglelist = NULL;
  Out.triangleattributelist = NULL;
  Out.trianglearealist = NULL;
  Out.neighborlist = NULL;
  Out.segmentlist = NULL;
  Out.segmentmarkerlist = NULL;
  Out.holelist = NULL;
  Out.regionlist = NULL;
  Out.edgelist = NULL;
  Out.edgemarkerlist = NULL;
  Out.normlist = NULL;
 
  opts.seekp(std::ios::beg);

//OutPut("MESHGEN_REF_QUALIT " << TDatabase::ParamDB->MESHGEN_REF_QUALITY << endl);

  opts<<'p'; // Constrained Delaunay Triangulation:
           // initial values - only points defined on the boundary of the domain;
           // triangulation near boundary may variate from Delaunay criterion
  opts<<"q"<<TDatabase::ParamDB->MESHGEN_REF_QUALITY;
    // Quality mesh generation with no angles smaller than given degrees;
    // an alternate minimum angle may be specified after the `q'.
  //opts<<"a0.00001";
  opts<<'e'; // Outputs a list of edges of the triangulation
  opts<<'z'; // Numbers if items starting from 0
  //opts<<"VVVV"; // Gives detailed information about what Triangle is doing
  opts<<'Q'; // Supress all explanation of what Triangle is doing, unless an error occurs
  opts<<'Y'; // Supress adding vertices on boundary edges
//   opts<<'j'; //Jettisons(discard) vertices that are not part of the final triangulation from the output .node file (including duplicate input vertices and vertices ``eaten'' by holes). 
//  opts<<"a0.04"; // Imposes a maximum triangle area.
  opts<<"a"<< area; // Imposes a maximum triangle area.
  opts << "nA" << ends;
  
  In.numberofpoints = N_MovVert[0]+N_MovVert[1]+N_MovVert[2];
  In.pointlist = new double[2*In.numberofpoints];
  In.pointmarkerlist = new int[In.numberofpoints];
  In.numberofpointattributes = 0;

  In.numberofsegments = In.numberofpoints;
  In.segmentlist = new int[2*In.numberofsegments];
  In.segmentmarkerlist = new int[In.numberofsegments];
  In.numberofholes = 0;
  In.holelist = NULL;
  In.numberofregions = 0;
  In.regionlist = NULL;

  In_Index = 0;
  CurrComp = 1;

  // points and segments on the solid boundary (marker=1)
  for(i=0;i<N_MovVert[0];i++) // without last point
   {
    In.pointlist[2*In_Index] = S_BX[i];
    In.pointlist[2*In_Index+1] = S_BY[i];
//     cout<<In_Index<< " x :" << S_BX[i]<< " -----------------y: " <<S_BY[i]<< endl;
    In.pointmarkerlist[In_Index] = CurrComp;
    In.segmentlist[2*In_Index] = In_Index;
    In.segmentlist[2*In_Index+1] = In_Index+1;
    In.segmentmarkerlist[In_Index] = CurrComp;
    In_Index++;
   }

   CurrComp++;  

//   Sx = new double[N_MovVert[2]+1];
//   Sy = new double[N_MovVert[2]+1];
 
//   for(i=0;i<N_MovVert[2];i++) // without last point  
//    MovBoundVert[2][i]->GetCoords(Sx[i], Sy[i]);   
  
//   if(TDatabase::ParamDB->P10==0)
//    {
//     Sx[N_MovVert[2]] = A_BX[0];
//     Sy[N_MovVert[2]] = A_BY[0]; 
//     EqDist_Pts(N_MovVert[2]+1, Sx, Sy);         
//    }

  // points and segments on the free boundary (marker=2)
  for(i=0;i<N_MovVert[2];i++) // without last point
   {
    MovBoundVert[2][i]->GetCoords(tx, ty);
    In.pointlist[2*In_Index] = tx;
    In.pointlist[2*In_Index+1] = ty;
//    cout<<In_Index<< " x :" << tx<< " -----------------y: " <<ty<< endl;
    In.pointmarkerlist[In_Index] = CurrComp;
    In.segmentlist[2*In_Index] = In_Index;
    In.segmentlist[2*In_Index+1] = In_Index+1;
    In.segmentmarkerlist[In_Index] = CurrComp;
    In_Index++;
   }
//    delete [] Sx; delete [] Sy; 

   CurrComp++;
  // points and segments on the solid boundary (marker=1)
  for(i=0;i<N_MovVert[1];i++) // without last point
   {
    In.pointlist[2*In_Index] = 0.;
    In.pointlist[2*In_Index+1] = A_BY[i];
//     cout<<In_Index<< " x :" << A_BX[i]<< " -----------------y: " <<A_BY[i]<< endl;
    In.pointmarkerlist[In_Index] = CurrComp;
    In.segmentlist[2*In_Index] = In_Index;
    In.segmentlist[2*In_Index+1] = In_Index+1;
    In.segmentmarkerlist[In_Index] = CurrComp;
    In_Index++;
   }
   
   In.segmentlist[2*(In.numberofsegments-1)+1] = 0; 

  if(Out.pointlist!=NULL) {
    free(Out.pointlist); Out.pointlist = NULL;}
  if(Out.pointattributelist!=NULL) {
    free(Out.pointattributelist); Out.pointattributelist = NULL;}
  if(Out.pointmarkerlist!=NULL) {
    free(Out.pointmarkerlist); Out.pointmarkerlist = NULL;}
  if(Out.trianglelist!=NULL) {
    free(Out.trianglelist); Out.trianglelist = NULL;}
  if(Out.triangleattributelist!=NULL) {
    free(Out.triangleattributelist); Out.triangleattributelist = NULL;}
  if(Out.trianglearealist!=NULL) {
    free(Out.trianglearealist); Out.trianglearealist = NULL;}
  if(Out.neighborlist!=NULL) {
    free(Out.neighborlist); Out.neighborlist = NULL;}
  if(Out.segmentlist!=NULL) {
    free(Out.segmentlist); Out.segmentlist = NULL;}
  if(Out.segmentmarkerlist!=NULL) {
    free(Out.segmentmarkerlist); Out.segmentmarkerlist = NULL;}
  if(Out.holelist!=NULL) {
    free(Out.holelist); Out.holelist = NULL;}
  if(Out.regionlist!=NULL) {
    free(Out.regionlist); Out.regionlist = NULL;}
  if(Out.edgelist!=NULL) {
    free(Out.edgelist); Out.edgelist = NULL;}
  if(Out.edgemarkerlist!=NULL) {
    free(Out.edgemarkerlist); Out.edgemarkerlist = NULL;}
  if(Out.normlist!=NULL) {
    free(Out.normlist); Out.normlist = NULL;}
    
  triangulate((char*)opts.str().c_str(), &In, &Out, (struct triangulateio *)NULL);    

  
  Old_Coll = Domain->GetCollection(It_Finest, 0);
  Old_N_Cells = Old_Coll->GetN_Cells();
  Domain->GetTreeInfo(Old_CellTree, Old_N_RootCells);  
  if(Old_N_Cells != Old_N_RootCells) 
   exit(-1);
  
  VertexDel = new TVertex*[3*Old_N_RootCells]; 
  CurrVertex = 0;
    
  for(i=0;i<Old_N_Cells;i++)
    {
      cell = Old_Coll->GetCell(i);
      N_Joints = cell->GetN_Joints();
      N_Vertices = cell->GetN_Vertices();
      for(j=0;j<N_Joints;j++)
        {
         if(CurrVertex==0)
          {
              VertexDel[CurrVertex] = cell->GetVertex(j);
              CurrVertex++;
           }
          else
           {
            ID = 0;
            for(k=0;k<CurrVertex;k++)
            if(VertexDel[k]==cell->GetVertex(j))
             {
              ID = 1; break;
             }
            if(ID!=1)
             {
              VertexDel[CurrVertex] = cell->GetVertex(j);
              CurrVertex++;
             }
           } // else if(CurrVertex==0)

          ID = 0;
          for(k=0;k<CurrVertex;k++)
          if(VertexDel[k]==cell->GetVertex((j+1)%N_Vertices))
           {
            ID = 1; break;
           }
            if(ID!=1)
           {
            VertexDel[CurrVertex] = cell->GetVertex((j+1)%N_Vertices);
            CurrVertex++;
           }
        } // for j
    } // for i
    
  N_RootCells = Out.numberoftriangles;   
  Coordinates = Out.pointlist;
  Triangles = Out.trianglelist;
  
  N_G = Out.numberofpoints;
  NewVertices = new TVertex*[N_G];  

  for (i=0;i<N_G;i++)
     NewVertices[i] = new TVertex(Coordinates[2*i], Coordinates[2*i+1]);  

  // set bounding box
  left = bottom = 1e8;
  right = top = -1e8;

   for(i=0;i<In.numberofpoints;i++)
    {
      if(left>In.pointlist[2*i]) left = In.pointlist[2*i];
      if(right<In.pointlist[2*i]) right = In.pointlist[2*i];
      if(top<In.pointlist[2*i+1]) top = In.pointlist[2*i+1];
      if(bottom>In.pointlist[2*i+1]) bottom = In.pointlist[2*i+1];
    }
  
  
  Domain->SetBoundBox(right-left,top-bottom);
  Domain->SetBoundBoxstart(left,bottom); 
  
  UpdateSlipBound->SetParams(Lx, Ly, Rx-Lx, Ry-Ly);
  UpdateAxialBound->SetParams(0., A_BY[0], 0., Ly-A_BY[0]);  

  // generate cells
  CellTree = new TBaseCell*[N_RootCells];

  for (i=0;i<N_RootCells;i++)
  {
    CellTree[i] = new TMacroCell(TDatabase::RefDescDB[Triangle], 0);

    CellTree[i]->SetVertex(0, NewVertices[Out.trianglelist[3*i    ]]);
    CellTree[i]->SetVertex(1, NewVertices[Out.trianglelist[3*i + 1]]);
    CellTree[i]->SetVertex(2, NewVertices[Out.trianglelist[3*i + 2]]);

      ((TMacroCell *) CellTree[i])->SetSubGridID(0);
  }

  Domain->SetTreeInfo(CellTree, N_RootCells);
  
  // initialize iterators
  TDatabase::IteratorDB[It_EQ]->SetParam(Domain);
  TDatabase::IteratorDB[It_LE]->SetParam(Domain);
  TDatabase::IteratorDB[It_Finest]->SetParam(Domain);
  TDatabase::IteratorDB[It_Between]->SetParam(Domain);
  TDatabase::IteratorDB[It_OCAF]->SetParam(Domain);
 
  // search neighbours
  N_G = Out.numberofpoints;
  PointNeighb = new int[N_G];  
  
  memset(PointNeighb, 0, N_G *SizeOfInt);

  for (i=0;i<3*N_RootCells;i++)
    PointNeighb[Triangles[i]]++;

  for (i=0;i<N_G;i++)
    if (PointNeighb[i] > maxEpV) maxEpV = PointNeighb[i];

  delete [] PointNeighb;

  PointNeighb = new int[++maxEpV * N_G];

  memset(PointNeighb, 0, maxEpV * N_G *SizeOfInt);
  
  for(i=0;i<3*N_RootCells;i++)
  {
    j = Triangles[i]*maxEpV;
    PointNeighb[j]++;
    PointNeighb[j + PointNeighb[j]] = i/3;
  }
  
 // generate new edges
  N_G = Out.numberofedges;
  for (i=0;i<N_G;i++)
  {
    a = Out.edgelist[2*i];
    b = Out.edgelist[2*i+1];
    Neib[0] = -1;
    Neib[1] = -1;
    CurrNeib = 0;

    len1 = PointNeighb[a*maxEpV];
    len2 = PointNeighb[b*maxEpV];
// find indexes of cells containing the current edge
   for (j=1;j<=len1;j++)
    {
      Neighb_tmp = PointNeighb[a*maxEpV + j];
       for (k=1;k<=len2;k++)
        if (Neighb_tmp == PointNeighb[b*maxEpV + k])
        {
          Neib[CurrNeib++] = Neighb_tmp;
          break;
        }
      if (CurrNeib == 2) break;
    }

   if(Out.edgemarkerlist[i]) // 0 for inner edges and Boundcomp+1 for Boundedge respect
    {
      CurrComp = Out.edgemarkerlist[i] - 1;
      if (CurrComp >= 100000) CurrComp -= 100000;


      if(Domain->GetBdPart(0)->GetBdComp(CurrComp)->GetTofXY(
            NewVertices[a]->GetX(), NewVertices[a]->GetY(), T_a) ||
          Domain->GetBdPart(0)->GetBdComp(CurrComp)->GetTofXY(
            NewVertices[b]->GetX(), NewVertices[b]->GetY(), T_b))
       {
          cerr<<"Error: could not set parameter values"<<endl;
          OutPut(NewVertices[a]<<endl);
          OutPut(NewVertices[b]<<endl);
          cout << " CurrComp " << CurrComp <<endl;
        //  exit(0);
       }

      if (CurrNeib == 2)    // 2 cells contain the current edge
        if(Domain->GetBdPart(0)->GetBdComp(CurrComp)->IsFreeBoundary())
          Joint = new TIsoInterfaceJoint(Domain->GetBdPart(0)->GetBdComp(CurrComp), T_a, T_b, 
                                         CellTree[Neib[0]], CellTree[Neib[1]]);
        else
          Joint = new TInterfaceJoint(Domain->GetBdPart(0)->GetBdComp(CurrComp),
                  T_a, T_b, CellTree[Neib[0]], CellTree[Neib[1]]);
      else
        if(Domain->GetBdPart(0)->GetBdComp(CurrComp)->IsFreeBoundary())
          Joint = new TIsoBoundEdge(Domain->GetBdPart(0)->GetBdComp(CurrComp), T_a, T_b);
        else
          Joint = new TBoundEdge(Domain->GetBdPart(0)->GetBdComp(CurrComp), T_a, T_b);
    }
   else // inner edge
    {
    if (CurrNeib != 2)
        cerr << "Error!!!!!!!! not enough neighbours!" << endl;

     Joint = new TJointEqN(CellTree[Neib[0]], CellTree[Neib[1]]);
    }
 // find the local index for the point 'a' on the cell
    for (j=0;j<3;j++)
      if (Triangles[3*Neib[0]+j] == a) break;

    // find the local index for the point 'b' on the cell
    for (k=0;k<3;k++)
      if (Triangles[3*Neib[0]+k] == b) break;
   k = k*10 + j;

    switch (k)
    {
      case  1:
      case 10:
        j = 0;
        break;
      case 12:
      case 21:
        j = 1;
        break;
      case  2:
      case 20:
        j = 2;
        break;
    }
   CellTree[Neib[0]]->SetJoint(j, Joint);

   if (Neib[1] != -1)
    {
      // find the local index for the point 'a' on the cell
      for (j=0;j<3;j++)
        if (Triangles[3*Neib[1]+j] == a) break;

      // find the local index for the point 'b' on the cell
      for (k=0;k<3;k++)
        if (Triangles[3*Neib[1]+k] == b) break;

      k = k*10 + j;
    switch (k) // j will contain the local index for the current
      {
        case  1:
        case 10:
          j = 0;
          break;
        case 12:
        case 21:
          j = 1;
          break;
        case  2:
        case 20:
          j = 2;
          break;
      }
    CellTree[Neib[1]]->SetJoint(j, Joint);
    }

  if (Joint->GetType() == InterfaceJoint ||
        Joint->GetType() == IsoInterfaceJoint)
      ((TInterfaceJoint *) Joint)->CheckOrientation();
  
  } //  for (i=0;i<N_G;i++)

  delete [] NewVertices;
  delete [] PointNeighb;
  delete [] In.pointlist;
  delete [] In.pointmarkerlist;
  delete [] In.segmentlist;
  delete [] In.segmentmarkerlist; 
  
  if(Out.pointlist!=NULL) {
    free(Out.pointlist); Out.pointlist = NULL;}
  if(Out.pointattributelist!=NULL) { 
    free(Out.pointattributelist); Out.pointattributelist = NULL;}
  if(Out.pointmarkerlist!=NULL) {
    free(Out.pointmarkerlist); Out.pointmarkerlist = NULL;}
  if(Out.trianglelist!=NULL) {
    free(Out.trianglelist); Out.trianglelist = NULL;}
  if(Out.triangleattributelist!=NULL) {
    free(Out.triangleattributelist); Out.triangleattributelist = NULL;}
  if(Out.trianglearealist!=NULL) {
    free(Out.trianglearealist); Out.trianglearealist = NULL;}
  if(Out.neighborlist!=NULL) {
    free(Out.neighborlist); Out.neighborlist = NULL;}
  if(Out.segmentlist!=NULL) {
    free(Out.segmentlist); Out.segmentlist = NULL;}
  if(Out.segmentmarkerlist!=NULL) {
    free(Out.segmentmarkerlist); Out.segmentmarkerlist = NULL;}
  if(Out.holelist!=NULL) {
    free(Out.holelist); Out.holelist = NULL;}
  if(Out.regionlist!=NULL) {
    free(Out.regionlist); Out.regionlist = NULL;}
  if(Out.edgelist!=NULL) {
    free(Out.edgelist); Out.edgelist = NULL;}
  if(Out.edgemarkerlist!=NULL) {
    free(Out.edgemarkerlist); Out.edgemarkerlist = NULL;}
  if(Out.normlist!=NULL) {
    free(Out.normlist); Out.normlist = NULL;} 
  
// ======================================================================
// Triangular for grid generation end
// ======================================================================
 
//       // write grid into an Postscript file
//       os.seekp(std::ios::beg);
//       os << "Domain" << ".ps" << ends;
//       Domain->PS(os.str().c_str(),It_Finest,0);
// //      exit(0);
  
//======================================================================
// construct all finite element spaces
//======================================================================
  coll=Domain->GetCollection(It_Finest, 0);
  N_Cells = coll->GetN_Cells();
 
#ifdef __SURFACT__  
  TDatabase::IteratorDB[It_EQ]->SetParam(SurfDomain);
  TDatabase::IteratorDB[It_LE]->SetParam(SurfDomain);
  TDatabase::IteratorDB[It_Finest]->SetParam(SurfDomain);
  TDatabase::IteratorDB[It_Between]->SetParam(SurfDomain);
  TDatabase::IteratorDB[It_OCAF]->SetParam(SurfDomain);
    
  Domain2DSurf_2Phase(coll, SurfDomain, N_List);

  IFace_Coll = SurfDomain->GetCollection(It_Finest, 0);
  N_IFaceCells= IFace_Coll->GetN_Cells();
  cout<< " N_IFaceCells " << N_IFaceCells <<endl; 
  
//   delete [] FE1D_List; deleted in fespace1D
  FE1D_List = new FE1D[N_IFaceCells];
  for(j=0;j<N_IFaceCells;j++)
   FE1D_List[j] = FE1D(TDatabase::ParamDB->ANSATZ_ORDER);

  TDatabase::IteratorDB[It_EQ]->SetParam(Domain);
  TDatabase::IteratorDB[It_LE]->SetParam(Domain);
  TDatabase::IteratorDB[It_Finest]->SetParam(Domain);
  TDatabase::IteratorDB[It_Between]->SetParam(Domain);
  TDatabase::IteratorDB[It_OCAF]->SetParam(Domain);  
#endif   

  ORDER = 0;
  VSP = TDatabase::ParamDB->VELOCITY_SPACE;

  if (abs(VSP) > 20)
   {ORDER = abs(VSP) - 20;}
  else if ( abs(VSP) > 10)
    {ORDER = abs(VSP) - 10;}
  else ORDER = abs(VSP);

//======================================================================
// construct all finite element spaces
//======================================================================
  // get velocity and pressure spacess
  GetVelocityAndPressureSpace(coll,BoundCondition,
                              mortarcoll, velocity_space,
                              pressure_space, &pressure_space_code,
                              TDatabase::ParamDB->VELOCITY_SPACE,
                              TDatabase::ParamDB->PRESSURE_SPACE);
  N_U = velocity_space->GetN_DegreesOfFreedom();
  N_P = pressure_space->GetN_DegreesOfFreedom();
  GlobalNumbers = velocity_space->GetGlobalNumbers();
  BeginIndex = velocity_space->GetBeginIndex();

  // mesh velocity space 
  delete FESpaces_All[2];
  FESpaces_All[2] = new TFESpace2D(coll, NameString, WString, GridBoundCondition, 1, NULL);
  N_G = FESpaces_All[2]->GetN_DegreesOfFreedom();
  
  // thermal space
#ifdef __ENERGY__    
  delete FESpaces_All[3];
  FESpaces_All[3] = new TFESpace2D(coll, NameString, TString, HeatBoundCondition,
                                  TDatabase::ParamDB->ANSATZ_ORDER, NULL); 
  N_ThermalDOF = FESpaces_All[3]->GetN_DegreesOfFreedom();
#endif  
  
#ifdef __SURFACT__  
   SurfactSpace  = new TFESpace2D(coll, NameString, CString, SurfactBoundCondition,
                                  TDatabase::ParamDB->ANSATZ_ORDER, NULL);
 
   N_surfactDOF = SurfactSpace->GetN_DegreesOfFreedom();
   OutPut("N_SurfactDOF    : "<< setw(10) << N_surfactDOF  << endl);
 
   delete IFaceFeSpaces[0];
   IFaceFeSpaces[0] = new TFESpace1D(IFace_Coll, IFaceSString, IFaceSString, FE1D_List);
 
   N_IsurfactDOF = IFaceFeSpaces[0]->GetN_DegreesOfFreedom();
   OutPut("N_SurfSurfactDOF    : "<< setw(10) << N_IsurfactDOF  << endl);   
   
   

   SurfSurfactSpace =  new TFESpace2D(coll, NameString, IFaceSString, SurfactBoundCondition,
                                     TDatabase::ParamDB->ANSATZ_ORDER, NULL);

   N_S =  SurfSurfactSpace->GetN_DegreesOfFreedom();     
#endif 
//===========================================================================  

  delete [] Bound_Joint[0];
  delete [] Bound_Joint[1];                      

  delete [] MovBoundVert[0];
  delete [] MovBoundVert[1];   
  delete [] MovBoundVert[2];  

  delete [] Free_Joint;    
  delete [] Free_Cells;  
  delete [] IsoCellEdgeNos[0];  
  delete [] IsoCellEdgeNos[1];  
    
  Ry = 0.;
  GetMovingBoundData(coll, N_MovVert, Bound_Joint, MovBoundVert, Free_Joint,
                     Free_Cells, IsoCellEdgeNos, Rx, Ry);
  
//======================================================================
// construct all finite element functions
//======================================================================
  N_Unknowns = 2*N_U + N_P;
  
  OutPut("dof velocity : "<< setw(10) << 2* N_U << endl);
  OutPut("dof pressure : "<< setw(10) << N_P << endl);
  OutPut("dof all      : "<<  setw(10) << N_Unknowns  << endl);  
  
  delete [] Rhs_All[0];
  Rhs_All[0] = new double[N_Unknowns];
  sol = new double[N_Unknowns];
  memset(sol, 0, N_Unknowns*SizeOfDouble);
   
  u =  new TFEVectFunct2D(velocity_space, UString, UString, sol, N_U, 2);
  p = new TFEFunction2D(pressure_space, PString,  PString,  sol+2*N_U, N_P);

  u->Interpolate(FEVectFuncts_All[0]);
//====================================================================== 
// impose no-penetration condition 
//====================================================================== 
  ValuesU2 = sol + N_U;
  
  for(i=0;i<N_Cells;i++)
   {
    cell = coll->GetCell(i);
    k = cell->GetN_Edges();
    
    for(l=0;l<k;l++)
     {
      Joint = cell->GetJoint(l);    
      
      if(Joint->GetType() == BoundaryEdge)
        {
         FeId = velocity_space->GetFE2D(i, cell);  
         FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FeId);
         JointDOF = FeDesc->GetJointDOF(l); 
         N_DOF = FeDesc->GetN_JointDOF();
         DOF = GlobalNumbers+BeginIndex[i];
    
         Solid_Joint = (TBoundEdge *)Joint;
         BoundComp = Solid_Joint->GetBoundComp();
         comp=BoundComp->GetID();

         if(comp==0) 
          {           
           for(j=0;j<N_DOF;j++)
            ValuesU2[DOF[JointDOF[j]]] = 0.;
          }  
         else if(comp==2)
          {
           for(j=0;j<N_DOF;j++)
            sol[DOF[JointDOF[j]]] = 0.;  
          }
 
        }// if(Joint->GetType() == Bo
     } // for(l=0;l<k;l++)
   } //  for(i=0;i<N_Cells;i++)
//======================================================================   
// no need to interpolate pressure for DirectSolvers !!!
//======================================================================    
  delete FESpaces_All[0];
  delete FESpaces_All[1]; 
 
  FESpaces_All[0] = velocity_space;
  FESpaces_All[1] = pressure_space;    

  delete [] Sol_All[0];
  Sol_All[0] = sol;
 
  delete FEVectFuncts_All[0]; 
  FEVectFuncts_All[0] = u;
 
  FEFunctions_All[0] = FEVectFuncts_All[0]->GetComponent(0);
  FEFunctions_All[1] = FEVectFuncts_All[0]->GetComponent(1); 

  delete FEFunctions_All[2];
  FEFunctions_All[2] = p;
  
//======================================================================
// grid space finite element functions
//====================================================================== 
  delete [] Sol_All[1];
  Sol_All[1] = new double[2*N_G];
  delete [] Rhs_All[1]; 
  Rhs_All[1] = new double[2*N_G];  
  
  memset(Sol_All[1], 0, 2*N_G*SizeOfDouble); 
  delete FEVectFuncts_All[1];
  FEVectFuncts_All[1]  = new TFEVectFunct2D(FESpaces_All[2], WString, WString, Sol_All[1], N_G, 2);
  
  FEFunctions_All[3] = FEVectFuncts_All[1]->GetComponent(0);
  FEFunctions_All[4] = FEVectFuncts_All[1]->GetComponent(1); 

//======================================================================
// thermal space finite element functions
//======================================================================    
#ifdef __ENERGY__  
  delete [] Rhs_All[2];
  Rhs_All[2] = new double[N_ThermalDOF]; 
  memset(Rhs_All[2], 0, N_ThermalDOF*SizeOfDouble);
    
  delete [] Sol_All[2];
  Sol_All[2] = new double[N_ThermalDOF];

  memset(Sol_All[2], 0, N_ThermalDOF*SizeOfDouble);
  
  // thermal fefunction
  delete FEFunctions_All[5];
  FEFunctions_All[5] = new TFEFunction2D(FESpaces_All[3], TString, TString, Sol_All[2], N_ThermalDOF);
#endif
  
//======================================================================
// surfactant space finite element functions
//======================================================================    
#ifdef __SURFACT__
  csol  = new double[N_surfactDOF];
  memset(csol, 0, N_surfactDOF*SizeOfDouble);
  New_c = new TFEFunction2D(SurfactSpace, CString,  CString, csol, N_surfactDOF);
  
#ifdef  __SOLUBLE__  
  New_c->Interpolate(FEFunctions_All[6]);  
#endif
  
  delete [] Sol_All[3]; 
  delete [] Rhs_All[3];
  delete FEFunctions_All[6]; 
  delete FESpaces_All[4];
   
  Sol_All[3] = csol;
  Rhs_All[3] = new double[N_surfactDOF];    
  FEFunctions_All[6] = New_c;
  FESpaces_All[4] = SurfactSpace;
  
  // surf surfact
  delete IFaceFeFunct[0];
  delete [] Sol_All[4];
  delete [] Rhs_All[4];
  
  Sol_All[4] = new double[N_IsurfactDOF]; 
  Rhs_All[4] = new double[N_IsurfactDOF]; 
  
  IFaceFeFunct[0] = new TFEFunction1D(IFaceFeSpaces[0], IFaceSString, IFaceSString, Sol_All[4], N_IsurfactDOF);
  
  Surfcsol  =  new double[N_S];
  SurfNew_c  = new TFEFunction2D(SurfSurfactSpace, IFaceSString,  IFaceSString, Surfcsol,  N_S); 
  SurfNew_c->Interpolate(FEFunctions_All[7]);   

  
  delete [] Sol_All[5];
  delete FEFunctions_All[7];
  delete FESpaces_All[5];     
  
  Sol_All[5] = Surfcsol;  
  FEFunctions_All[7] = SurfNew_c;  
  FESpaces_All[5] = SurfSurfactSpace;

  memset(Rhs_All[3], 0, N_surfactDOF*SizeOfDouble);  
#endif 
// ======================================================================
// allocate memory for all matrices
// ======================================================================  
  delete Structure_All[0];  delete Structure_All[1];
  
  Structure_All[0] = new TStructure2D(FESpaces_All[1], FESpaces_All[0]); // B
  Structure_All[1] = new TStructure2D(FESpaces_All[0], FESpaces_All[1]); // BT

  delete SquareStructure_All[0];  delete SquareStructure_All[1]; 

  
  //velo 
  SquareStructure_All[0] = new TSquareStructure2D(FESpaces_All[0]);  
  SquareStructure_All[0]->Sort();  
  
  // grid 
  SquareStructure_All[1] = new TSquareStructure2D(FESpaces_All[2]); 
  SquareStructure_All[1]->Sort();  
   
  //thermal
#ifdef __ENERGY__    
  delete SquareStructure_All[2];
  SquareStructure_All[2] = new TSquareStructure2D(FESpaces_All[3]);
  SquareStructure_All[2]->Sort(); 
#endif
  
#ifdef __SURFACT__
  //surfact
  delete SquareStructure_All[3];
  SquareStructure_All[3] = new TSquareStructure2D(FESpaces_All[4]);
  SquareStructure_All[3]->Sort();
 
  /* interface surfactant matrices */
  delete IFaceStruct[0];
  IFaceStruct[0] = new TSquareStructure1D(IFaceFeSpaces[0]);
  IFaceStruct[0]->Sort();
#endif      
  
  // u
  for(i=0; i<10; i++)
   {
    delete SqMat_All[i];
    SqMat_All[i] = new TSquareMatrix2D(SquareStructure_All[0]);
   }

  // B
  for(i=0; i<2; i++)
   {
    delete Mat_All[i];
    Mat_All[i] = new TMatrix2D(Structure_All[0]);   
   }  

  // BT
  for(i=2; i<4; i++)
   {
    delete Mat_All[i];
    Mat_All[i] = new TMatrix2D(Structure_All[1]);   
   }  

  // for mesh
  for(i=10; i<14; i++)
   {
    delete SqMat_All[i];
    SqMat_All[i] = new TSquareMatrix2D(SquareStructure_All[1]);  
   }  

  // for heat
#ifdef __ENERGY__    
  for(i=14; i<16; i++)
   {
    delete SqMat_All[i];
    SqMat_All[i] = new TSquareMatrix2D(SquareStructure_All[2]);  
   }  
#endif   

#ifdef __SURFACT__  
  //bulk surfact
  for(i=16; i<18; i++)
   {
    delete SqMat_All[i];
    SqMat_All[i] = new TSquareMatrix2D(SquareStructure_All[3]);   // C_M,  C_A
   }  
 
 /* interface surfactant matrices */
  for(i=0; i<2; i++)
   {
    delete SqMat_IFace[i];
    SqMat_IFace[i] = new TSquareMatrix1D(IFaceStruct[0]); // Gamma_A,Gamma_M   
   } 
#endif   

// ======================================================================
// delete old mesh
// ======================================================================  
    
  for(i=0;i<CurrVertex;i++)
    delete VertexDel[i];

  delete [] VertexDel;
  OutPut(CurrVertex<<" vertices were deleted"<<endl);
  
   // remove all existing cells and joints
  for(i=0;i<Old_N_RootCells;i++)
     delete (TGridCell*)Old_CellTree[i];
  OutPut(Old_N_RootCells<<" cells were deleted"<<endl);
   delete [] Old_CellTree;

  delete Old_Coll;

  delete [] S_BX;
  delete [] S_BY;
//   delete [] A_BX;
  delete [] A_BY; 
 
}// RemeshAxial3D_ImpDrop



 void GetSurfactMass(TFEFunction2D *fefunction, TFEFunction1D *fefunct_low,
                   int *Cell_array, int *Joint_array, double *errors)
{
  int i, j, k, l, m, n, N, N_Cells_low, N_LocalUsedElements, local_i, local_dof, ORDER;
  int N_BaseFunct,  N_Points, N_JointDOF, begin, end, *N_BaseFuncts, IJoint;
  TCollection *Coll, *Coll_low;
  TBaseCell *Me, *Me_low;
  FE2D FEId;
  TFE2D *ele;
  FE1D FEId_low;
  TFE1D *Element;
  TFESpace2D *fespace;
  TFESpace1D *fespace_low;
  BaseFunct2D LocBF[N_BaseFuncts2D];
  BaseFunct2D *BaseFuncts;
  int *DOF, *JointDOF, *BeginIndex, *GlobalNumbers, N_BaseFunct_low;
  TFEDesc2D *FeDesc;
  TFEDesc1D *FeDesc_low;
  BF2DRefElements RefElement;
  RefTrans2D RefTrans;
  TRefTrans2D *F_K;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1;
  int N_LinePoints;
  double *LineWeights, t0, t1, normn, n0, n1, X_P, Y_P;
  double *u, **uref, **uxiref, **uetaref, **uzetaref;
  double *Weights, *p1, *p2, *zeta, LocL2U, LocH1U,  Surf[100], Surf_DOF[10], Surf_val[10];
  double r_axial;
  double  X_B[MaxN_QuadPoints_2D], Y_B[MaxN_QuadPoints_2D], Exact_Surf[10], v, Mult;
  double uorig[MaxN_BaseFunctions2D], uxorig[MaxN_BaseFunctions2D];
  double uyorig[MaxN_BaseFunctions2D];
  double test00, test10, test01;
  double U, ux, uy, d1, d2;
//   ngrad, h_K =0., h_K_min=1e8, h_K_max=-1e8;

  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFuncts = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  fespace=fefunction->GetFESpace2D();
  Coll = fespace->GetCollection();
  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();
  u = fefunction->GetValues();

  fespace_low=fefunct_low->GetFESpace1D();
  Coll_low = fespace_low->GetCollection();
  N_Cells_low = Coll_low->GetN_Cells();

  errors[0] =0.; // surfactant mass
  errors[1] =0.; // surface area

// ########################################################################
// loop over all surf cells
// ########################################################################
  for(i=0;i<N_Cells_low;i++)
  {
//     h_K =0.;
    Me_low = Coll_low->GetCell(i);
    FEId_low = fespace_low->GetFE1D(i, Me_low);
    Element = TFEDatabase2D::GetFE1D(FEId_low);
    N_BaseFunct_low = Element->GetN_DOF();

    N = Cell_array[i];
    Me = Coll->GetCell(N);
    IJoint = Joint_array[i];

    FEId = fespace->GetFE2D(N, Me);
    ele = TFEDatabase2D::GetFE2D(FEId);

    RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
    FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
    N_JointDOF =  FeDesc->GetN_JointDOF();
    JointDOF = FeDesc->GetJointDOF(IJoint);
    N_BaseFunct = FeDesc->GetN_DOF();
    DOF = GlobalNumbers + BeginIndex[N];

    l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
    ORDER = TFEDatabase2D::GetAccuracyFromFE2D(FEId);
    LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
    qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);

    switch(RefElement)
     {
      case BFUnitSquare:
          RefTrans = QuadIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TQuadIsoparametric *)F_K)->SetCell(Me);
          ((TQuadIsoparametric *)F_K)->GetOrigBoundFromRef(IJoint, N_LinePoints, zeta, X_B, Y_B);
      break;

      case BFUnitTriangle:
          RefTrans = TriaIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TTriaIsoparametric *)F_K)->SetApproximationOrder(ORDER);	  
          ((TTriaIsoparametric *)F_K)->SetCell(Me);
          ((TTriaIsoparametric *)F_K)->GetOrigBoundFromRef(IJoint, N_LinePoints, zeta, X_B, Y_B);

      break;
      
      default:
       OutPut("RefElement N_Unknowns GetSurfactMass " << endl);
       exit(1);      
      
    } // endswitch

    TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(LineQuadFormula);


      uref = TFEDatabase2D::GetJointValues2D(BaseFuncts[FEId],
                        LineQuadFormula, IJoint);
      uxiref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId],
                        LineQuadFormula, IJoint, D10);
      uetaref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId],
                        LineQuadFormula, IJoint, D01);


    for(k=0;k<N_LinePoints;k++)
     {

         switch(RefElement)
          {
            case BFUnitSquare:
              ((TQuadIsoparametric *)F_K)->SetCell(Me);
              ((TQuadIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
              ((TQuadIsoparametric *) F_K)->GetTangent(IJoint, zeta[k], t0, t1);
            break;

            case BFUnitTriangle:
              ((TTriaIsoparametric *)F_K)->SetApproximationOrder(ORDER);
              ((TTriaIsoparametric *)F_K)->SetCell(Me); 
              ((TTriaIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
              ((TTriaIsoparametric *) F_K)->GetTangent(IJoint, zeta[k], t0, t1);
            break;
          } // endswitch


      normn = sqrt(t0*t0+t1*t1);
      n0 =  t1/normn;
      n1 = -t0/normn;

//  //  Projection of all points to the free boundary in normal direction  !! only for circle!!!
//       v =  sqrt(X_B[k]* X_B[k] + Y_B[k]*Y_B[k] );// /r;
      X_P =  X_B[k]  ;
      Y_P =  Y_B[k]  ;
      r_axial = fabs(X_P);

//       get solution 
     U=0.;
     for(l=0;l<N_BaseFunct_low;l++)
      {
       local_dof   = JointDOF[l];
       test00 = uorig[local_dof];

       m = DOF[local_dof];
       U  += u[m]*test00;
      }
// cout << "ux " << ux <<endl;
//       ExactS(X_P, Y_P, Exact_Surf);
//       if(Y_B[k]==0)
// 	cout<< "x " << X_P << " U " << U <<endl;
	
      Mult = LineWeights[k]*normn*r_axial;
//       h_K +=normn;
//       errors[0] +=Mult*(Exact_Surf[0]-U)*(Exact_Surf[0]-U);
      errors[0] +=U*Mult;
      errors[1] +=Mult;
     } //  for(k=0;k<N_LinePoints;k++)
  } // for(i=0;i<N
//     OutPut("h_K_min and h_K_max of free surface: "<< h_K_min << " " << h_K_max<<endl; );

   errors[0] *=2.*Pi;
   errors[1] *=2.*Pi;

//    OutPut( "Time, Surfactant Mass " <<TDatabase::TimeDB->CURRENTTIME<< " " <<errors[0]<< " "<<endl);
//    OutPut( "Time, Surface area " <<TDatabase::TimeDB->CURRENTTIME<< " " <<errors[1]<< " "<<endl);

// exit(0);
}



void SurfSurfactCorrectMass(TFEFunction2D *fefunction, TFEFunction1D *fefunct_low,
                            int *Cell_array, int *Joint_array, double OldMass)
{
  int i, Length;
  
  double sum=0., *Values1D, *OrigSol, *w, params[3], WeightedSurArea, MeanDiff;
  
  Length = fefunct_low->GetLength();
  Values1D = fefunct_low->GetValues();
  
  for(i=0;i<Length;i++)
   {
    sum +=Values1D[i];
   }
  
  // if the mean value or the value is zero, do nothing
  if(fabs(sum)>1e-12 )
   { 
    // add the MeanMassDiff to the sol with weights to avoid negative values
    w = new double[Length];
    OrigSol = new double[Length];
   
    // store the orignal values of the fefunction
    memcpy(OrigSol, Values1D, Length*SizeOfDouble);
   
   //weights, sol with 0 values will not be altered, otherwise sol vale become negative due to correction 
   for(i=0;i<Length;i++)
    { 
      w[i] = Values1D[i]/sum; 
//        w[i] = 1.;
    }   
   
    // calculate the weitgted surfce area, if w=1, this step is not needed, but sol may take negative
    memcpy(Values1D, w, Length*SizeOfDouble);
    MapSurfToDomain(fefunct_low, fefunction, Cell_array, Joint_array);
    GetSurfactMass(fefunction, fefunct_low, Cell_array, Joint_array, params);
    WeightedSurArea = params[0];
    
    // copy back the orig values to fefunction
    memcpy(Values1D, OrigSol, Length*SizeOfDouble);   
    MapSurfToDomain(fefunct_low, fefunction, Cell_array, Joint_array);
    GetSurfactMass(fefunction, fefunct_low, Cell_array, Joint_array, params);   
    MeanDiff = (OldMass - params[0])/WeightedSurArea;
     
//    cout << "Surface MeanDiff " << MeanDiff << endl;
    
    for(i=0;i<Length;i++)
     {Values1D[i] += (MeanDiff*w[i]); }  
   
    delete [] w;
    delete [] OrigSol;
   } // if(fabs(sum)>1e-12 )
} //SurfSurfactCorrectMass


void Surfact2D_InterfaceInt(TSquareMatrix2D *A, double *rhs, BoundCondFunct2D SurfactBoundCondition, 
                           TFEFunction2D *FeFunct, TFEFunction1D *SurfaceFeFunct, 
                           TFESpace2D *GlobalFE_Space, int *Cell_array, int *Joint_array, double *XmaxVal)
{
  int i, j, k, l, N_Cells, N_Vertices, N_Edges, ORDER;
  int *BeginIndex, *GlobalNumbers, *DOF, *DOF_Surface, TestDOF, AnsatzDOF;
  int *BeginIndex_Surface, *GlobalNumbers_Surface, dof_surface;
  int JointNumbers[MAXN_JOINTS], N_IsoJoints, N_JointDOF;
  int *KCol, *RowPtr, *JointDOF, N_DOF, N, IJoint, local_dof, m;
  int N_LinePoints, ActiveBound, N_BaseFunct_Surface, Global_N;
  int N_BaseFunct, *N_BaseFuncts, index1, Surf_cell_No, N_Cells_Surface;

  double **uref, **uxiref, **uetaref, val, rhsval, *V, *Surfactant;
  double uorig[MaxN_BaseFunctions2D], uxorig[MaxN_BaseFunctions2D];
  double uyorig[MaxN_BaseFunctions2D], U, test00, C;
  double *LineWeights, *zeta, X_B[100], Y_B[100];
  double t0, t1, r_axial, normn, surf, surfx, surfy, surf_flux, CMass, GammaMass;
//   double C_Infty, Gamma_Infty, Bi_a, Bi_d;
  double Bi, Da, beta;
   
  BoundCond Cond0, Cond1;
  FE2D FEId;
  FE1D FEId_Surface;
  TFE2D *ele, *ele_surf;
  TFE1D *Element;
  RefTrans2D RefTrans;
  TRefTrans2D *F_K;
  BF2DRefElements RefElement;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1;
  TBaseCell *Me, *Me_Surface;
  TCollection *GlobalColl, *Coll, *Coll_Surface;
  TFESpace2D *fespace;
  TFESpace1D *fespace_surface;
  TJoint *joint;
  TFEDesc2D *FeDesc;
  TFEDesc1D *FeDesc_low;

  GlobalColl = GlobalFE_Space->GetCollection();

  fespace = A->GetFESpace();    // outer phase
  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();
  BaseFunct2D *BaseFuncts;


  ActiveBound = fespace->GetActiveBound();
  GlobalNumbers = fespace->GetGlobalNumbers();
  BeginIndex = fespace->GetBeginIndex();
  Surfactant = FeFunct->GetValues();

  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFuncts = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  fespace_surface=SurfaceFeFunct->GetFESpace1D();
  Coll_Surface = fespace_surface->GetCollection();
  N_Cells_Surface = Coll_Surface->GetN_Cells();

  GlobalNumbers_Surface = fespace_surface->GetGlobalNumbers();
  BeginIndex_Surface = fespace_surface->GetBeginIndex();
  V = SurfaceFeFunct->GetValues();

//   Gamma_Infty = TDatabase::ParamDB->REACTOR_P13;
//   C_Infty = TDatabase::ParamDB->REACTOR_P14;
//   Bi_a = TDatabase::ParamDB->REACTOR_P15;
//   Bi_d = TDatabase::ParamDB->REACTOR_P16;
  
  Bi = TDatabase::ParamDB->REACTOR_P13; 
  Da = TDatabase::ParamDB->REACTOR_P14;
  beta = TDatabase::ParamDB->REACTOR_P15;
    
  XmaxVal[0] = -1.e10;
  XmaxVal[1] =  0.;
  surf_flux =  0.;

 for(i=0;i<N_Cells_Surface;i++)
  {
    Me_Surface = Coll_Surface->GetCell(i);
    FEId_Surface = fespace_surface->GetFE1D(i, Me_Surface);
    Element = TFEDatabase2D::GetFE1D(FEId_Surface);
    N_BaseFunct_Surface = Element->GetN_DOF();
    DOF_Surface = GlobalNumbers_Surface + BeginIndex_Surface[i];

    Global_N = Cell_array[i];     // outerphase
    Me = GlobalColl->GetCell(Global_N);
//     N = Me->GetLocalCellNo();
    N=Global_N; // free surf flows N=Global_N;
    IJoint = Joint_array[i];

    FEId = fespace->GetFE2D(N, Me);
    ele = TFEDatabase2D::GetFE2D(FEId);

    RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
    FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
    N_JointDOF =  FeDesc->GetN_JointDOF();
    JointDOF = FeDesc->GetJointDOF(IJoint);
    N_BaseFunct = FeDesc->GetN_DOF();
    DOF = GlobalNumbers + BeginIndex[N];

    l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
    ORDER = TFEDatabase2D::GetAccuracyFromFE2D(FEId);
    LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
    qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);

    switch(RefElement)
     {
      case BFUnitSquare:
          RefTrans = QuadIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TQuadIsoparametric *)F_K)->SetCell(Me);
          ((TQuadIsoparametric *)F_K)->GetOrigBoundFromRef(IJoint, N_LinePoints, zeta, X_B, Y_B);
      break;

      case BFUnitTriangle:
          RefTrans = TriaIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TTriaIsoparametric *)F_K)->SetApproximationOrder(ORDER);  
          ((TTriaIsoparametric *)F_K)->SetCell(Me);
          ((TTriaIsoparametric *)F_K)->GetOrigBoundFromRef(IJoint, N_LinePoints, zeta, X_B, Y_B);
      break;
    } // endswitch

    TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)
                  ->MakeRefElementData(LineQuadFormula);


    uref = TFEDatabase2D::GetJointValues2D(BaseFuncts[FEId],
                        LineQuadFormula, IJoint);
    uxiref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId],
                        LineQuadFormula, IJoint, D10);
    uetaref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId],
                        LineQuadFormula, IJoint, D01);

    for(k=0;k<N_LinePoints;k++)
     {

         switch(RefElement)
          {
            case BFUnitSquare:
              ((TQuadIsoparametric *)F_K)->SetCell(Me);
              ((TQuadIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
              ((TQuadIsoparametric *) F_K)->GetTangent(IJoint, zeta[k], t0, t1);
            break;

            case BFUnitTriangle:
              ((TTriaIsoparametric *)F_K)->SetApproximationOrder(ORDER);
              ((TTriaIsoparametric *)F_K)->SetCell(Me); 
              ((TTriaIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
              ((TTriaIsoparametric *) F_K)->GetTangent(IJoint, zeta[k], t0, t1);
            break;
          } // endswitch

     normn = sqrt(t0*t0+t1*t1);
     r_axial = fabs(X_B[k]);
     
         if(X_B[k]<=0)
          {
           cout <<"X_B[k] negative in SurfInt  change Quad rule " <<  X_B[k] <<endl;
//         exit(0);
          }

     // get the outer phase surfactant value at this point
     C=0.;
     for(l=0;l<N_BaseFunct;l++)
       C  += Surfactant[DOF[l]]*uorig[l];


     // get surfactant solution  on the interface
     U=0.;
#ifndef __MASSTRANSTEST__
     for(l=0;l<N_BaseFunct_Surface;l++)
      {
       local_dof   = JointDOF[l];
       test00 = uorig[local_dof];
       m = DOF_Surface[l];
       U  += V[m]*test00;
      }
#endif
     
     //U=0.; // for bulk concentration test problem
//     if( XmaxVal[0]< r_axial )
//      {
//       XmaxVal[0] = r_axial;
//       XmaxVal[1] = U;
//      }
     if( XmaxVal[1]< U )
      {
       XmaxVal[0] = r_axial;
       XmaxVal[1] = U;
      }
     rhsval =  -beta*C*(1. - U) + Bi*Da*U;
     rhsval *= normn * LineWeights[k]*r_axial;
     
//   OutPut("Flux rhsval " <<  rhsval <<endl);
  
//         // update rhs for all test functions
     for(l=0;l<N_BaseFunct;l++)
       {
         if((index1 = DOF[l])<ActiveBound)
          {
           rhs[index1] += rhsval*uorig[l];
           surf_flux += rhsval*uorig[l];
//            cout <<"surf_flux " << surf_flux << " rhs "<<  rhsval*uorig[l] <<endl;
          }
       }//  for(l=0;l<N_BaseFunct;l++)
//          cout  <<endl;
    }  // for(k=0;k<N_LinePoints;k
  } // for(i=0;i<N_Cells;

//   OutPut("Flux at the interface " << 2*Pi*surf_flux <<endl);
//   cout << " GammaMass on Interface; " << 2*Pi*GammaMass <<endl;
//   cout << " CMass on Interface; " << 2*Pi*CMass <<endl;
//    cout << " xpos; " << xpos << " XmaxVal; " << XmaxVal << endl;

// exit(0);
}


void Surfact2D_InterfaceInt_Implicit(TSquareMatrix2D *A, double *rhs,
                           BoundCondFunct2D SurfactBoundCondition,
                           TFEFunction2D *FeFunct, TFEFunction1D *SurfaceFeFunct, 
                           TFESpace2D *GlobalFE_Space, int *Cell_array, int *Joint_array, double *XmaxVal)
{
  int i, j, k, l, N_Cells, N_Vertices, N_Edges, ORDER;
  int *BeginIndex, *GlobalNumbers, *DOF, *DOF_Surface, TestDOF, AnsatzDOF;
  int *BeginIndex_Surface, *GlobalNumbers_Surface, dof_surface;
  int JointNumbers[MAXN_JOINTS], N_IsoJoints, N_JointDOF;
  int *KCol, *RowPtr, *JointDOF, N_DOF, N, IJoint, local_dof, m;
  int N_LinePoints, ActiveBound, N_BaseFunct_Surface, Global_N, index2;
  int N_BaseFunct, *N_BaseFuncts, index1, Surf_cell_No, N_Cells_Surface;

  double **uref, **uxiref, **uetaref, val, rhsval, *V, *Surfactant;
  double uorig[MaxN_BaseFunctions2D], uxorig[MaxN_BaseFunctions2D];
  double uyorig[MaxN_BaseFunctions2D], U, test00, C;
  double *LineWeights, *zeta, X_B[100], Y_B[100];
  double t0, t1, r_axial, normn, surf, surfx, surfy, surf_flux, CMass, GammaMass;
//   double C_Infty, Gamma_Infty, Bi_a, Bi_d, *ValuesA;
  double Bi, Da, beta, *ValuesA;
  
  bool IsIsoparametric;
   
  BoundCond Cond0, Cond1;
  FE2D FEId;
  FE1D FEId_Surface;
  TFE2D *ele, *ele_surf;
  TFE1D *Element;
  RefTrans2D RefTrans;
  TRefTrans2D *F_K;
  BF2DRefElements RefElement;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1;
  TBaseCell *Me, *Me_Surface;
  TCollection *GlobalColl, *Coll, *Coll_Surface;
  TFESpace2D *fespace;
  TFESpace1D *fespace_surface;
  TJoint *joint;
  TFEDesc2D *FeDesc;
  TFEDesc1D *FeDesc_low;

  GlobalColl = GlobalFE_Space->GetCollection();

  fespace = A->GetFESpace();    // outer phase
  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();
  BaseFunct2D *BaseFuncts;

  ActiveBound = fespace->GetActiveBound();
  GlobalNumbers = fespace->GetGlobalNumbers();
  BeginIndex = fespace->GetBeginIndex();
  Surfactant = FeFunct->GetValues();

  ValuesA  = A->GetEntries();
  RowPtr = A->GetRowPtr();
  KCol = A->GetKCol();

 
  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFuncts = TFEDatabase2D::GetN_BaseFunctFromFE2D();
  fespace_surface=SurfaceFeFunct->GetFESpace1D();
  Coll_Surface = (SurfaceFeFunct->GetFESpace1D())->GetCollection();
  N_Cells_Surface = Coll_Surface->GetN_Cells();
 
  GlobalNumbers_Surface = fespace_surface->GetGlobalNumbers();
  BeginIndex_Surface = fespace_surface->GetBeginIndex();
  V = SurfaceFeFunct->GetValues();

  Bi = TDatabase::ParamDB->REACTOR_P13; 
  Da = TDatabase::ParamDB->REACTOR_P14;
  beta = TDatabase::ParamDB->REACTOR_P15;
  
  XmaxVal[0] = -1.e10;
  XmaxVal[1] =  0.;

  
 for(i=0;i<N_Cells_Surface;i++)
  {
    Me_Surface = Coll_Surface->GetCell(i);
    FEId_Surface = fespace_surface->GetFE1D(i, Me_Surface);
    Element = TFEDatabase2D::GetFE1D(FEId_Surface);
    N_BaseFunct_Surface = Element->GetN_DOF();
    DOF_Surface = GlobalNumbers_Surface + BeginIndex_Surface[i];
    Global_N = Cell_array[i];     // outerphase
    Me = GlobalColl->GetCell(Global_N);
//     N = Me->GetLocalCellNo();
    N=Global_N; // free surf flows N=Global_N;
    IJoint = Joint_array[i];

    FEId = fespace->GetFE2D(N, Me);
    ele = TFEDatabase2D::GetFE2D(FEId);
    RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
    FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
    N_JointDOF =  FeDesc->GetN_JointDOF();
    JointDOF = FeDesc->GetJointDOF(IJoint);
    
    N_BaseFunct = FeDesc->GetN_DOF();
    DOF = GlobalNumbers + BeginIndex[N];
    l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
    ORDER = TFEDatabase2D::GetAccuracyFromFE2D(FEId);
    LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
    qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);

//     IsIsoparametric = FALSE;
//     if (TDatabase::ParamDB->USE_ISOPARAMETRIC)
//     {
//       
//     }
//    exit(0);
   
    switch(RefElement)
     {
      case BFUnitSquare:
          RefTrans = QuadIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TQuadIsoparametric *)F_K)->SetCell(Me);
          ((TQuadIsoparametric *)F_K)->GetOrigBoundFromRef(IJoint,
	                            N_LinePoints, zeta, X_B, Y_B);
      break;

      case BFUnitTriangle:
          RefTrans = TriaIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TTriaIsoparametric *)F_K)->SetApproximationOrder(ORDER);  
          ((TTriaIsoparametric *)F_K)->SetCell(Me);
          ((TTriaIsoparametric *)F_K)->GetOrigBoundFromRef(IJoint, N_LinePoints, zeta, X_B, Y_B);
      break;
    } // endswitch

    TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)
                  ->MakeRefElementData(LineQuadFormula);

    uref = TFEDatabase2D::GetJointValues2D(BaseFuncts[FEId],
                        LineQuadFormula, IJoint);
    uxiref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId],
                        LineQuadFormula, IJoint, D10);
    uetaref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId],
                        LineQuadFormula, IJoint, D01);

    for(k=0;k<N_LinePoints;k++)
     {

         switch(RefElement)
          {
            case BFUnitSquare:
              ((TQuadIsoparametric *)F_K)->SetCell(Me);
              ((TQuadIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
              ((TQuadIsoparametric *) F_K)->GetTangent(IJoint, zeta[k], t0, t1);
            break;

            case BFUnitTriangle:
              ((TTriaIsoparametric *)F_K)->SetApproximationOrder(ORDER);
              ((TTriaIsoparametric *)F_K)->SetCell(Me); 
              ((TTriaIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
              ((TTriaIsoparametric *) F_K)->GetTangent(IJoint, zeta[k], t0, t1);
            break;
          } // endswitch

      normn = sqrt(t0*t0+t1*t1);
      r_axial = fabs(X_B[k]);
     
         if(X_B[k]<=0)
          {
           cout <<"X_B[k] negative in SurfInt  change Quad rule " <<  X_B[k] <<endl;
//         exit(0);
          }


     // get surfactant solution  on the interface
     U=0.;
#ifndef __MASSTRANSTEST__
     for(l=0;l<N_BaseFunct_Surface;l++)
      {
       local_dof   = JointDOF[l];
       test00 = uorig[local_dof];
       m = DOF_Surface[l];
       U  += V[m]*test00;
      }
#endif

//      if( XmaxVal[0]< r_axial )
//       {
//        XmaxVal[0] = r_axial;
//        XmaxVal[1] = U;
//       }

     if( XmaxVal[1]< U )
      {
       XmaxVal[0] = r_axial;
       XmaxVal[1] = U;
      }
      
      rhsval = Bi*Da*U;
      rhsval *= normn*LineWeights[k]*r_axial;

//    update the matrix for all test functions
      for(l=0;l<N_BaseFunct;l++)
       {
        TestDOF = DOF[l];

        if(TestDOF<ActiveBound)
          rhs[TestDOF] += rhsval*uorig[l];

        index2 = RowPtr[TestDOF+1];

        for(m=0;m<N_BaseFunct;m++)
         {
          AnsatzDOF = DOF[m];
          // cout << AnsatzDOF << " -- " << TestDOF << endl;
          index1 = RowPtr[TestDOF];

          if(index1+1 == index2) continue;
           while(KCol[index1] != AnsatzDOF) index1++;

          val =  beta*(1. - U)*uorig[m]*uorig[l];
          val *= normn*LineWeights[k]*r_axial;

          ValuesA[index1] += val;
         } // for(m=0;m<N_BaseFunct;m++)
       }//  for(l=0;l<N_BaseFunct;l++)
    }  // for(k=0;k<N_LinePoints;k

  } // for(i=0;i<N_Cells;

//   cout<< " Flux at the interface " << 2*Pi*surf_flux <<endl;
//   cout << " GammaMass on Interface; " << 2*Pi*GammaMass <<endl;
//   cout << " CMass on Interface; " << 2*Pi*CMass <<endl;
//    cout << " XmaxVal; " << XmaxVal[1] << endl;
// 
// exit(0);
}

void  AssembleSurf1D_Mmat(TFESpace2D *fespace, TFESpace1D *fespaces_low,
                     TSquareMatrix1D *sqmatrices_low, int *Cell_array, int *Joint_array)
{
  int i, j, k, l, m, n, N_Cells_low, N, N_LocalUsedElements, local_i, local_j, ORDER;
  int N_BaseFunct, N_BaseFunct_low,  N_Points, begin, end, *N_BaseFuncts;
  int *BeginIndex_low, *GlobalNumbers_low, *DOF_LOW, TestDOF, AnsatzDOF, IJoint ;
  int *GlobalNumbers_Outer, *BeginIndex_Outer, N_BaseFunct_Outer;
  int LocN_BF[N_BaseFuncts2D], N_LinePoints, *KCol, *RowPtr, *JointDOF_Outer, N_Outer;
  int *DOF_Outer, N_JointDOF_Outer;

  double x0, y0, x1, y1, t0, t1, n0, n1, normn;
  double AbsDetjk[MaxN_QuadPoints_2D], Mult;
  double **uref_Outer, **uxiref_Outer, **uetaref_Outer,  C, CX, CY;
  double *LineWeights, *zeta, *ValuesA, *ValuesM;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double uorig_Outer[MaxN_BaseFunctions2D], uxorig_Outer[MaxN_BaseFunctions2D];
  double uyorig_Outer[MaxN_BaseFunctions2D];
  double  val, X_B[100], Y_B[100], r_axial, d1, d2, e1, e2;
  double LocMatrixM[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double ansatz00, test00;
 
  BaseFunct2D LocBF[N_BaseFuncts2D];
  BaseFunct2D *BaseFuncts;
  BF2DRefElements RefElement;
  RefTrans2D RefTrans;
  TRefTrans2D *F_K;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1;
  TFEDesc2D *FeDesc_Outer;
  TFEDesc1D *FeDesc_low;
  TCollection *Coll, *Coll_low;
  TBaseCell *Me, *Me_low;
  FE2D FEId_Outer;
  FE1D FEId_low;
  TFE1D *Element;
  TFE2D *ele;
 
// // ########################################################################
// // store information in local arrays
// // ########################################################################
  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFuncts = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  Coll = fespace->GetCollection(); // all spaces use same Coll
  GlobalNumbers_Outer = fespace->GetGlobalNumbers();
  BeginIndex_Outer = fespace->GetBeginIndex();

  Coll_low = fespaces_low->GetCollection(); // all low spaces use same Coll
  N_Cells_low = Coll_low->GetN_Cells();
  BeginIndex_low =  fespaces_low->GetBeginIndex();
  GlobalNumbers_low =  fespaces_low->GetGlobalNumbers();

  RowPtr = sqmatrices_low->GetRowPtr();
  KCol = sqmatrices_low->GetKCol();
  ValuesM = sqmatrices_low->GetEntries();

// ########################################################################
// loop over all low space cells
// ########################################################################
  for(i=0;i<N_Cells_low;i++)
  {
    N = Cell_array[i];
    Me = Coll->GetCell(N);
    IJoint = Joint_array[i];

    N_Outer = N; // free surf flows N=N_Outer;
    FEId_Outer = fespace->GetFE2D(N_Outer, Me);  // FEID of surfactant space in the outer domain
    RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId_Outer);    
    FeDesc_Outer = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId_Outer);
    N_JointDOF_Outer = FeDesc_Outer->GetN_JointDOF();
    JointDOF_Outer = FeDesc_Outer->GetJointDOF(IJoint);   
    N_BaseFunct_Outer = FeDesc_Outer->GetN_DOF();
    DOF_Outer = GlobalNumbers_Outer + BeginIndex_Outer[N_Outer];

    DOF_LOW = GlobalNumbers_low + BeginIndex_low[i];
    Me_low = Coll_low->GetCell(i);
    FEId_low = fespaces_low->GetFE1D(i, Me_low);
    Element = TFEDatabase2D::GetFE1D(FEId_low);
    N_BaseFunct_low = Element->GetN_DOF();
    
    if(N_JointDOF_Outer != N_BaseFunct_low )
     {
      cout<< " N_JointDOF_Outer != N_BaseFunct_low " <<endl;
      exit(0);
    }

    memset(LocMatrixM, 0, N_BaseFunct_low*N_BaseFunct_low*SizeOfDouble);
    
    l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId_Outer);
    LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
    qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);
    TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId_Outer)->MakeRefElementData(LineQuadFormula);

    ORDER = TFEDatabase2D::GetAccuracyFromFE2D(FEId_Outer);

    switch(RefElement)
      {
        case BFUnitSquare:
          RefTrans = QuadIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TQuadIsoparametric *)F_K)->SetCell(Me);
          ((TQuadIsoparametric *)F_K)->GetOrigBoundFromRef(IJoint, N_LinePoints, zeta, X_B, Y_B);
        break;

        case BFUnitTriangle:
          RefTrans = TriaIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TTriaIsoparametric *)F_K)->SetApproximationOrder(ORDER);  
          ((TTriaIsoparametric *)F_K)->SetCell(Me);
          ((TTriaIsoparametric *)F_K)->GetOrigBoundFromRef(IJoint, N_LinePoints, zeta, X_B, Y_B);
        break;
      } // endswitch
      uref_Outer = TFEDatabase2D::GetJointValues2D(BaseFuncts[FEId_Outer], LineQuadFormula, IJoint);
      uxiref_Outer = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId_Outer], LineQuadFormula, IJoint, D10);
      uetaref_Outer = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId_Outer], LineQuadFormula, IJoint, D01);
      for(k=0;k<N_LinePoints;k++)
        {
          F_K->GetTangent(IJoint, zeta[k], t0, t1);
          normn = sqrt(t0*t0+t1*t1);
          n0 =  t1/normn;
          n1 = -t0/normn;

          switch(RefElement)
           {
            case BFUnitSquare:
              ((TQuadIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct_Outer, uref_Outer[k], uxiref_Outer[k], uetaref_Outer[k],
                        uorig_Outer, uxorig_Outer, uyorig_Outer);
            break;

            case BFUnitTriangle:
              ((TTriaIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct_Outer, uref_Outer[k], uxiref_Outer[k], uetaref_Outer[k],
                        uorig_Outer, uxorig_Outer, uyorig_Outer);


            break;
           } // endswitch

         r_axial = fabs(X_B[k]); // r value in the axial symmetric integral
         if(X_B[k]<=0)
          {
           cout <<"X_B[k] negative in Assemble1D change Quad rule " <<  X_B[k] <<endl;
//         exit(0);
          }
// //           cout << " x " << r_axial<< " y " << Y_B[k]<< endl;

          Mult = sqrt(t0*t0+t1*t1)*(LineWeights[k])*r_axial;

          for(l=0;l<N_BaseFunct_low;l++)
           {
            local_j   = JointDOF_Outer[l];
            test00  = uorig_Outer[local_j];
   
            for(m=0;m<N_BaseFunct_low;m++)
             {
              local_i   = JointDOF_Outer[m];
              ansatz00 = uorig_Outer[local_i];
	      
	      // mass matrix
              val  = test00*ansatz00;
              val *= Mult;
              LocMatrixM[l*N_BaseFunct_low+m] += val;
            }
          } //  for(l=0;l<N_Joint      
        } //  for(k=0;k<N_

//       for(l=0;l<N_BaseFunct_low;l++)
//        for(m=0;m<N_BaseFunct_low;m++)
//         cout << " LocMatrixM " << LocMatrixM[l*N_BaseFunct_low+m]<<endl;

//   add to global matrices
    for(l=0;l<N_BaseFunct_low;l++)
     {
      TestDOF = DOF_LOW[l];


      begin = RowPtr[TestDOF];
      end = RowPtr[TestDOF+1];
      for(n=begin;n<end;n++)
       {
       for(m=0;m<N_BaseFunct_low;m++)
        {
         if(KCol[n] == DOF_LOW[m])
          {
           ValuesM[n] +=LocMatrixM[l*N_BaseFunct_low+m];
           break;
          }
        } // for(m=0;m<N_BaseFunct_low
      } // for(n=begin;n<end;n++)
     } // for(l=0;l<N_BaseFunct_low
    } //  for(i=0;i<N_Cells_low
 
 
// // exit(0);
}


void  AssembleSurf1D_SolubleSurfact(int n_fespaces, TFESpace2D **fespaces, TFEFunction2D **fefunctions,
                     int N_FESpaces_low, TFESpace1D **fespaces_low, TFEFunction1D *SurfaceFeFunct, int N_SquareMatrices,
                     TSquareMatrix1D **sqmatrices_low, int N_Rhs, double **RHSs, 
                     TFESpace1D **ferhs_low, int *Cell_array, int *Joint_array, double *C_Outer)
{
  int i, j, k, l, m, n, N_Cells_low, N, N_LocalUsedElements, local_i, local_j, ORDER;
  int N_BaseFunct, N_BaseFunct_low,  N_Points, begin, end, *N_BaseFuncts;
  int *BeginIndex_low, *GlobalNumbers_low, *DOF, *DOF_LOW, TestDOF, AnsatzDOF, IJoint ;
  int *BeginIndex, *GlobalNumbers, *GlobalNumbers_Outer, *BeginIndex_Outer, N_BaseFunct_Outer;
  int LocN_BF[N_BaseFuncts2D], N_LinePoints, *KCol, *RowPtr, *JointDOF_Outer, N_Outer;
  int *DOF_Outer, N_JointDOF_Outer;

  double x0, y0, x1, y1, t0, t1, n0, n1, normn;
  double AbsDetjk[MaxN_QuadPoints_2D], Mult;
  double *weights, *xi, *eta;
  double **uref, **uxiref, **uetaref;
  double **uref_Outer, **uxiref_Outer, **uetaref_Outer,  C, CX, CY;
  double *LineWeights, *zeta, *ValuesA, *ValuesM;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double uorig[MaxN_BaseFunctions2D], uxorig[MaxN_BaseFunctions2D], uyorig[MaxN_BaseFunctions2D];
  double uorig_Outer[MaxN_BaseFunctions2D], uxorig_Outer[MaxN_BaseFunctions2D];
  double uyorig_Outer[MaxN_BaseFunctions2D];
  double c0, r2, c1, c2, Pe_s_solid;

  double val, rhsval, theta, ngrad_ansatz, ngrad_test, TangDivU;
  double  X_B[100], Y_B[100], r_axial, d1, d2, e1, e2;
  double LocMatrixA[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixM[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocRhs[MaxN_BaseFunctions2D];
  double ansatz00, ansatz10, ansatz01 ;
  double test00, test10, test01, *u1, *u2, u1x, u2x, u1y, u2y, U1;
  double *RHS, *V, U, AddedMass, CMass, GammaMass;
  double Bi, Da, beta, Pe_s, NGrad_C;
  
  BaseFunct2D LocBF[N_BaseFuncts2D];
  BaseFunct2D *BaseFuncts;
  boolean *SecondDer;
  BF2DRefElements RefElement;
  RefTrans2D RefTrans;
  TRefTrans2D *F_K;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1;
  TFEDesc2D *FeDesc, *FeDesc_Outer;
  TFEDesc1D *FeDesc_low;
  TCollection *Coll, *Coll_low;
  TBaseCell *Me, *Me_low;
  FE2D FEId, FEId_Outer;
  FE1D FEId_low;
  TFE1D *Element;
  TFE2D *ele;

  SecondDer = new boolean[n_fespaces];
// ########################################################################
// store information in local arrays
// ########################################################################
  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFuncts = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  Coll = fespaces[0]->GetCollection(); // all spaces use same Coll
  GlobalNumbers = fespaces[0]->GetGlobalNumbers();
  BeginIndex = fespaces[0]->GetBeginIndex();
  u1 = fefunctions[0]->GetValues();
  u2 = fefunctions[1]->GetValues();

  GlobalNumbers_Outer = fespaces[1]->GetGlobalNumbers();
  BeginIndex_Outer = fespaces[1]->GetBeginIndex();
  V = SurfaceFeFunct->GetValues();

  Coll_low = fespaces_low[0]->GetCollection(); // all low spaces use same Coll
  N_Cells_low = Coll_low->GetN_Cells();
  BeginIndex_low =  fespaces_low[0]->GetBeginIndex();
  GlobalNumbers_low =  fespaces_low[0]->GetGlobalNumbers();

  RowPtr = sqmatrices_low[0]->GetRowPtr();
  KCol = sqmatrices_low[0]->GetKCol();

  ValuesA = sqmatrices_low[0]->GetEntries();
  ValuesM = sqmatrices_low[1]->GetEntries();
  RHS = RHSs[0];

  N_LocalUsedElements = n_fespaces;
  for(j=0;j<n_fespaces;j++)
    SecondDer[j]=FALSE;
  
  Bi = TDatabase::ParamDB->REACTOR_P13; 
  Da = TDatabase::ParamDB->REACTOR_P14;
  beta = TDatabase::ParamDB->REACTOR_P15;
  
  Pe_s = TDatabase::ParamDB->REACTOR_P17; // Peclet number
  if(Pe_s==0.)
  { c1=0.; }
  else
  { c1 =  1./Pe_s; }
  
  Pe_s_solid = TDatabase::ParamDB->REACTOR_P22; // Peclet number
  if(Pe_s_solid==0.)
  { c2=0.; }
  else
  { c2 = 1./Pe_s_solid; }
  
//   AddedMass = 0.;
//   cout << " N_Cells_low " <<  N_Cells_low << endl;
//   exit(0);
  
// ########################################################################
// loop over all low space cells
// ########################################################################
  for(i=0;i<N_Cells_low;i++)
  {
    N = Cell_array[i];
    Me = Coll->GetCell(N);
    IJoint = Joint_array[i];

    // on liquid-solid interface, Pe is different
    if((Me->GetJoint(IJoint))->GetType() == BoundaryEdge)
    { c0=c2; }
    else
    { c0=c1; }
    
    FEId = fespaces[0]->GetFE2D(N, Me);  // FEID of velocity space in the outer domain
    ele = TFEDatabase2D::GetFE2D(FEId);
    RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
    FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
    N_BaseFunct = FeDesc->GetN_DOF();
    DOF = GlobalNumbers + BeginIndex[N];

    N_Outer = N; // free surf flows N=N_Outer;
    FEId_Outer = fespaces[1]->GetFE2D(N_Outer, Me);  // FEID of surfactant space in the outer domain
    FeDesc_Outer = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId_Outer);
    N_JointDOF_Outer = FeDesc_Outer->GetN_JointDOF();
    JointDOF_Outer = FeDesc_Outer->GetJointDOF(IJoint);   
    N_BaseFunct_Outer = FeDesc_Outer->GetN_DOF();
    DOF_Outer = GlobalNumbers_Outer + BeginIndex_Outer[N_Outer];

    DOF_LOW = GlobalNumbers_low + BeginIndex_low[i];
    Me_low = Coll_low->GetCell(i);
    FEId_low = fespaces_low[0]->GetFE1D(i, Me_low);
    Element = TFEDatabase2D::GetFE1D(FEId_low);
    N_BaseFunct_low = Element->GetN_DOF();
    
    if(N_JointDOF_Outer != N_BaseFunct_low )
     {
      cout<< " N_JointDOF_Outer != N_BaseFunct_low " <<endl;
      exit(0);
    }

    memset(LocMatrixA, 0, N_BaseFunct_low*N_BaseFunct_low*SizeOfDouble);
    memset(LocMatrixM, 0, N_BaseFunct_low*N_BaseFunct_low*SizeOfDouble);
    memset(LocRhs, 0, N_BaseFunct_low*SizeOfDouble);

    l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
    LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
    qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);

    TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(LineQuadFormula);
    TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId_Outer)->MakeRefElementData(LineQuadFormula);

    ORDER = TFEDatabase2D::GetAccuracyFromFE2D(FEId);
    switch(RefElement)
      {
        case BFUnitSquare:
          RefTrans = QuadIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TQuadIsoparametric *)F_K)->SetCell(Me);
          ((TQuadIsoparametric *)F_K)->GetOrigBoundFromRef(IJoint, N_LinePoints, zeta, X_B, Y_B);
        break;

        case BFUnitTriangle:
          RefTrans = TriaIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TTriaIsoparametric *)F_K)->SetApproximationOrder(ORDER);  
          ((TTriaIsoparametric *)F_K)->SetCell(Me);
          ((TTriaIsoparametric *)F_K)->GetOrigBoundFromRef(IJoint, N_LinePoints, zeta, X_B, Y_B);
        break;
      } // endswitch


      uref_Outer = TFEDatabase2D::GetJointValues2D(BaseFuncts[FEId_Outer],
                     LineQuadFormula, IJoint);
      uxiref_Outer = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId_Outer],
                        LineQuadFormula, IJoint, D10);
      uetaref_Outer = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId_Outer],
                        LineQuadFormula, IJoint, D01);

      uref = TFEDatabase2D::GetJointValues2D(BaseFuncts[FEId],
                        LineQuadFormula, IJoint);
      uxiref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId],
                        LineQuadFormula, IJoint, D10);
      uetaref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId],
                        LineQuadFormula, IJoint, D01);

      for(k=0;k<N_LinePoints;k++)
        {
          F_K->GetTangent(IJoint, zeta[k], t0, t1);
          normn = sqrt(t0*t0+t1*t1);
          n0 =  t1/normn;
          n1 = -t0/normn;

          switch(RefElement)
          {
            case BFUnitSquare:
              ((TQuadIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
              ((TQuadIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct_Outer, uref_Outer[k], uxiref_Outer[k], uetaref_Outer[k],
                        uorig_Outer, uxorig_Outer, uyorig_Outer);
            break;

            case BFUnitTriangle:
              ((TTriaIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
              ((TTriaIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct_Outer, uref_Outer[k], uxiref_Outer[k], uetaref_Outer[k],
                        uorig_Outer, uxorig_Outer, uyorig_Outer);


            break;
          } // endswitch

         r_axial = fabs(X_B[k]); // r value in the axial symmetric integral
         if(X_B[k]<=0)
          {
           cout <<"X_B[k] negative in Assemble1D change Quad rule " <<  X_B[k] <<endl;
//         exit(0);
          }
//           cout << " x " << r_axial<< " y " << Y_B[k]<< endl;

//       value of C at this integral point
         C = 0.;
#ifndef __MASSTRANSTEST__
         for(l=0;l<N_BaseFunct_Outer;l++)
           C += C_Outer[DOF_Outer[l]]*uorig_Outer[l];
#endif
//           cout <<" C:  " << C <<endl;

//       get surfactant solution  on the interface
        U=0.;
#ifndef __MASSTRANSTEST__
        for(l=0;l<N_BaseFunct_low;l++)
         {
          local_j   = JointDOF_Outer[l];
          test00 = uorig[local_j];
          m = DOF_LOW[l];
          U  += V[m]*test00;
         }
#endif
         //         cout <<" U:  out " << U <<endl;

          //       get velocity gradients
          U1 = 0.;  u1x=0.; u2x=0.; u1y=0.; u2y=0.;
// #ifndef __MASSTRANSTEST__
          for(l=0;l<N_BaseFunct;l++)
            {
             m = DOF[l];
             U1  += u1[m]*uorig[l];
             u1x += u1[m]*uxorig[l];
             u1y += u1[m]*uyorig[l];
             u2x += u2[m]*uxorig[l];
             u2y += u2[m]*uyorig[l];
            }
// #endif
          TangDivU =  u1x - (u1x*n0 + u1y*n1)*n0  + U1/r_axial
                    + u2y - (u2x*n0 + u2y*n1)*n1;

//           cout <<" u1x:  " << u1x <<" u1y:  " << u1y <<endl;
//           cout <<" u2x:  " << u2x <<" u2y:  " << u2y <<endl;
//           cout <<" TangDivU:  " << TangDivU <<endl;

          Mult = sqrt(t0*t0+t1*t1)*(LineWeights[k])*r_axial;
#ifndef __MASSTRANSTEST__
          rhsval = (beta/Da)*C*(1. - U)  -  Bi*U;
          rhsval *= Mult;  
          U *= Mult;
          C *= Mult;  
#endif
          for(l=0;l<N_BaseFunct_low;l++)
           {
            local_j   = JointDOF_Outer[l];

            test00  = uorig_Outer[local_j];
            test10  = uxorig_Outer[local_j];
            test01  = uyorig_Outer[local_j];

            ngrad_test= n0*test10 + n1*test01;
            d1 = test10 - ngrad_test*n0;
            d2 = test01 - ngrad_test*n1;
#ifndef __MASSTRANSTEST__
//          rhs
            LocRhs[l] += rhsval*test00; // soluble surfactant relation explicit (in C) form
#endif          
            for(m=0;m<N_BaseFunct_low;m++)
             {
              local_i   = JointDOF_Outer[m];

              ansatz00 = uorig_Outer[local_i];
              ansatz10 = uxorig_Outer[local_i];
              ansatz01 = uyorig_Outer[local_i];

//              cout << local_i << " -- " << local_j << endl;
              ngrad_ansatz= n0*ansatz10 + n1*ansatz01;
              e1 = ansatz10 - ngrad_ansatz*n0;
              e2 = ansatz01 - ngrad_ansatz*n1;

//              cout << " Tgrad . n  " << e1*n0 + e2*n1 << endl;
              val = c0*(d1*e1 + d2*e2);
// #ifndef __MASSTRANSTEST__
              val +=TangDivU*test00*ansatz00;
// #endif
              val *= Mult;
              LocMatrixA[l*N_BaseFunct_low+m] += val;

//             cout << local_i << " -- " << local_j<< " -- " << val << endl;

              val  = test00*ansatz00;
              val *= Mult;
              LocMatrixM[l*N_BaseFunct_low+m] += val;
            }
          } //  for(l=0;l<N_Joint
          
//         AddedMass +=Mult;         
        } //  for(k=0;k<N_

//       for(l=0;l<N_BaseFunct_low;l++)
//        for(m=0;m<N_BaseFunct_low;m++)
//         cout << " LocMatrixA " << LocMatrixA[l*N_BaseFunct_low+m]<<endl;

//   add to global matrices
    for(l=0;l<N_BaseFunct_low;l++)
     {
      TestDOF = DOF_LOW[l];
      RHS[TestDOF] += LocRhs[l]; // soluble surfactant relation

      begin = RowPtr[TestDOF];
      end = RowPtr[TestDOF+1];
      for(n=begin;n<end;n++)
       {
       for(m=0;m<N_BaseFunct_low;m++)
        {
         if(KCol[n] == DOF_LOW[m])
          {
           ValuesA[n] +=LocMatrixA[l*N_BaseFunct_low+m];
           ValuesM[n] +=LocMatrixM[l*N_BaseFunct_low+m];
           break;
          }
        } // for(m=0;m<N_BaseFunct_low
      } // for(n=begin;n<end;n++)
     } // for(l=0;l<N_BaseFunct_low
    } //  for(i=0;i<N_Cells_low

//   cout << " AddedMass on Interface; " << 2*Pi*AddedMass <<endl;
//   cout << " GammaMass at the Interface; " << 2*Pi*GammaMass <<endl;
//   cout << " CMass at the Interface; " << 2*Pi*CMass <<endl;
 delete [] SecondDer;
// exit(0);
}

void  AssembleSurf1D_SolubleSurfact_Implicit(int n_fespaces, TFESpace2D **fespaces, 
                     TFEFunction2D **fefunctions, int N_FESpaces_low, TFESpace1D **fespaces_low,
                     double *old_c_Gamma, int N_SquareMatrices,
                     TSquareMatrix1D **sqmatrices_low, int N_Rhs, double **RHSs, 
                     TFESpace1D **ferhs_low, int *Cell_array, int *Joint_array, double *C_Outer, double dt)
{
  int i, j, k, l, m, n, N_Cells_low, N, N_LocalUsedElements, local_i, local_j, ORDER;
  int N_BaseFunct, N_BaseFunct_low,  N_Points, begin, end, *N_BaseFuncts;
  int *BeginIndex_low, *GlobalNumbers_low, *DOF, *DOF_LOW, TestDOF, AnsatzDOF, IJoint ;
  int *BeginIndex, *GlobalNumbers, *GlobalNumbers_Outer, *BeginIndex_Outer, N_BaseFunct_Outer;
  int LocN_BF[N_BaseFuncts2D], N_LinePoints, *KCol, *RowPtr, *JointDOF_Outer, *JointDOF, N_Outer;
  int *DOF_Outer, count=0, N_JointDOF_Outer, N_JointDOF;

  double x0, y0, x1, y1, t0, t1, n0, n1, normn;
  double AbsDetjk[MaxN_QuadPoints_2D], Mult;
  double *weights, *xi, *eta;
  double **uref, **uxiref, **uetaref;
  double **uref_Outer, **uxiref_Outer, **uetaref_Outer,  C, CX, CY;
  double *LineWeights, *zeta, *ValuesA, *ValuesM;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double uorig[MaxN_BaseFunctions2D], uxorig[MaxN_BaseFunctions2D], uyorig[MaxN_BaseFunctions2D];
  double uorig_Outer[MaxN_BaseFunctions2D], uxorig_Outer[MaxN_BaseFunctions2D];
  double uyorig_Outer[MaxN_BaseFunctions2D];
  double c0, r2, c1, c2, Pe_s_solid;
  double Pr = TDatabase::ParamDB->PR_NR;
  double val, rhsval, theta, ngrad_ansatz, ngrad_test, TangDivU;
  double  X_B[100], Y_B[100], r_axial, d1, d2, e1, e2;
  double LocMatrixA[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixM[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocRhs[MaxN_BaseFunctions2D];
  double ansatz00, ansatz10, ansatz01 ;
  double test00, test10, test01, *u1, *u2, u1x, u2x, u1y, u2y, U1, U2;
  double *RHS, U, AddedMass, CMass, GammaMass;
  double Bi, Da, beta, Pe_s, NGrad_C, maxTangDivU=0;
  
  bool IsIsoparametric;
  
  BaseFunct2D LocBF[N_BaseFuncts2D];
  BaseFunct2D *BaseFuncts;
//   boolean *SecondDer;
  BF2DRefElements RefElement;
  RefTrans2D RefTrans;
  TRefTrans2D *F_K;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1;
  TFEDesc2D *FeDesc, *FeDesc_Outer;
  TFEDesc1D *FeDesc_low;
  TCollection *Coll, *Coll_low;
  TBaseCell *Me, *Me_low;
  FE2D FEId, FEId_Outer;
  FE1D FEId_low;
  TFE1D *Element;
//   TFE2D *ele;
  TJoint *joint;
  
//   SecondDer = new boolean[n_fespaces];
// ########################################################################
// store information in local arrays
// ########################################################################
  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFuncts = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  Coll = fespaces[0]->GetCollection(); // all spaces use same Coll
  GlobalNumbers = fespaces[0]->GetGlobalNumbers(); // velo space
  BeginIndex = fespaces[0]->GetBeginIndex();
  u1 = fefunctions[0]->GetValues();
  u2 = fefunctions[1]->GetValues();

  GlobalNumbers_Outer = fespaces[1]->GetGlobalNumbers(); // surfactant space
  BeginIndex_Outer = fespaces[1]->GetBeginIndex();

//   V = SurfaceFeFunct->GetValues();

  Coll_low = fespaces_low[0]->GetCollection(); // all low spaces use same Coll
  N_Cells_low = Coll_low->GetN_Cells();
  BeginIndex_low =  fespaces_low[0]->GetBeginIndex();
  GlobalNumbers_low =  fespaces_low[0]->GetGlobalNumbers();

  RowPtr = sqmatrices_low[0]->GetRowPtr();
  KCol = sqmatrices_low[0]->GetKCol();

  ValuesA = sqmatrices_low[0]->GetEntries();
  ValuesM = sqmatrices_low[1]->GetEntries();
  RHS = RHSs[0];

  N_LocalUsedElements = n_fespaces;
//   for(j=0;j<n_fespaces;j++)
//     SecondDer[j]=FALSE;
  
  Bi = TDatabase::ParamDB->REACTOR_P13; 
  Da = TDatabase::ParamDB->REACTOR_P14;
  beta = TDatabase::ParamDB->REACTOR_P15;
  
 
  Pe_s = TDatabase::ParamDB->REACTOR_P17; // Peclet number
  if(Pe_s==0.)
  { c1=0.; }
  else
  { c1 =  1./Pe_s; }
  
  Pe_s_solid = TDatabase::ParamDB->REACTOR_P22; // Peclet number
  if(Pe_s_solid==0.)
  { c2=0.; }
  else
  { c2 = 1./Pe_s_solid; }
//   AddedMass = 0.;
//   CMass =0.;
// ########################################################################
// loop over all low space cells
// ########################################################################
  for(i=0;i<N_Cells_low;i++)
  {
    N = Cell_array[i];
    Me = Coll->GetCell(N);
    IJoint = Joint_array[i];

    // on liquid-solid interface, Pe is different
    if((Me->GetJoint(IJoint))->GetType() == BoundaryEdge)
    { c0=c2; }
    else
    { c0=c1; }
    
    
    FEId = fespaces[0]->GetFE2D(N, Me);  // FEID of velocity space in the outer domain
//     ele = TFEDatabase2D::GetFE2D(FEId);

    FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
    N_JointDOF = FeDesc->GetN_JointDOF();
    JointDOF = FeDesc->GetJointDOF(IJoint);    
    
    N_BaseFunct = FeDesc->GetN_DOF();
    DOF = GlobalNumbers + BeginIndex[N];

//     N_Outer = Me->GetLocalCellNo();
    N_Outer = N; // free surf flows N=N_Outer;    
    FEId_Outer = fespaces[1]->GetFE2D(N_Outer, Me);  // FEID of surfactant space in the outer domain
    RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId_Outer);
    FeDesc_Outer = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId_Outer);
    N_JointDOF_Outer = FeDesc_Outer->GetN_JointDOF();
    JointDOF_Outer = FeDesc_Outer->GetJointDOF(IJoint);    
    N_BaseFunct_Outer = FeDesc_Outer->GetN_DOF();
    DOF_Outer = GlobalNumbers_Outer + BeginIndex_Outer[N_Outer];

    DOF_LOW = GlobalNumbers_low + BeginIndex_low[i];
    Me_low = Coll_low->GetCell(i);
    FEId_low = fespaces_low[0]->GetFE1D(i, Me_low);
    Element = TFEDatabase2D::GetFE1D(FEId_low);
    N_BaseFunct_low = Element->GetN_DOF();
    if(N_JointDOF_Outer != N_BaseFunct_low )
     {
      cout<< " N_JointDOF_Outer != N_BaseFunct_low " <<endl;
      exit(0);
    }

    memset(LocMatrixA, 0, N_BaseFunct_low*N_BaseFunct_low*SizeOfDouble);
    memset(LocMatrixM, 0, N_BaseFunct_low*N_BaseFunct_low*SizeOfDouble);
    memset(LocRhs, 0, N_BaseFunct_low*SizeOfDouble);

    l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId_Outer);
    LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
    qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);

    TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(LineQuadFormula);
    TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId_Outer)->MakeRefElementData(LineQuadFormula);

    ORDER = TFEDatabase2D::GetAccuracyFromFE2D(FEId_Outer);
    
    IsIsoparametric = FALSE;

    if (TDatabase::ParamDB->USE_ISOPARAMETRIC)
     {
      joint = Me->GetJoint(IJoint);

      if(joint->GetType() == IsoBoundEdge   || joint->GetType() == IsoInterfaceJoint
           || joint->GetType() == InterfaceJoint)
      {
       IsIsoparametric=TRUE; 
      }
     }

    switch(RefElement)
      {
        case BFUnitSquare:
          RefTrans = QuadIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TQuadIsoparametric *)F_K)->SetCell(Me);
          ((TQuadIsoparametric *)F_K)->GetOrigBoundFromRef(IJoint, N_LinePoints, zeta, X_B, Y_B);
        break;

        case BFUnitTriangle:
         if(IsIsoparametric)
	 {
          RefTrans = TriaIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TTriaIsoparametric *)F_K)->SetApproximationOrder(ORDER);  
          ((TTriaIsoparametric *)F_K)->SetCell(Me);
          ((TTriaIsoparametric *)F_K)->GetOrigBoundFromRef(IJoint, N_LinePoints, zeta, X_B, Y_B);
	 }
	 else
	 {
          RefTrans = TriaAffin; 
// 	  cout <<"TTriaAffin setting Cell IJoint " << IJoint << endl;
// 	  
// 	  
// 	   Me->GetVertex(0)->GetCoords(x0, y0);
// 	   cout<< "x " << x0 << " y " << y0 <<endl;
// 	   Me->GetVertex(1)->GetCoords(x0, y0);
// 	   cout<< "x " << x0 << " y " << y0 <<endl;
// 	   Me->GetVertex(2)->GetCoords(x0, y0);
// 	   cout<< "x " << x0 << " y " << y0 <<endl;
	   
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TTriaAffin *)F_K)->SetCell(Me); 
          ((TTriaAffin *)F_K)->GetOrigBoundFromRef(IJoint, N_LinePoints, zeta, X_B, Y_B); 
         }
 
        break;
      } // endswitch

      uref_Outer = TFEDatabase2D::GetJointValues2D(BaseFuncts[FEId_Outer], LineQuadFormula, IJoint);
      uxiref_Outer = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId_Outer], LineQuadFormula, IJoint, D10);
      uetaref_Outer = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId_Outer], LineQuadFormula, IJoint, D01);

      uref = TFEDatabase2D::GetJointValues2D(BaseFuncts[FEId], LineQuadFormula, IJoint);
      uxiref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId], LineQuadFormula, IJoint, D10);
      uetaref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId], LineQuadFormula, IJoint, D01);

      for(k=0;k<N_LinePoints;k++)
        {
          F_K->GetTangent(IJoint, zeta[k], t0, t1);
          normn = sqrt(t0*t0+t1*t1);
          n0 =  t1/normn;
          n1 = -t0/normn;

          switch(RefElement)
          {
            case BFUnitSquare:
              ((TQuadIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k], uorig, uxorig, uyorig);
              ((TQuadIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct_Outer, uref_Outer[k], uxiref_Outer[k], uetaref_Outer[k],
                        uorig_Outer, uxorig_Outer, uyorig_Outer);
            break;

            case BFUnitTriangle:

             if(IsIsoparametric)
              {      
               ((TTriaIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
               ((TTriaIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct_Outer, uref_Outer[k], uxiref_Outer[k], uetaref_Outer[k],
                        uorig_Outer, uxorig_Outer, uyorig_Outer);
              }
             else              
	     {
               ((TTriaAffin *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);     
	       
               ((TTriaAffin *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct_Outer, uref_Outer[k], uxiref_Outer[k], uetaref_Outer[k],
                        uorig_Outer, uxorig_Outer, uyorig_Outer);
/*
              cout<< " TriaAffin     " <<endl;   
	      exit(0);*/
             }
            break;
          } // endswitch

         r_axial = fabs(X_B[k]); // r value in the axial symmetric integral
         if(X_B[k]<=0)
          {
           cout <<"X_B[k] negative in Assemble1D change Quad rule " <<  X_B[k] <<endl;
//         exit(0);
          }
//           cout << " x " << r_axial<< " y " << Y_B[k]<< endl;

         // value of C at this integral point
         C = 0.;
         for(l=0;l<N_BaseFunct_Outer;l++)
           C += C_Outer[DOF_Outer[l]]*uorig_Outer[l];
 
//           cout <<" C:  " << C <<endl;
	 
// //       get surfactant solution  on the interface
//          U=0.;
//         for(l=0;l<N_BaseFunct_low;l++)
//          {
//           local_j   = JointDOF_Outer[l];
//           test00 = uorig[local_j];
//           m = DOF_LOW[l];
//           U  += old_c_Gamma[m]*test00;
//          }

          // get velocity gradients
          U1 = 0.;  u1x=0.; u2x=0.; u1y=0.; u2y=0.; U2=0;
          for(l=0;l<N_BaseFunct;l++)
            {
	     m = DOF[l];
             U1  += u1[m]*uorig[l];
	     U2  += u2[m]*uorig[l];
             u1x += u1[m]*uxorig[l];
             u1y += u1[m]*uyorig[l];
             u2x += u2[m]*uxorig[l];
             u2y += u2[m]*uyorig[l];

//            if(!IsIsoparametric)
// 	      cout <<" u2[m]:  " << u2[m] <<" uxorig[l]:  " << uxorig[l] << endl;    
            }

          TangDivU =  u1x - (u1x*n0 + u1y*n1)*n0  + U1/r_axial + u2y - (u2x*n0 + u2y*n1)*n1;
 
//         if(!SolidJoint)
//             if(maxTangDivU>(TangDivU) ) maxTangDivU= TangDivU; 
//            if(!IsIsoparametric)
// 	   {
// 	     cout   <<endl;
// 	     cout << " cell " << i << " x " << r_axial<< " y " << Y_B[k];
//              cout <<" U1:  " << U1<<" U2:  " << U2 << endl;    
//              cout <<" u1x:  " << u1x <<" u1y:  " << u1y <<endl;
//              cout <<" u2x:  " << u2x <<" u2y:  " << u2y <<endl;
//              cout <<u1x + U1/r_axial <<  " TangDivU:  " << TangDivU <<endl;
// 	     cout   <<endl;
// // 	     exit(0);
// 	   }
//          if(maxTangDivU<-100)
// 	  {  
//            cout << " maxTangDivU on Interface; " << maxTangDivU <<endl;
//             exit(0);   
//            
// 	  }
	  
          Mult = sqrt(t0*t0+t1*t1)*LineWeights[k]*r_axial;
          rhsval =   (beta/Da)*C;
//           AddedMass +=rhsval;
//           CMass +=Mult*C;
	  
// 	  if(!IsIsoparametric)
//           cout<<  " x " << r_axial<< " y " << Y_B[k] << " U*TangDivU " << U*TangDivU <<endl; 
	  
	  
          for(l=0;l<N_BaseFunct_low;l++)
           {
            local_j   = JointDOF_Outer[l];

            test00  = uorig_Outer[local_j];
            test10  = uxorig_Outer[local_j];
            test01  = uyorig_Outer[local_j];

            ngrad_test= n0*test10 + n1*test01;
            d1 = test10 - ngrad_test*n0;
            d2 = test01 - ngrad_test*n1;
	    
// 	  if(SolidJoint)
//             cout<< " test   " << test00 ;    
	    
//          rhs
// 	    rhsval += dt*TangDivU*U;
            LocRhs[l] += rhsval*Mult*test00; 

            for(m=0;m<N_BaseFunct_low;m++)
             {
              local_i   = JointDOF_Outer[m];

              ansatz00 = uorig_Outer[local_i];
              ansatz10 = uxorig_Outer[local_i];
              ansatz01 = uyorig_Outer[local_i];

//              cout << local_i << " -- " << local_j << endl;
              ngrad_ansatz= n0*ansatz10 + n1*ansatz01;
              e1 = ansatz10 - ngrad_ansatz*n0;
              e2 = ansatz01 - ngrad_ansatz*n1;

#ifndef __MASSTRANSTEST__
//              cout << " Tgrad . n  " << e1*n0 + e2*n1 << endl;
              val = c0*(d1*e1 + d2*e2);
              val += ((beta/Da)*C  +  Bi)*ansatz00*test00;
              val += TangDivU*test00*ansatz00;
              val *=  Mult;
              LocMatrixA[l*N_BaseFunct_low+m] += val;
// 	        cout << " LocMatrixA   val  " << val  ;
#endif   
              val  = test00*ansatz00;
              val *=  Mult ;
              LocMatrixM[l*N_BaseFunct_low+m] += val; 
// 	      cout << " LocMatrixM   val  " << val << endl;     
            }
          } //  for(l=0;l<N_Joint
// 	  if(SolidJoint)          
//           cout<<endl;
        } //  for(k=0;k<N_

//       for(l=0;l<N_BaseFunct_low;l++)
//        for(m=0;m<N_BaseFunct_low;m++)
//         cout << " LocMatrixA " << LocMatrixA[l*N_BaseFunct_low+m]<<endl;

//   add to global matrices
    for(l=0;l<N_BaseFunct_low;l++)
     {
      TestDOF = DOF_LOW[l];
      RHS[TestDOF] += LocRhs[l]; // soluble surfactant relation

      begin = RowPtr[TestDOF];
      end = RowPtr[TestDOF+1];
      for(n=begin;n<end;n++)
       {
       for(m=0;m<N_BaseFunct_low;m++)
        {
         if(KCol[n] == DOF_LOW[m])
          {
           ValuesA[n] +=LocMatrixA[l*N_BaseFunct_low+m];
           ValuesM[n] +=LocMatrixM[l*N_BaseFunct_low+m];
 
//            cout << TestDOF  << ", " << DOF_LOW[m] <<" M: " 
//            << LocMatrixM[l*N_BaseFunct_low+m] <<endl;   
           break;
          }
        } // for(m=0;m<N_BaseFunct_low
      } // for(n=begin;n<end;n++)
     } // for(l=0;l<N_BaseFunct_low
    } //  for(i=0;i<N_Cells_low
    
//   cout << " maxTangDivU on Interface; " << maxTangDivU <<endl;
//   cout << " AddedMass on Interface; " << 2*Pi*AddedMass <<endl;
//   cout << " GammaMass at the Interface; " << 2*Pi*GammaMass <<endl;
//   cout << " CMass at the Interface; " << 2*Pi*CMass <<endl;
//  delete [] SecondDer;
// exit(0);
}

// ======================================================================
// utilities for main program
// ======================================================================

int main(int argc, char* argv[])
{
  TDomain *Domain;
  TDatabase *Database = new TDatabase();
  TFEDatabase2D *FEDatabase = new TFEDatabase2D();
  TCollection *coll, *NSE_coll, *mortarcoll = NULL;
  TBaseCell *Me, *cell, **Free_Cells, **NSE_Cells;
  TGridCell **DelCell;
  TFESpace2D *velocity_space, *pressure_space, *streamfunction_space, *convolution_space, *fesps;
  TFESpace2D  *velocity_space_output, *pressure_space_output;
  TFESpace2D *Grid_space, *vorticity_space, *thermal_space,*grid_space;

  TOutput2D *Output;
  TFEVectFunct2D *RefGridPos, *AuxGridPos, *GridPos, *ReparamPos;
  TFEFunction2D *fefct[4];
  TFESpace2D *fesp[3], *ferhs_T[3], *ferhs[2];
  TAuxParam2D *aux;
  TDiscreteForm2D *DiscreteFormMatrixT_MRhs;
  TDiscreteForm2D *DiscreteForm;
  TMatrix2D *MATRICES[4];
  TSquareMatrix2D *SQMATRICES[8], *SQMATRICES_GRID[4];
  TDiscreteForm2D *DiscreteFormGalerkin;
  TDiscreteForm2D *DiscreteFormNLGalerkin;
  TDiscreteForm2D *DiscreteFormGrid, *DiscreteFormHeat;
  TSquareMatrix **sqmatrices = (TSquareMatrix **)SQMATRICES;
  TMatrix **matrices = (TMatrix **)MATRICES;
  MatVecProc *MatVect;
  DefectProc *Defect;

  TFESpace2D **FESpaces_All = new TFESpace2D *[6];      
  TFEFunction2D **FEFunctions_All = new TFEFunction2D *[8];    
  TFEVectFunct2D **FEVectFuncts_All = new TFEVectFunct2D*[2];
  TStructure2D **Structure_All = new TStructure2D *[2];
  TSquareStructure2D **SquareStructure_All = new TSquareStructure2D *[4];
  TSquareMatrix2D **SqMat_All = new TSquareMatrix2D *[18];
  TMatrix2D **Mat_All = new TMatrix2D *[4];
  
#if defined(__SURFACT__) || defined(__SURFACTCONVTEST__)
  // variables for surfactant
  TDomain *IFaceDomain = new TDomain();  
  
  int **N_List = new int*[4], N_IFaceCells, N_FESpaces_low, N_IActive;   
  int N_S, N_surfactDOF, N_surfactActive, N_surfactNonActive,N_IsurfactDOF;
  int Max_It_scalar, N_BData_S=1, N_BData_F=1, surf_couple_var, N_E, N_Surf_M_MatValues;
  
  double *C_defect, *Csol_old, *oldsol, *C_B, *CRhs_old, *Surf_M_MatValues;
  double *I_defect, *Isol_old, *I_B, *IRhs_old, *Csol_nonlinearstep;
  double  *CRHSs[1], GammaXmaxVal[2];
  double Surf_Mass[2], Initial_SurfactMass, Initial_IFaceSurfactMass;
  double residual_scalar, oldresidual_scalar, oldresidual_surfscalar, limit_scalar, *SRHSs[1];
  double errors[7], Terrors[2],  T_inftyL2, T_inftyH1, max_r;  
  double l_2_l_2u, T_inftyL2_time, olderror;
  
  
  FE1D *FE1D_List;
  TCollection *IFace_Coll;
  TFESpace1D *IFaceSurfact_space;
  TFESpace1D **IFaceFeSpaces = new TFESpace1D*[2];
  TSquareMatrix1D *SQMATRICES_IFace[2];
  TFESpace1D *IFacefesp[1], *IFaceferhs[1];  
  TFESpace2D *surfact_space;
  TSquareMatrix1D **SqMat_IFace = new TSquareMatrix1D*[2];
  TSquareStructure1D **IFaceStruct = new TSquareStructure1D*[1];
  TFEFunction1D **IFaceFeFunct = new TFEFunction1D*[1];  
  TSquareMatrix2D *SQMATRICES_SURFACT[4];  
  TDiscreteForm2D *DiscreteFormSurfact, *DiscreteFormSurfact_SUPG;  
  BoundCondFunct2D *SurfactBoundaryConditions[1];
  BoundValueFunct2D *SurfactBoundValues[1];  
  
  TSystemTCD2D *SystemMatrix;    
#endif  


    
  double total_time,*Coordinates;
  double  t, teta, dt,x,y,gamma, tx,ty,sx,sy, R_Theta[3];
  double left, right, top, bottom,T_a, T_b;
  double x0, y0,x1,y1,hi, residual, impuls_residual, oldresidual, solver_time;
  double *oldsol_T, sphericity, radius, KE;
  double end_time, t1, t2, t4, t3;
  double *B, *defect;
  double *RHSs[3], *refpos, *auxpos, *pos, *ReparamDisp, *ReparamMeshVelo;
  double  TX[2], TY[2], solver_time_curr;
  double SLPX, SLPY, *Entries[4], tau, oldtau, limit, *sol_output, Params[10], InitVolume, CurrVolume;  
  double Lx, Ly, Rx, Ry,  x2, y2, x3, y3, x4, y4, fh, fhlimit, fhtot, fhmin, fhmax;
  double *Angle = new double[2];  
  double **Sol_All = new double *[6], *tmp_Gsol, *tmp_Gd;
  double **Rhs_All = new double *[6], MaxWetD=0., T_MaxWetD=0., h_interface;
  double MaxApexHeight=0, ApexHeight=0, old_c_mass, old_c_Gamma_mass;
  
  int ret,N_RootCells,N_Cells,N_Joints, N_Vertices,N_G, N_NSE_Cells, N_NonNSE_Cells;
  int N_SlipBound_Vert,  N_FreeBound_Vert,  N_AxialBound_Vert,N_Interf_Vertices;
  int In_Index,CurrComp,CurrVertex, img=1, RemeshImg=1;
  int ID,CurrNeib,Neib[2], N_SquareMatrices, N_RectMatrices;
  int a,b,i,X,Y,j,k,l,len1, len2,Neighb_tmp,Indextemp;
  int *PartMarker, *Triangles,comp, *NSE_GlobalCllNo;
  int *PointNeighb,maxEpV = 0, Max_It, N_U_output, N_P_output;
  int  N_U, N_P,N_Unknowns,N_pressureDOF,N_Rhs,N_FESpaces;
  int  N_Hori1, N_Hori2,N_Verti,N_Boundary_Vert,N_thermalDOF,N_thermalActive,N_thermalNonActive;
  int velocity_space_code, pressure_space_code;
  int m,m0,m1,m2,m3,m4,m5,m6, N_Active, N_LinIterCurr, N_LinIter;
  int ORDER,VSP, N_GActive, N_GBoundaryNodes, N_SubSteps, very_first_time=1 ;
  int **IsoCellEdgeNos, *GridKCol, *GridRowPtr, *RowPtr, last_sq;
  int  *JointDOF, N_DOF, dof, *DOF, *GlobalNumbers, *BeginIndex, N_ReParam=0, N_Remesh=0;
  
  
  char *PRM, *GEO, *PsBaseName, *VtkBaseName;
  char ReadinDat[] = "readin.dat";
  char TString[] = "T";
  char NameString[]  = "name";
  char UString[] = "u";
  char PString[] = "p";
  char WString[] = "w";
  const char vtkdir[] = "VTK";
  const char BDdir[] = "BDData";
  char CString[] = "C";
  char IFaceSString[] = "C_I";
  char Surf_F[] = "Gamma_F"; 
  char Surf_S[] = "Gamma_S"; 
  
  bool remeshed=FALSE, reparam = FALSE, SurfSurfactReparam=FALSE, AdaptTimeStep=TRUE, SurfactMassCorrection=FALSE;

  TBoundPart *BoundPart;
  TBoundComp *BoundComp;
  TBdLine *UpdateSlipBound, *UpdateAxialBound;
  TBdCircle *UpdateFreeBound;
  TBaseCell **CellTree;
  TVertex *vert, **VertexDel,**NewVertices;
  TJoint *Joint;
  TBoundEdge *Solid_Joint;
  TBoundEdge ***Bound_Joint = new TBoundEdge**[2];
  TIsoBoundEdge **Free_Joint, *IsoJoint;
  TVertex ***MovBoundVert = new TVertex**[3];
  TVertex *temp_Mov;
  TBoundEdge *tempSlip_Joint;
  FE2D FeId;
  TFEDesc2D *FeDesc;

  BoundCondFunct2D *BoundaryConditions[2];
  BoundValueFunct2D *BoundValues[2];

  BoundCondFunct2D *GridBoundaryConditions[1];
  BoundValueFunct2D *GridBoundValues[1];

#ifdef __ENERGY__
  BoundCondFunct2D *HeatBoundaryConditions[1];
  BoundValueFunct2D *HeatBoundValues[1];
#endif

  struct triangulateio In;

//   std::ostringstream opts;
  std::ostringstream os;

  os << " ";
//   opts << " ";
//======================================================================
// read parameter file
//======================================================================
//   total_time = GetTime();
//   if(argc>=2)
//     ret=Domain->ReadParam(argv[1]);
//   else
//     ret=Domain->ReadParam(ReadinDat);

  // ======================================================================
  // set the database values and generate mesh
  // ======================================================================    
 /** set variables' value in TDatabase using argv[1] (*.dat file) */
  Domain = new TDomain(argv[1]);  

//   if(ret==-1)
//   {
//     exit(-1);
//   }

  OpenFiles();
  OutFile.setf(std::ios::scientific);

  Database->CheckParameterConsistencyNSE();
  Database->WriteParamDB(argv[0]);
  Database->WriteTimeDB();
  ExampleFile();
//======================================================================
// copy read parameters into local variables
//======================================================================
  if(TDatabase::ParamDB->DISCTYPE==2)
  {
    OutPut("SDFEM does not work!" << endl);
    Error("SDFEM does not work!" << endl);
    exit(4711);
  }
  if( TDatabase::ParamDB->DISCTYPE==5 )
  {
    OutPut("DISCTYPE 5 NOT IMPLEMENTED!" << endl);
    Error("DISCTYPE 5 NOT IMPLEMENTED!" << endl);
    exit(4711);
  }

  PRM = TDatabase::ParamDB->BNDFILE;
  GEO = TDatabase::ParamDB->GEOFILE;
  
  PsBaseName = TDatabase::ParamDB->BASENAME;
  VtkBaseName = TDatabase::ParamDB->VTKBASENAME;

  mkdir(vtkdir, 0777);
  mkdir(BDdir, 0777);
  
  
#define __AXIAL3D__ 

//======================================================================
// Triangular for grid generation begin
//======================================================================
  BoundPart = Domain->GetBdPart(0);
  UpdateSlipBound = (TBdLine*)BoundPart->GetBdComp(0);
  UpdateFreeBound = (TBdCircle*)BoundPart->GetBdComp(1);
  UpdateAxialBound = (TBdLine*)BoundPart->GetBdComp(2);

  double *Sx, *Sy, area = TDatabase::ParamDB->Area;
  double  phi, r = TDatabase::ParamDB->P4, mod_r;
  double  refX;
  int N_refX;
  int *N_MovVert = new int[3];

  N_FreeBound_Vert = (int)TDatabase::ParamDB->P6;    //Freesurf except end point
  N_AxialBound_Vert = 20;
  
  x = (Pi/2.)/TDatabase::ParamDB->P6;
  N_AxialBound_Vert = (int)(1.0/x);
  N_SlipBound_Vert = N_AxialBound_Vert; 
  
#ifndef __SURFACTCONVTEST__
    N_SlipBound_Vert = 2;     // Initially only two points on solid bound (except end point)
#endif  
  
  N_Interf_Vertices = N_FreeBound_Vert+N_SlipBound_Vert+N_AxialBound_Vert;
  In.numberofpoints = N_Interf_Vertices;
  In.pointlist = new double[2*In.numberofpoints];
  In.pointmarkerlist = new int[In.numberofpoints];
  In.numberofpointattributes = 0;

  In.numberofsegments = In.numberofpoints;
  In.segmentlist = new int[2*In.numberofsegments];
  In.segmentmarkerlist = new int[In.numberofsegments];
  In.numberofholes = 0;
  In.holelist = NULL;
  In.numberofregions = 0;
  In.regionlist = NULL;

  In_Index = 0;
  CurrComp = 1;
  
          a=0;
          b=r;
  teta = -Pi/2.; // end degree value of freesurface
  dt = (Pi/2. - teta)/(double)(N_SlipBound_Vert+N_FreeBound_Vert);  
#ifdef __SURFACTCONVTEST__
  dt = r/(double)(N_SlipBound_Vert);  
#endif
  
  t = teta;
  ty = r;      
  for(i=0;i<N_SlipBound_Vert;i++) // without last point
   {
#ifdef __EXPERIMENTAL__ 
    x = r*cos(t);
    y = r*sin(t); 

    phi = atan2(y, x);
    mod_r = 1.0 + 0.29*(sqrt(5./Pi)*(3.*cos(phi+Pi/2.)*cos(phi+Pi/2.) -1. )/4.); 
    x = a+mod_r*cos(t);
    y = b+mod_r*sin(t);
    
    if(i==0)   ty = mod_r;      
#else
   mod_r = r;
   ty = r;    
#endif
   
#ifdef __SURFACTCONVTEST__
    x = dt*(double)i;  
#else   
    x = a + mod_r*cos(t);
#endif    
    y=0;
 
    if(i==0) x = 0.;    
    
    In.pointlist[2*In_Index] = x;
    In.pointlist[2*In_Index+1] = y;
       cout<<i<< " x : "<< x << " y : "<< y<<endl;
    In.pointmarkerlist[In_Index] = CurrComp;
    In.segmentlist[2*In_Index] = In_Index;
    In.segmentlist[2*In_Index+1] = In_Index+1;
    In.segmentmarkerlist[In_Index] = CurrComp;
    In_Index++;  
#ifndef __SURFACTCONVTEST__    
    t = teta + double(i+1.)*dt;
#endif    
   }
  
   CurrComp++;
    cout<<endl; 
//      exit(0);
     
   Sx = new double[N_FreeBound_Vert+1];
   Sy = new double[N_FreeBound_Vert+1];
 
#ifdef __SURFACTCONVTEST__   
    a = 0.;
    ty = 0.;
    t=0;
    teta=0.;
    dt = (Pi/2.)/(double)N_FreeBound_Vert;      
#endif   
   
   for(i=0;i<N_FreeBound_Vert;i++) // without last point
    {
#ifdef __EXPERIMENTAL__   
     x = r*cos(t);
     y = r*sin(t);  

     phi = atan2(y, x); 
     mod_r = 1.0 + 0.29*(sqrt(5./Pi)*(3.*cos(phi+Pi/2.)*cos(phi+Pi/2.) -1. )/4.);           
#else
   mod_r = r;    
#endif   
   
      x = a + mod_r*cos(t);
      y = ty + mod_r*sin(t); 
      cout<< " i " <<"  ty " << x << " y " << y <<endl;
// #ifndef __SURFACTDIFFTEST__  
     if(i==0) y =0.;        
// #endif   

     if (fabs(x)<1.e-10) x = 0;


     Sx[i] = x;
     Sy[i] = y;   
 #ifdef __SURFACTCONVTEST__   
     t = (double)(i+1)*dt;      
#else     
     t = teta + (double)(N_SlipBound_Vert+i+1)*dt;    
#endif     
    }// for(i=0;i<N_FreeBound_Ver
//      exit(0);
   //endpoint 
     Sx[N_FreeBound_Vert] = 0;
     Sy[N_FreeBound_Vert] = 2*ty;       
   
#ifdef __EXPERIMENTAL__       
     EqDist_Pts(N_FreeBound_Vert+1, Sx, Sy); 
#endif       
     
   for(i=0;i<N_FreeBound_Vert;i++) // without last point
    {
      x = Sx[i];
      y = Sy[i];
      
     if (fabs(x)<1.e-10) x = 0;
     
// #ifdef __SURFACTDIFFTEST__   
//      
//      SLPX = x;
//      
//     if(i==N_FreeBound_Vert-2) 
//        h_interface = sqrt( (In.pointlist[2*(In_Index-1)] - x)*(In.pointlist[2*(In_Index-1)] - x) 
//                + (In.pointlist[2*(In_Index-1) +1] - y)*(In.pointlist[2*(In_Index-1)+1] - y));   
// #else     
     if(i==0) 
       {
        y =0.;   SLPX = x;
        Indextemp = In_Index;    
       }
      else if(i==N_FreeBound_Vert-2) 
      {
       h_interface = sqrt( (In.pointlist[2*(In_Index-1)] - x)*(In.pointlist[2*(In_Index-1)] - x) 
               + (In.pointlist[2*(In_Index-1) +1] - y)*(In.pointlist[2*(In_Index-1)+1] - y));
      }
// #endif

      In.pointlist[2*In_Index] = x;
      In.pointlist[2*In_Index+1] = y;   
//        cout<<i<< " x : "<< x << " y : "<< y<<endl;
      In.pointmarkerlist[In_Index] = CurrComp;
      In.segmentlist[2*In_Index] = In_Index;
      In.segmentlist[2*In_Index+1] = In_Index+1;
      In.segmentmarkerlist[In_Index] = CurrComp;
      In_Index++;
    }
//    exit(0);
  cout<<endl; 
   CurrComp++;
   dt=-2.*ty/N_AxialBound_Vert;
   t = -teta;
   
#ifdef __EXPERIMENTAL__  
    phi = Pi/2;     
    mod_r = 1.0 + 0.29*(sqrt(5./Pi)*(3.*cos(phi+Pi/2.)*cos(phi+Pi/2.) -1. )/4.); 
#else
    mod_r = r;       
#endif    

    y = ty+mod_r*sin(t);
   
#ifdef __EXPERIMENTAL__  
      y =  mod_r*sin(t);
#endif   
      
#ifdef __SURFACTCONVTEST__
      y = 1.;    
#endif      
      
   dt= -y / (double)(N_AxialBound_Vert);
   x = 0.;
   t = y;   

   
   UpdateAxialBound->SetParams(x, y, 0, -y);  
    for(i=0;i<N_AxialBound_Vert;i++) // without last point
     {
      if (fabs(y)<1e-12) y = 0.;
      In.pointlist[2*In_Index] = x;
      In.pointlist[2*In_Index+1] = y;
      cout<<" x : "<< x << " y : "<< y<<endl;
      In.pointmarkerlist[In_Index] = CurrComp;
      In.segmentlist[2*In_Index] = In_Index;
      In.segmentlist[2*In_Index+1] = In_Index+1;
      In.segmentmarkerlist[In_Index] = CurrComp;
      In_Index++;
      y = t + (double)(i+1)*dt;
    }   

  In.segmentlist[2*(In_Index-1)+1] = 0;
 
 delete [] Sx;
 delete [] Sy; 
//   exit(0);

 UpdateSlipBound->SetParams(In.pointlist[0],In.pointlist[1],
                          In.pointlist[2*N_SlipBound_Vert]-In.pointlist[0],
                          In.pointlist[2*N_SlipBound_Vert+1]-In.pointlist[1]); 
 
#ifdef __SURFACTCONVTEST__    
// Free boundary xmid, ymid, radius_a, radius_b, start angle, end angle
 UpdateFreeBound ->SetParams(0.0, 0.0, 1.0, 1.0, 0., Pi/2.);  
#endif  
 
  /** generate the mesh for the input object IN */ 
  Domain->TriMeshGen(&In);
  
//======================================================================
// Triangular for grid generation end
//======================================================================
#ifdef __SURFACTCONVTEST__  
  // refine grid up to the coarsest level
  for(i=0;i<TDatabase::ParamDB->UNIFORM_STEPS;i++)
    Domain->RegRefineAll();
#endif

//      exit(0);
//       write grid into an Postscript file
      os.seekp(std::ios::beg);
      os << "Domain" << ".ps" << ends;
      Domain->PS(os.str().c_str(),It_Finest,0);
//      exit(0);

   //Initialize DiscreteForms         
  InitializeDiscreteForms_Moving(DiscreteFormGalerkin, DiscreteFormNLGalerkin,
                                 DiscreteFormGrid, LinCoeffs, GridCoeffs);
 
  BoundaryConditions[0] = BoundCondition;
  BoundaryConditions[1] = BoundCondition;

  BoundValues[0] = U1BoundValue;
  BoundValues[1] = U2BoundValue;

  GridBoundaryConditions[0] = GridBoundCondition;
  GridBoundValues[0] = GridBoundValue;

#ifdef __ENERGY__
  HeatBoundaryConditions[0] = HeatBoundCondition;
  HeatBoundValues[0] = TBoundValue;
#endif

#ifdef __SURFACT__
  InitializeDiscreteForms_Moving(DiscreteFormSurfact, DiscreteFormSurfact_SUPG, SurfactCoeffs);   
  SurfactBoundaryConditions[0] = SurfactBoundCondition;
  SurfactBoundValues[0] = SurfactBoundValue;
#endif
//======================================================================
// construct all finite element spaces
//======================================================================
  coll=Domain->GetCollection(It_Finest, 0);
  N_Cells = coll->GetN_Cells();
  
#ifdef __SURFACT__
//  list of outer phase cells containing interafce
  N_List[0] = new int[N_Cells];   // Cell_No
  N_List[1] = new int[N_Cells];   //Joint_No
  Domain2DSurf_2Phase(coll, IFaceDomain, N_List);  
  
  IFace_Coll = IFaceDomain->GetCollection(It_Finest, 0);
  N_IFaceCells= IFace_Coll->GetN_Cells();
  
  OutPut("N_Cells      : "<< setw(10) << N_Cells  << endl);
  OutPut("N_IFaceCells : "<<  setw(10) <<  N_IFaceCells  << endl); 
    
  FE1D_List = new FE1D[N_Cells];
  for(j=0;j<N_Cells;j++)
  FE1D_List[j] = FE1D(TDatabase::ParamDB->ANSATZ_ORDER);

  TDatabase::IteratorDB[It_EQ]->SetParam(Domain);
  TDatabase::IteratorDB[It_LE]->SetParam(Domain);
  TDatabase::IteratorDB[It_Finest]->SetParam(Domain);
  TDatabase::IteratorDB[It_Between]->SetParam(Domain);
  TDatabase::IteratorDB[It_OCAF]->SetParam(Domain);    
# endif
 
  ORDER = 0;
  VSP = TDatabase::ParamDB->VELOCITY_SPACE;

  if(abs(VSP) > 20)
   {ORDER = abs(VSP) - 20;}
  else if (abs(VSP) > 10)
    {ORDER = abs(VSP) - 10;}
  else ORDER = abs(VSP);

//   NSE_coll = new TCollection(N_NSE_Cells, NSE_Cells);
//======================================================================
// construct all finite element spaces
//======================================================================
  // get velocity and pressure spacess
  GetVelocityAndPressureSpace(coll,BoundCondition,
                              mortarcoll, velocity_space,
                              pressure_space, &pressure_space_code,
                              TDatabase::ParamDB->VELOCITY_SPACE,
                              TDatabase::ParamDB->PRESSURE_SPACE);
  velocity_space_code =   TDatabase::ParamDB->VELOCITY_SPACE;
  FESpaces_All[0] = velocity_space;
  FESpaces_All[1] =  pressure_space;  
  N_Active =  FESpaces_All[0]->GetActiveBound();
  N_U = FESpaces_All[0]->GetN_DegreesOfFreedom();
  N_P = FESpaces_All[1]->GetN_DegreesOfFreedom();

  GlobalNumbers = FESpaces_All[0]->GetGlobalNumbers();
  BeginIndex = FESpaces_All[0]->GetBeginIndex();
  
// mesh velocity space 
   grid_space = new TFESpace2D(coll, NameString, TString, GridBoundCondition, 1, NULL);
   FESpaces_All[2] =  grid_space;     
   N_G = FESpaces_All[2]->GetN_DegreesOfFreedom();
   N_GActive = FESpaces_All[2]->GetActiveBound();
   N_GBoundaryNodes = N_G - N_GActive;
   
   OutPut("N_G    : "<< setw(10) << N_G  << endl);
   OutPut("N_GActive    : "<< setw(10) << N_GActive  << endl);
   
// thermal space
#ifdef __ENERGY__
   thermal_space = new TFESpace2D(coll, NameString, TString, HeatBoundCondition,
                                  TDatabase::ParamDB->ANSATZ_ORDER, NULL);
   FESpaces_All[3] =  thermal_space;  
   N_thermalDOF = FESpaces_All[3]->GetN_DegreesOfFreedom();
   N_thermalActive = FESpaces_All[3]->GetActiveBound();
   N_thermalNonActive = N_thermalDOF - N_thermalActive;
   OutPut("N_thermalDOF    : "<< setw(10) << N_thermalDOF  << endl);
#endif

   // surfact space
#ifdef __SURFACT__
   surfact_space = new TFESpace2D(coll, NameString, CString, SurfactBoundCondition,
                                  TDatabase::ParamDB->ANSATZ_ORDER, NULL);
   FESpaces_All[4] =  surfact_space;  
   N_surfactDOF = FESpaces_All[4]->GetN_DegreesOfFreedom();
   N_surfactActive = FESpaces_All[4]->GetActiveBound();
   N_surfactNonActive = N_surfactDOF - N_surfactActive;
   OutPut("N_SurfactDOF    : "<< setw(10) << N_surfactDOF  << endl);
//     N_CActive = FESpaces_All[4]->GetActiveBound(); 
   IFaceSurfact_space = new TFESpace1D(IFace_Coll , IFaceSString, IFaceSString, FE1D_List);

   IFaceFeSpaces[0] = IFaceSurfact_space;
   N_IsurfactDOF = IFaceFeSpaces[0]->GetN_DegreesOfFreedom();
   N_IActive = IFaceFeSpaces[0]->GetActiveBound();
   OutPut("N_SurfSurfactDOF    : "<< setw(10) << N_IsurfactDOF  << endl);   
   
   FESpaces_All[5] =  new TFESpace2D(coll, NameString, IFaceSString, SurfactBoundCondition,
                                     TDatabase::ParamDB->ANSATZ_ORDER, NULL);

   N_S =  FESpaces_All[5]->GetN_DegreesOfFreedom();   
#endif   
   
//======================================================================
// construct all finite element functions
//======================================================================
  N_Unknowns = 2*N_U + N_P;
  Sol_All[0] = new double[N_Unknowns];  
  Rhs_All[0] = new double[N_Unknowns];
  
  OutPut("dof velocity : "<< setw(10) << 2* N_U << endl);
  OutPut("dof pressure : "<< setw(10) << N_P << endl);
  OutPut("dof all      : "<<  setw(10) << N_Unknowns  << endl);  

  B = new double[N_Unknowns];
  defect = new double[N_Unknowns];
  
  memset(Sol_All[0], 0, N_Unknowns*SizeOfDouble);
  memset(Rhs_All[0], 0, N_Unknowns*SizeOfDouble);

  //velo vect
  FEVectFuncts_All[0] =  new TFEVectFunct2D(FESpaces_All[0], UString, UString, Sol_All[0], N_U, 2);
  FEFunctions_All[0] = FEVectFuncts_All[0]->GetComponent(0);
  FEFunctions_All[1] = FEVectFuncts_All[0]->GetComponent(1); 
  
#ifndef __SURFACTCONVTEST__      
// #ifndef __SURFACTDIFFTEST__   
  FEFunctions_All[0]->Interpolate(InitialU1);
  FEFunctions_All[1]->Interpolate(InitialU2);
// #endif
#endif
  
  //pressure FeFunction
  FEFunctions_All[2] = new TFEFunction2D(FESpaces_All[1], PString,  PString,  Sol_All[0]+2*N_U, N_P);

  tmp_Gd = new double[2*N_G];
  tmp_Gsol = new double[2*N_G];
  Sol_All[1] = new double[2*N_G];
  Rhs_All[1] = new double[2*N_G];   
  
  memset(Sol_All[1], 0, 2*N_G*SizeOfDouble);   
  FEVectFuncts_All[1]  = new TFEVectFunct2D(FESpaces_All[2], WString, WString, Sol_All[1], N_G, 2);
  
  FEFunctions_All[3] = FEVectFuncts_All[1]->GetComponent(0);
  FEFunctions_All[4] = FEVectFuncts_All[1]->GetComponent(1);
//======================================================================
// thermal space finite element functions
//======================================================================
#ifdef __ENERGY__
  Sol_All[2] = new double[N_thermalDOF];
  Rhs_All[2] = new double[N_thermalDOF];
  oldsol_T = new double[N_thermalDOF];
  
  memset(Sol_All[2], 0, N_thermalDOF*SizeOfDouble);
  memset(oldsol_T, 0, N_thermalDOF*SizeOfDouble);
  memset(Rhs_All[2], 0, N_thermalDOF*SizeOfDouble);

  // thermal fefunction
  FEFunctions_All[5] = new TFEFunction2D(FESpaces_All[3], TString, TString, Sol_All[2], N_thermalDOF);
#endif
  
//======================================================================
// surfact space finite element functions
//======================================================================
#ifdef __SURFACT__
  Sol_All[3] = new double[N_surfactDOF];
  Rhs_All[3] = new double[N_surfactDOF];
  Csol_old = new double[N_surfactDOF];
  CRhs_old = new double[N_surfactDOF];
  C_B = new double[N_surfactDOF];  
  C_defect = new double[N_surfactDOF];  
  Csol_nonlinearstep = new double[N_surfactDOF];
   
  memset(Sol_All[3], 0, N_surfactDOF*SizeOfDouble);
  memset(Csol_old, 0, N_surfactDOF*SizeOfDouble);
  memset(Rhs_All[3], 0, N_surfactDOF*SizeOfDouble);
  memset(CRhs_old, 0, N_surfactDOF*SizeOfDouble);
  // surfact fefunction
  FEFunctions_All[6] = new TFEFunction2D(FESpaces_All[4], CString, CString, Sol_All[3], N_surfactDOF);
  
// #ifndef __SURFACTDIFFTEST__ 
#ifdef  __SOLUBLE__
  FEFunctions_All[6]->Interpolate(InitialSuract);
#endif 
// #endif  
  
   Sol_All[4] = new double[N_IsurfactDOF];
   Rhs_All[4] = new double[N_IsurfactDOF];

   I_defect = new double[N_IsurfactDOF];
   Isol_old = new double[N_IsurfactDOF];
   I_B = new double[N_IsurfactDOF];
   IRhs_old = new double[N_IsurfactDOF];

   
   memset(Sol_All[4], 0, N_IsurfactDOF*SizeOfDouble);
   memset(Rhs_All[4], 0, N_IsurfactDOF*SizeOfDouble);

   memset(I_defect, 0, N_IsurfactDOF*SizeOfDouble);
   memset(Isol_old, 0, N_IsurfactDOF*SizeOfDouble);
   memset(I_B, 0, N_IsurfactDOF*SizeOfDouble);
   memset(IRhs_old, 0, N_IsurfactDOF*SizeOfDouble);
   
   IFaceFeFunct[0] = new TFEFunction1D(IFaceFeSpaces[0], IFaceSString, IFaceSString, Sol_All[4], N_IsurfactDOF);
//    IFaceFeFunct[0]->Interpolate(InitialS);   
 
   
   //for output 
   Sol_All[5] =  new double[N_S];
   FEFunctions_All[7]  = new TFEFunction2D(FESpaces_All[5], IFaceSString,  IFaceSString, Sol_All[5],  N_S);  
   
   FEFunctions_All[7]->Interpolate(InitialS); 
   MapDomainToSurf(FEFunctions_All[7], IFaceFeFunct[0], N_List[0], N_List[1]);   
   memset(Sol_All[5], 0, N_S*SizeOfDouble);   
 #endif 

  
//======================================================================
// allocate memory for all matrices
//======================================================================
  Structure_All[0] = new TStructure2D(FESpaces_All[1], FESpaces_All[0]);  // B
  Structure_All[1] = new TStructure2D(FESpaces_All[0], FESpaces_All[1]); // BT
  
  //velo
  SquareStructure_All[0] = new TSquareStructure2D(FESpaces_All[0]);  
  SquareStructure_All[0]->Sort();

  // grid 
  SquareStructure_All[1] = new TSquareStructure2D(FESpaces_All[2]); 
  SquareStructure_All[1]->Sort();
  
  //thermal
#ifdef __ENERGY__
  SquareStructure_All[2] = new TSquareStructure2D(FESpaces_All[3]);
  SquareStructure_All[2]->Sort();
#endif
 
#ifdef __SURFACT__
  //surfact
  SquareStructure_All[3] = new TSquareStructure2D(FESpaces_All[4]);
  SquareStructure_All[3]->Sort();
 
  /* interface surfactant matrices */
  IFaceStruct[0] = new TSquareStructure1D(IFaceFeSpaces[0]);
  IFaceStruct[0]->Sort();
#endif     
  
  // for NSE
  MatVect = MatVect_NSE4;
  Defect = Defect_NSE4;

  SqMat_All[0] = new TSquareMatrix2D(SquareStructure_All[0]); // M11
  SqMat_All[1] = new TSquareMatrix2D(SquareStructure_All[0]); // M12
  SqMat_All[2] = new TSquareMatrix2D(SquareStructure_All[0]); // M21
  SqMat_All[3] = new TSquareMatrix2D(SquareStructure_All[0]); // M22
  
  SqMat_All[4] = new TSquareMatrix2D(SquareStructure_All[0]); // A11
  SqMat_All[5] = new TSquareMatrix2D(SquareStructure_All[0]); // A12
  SqMat_All[6] = new TSquareMatrix2D(SquareStructure_All[0]); // A21
  SqMat_All[7] = new TSquareMatrix2D(SquareStructure_All[0]); // A22
  SqMat_All[8] = new TSquareMatrix2D(SquareStructure_All[0]); // F11
  SqMat_All[9] = new TSquareMatrix2D(SquareStructure_All[0]); // F22

  Mat_All[0] = new TMatrix2D(Structure_All[0]); // B1
  Mat_All[1] = new TMatrix2D(Structure_All[0]); // B2
  Mat_All[2] = new TMatrix2D(Structure_All[1]); // B1T
  Mat_All[3] = new TMatrix2D(Structure_All[1]); // B2T

  // for mesh
  SqMat_All[10] = new TSquareMatrix2D(SquareStructure_All[1]); // G11
  SqMat_All[11] = new TSquareMatrix2D(SquareStructure_All[1]); // G12
  SqMat_All[12] = new TSquareMatrix2D(SquareStructure_All[1]); // G21
  SqMat_All[13] = new TSquareMatrix2D(SquareStructure_All[1]); // G22

  // for heat
#ifdef __ENERGY__  
  SqMat_All[14]  = new TSquareMatrix2D(SquareStructure_All[2]); // T_M
  SqMat_All[15] = new TSquareMatrix2D(SquareStructure_All[2]); // T_A
#endif     
  
#ifdef __SURFACT__  
  //bulk surfact
  SqMat_All[16]  = new TSquareMatrix2D(SquareStructure_All[3]); // C_M
  SqMat_All[17] = new TSquareMatrix2D(SquareStructure_All[3]); // C_A
  
 /* interface surfactant matrices */
  SqMat_IFace[0] = new TSquareMatrix1D(IFaceStruct[0]); // Gamma_A
  SqMat_IFace[1] = new TSquareMatrix1D(IFaceStruct[0]); // Gamma_M  
#endif   
  
  
  IsoCellEdgeNos = new int *[2];
 
 //  ====================================================================================   
  y = 0.;
  GetMovingBoundData(coll, N_MovVert, Bound_Joint, MovBoundVert, Free_Joint,
                     Free_Cells, IsoCellEdgeNos, SLPX, y);
                                          
//  ====================================================================================  
// assemble matrix for grid moving - begin
//  ====================================================================================  
    fesp[0] = FESpaces_All[2];
    SQMATRICES_GRID[0] = SqMat_All[10];
    SQMATRICES_GRID[0]->Reset();
    SQMATRICES_GRID[1] = SqMat_All[11];
    SQMATRICES_GRID[1]->Reset();
    SQMATRICES_GRID[2] = SqMat_All[12];
    SQMATRICES_GRID[2]->Reset();
    SQMATRICES_GRID[3] = SqMat_All[13];
    SQMATRICES_GRID[3]->Reset();
    aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);  
       
    Assemble2D(1, fesp,
             4, SQMATRICES_GRID,
             0, NULL,
             0, NULL, NULL,
             DiscreteFormGrid,
             GridBoundaryConditions,
             GridBoundValues,
             aux);
    delete aux;   
   
     Entries[0] = SqMat_All[10]->GetEntries();
     Entries[1] = SqMat_All[11]->GetEntries();
     Entries[2] = SqMat_All[12]->GetEntries();
     Entries[3] = SqMat_All[13]->GetEntries();

     GridKCol = SquareStructure_All[1]->GetKCol();
     GridRowPtr = SquareStructure_All[1]->GetRowPtr();

  // for Dirichlet rows in off-diagonal matrices
  memset(Entries[1] + GridRowPtr[N_GActive], 0, (GridRowPtr[N_G] - GridRowPtr[N_GActive])*SizeOfDouble);
  memset(Entries[2] + GridRowPtr[N_GActive], 0, (GridRowPtr[N_G] - GridRowPtr[N_GActive])*SizeOfDouble);

  refpos = new double[2*N_G];
  auxpos = new double[2*N_G];
  pos = new double[2*N_G];
  ReparamMeshVelo = new double[2*N_G];  
  ReparamDisp = new double[2*N_G]; 
  ReparamPos = new TFEVectFunct2D(FESpaces_All[2], WString, WString, ReparamDisp, N_G, 2);  
  RefGridPos = new TFEVectFunct2D(FESpaces_All[2], WString, WString, refpos, N_G, 2);
  AuxGridPos = new TFEVectFunct2D(FESpaces_All[2], WString, WString, auxpos, N_G, 2);
  GridPos = new TFEVectFunct2D(FESpaces_All[2], WString, WString, pos, N_G, 2);

  
  RefGridPos->GridToData();
  AuxGridPos->GridToData();
  GridPos->GridToData();

  // prepare output (maxn_fespaces,  maxn_scalar,  maxn_vect, maxn_parameters, domain)
  Output = new TOutput2D(1, 3, 1, 2, Domain);
  Output->AddFEVectFunct(FEVectFuncts_All[0]);
  Output->AddFEFunction(FEFunctions_All[2]);   
#ifdef __ENERGY__      
  Output->AddFEFunction(FEFunctions_All[5]);
#endif      
  
#ifdef __SURFACT__      
  Output->AddFEFunction(FEFunctions_All[6]);
  
  Output->AddFEFunction(FEFunctions_All[7]);
#endif     
//   Output->AddFEVectFunct(FEVectFuncts_All[1]);  
      
  os.seekp(std::ios::beg);
  Output->AddParameter(TDatabase::TimeDB->CURRENTTIME,os.str().c_str());      
 

// #ifdef __SURFACTDIFFTEST__    
#ifndef __SURFACTCONVTEST__     
   MapSurfToDomain(IFaceFeFunct[0], FEFunctions_All[7], N_List[0], N_List[1]);
   GetSurfErrors(FEFunctions_All[7], IFaceFeFunct[0], errors, N_List[0], N_List[1]);
      
   OutPut( "L2: " << errors[0] << endl);
   OutPut( "H1-semi: " << errors[1] << endl);
#else
   TDatabase::ParamDB->REACTOR_P13 = 0; // Bi =
   TDatabase::ParamDB->REACTOR_P15 = 0; // beta =
   TDatabase::ParamDB->REACTOR_P7 = 1; // fixpt itr =  
   TDatabase::ParamDB->REACTOR_P14 = 1; // Da = 
   
  TDatabase::TimeDB->TIMESTEPLENGTH = pow(h_interface, TDatabase::ParamDB->REACTOR_P30);   
#endif    
   olderror=0.;
   Terrors[0] = 0;
   Terrors[1] = 0;
   T_inftyL2 = 0;
   T_inftyH1 = 0;
   
// #ifdef __SURFACTDIFFTEST__     
// 
// #endif
//   OutPut( "tau = " << tau <<endl);
// #endif
   Get_KE(FEVectFuncts_All[0], Params);
   InitVolume = CurrVolume = Params[0];
   KE = Params[1];
#ifdef __SURFACT__            
       // interface surfactant
      MapSurfToDomain(IFaceFeFunct[0], FEFunctions_All[7], N_List[0], N_List[1]);
     
      GetSurfactMass(FEFunctions_All[7], IFaceFeFunct[0], N_List[0], N_List[1], Surf_Mass);
      radius = pow((3.*CurrVolume), (1./3.)); // 2Pi get acanelled to to volume calculation in Get_KE
//       OutPut("CurrVolume : " << CurrVolume << " radius : " << radius <<endl);
      sphericity =  2.*Pi*radius*radius/(Surf_Mass[1]);  // hemisphire, so 2. * only   
 
//       Get_FeFunction2DMass(FEFunctions_All[6], Params);   
      FEFunctions_All[6]->GetMassAndMean(Params);     
      Initial_IFaceSurfactMass = TDatabase::ParamDB->REACTOR_P14*Surf_Mass[0];
      Initial_SurfactMass = Params[0];    
      
      OutPut( "Time, Surfactant_Mass, Da*InterfaceSurfactant_Mass, InterfaceSurfactant_Conc " <<TDatabase::TimeDB->CURRENTTIME<<
              " " <<Params[0]<< " "<<TDatabase::ParamDB->REACTOR_P14*Surf_Mass[0]<< " "<<Surf_Mass[0]/Surf_Mass[1]<< " "<<endl);
      OutPut( "Time, Surfactant_Mass_dif, InterfaceSurfactant_Mass_diff " <<TDatabase::TimeDB->CURRENTTIME<<
              " " <<Params[0] - Initial_SurfactMass<< " "<<TDatabase::ParamDB->REACTOR_P14*Surf_Mass[0]-Initial_IFaceSurfactMass<< endl);
      OutPut( "Time, Total Mass, Total Mass diff, RelativeTotal Mass diff, sphericity, KE " <<TDatabase::TimeDB->CURRENTTIME<<
          " " <<Params[0]+TDatabase::ParamDB->REACTOR_P14*Surf_Mass[0]<< " " << ((Params[0]+TDatabase::ParamDB->REACTOR_P14*Surf_Mass[0]) - 
            (Initial_IFaceSurfactMass + Initial_SurfactMass)) << " " << ((Params[0]+TDatabase::ParamDB->REACTOR_P14*Surf_Mass[0]) - 
            (Initial_IFaceSurfactMass + Initial_SurfactMass))
            /(Initial_IFaceSurfactMass + Initial_SurfactMass) <<" " << sphericity << " " << KE<<endl);   
      
      
      PrintSurfSurfactant(N_MovVert[2], MovBoundVert[2], FEFunctions_All[7], N_BData_F, Surf_F);   
      PrintSurfSurfactant(N_MovVert[0], MovBoundVert[0], FEFunctions_All[7], N_BData_S, Surf_S);         
#endif        
      
   if(TDatabase::ParamDB->WRITE_VTK)
     {       
      os.seekp(std::ios::beg);
        if(img<10) os << "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
         else if(img<100) os <<"VTK/"<< VtkBaseName<<".000"<<img<<".vtk" << ends;
         else if(img<1000) os <<"VTK/"<< VtkBaseName<<".00"<<img<<".vtk" << ends;
         else if(img<10000) os << "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
         else  os << "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
      Output->WriteVtk(os.str().c_str());
      
//        os.seekp(std::ios::beg);
//        if(img<10) os << "VTK/"<<VtkBaseName<<".0000"<<img<<".plt" << ends;
//        else if(img<100) os << "VTK/"<<VtkBaseName<<".000"<<img<<".plt" << ends;
//        else if(img<1000) os << "VTK/"<<VtkBaseName<<".00"<<img<<".plt" << ends;
//        else if(img<10000) os << "VTK/"<<VtkBaseName<<".0"<<img<<".plt" << ends;
//        else  os << "VTK/"<<VtkBaseName<<"."<<img<<".plt" << ends; 
//        Output->WriteBinaryPlt(os.str().c_str());      
//       
      img++;
     }
      
// #ifndef __SURFACTDIFFTEST__       
   MovBoundVert[0][0]->GetCoords(Lx, Ly);
   MovBoundVert[2][0]->GetCoords(Rx, Ry);
   OutPut(setw(20)<<"T, Wett Len d : " << TDatabase::TimeDB->CURRENTTIME<<"   "<< Rx-Lx <<endl);
// #endif
   OutPut(setw(20)<<"T, Volume : " << TDatabase::TimeDB->CURRENTTIME<<"   "<< CurrVolume<<endl);
   OutPut(setw(20)<<"T, Volume Diff : "<< TDatabase::TimeDB->CURRENTTIME<<"   "<< CurrVolume - InitVolume << endl);
   
   Getcellangle(FESpaces_All[2], Angle);   
   OutPut( "MinAngle : "<< Angle[0]<< "  MaxAngle : "<<Angle[1]<< endl);
 
   
   TDatabase::ParamDB->P10 = 1; // free surf reparam
   TDatabase::ParamDB->P5 = 0;  // move boundary with velo
   
#ifdef __SURFACTCONVTEST__    
   olderror = sqrt(Ddot(N_IsurfactDOF,Sol_All[4],Sol_All[4]));
#endif  
// exit(0);
      OutPut( "h_interface : "<< h_interface << endl);
//    TDatabase::TimeDB->TIMESTEPLENGTH = TDatabase::TimeDB->TIMESTEPLENGTH/10.;
//    exit(0);
      
   N_Surf_M_MatValues = SqMat_IFace[1]->GetN_Entries();
   Surf_M_MatValues = new double[N_Surf_M_MatValues];

//======================================================================
// start of time cycle
//======================================================================
  end_time = TDatabase::TimeDB->ENDTIME;
  N_SubSteps = GetN_SubSteps(); 
  TDatabase::TimeDB->CURRENTTIME = TDatabase::TimeDB->STARTTIME;
  m=0;
  limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE;
  Max_It = TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE; 
  solver_time = 0.0;
  N_LinIter = 0; 
  t3 = GetTime();
  total_time = t3 - total_time;
  
#ifdef __SURFACT__       
  surf_couple_var = (int)TDatabase::ParamDB->REACTOR_P7;
  Max_It_scalar = (int)TDatabase::ParamDB->REACTOR_P9;
  if(surf_couple_var==1) Max_It_scalar = 1; // explicit coupling 
  limit_scalar = TDatabase::ParamDB->REACTOR_P8;  
#endif     
  
  while(TDatabase::TimeDB->CURRENTTIME< end_time)
  {
    // time cycle
    m++;
    TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;

      for(l=0;l<N_SubSteps;l++)   // sub steps of fractional step theta
      {    
        SetTimeDiscParameters(1);

        if (m==1)
        {
          OutPut("Theta1: " << TDatabase::TimeDB->THETA1<< endl);
          OutPut("Theta2: " << TDatabase::TimeDB->THETA2<< endl);
          OutPut("Theta3: " << TDatabase::TimeDB->THETA3<< endl);
          OutPut("Theta4: " << TDatabase::TimeDB->THETA4<< endl);
        }

        tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;

// 	cout << " tau " << tau <<  endl;

        TDatabase::TimeDB->CURRENTTIME += tau;

        if (very_first_time)
            oldtau=tau;
// #ifndef __SURFACTCONVTEST__  
// #ifndef __SURFACTDIFFTEST__  
        // working rhs array for NSE
        memset(B, 0, N_Unknowns*SizeOfDouble);

        GridVelo_imping(Entries, tmp_Gsol, tmp_Gd, Rhs_All[1],
                        GridKCol, GridRowPtr,
                        GridPos, AuxGridPos,
                        FEVectFuncts_All[0], tau,
                        FEVectFuncts_All[1], MovBoundVert, N_MovVert,
                        Free_Cells, IsoCellEdgeNos, reparam, RefGridPos);
	


// ============================================================================================
//  Assemble NSE
// ============================================================================================     
        DiscreteForm = DiscreteFormGalerkin;

        SQMATRICES[0] = SqMat_All[4];
        SQMATRICES[1] = SqMat_All[5];
        SQMATRICES[2] = SqMat_All[6];
        SQMATRICES[3] = SqMat_All[7];
        SQMATRICES[4] = SqMat_All[0];
        SQMATRICES[5] = SqMat_All[3];

        SQMATRICES[6] = SqMat_All[1];
        SQMATRICES[7] = SqMat_All[2];

        MATRICES[0] = Mat_All[0];
        MATRICES[1] = Mat_All[1];
        MATRICES[2] = Mat_All[2];
        MATRICES[3] = Mat_All[3];

        SQMATRICES[0]->Reset();
        SQMATRICES[1]->Reset();
        SQMATRICES[2]->Reset();
        SQMATRICES[3]->Reset();
        SQMATRICES[4]->Reset();
        SQMATRICES[5]->Reset();
        SQMATRICES[6]->Reset();
        SQMATRICES[7]->Reset();

        MATRICES[0]->Reset();
        MATRICES[1]->Reset();
        MATRICES[2]->Reset();
        MATRICES[3]->Reset();
	
        N_SquareMatrices = 6;
        N_RectMatrices = 4;

       // parameters which are the same for all NSTYPEs
        N_Rhs = 2;
        N_FESpaces = 3;

        fesp[0] = FESpaces_All[0];
        fesp[1] = FESpaces_All[1];
        fesp[2] = FESpaces_All[2];

        fefct[0] = FEFunctions_All[0];
        fefct[1] = FEFunctions_All[1];
        fefct[2] = FEFunctions_All[3];
        fefct[3] = FEFunctions_All[4];
 
        ferhs[0] = FESpaces_All[0];
        ferhs[1] = FESpaces_All[0];	

        RHSs[0] = Rhs_All[0];
        RHSs[1] = Rhs_All[0] + N_U;
        RHSs[2] = Rhs_All[0] + 2*N_U;

        memset(Rhs_All[0], 0, N_Unknowns*SizeOfDouble);

       // 4 parameters are needed for assembling (u1_old, u2_old)
        aux =  new TAuxParam2D(MovingTNSN_FESpaces_Axial3D, MovingTNSN_Fct_Axial3D,
                               MovingTNSN_ParamFct_Axial3D,
                               MovingTNSN_FEValues_Axial3D,
                               fesp, fefct,
                               MovingTNSFct_Axial3D,
                               MovingTNSFEFctIndex_Axial3D,
                               MovingTNSFEMultiIndex_Axial3D,
                               MovingTNSN_Params_Axial3D, MovingTNSBeginParam_Axial3D);                        
                        
       //======================================================================
       // assembling of matrices for each level
       // A_11 , (A_12), (A_21), (A_22), M_11, (M_22)
       //======================================================================
       Assemble2D(N_FESpaces, fesp,
                  N_SquareMatrices, SQMATRICES,
                  N_RectMatrices, MATRICES,
                  N_Rhs, RHSs, ferhs,
                  DiscreteForm,
                  BoundaryConditions,
                  BoundValues,
                  aux);

      delete aux;
      
//            cout << " Ddot " <<  Ddot(N_Unknowns, Rhs_All[0], Rhs_All[0]) << endl;
//        exit(0);   

     SqMat_All[8]->Reset(); // Matrix entries for freesurf int;
     SqMat_All[9]->Reset(); // no need to calculate in nonlinear steps

#ifdef __SURFACT__  
     // interface surfactant
     MapSurfToDomain(IFaceFeFunct[0], FEFunctions_All[7], N_List[0], N_List[1]);  
     
     FreeSurf_axial3D_new(SqMat_All[8], SqMat_All[9],  RHSs[0], RHSs[1], BoundCondition, tau,
                          FEFunctions_All[0]->GetValues(), FEFunctions_All[7], Params);     
#else     
     FreeSurf_axial3D_new(SqMat_All[8], SqMat_All[9],  RHSs[0], RHSs[1], BoundCondition, tau,
                          FEFunctions_All[0]->GetValues(), NULL, Params);
#endif  
     
     // Adding freesurf entries to A11 and A22
     MatAdd(SqMat_All[4], SqMat_All[8], 1);
     MatAdd(SqMat_All[7], SqMat_All[9], 1);

     // set rows of Dirichlet dof in off diagonal matrix blocks to zero    
     // get row in off diagonal matrix where the Dirichlet nodes start
     RowPtr = SqMat_All[5]->GetRowPtr();
     // compute number of entries starting from this row to the end of the matrix
     j = RowPtr[N_Active];
     k = RowPtr[N_U]-j;
     // get number of active dof
     // set these entries to zero
     memset(SqMat_All[5]->GetEntries()+j, 0, SizeOfDouble*k);
     memset(SqMat_All[6]->GetEntries()+j, 0, SizeOfDouble*k);

     // slip type bc detected, manipulation of matrices is necessary
     // this is done only at the very beginning
     // the matrices A_12, A_12, M_11, M_12, M_21, M_22, B1T, B2T
     //     stay unchanged during the complete solution process
     // the matrices A_11, A_22 are manipulated after their new
     //     assembling during the nonlinear iteration

       if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
        {
          if (TDatabase::ParamDB->NSTYPE <4)
          {
            OutPut("For slip with friction bc NSTYPE 4 is ");
            OutPut("necessary !!!!! " << endl);
            exit(4711);
          }

          // prepare everything for the assembling of slip with friction bc
          // on all levels
          N_FESpaces = 1;
          N_SquareMatrices = 8;
          N_RectMatrices = 2;
          N_Rhs = 2;
          DiscreteForm = NULL;

          SQMATRICES[0] = SqMat_All[4];
          SQMATRICES[1] = SqMat_All[7];
          SQMATRICES[2] = SqMat_All[5];
          SQMATRICES[3] = SqMat_All[6];
          SQMATRICES[4] = SqMat_All[0];
          SQMATRICES[5] = SqMat_All[3];
          SQMATRICES[6] = SqMat_All[1];
          SQMATRICES[7] = SqMat_All[2];

          MATRICES[0] = Mat_All[2];
          MATRICES[1] = Mat_All[3];

          fesp[0] = FESpaces_All[0];
          ferhs[0] = FESpaces_All[0];
          ferhs[1] = FESpaces_All[0];

          RHSs[0] = Rhs_All[0];
          RHSs[1] = Rhs_All[0]+N_U;

          aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);

          Assemble2DSlipBC(N_FESpaces, fesp,
                           N_SquareMatrices, SQMATRICES,
                           N_RectMatrices, MATRICES,
                           N_Rhs, RHSs, ferhs,
                           DiscreteForm,
                           BoundaryConditions,
                           BoundValues,
                           aux, FEFunctions_All[0], FEFunctions_All[1]);

       delete aux;

      } // if (TDatabase::ParamDB->INTERNA

     //    scale the pressure matrices
     Dscal(Mat_All[2]->GetN_Entries(), tau, Mat_All[2]->GetEntries());
     Dscal(Mat_All[3]->GetN_Entries(), tau, Mat_All[3]->GetEntries());
     Dscal(Mat_All[0]->GetN_Entries(), tau, Mat_All[0]->GetEntries());
     Dscal(Mat_All[1]->GetN_Entries(), tau, Mat_All[1]->GetEntries());
     
     // update rhs
     Daxpy(N_Active, tau, Rhs_All[0], B);
     Daxpy(N_Active, tau, Rhs_All[0]+N_U, B+N_U);

     // update rhs by Laplacian and convective term initialy by current time step
     // scaled by current sub time step length and theta2
     // currently : M := M + gamma A
     // M = M + (-gamma - tau*TDatabase::TimeDB->THETA2) A
     MatAdd(SqMat_All[0], SqMat_All[4], - tau*TDatabase::TimeDB->THETA2);
     MatAdd(SqMat_All[1], SqMat_All[5], - tau*TDatabase::TimeDB->THETA2);
     MatAdd(SqMat_All[2], SqMat_All[6], - tau*TDatabase::TimeDB->THETA2);
     MatAdd(SqMat_All[3], SqMat_All[7], - tau*TDatabase::TimeDB->THETA2);

     // set current factor of steady state matrix
     gamma = -tau*TDatabase::TimeDB->THETA2;		     

     // defect = M * Sol
     // B:= B + defect (rhs)     
     MatVectActive(SqMat_All[0], Sol_All[0], defect);
     Daxpy(N_Active, 1, defect, B);
     MatVectActive(SqMat_All[1], Sol_All[0]+N_U, defect);
     Daxpy(N_Active, 1, defect, B);
     MatVectActive(SqMat_All[2], Sol_All[0], defect+N_U);
     Daxpy(N_Active, 1, defect+N_U, B+N_U);
     MatVectActive(SqMat_All[3], Sol_All[0]+N_U, defect+N_U);
     Daxpy(N_Active, 1, defect+N_U, B+N_U);     
     
    // set Dirichlet values
    // RHSs[0] still available from assembling
    memcpy(B+N_Active, RHSs[0]+N_Active, (N_U-N_Active)*SizeOfDouble);
    memcpy(B+N_Active+N_U, RHSs[1]+N_Active,(N_U-N_Active)*SizeOfDouble);

    // copy Dirichlet values from rhs into Sol[0][mg_level-1]
    memcpy(Sol_All[0]+N_Active, RHSs[0]+N_Active, (N_U-N_Active)*SizeOfDouble);
    memcpy(Sol_All[0]+N_Active+N_U, RHSs[1]+N_Active, (N_U-N_Active)*SizeOfDouble);

    //=====================================================================
    // the stiffness matrix is stored on M11, (M12, M21, M22)
    // assembling of system matrix
    //========================================================================
    // M = M + (-gamma + tau*TDatabase::TimeDB->THETA1) A
     MatAdd(SqMat_All[0], SqMat_All[4], -gamma + tau*TDatabase::TimeDB->THETA1);
     MatAdd(SqMat_All[1], SqMat_All[5], -gamma + tau*TDatabase::TimeDB->THETA1);
     MatAdd(SqMat_All[2], SqMat_All[6], -gamma + tau*TDatabase::TimeDB->THETA1);
     MatAdd(SqMat_All[3], SqMat_All[7], -gamma + tau*TDatabase::TimeDB->THETA1);
          
     // set current factor of steady state matrix
     gamma = tau*TDatabase::TimeDB->THETA1;
     
     OutPut(endl << "CURRENT TIME: ");
     OutPut(TDatabase::TimeDB->CURRENTTIME << endl);

     
     //======================================================================
     // nonlinear loop
     //======================================================================   
     N_LinIterCurr = 0;
     solver_time_curr = 0;     
     
     for(j=0;j<Max_It;j++)
      {
       memset(defect, 0, N_Unknowns*SizeOfDouble);

       SQMATRICES[0] = SqMat_All[0];
       SQMATRICES[1] = SqMat_All[1];
       SQMATRICES[2] = SqMat_All[2];
       SQMATRICES[3] = SqMat_All[3];
       MATRICES[0] = Mat_All[0];
       MATRICES[1] = Mat_All[1];
       MATRICES[2] = Mat_All[2];
       MATRICES[3] = Mat_All[3];      
       
      // compute defect
      Defect(sqmatrices, matrices, Sol_All[0], B, defect);

      residual =  Ddot(N_Unknowns, defect, defect);
      impuls_residual = Ddot(2*N_U, defect, defect);
      OutPut("TNSE nonlinear step " << setw(3) << j);
      OutPut(setw(14) << impuls_residual);
      OutPut(setw(14) << Ddot(N_P,defect+2*N_U,defect+2*N_U));
      OutPut(setw(14) << sqrt(residual));

      
      if(j>0)
       {
        OutPut(setw(14) << sqrt(residual)/oldresidual << endl);
       }
      else
       {
        OutPut(endl);
       }
       
      oldresidual = sqrt(residual);

      if ((((sqrt(residual)<=limit)||(j==Max_It-1)))  && (j>=TDatabase::ParamDB->SC_MINIT))
       {
        if (j==Max_It-1)
        j++;
//         OutPut("ITE : " << setw(3) << j);
//         OutPut(" (" << setw(3) << N_LinIterCurr << "/");
//         OutPut(setw(3) << N_LinIter << " LINITE)");
//         OutPut("  TIME FOR SOLVER : " << solver_time_curr << "/" << solver_time << "s");
//         OutPut("  RES : " <<  sqrt(residual) << endl);
        // count total running time
//         t4 =  GetTime();
//         total_time += t4 - t3;
//         t3 = t4;
//         OutPut(TDatabase::TimeDB->CURRENTTIME << " total elapsed running time "<< total_time << endl);
        break;
       }
//        cout<< "test " << endl;
//           exit(0);   
       //======================================================================
       // solve linear system
       //======================================================================
//         t1 = GetTime();
        DirectSolver(SQMATRICES[0], SQMATRICES[1], SQMATRICES[2], SQMATRICES[3],
                     MATRICES[2], MATRICES[3], MATRICES[0], MATRICES[1],
                     B, Sol_All[0]);
//         t2 = GetTime();
//         solver_time_curr = t2-t1;
//         solver_time += solver_time_curr;
//  cout<< "test " << endl;
       //======================================================================
       // end solve linear system
       //======================================================================
       // restore mass matrices by subtracting the A-matrices
       MatAdd(SqMat_All[0], SqMat_All[4], -gamma);
       MatAdd(SqMat_All[3], SqMat_All[7], -gamma);

       //======================================================================
       // assemble new matrix due to nonlinearity
       //======================================================================
        GridVelo_imping(Entries, tmp_Gsol, tmp_Gd, Rhs_All[1],
                        GridKCol, GridRowPtr,
                        GridPos, AuxGridPos,
                        FEVectFuncts_All[0], tau,
                        FEVectFuncts_All[1], MovBoundVert, N_MovVert,
                        Free_Cells, IsoCellEdgeNos, reparam, RefGridPos);
 
       DiscreteForm = DiscreteFormNLGalerkin;	 
       N_RectMatrices = 0;
       N_Rhs = 0;
       N_FESpaces = 3;

       SQMATRICES[0] = SqMat_All[4];
       SQMATRICES[1] = SqMat_All[7];
       SQMATRICES[0]->Reset();
       SQMATRICES[1]->Reset();

       N_SquareMatrices = 2;
       last_sq = 1;
       
       fesp[0] = FESpaces_All[0];
       fesp[1] = FESpaces_All[1];
       fesp[2] = FESpaces_All[2];

       fefct[0] = FEFunctions_All[0];
       fefct[1] = FEFunctions_All[1];
       fefct[2] = FEFunctions_All[3];
       fefct[3] = FEFunctions_All[4];
 
       //======================================================================
       // assembling of matrices for each level due to nonlinearity
       // A_11, (A_22)
       // no assembling of rhs
       //======================================================================
        aux =  new TAuxParam2D(MovingTNSN_FESpaces_Axial3D, MovingTNSN_Fct_Axial3D,
                               MovingTNSN_ParamFct_Axial3D,
                               MovingTNSN_FEValues_Axial3D,
                               fesp, fefct,
                               MovingTNSFct_Axial3D,
                               MovingTNSFEFctIndex_Axial3D,
                               MovingTNSFEMultiIndex_Axial3D,
                               MovingTNSN_Params_Axial3D, MovingTNSBeginParam_Axial3D);

         Assemble2D(N_FESpaces, fesp,
                    N_SquareMatrices, SQMATRICES,
                    N_RectMatrices, MATRICES,
                    N_Rhs, RHSs, ferhs,
                    DiscreteForm,
                    BoundaryConditions,
                    BoundValues,
                    aux); 
    
       // Adding freesurf entries to A11 and A22
       MatAdd(SqMat_All[4], SqMat_All[8], 1.);
       MatAdd(SqMat_All[7], SqMat_All[9], 1.); 
       

       if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
        {
          // prepare everything for the assembling of slip with friction bc
          // on all levels
          N_FESpaces = 1;
          N_SquareMatrices = 2;
          N_RectMatrices = 0;
          N_Rhs = 2;
          DiscreteForm = NULL;

          SQMATRICES[0] = SqMat_All[4];
          SQMATRICES[1] = SqMat_All[7];

          fesp[0] = FESpaces_All[0];
          ferhs[0] = FESpaces_All[0];
          ferhs[1] = FESpaces_All[0];

          RHSs[0] = Rhs_All[0];
          RHSs[1] = Rhs_All[0]+N_U;

          aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);

          Assemble2DSlipBC(N_FESpaces, fesp,
                           N_SquareMatrices, SQMATRICES,
                           N_RectMatrices, MATRICES,
                           N_Rhs, RHSs, ferhs,
                           DiscreteForm,
                           BoundaryConditions,
                           BoundValues,
                           aux, FEFunctions_All[0], FEFunctions_All[1]);

       delete aux;

      } // if (TDatabase::ParamDB->INTERNA
       
     //======================================================================
     // end of assemble new matrix due to nonlinearity
     //======================================================================       
     // build stiffness matrix for next nonlinear iteration step
     // stiffness matrix (left upper block) is stored on
     // M11, (M12, M21, M22)
     // M = M +  tau*TDatabase::TimeDB->THETA1 A      
       
     MatAdd(SqMat_All[0], SqMat_All[4], tau*TDatabase::TimeDB->THETA1);
     MatAdd(SqMat_All[3], SqMat_All[7], tau*TDatabase::TimeDB->THETA1);       
     } //   for(j=0;j<Max_It;j++)  
     
// #endif  // __SURFACTDIFFTEST__

// #endif // __SURFACTCONVTEST__  

   //======================================================================
   // end NSE nonlinear iteration     
   
   //======================================================================     
//    // get M mat values before moving the mesh
//     SqMat_IFace[1]->Reset();
// 
//       AssembleSurf1D_Mmat(FESpaces_All[4], IFaceFeSpaces[0], SqMat_IFace[1],  N_List[0], N_List[1]);
//       
// //       cout << " MatVal Ddot Before : " <<  setprecision(12) <<Ddot(N_Surf_M_MatValues,SqMat_IFace[1]->GetEntries(),SqMat_IFace[1]->GetEntries()) << endl;
//       memcpy(Surf_M_MatValues,  SqMat_IFace[1]->GetEntries(),  N_Surf_M_MatValues*SizeOfDouble);  
// //    exit(0);
   
   //======================================================================    
   // move the grid
#ifdef __SURFACT__      
//  ========================================================================================== 
// //    test
//       MapSurfToDomain(IFaceFeFunct[0], FEFunctions_All[7], N_List[0], N_List[1]);
//       GetSurfactMass(FEFunctions_All[7], IFaceFeFunct[0], N_List[0], N_List[1], Surf_Mass);
//       FEFunctions_All[6]->GetMassAndMean(Params);      
//   
//       old_c_mass = Surf_Mass[1];
//       old_c_Gamma_mass = Surf_Mass[0];
//       OutPut("solution A" << sqrt(Ddot(N_IsurfactDOF,Sol_All[4],Sol_All[4])) << endl);  
//  ==========================================================================================    
   // reparam includes surf surfact interpolation
   if(reparam)
    {
     MapSurfToDomain(IFaceFeFunct[0], FEFunctions_All[7], N_List[0], N_List[1]);
     SurfSurfactReparam = TRUE;
     
// #ifdef __SOLUBLE__  
     // surfactant Mass correction 
//      GetSurfactMass(FEFunctions_All[7], IFaceFeFunct[0], N_List[0], N_List[1], Surf_Mass);
//      FEFunctions_All[6]->GetMassAndMean(Params);      
     //save the old surfactant mass for correction
//      old_c_mass = Params[0];
//      old_c_Gamma_mass = Surf_Mass[0];   
//      SurfactMassCorrection = TRUE;
// #endif 
    }
#endif

   MoveGrid_imping(Entries, tmp_Gsol, tmp_Gd, Rhs_All[1],
                  GridKCol, GridRowPtr,
                  GridPos, FEVectFuncts_All[0], tau, FEVectFuncts_All[1],
                  AuxGridPos, MovBoundVert, N_MovVert,
                  Free_Cells, IsoCellEdgeNos,  
                  FEFunctions_All[6], FEFunctions_All[7],
                  reparam, N_ReParam);
   
#ifdef __SURFACT__    
    if(SurfSurfactReparam)
     {
      MapDomainToSurf(FEFunctions_All[7], IFaceFeFunct[0], N_List[0], N_List[1]);  
      SurfSurfactReparam = FALSE;  
      
      
// #ifdef __SOLUBLE__   
//      if(SurfactMassCorrection)
//       {
       // bulk and interface surfactant mass correction
//        FEFunctions_All[6]->CorrectMass(old_c_mass);     
//        SurfSurfactCorrectMass(FEFunctions_All[7], IFaceFeFunct[0], N_List[0], N_List[1], old_c_Gamma_mass); 
//        SurfactMassCorrection = FALSE;
//       }
// #endif     
     }
#endif    

  //=====================================================================
  // Updating the Quard points on the solid surfaces
  //======================================================================         
  GridPos->GridToData();
  RefGridPos->GridToData();
  
  // Updating solid boundary 
  MovBoundVert[0][0]->GetCoords(x, y);
  MovBoundVert[2][0]->GetCoords(SLPX, SLPY);  
  UpdateSlipBound->SetParams(x, y, SLPX-x, SLPY-y);
  for(k=0;k<N_MovVert[0];k++)
   if(k==N_MovVert[0]-1)
    { Bound_Joint[0][k]->UpdateParameters(MovBoundVert[0][k], MovBoundVert[2][0]); }
   else
    { Bound_Joint[0][k]->UpdateParameters(MovBoundVert[0][k], MovBoundVert[0][k+1]); }
         
  // Updating axis boundary
   MovBoundVert[1][0]->GetCoords(x, y);
   MovBoundVert[0][0]->GetCoords(SLPX, SLPY);
   UpdateAxialBound->SetParams(x, y, SLPX-x, SLPY-y);
   for(k=0;k<N_MovVert[1];k++)
    if(k==N_MovVert[1]-1)
     { Bound_Joint[1][k]->UpdateParameters(MovBoundVert[1][k], MovBoundVert[0][0]); }
    else
     { Bound_Joint[1][k]->UpdateParameters(MovBoundVert[1][k], MovBoundVert[1][k+1]); }

  //======================================================================          
  // end Updating the Quard points on the solid surface
  //======================================================================   
  
#ifdef __SURFACT__ 
   //======================================================================
   // soluble surfactants --- begin
   //======================================================================
//  OutPut("W old" << sqrt(Ddot(2*N_G,Sol_All[1],Sol_All[1])) << endl);       
   // Grid velocity calculated during move mesh
     GridVelo_imping(Entries, tmp_Gsol, tmp_Gd, Rhs_All[1],
                     GridKCol, GridRowPtr,
                     GridPos, AuxGridPos,
                     FEVectFuncts_All[0], tau,
                     FEVectFuncts_All[1], MovBoundVert, N_MovVert,
                     Free_Cells, IsoCellEdgeNos, reparam, RefGridPos);   
 
//  OutPut("w new" << sqrt(Ddot(2*N_G,Sol_All[1],Sol_All[1])) << endl);       
   // save the old time step solution
   memcpy(Csol_old, Sol_All[3], N_surfactDOF*SizeOfDouble);
   memcpy(Isol_old, Sol_All[4], N_IsurfactDOF*SizeOfDouble);  
      
   for(j=0;j<Max_It_scalar;j++)
   {
#ifdef __SOLUBLE__     
    memcpy(Csol_nonlinearstep,  Sol_All[3], N_surfactDOF*SizeOfDouble);
    
    //======================================================================
    // surfactant in bulk phase - begin
    //======================================================================    
    // assembling matrices
    SQMATRICES_SURFACT[0] = SqMat_All[17]; // A
    SQMATRICES_SURFACT[0]->Reset();
    SQMATRICES_SURFACT[1] = SqMat_All[16]; // M
    SQMATRICES_SURFACT[1]->Reset();
 
    fesp[0] = FESpaces_All[4]; // surfactant space
    fesp[1] = FESpaces_All[0];  // velocity space
    fesp[2] = FESpaces_All[2];  // mesh velocity space

    ferhs[0] = FESpaces_All[4]; // surfactant space for rhs
 
    fefct[0] = FEFunctions_All[0]; // u1
    fefct[1] = FEFunctions_All[1]; // u2
    fefct[2] = FEFunctions_All[3]; // w1
    fefct[3] = FEFunctions_All[4]; // w2
    
    CRHSs[0] =   Rhs_All[3];

    memset(CRHSs[0], 0, N_surfactDOF*SizeOfDouble);

     // (u1-w1, u2-w2)  parameters are needed for assembling
     // fesp is taken from fefct in aux
    aux =  new TAuxParam2D(MovingTNSN_FESpaces_Axial3D, MovingTNSN_Fct_Axial3D,
                           MovingTNSN_ParamFct_Axial3D,
                           MovingTNSN_FEValues_Axial3D,
                           fesp+1, fefct,
                           MovingTNSFct_Axial3D,
                           MovingTNSFEFctIndex_Axial3D,
                           MovingTNSFEMultiIndex_Axial3D,
                           MovingTNSN_Params_Axial3D, MovingTNSBeginParam_Axial3D);

     Assemble2D(3, fesp,
                2, SQMATRICES_SURFACT,
                0, NULL,
                1, CRHSs, ferhs,
                DiscreteFormSurfact,
                SurfactBoundaryConditions,
                SurfactBoundValues,
                aux);
     delete aux;    
        
     switch(surf_couple_var)
     {
      case 1: // explicit j=0 only, no iteration
       Surfact2D_InterfaceInt(SqMat_All[17], CRHSs[0], SurfactBoundCondition, FEFunctions_All[6],
                          IFaceFeFunct[0], FESpaces_All[4], N_List[0], N_List[1], GammaXmaxVal);
//        memcpy(CRhs_old, CRHSs[0], N_surfactDOF*SizeOfDouble);
      break;

      case 2: // implicit but explicit in fixed pt iteration step
      case 3:
       Surfact2D_InterfaceInt(SqMat_All[17], CRHSs[0], SurfactBoundCondition, FEFunctions_All[6],
                          IFaceFeFunct[0], FESpaces_All[4], N_List[0], N_List[1], GammaXmaxVal);
      break;

      case 4: // implicit but explicit in fixed pt iterartion step
      case 5: // fully implicit
       Surfact2D_InterfaceInt_Implicit(SqMat_All[17], CRHSs[0], SurfactBoundCondition, FEFunctions_All[6],
                                       IFaceFeFunct[0], FESpaces_All[4], N_List[0], N_List[1], GammaXmaxVal);
      break;

      default:
       OutPut("error in selecting linerizatin type for coupled surfaact eqns " << endl);
       exit(1);
     }


   // working rhs for surfactant
   memset(C_B, 0, N_surfactDOF*SizeOfDouble);     

     //if(TDatabase::TimeDB->THETA3 !=0)
    if(j==0)
     {
      // rhs from old sol  
      Daxpy(N_surfactActive, tau*TDatabase::TimeDB->THETA3, CRHSs[0], C_B);       
      //for given c^Gamma_{old}
      if(TDatabase::TimeDB->THETA2>0.)
       MatAdd(SqMat_All[16], SqMat_All[17], -tau*TDatabase::TimeDB->THETA2);
      
      gamma = -tau*TDatabase::TimeDB->THETA2;   
   
      memset(C_defect, 0, N_surfactDOF*SizeOfDouble);
      MatVectActive(SqMat_All[16], Csol_old, C_defect); 
      Daxpy(N_surfactActive, 1, C_defect, C_B); 
      
      // save the old rhs for further iterations
      memset(CRhs_old, 0, N_surfactDOF*SizeOfDouble); 
      memcpy(CRhs_old, C_B, N_surfactActive*SizeOfDouble);   
     }
    else
     {
      //copy rhs from j=0
      memcpy(C_B, CRhs_old, N_surfactActive*SizeOfDouble);
      gamma = 0.;
     }
    
   // copy the currtent iterative rhs,  
   Daxpy(N_surfactActive, tau*TDatabase::TimeDB->THETA4, CRHSs[0], C_B);
 
   // set Dirichlet values
   if( (N_surfactDOF - N_surfactActive)>0)
    {
     memcpy(C_B+N_surfactActive, CRHSs[0]+N_surfactActive, (N_surfactDOF - N_surfactActive)*SizeOfDouble);
     memcpy(Sol_All[3]+N_surfactActive, CRHSs[0]+N_surfactActive, (N_surfactDOF - N_surfactActive)*SizeOfDouble);
    }
    
   // assemble system matrix
   MatAdd(SqMat_All[16], SqMat_All[17], -gamma + tau*TDatabase::TimeDB->THETA1);

   // check the convergence of the fixed point linerization
   memset(C_defect, 0, N_surfactDOF*SizeOfDouble);
   ScalarDefect(SqMat_All[16], Sol_All[3], C_B, C_defect, residual_scalar);

   OutPut("Scalar nonlinear step " << setw(3) << j);
   OutPut(setw(14) << residual_scalar); // sqrt of residual_scalar is alread done in ScalarDefect
   
 
   if (j>0)
    { 
      if(fabs(oldresidual_scalar)>0)
      OutPut(setw(14) <<  residual_scalar/oldresidual_scalar ); 
      OutPut(endl);
    }
   else
    { OutPut(endl); }

   oldresidual_scalar = residual_scalar;
   if( ((residual_scalar<=limit_scalar)||(j==Max_It_scalar-1))  && (j>=TDatabase::ParamDB->SC_MINIT))
    { 
     break;
    }

    //solve the system
    DirectSolver(SqMat_All[16], C_B, Sol_All[3]);
   
#endif //__INSOLUBLE__
    //======================================================================
    //  surfactant in bulk phase - end
    //  surfactant on interphase phase - begin
    //======================================================================  
 
    //  assembling surf surfactant matrices
    N_SquareMatrices = 2;
    N_FESpaces = 2;
    N_Rhs =1;
    N_FESpaces_low=1;  
    
    SQMATRICES_IFace[0] = SqMat_IFace[0];
    SQMATRICES_IFace[0]->Reset();
    SQMATRICES_IFace[1] = SqMat_IFace[1];
    SQMATRICES_IFace[1]->Reset();

    fesp[0] = FESpaces_All[2]; // mesh velospace space 
    fesp[1] = FESpaces_All[4];  // bulk surfactant space 

    IFacefesp[0] = IFaceFeSpaces[0]; // Interface surfactant space
    IFaceferhs[0] = IFaceFeSpaces[0]; // Interface surfactant space for rhs

    //mesh velocity, even though interface moves in a Lagrangian manner due to reparam w is diff from u on interfaces
    fefct[0] = FEFunctions_All[3]; // wr
    fefct[1] = FEFunctions_All[4]; // wz

    SRHSs[0] =  Rhs_All[4];
    memset(SRHSs[0], 0, N_IsurfactDOF*SizeOfDouble);
  
     // working array for srhs is I_B, initialize I_B
     memset(I_B, 0,  N_IsurfactDOF*SizeOfDouble);    
    
     if(j==0)
      {
       switch(surf_couple_var)
        {
         case 1: // all source terms on rhs, using old solution
         case 2:  
         case 3:  
          AssembleSurf1D_SolubleSurfact(N_FESpaces, fesp, fefct, N_FESpaces_low,
                   IFacefesp, IFaceFeFunct[0], N_SquareMatrices, SQMATRICES_IFace, N_Rhs, SRHSs, 
                   IFaceferhs, N_List[0], N_List[1], Csol_old);
        break;

        case 4: // c^\Gamma terms in source are in LHS, but the solution  
        case 5: // c^\Gamma terms in source are in LHS, solution (c_old) 
         AssembleSurf1D_SolubleSurfact_Implicit(N_FESpaces, fesp, fefct, N_FESpaces_low,
                     IFacefesp, Isol_old, N_SquareMatrices, SQMATRICES_IFace, N_Rhs, SRHSs, 
                     IFaceferhs, N_List[0], N_List[1], Csol_old, tau); 
        break;
        default:
         OutPut("error in selecting linerizatin type for coupled surfaact eqns " << endl);
         exit(1);
        }     
        
      // copy the rhs, for given c^{old}
      Daxpy(N_IActive, tau*TDatabase::TimeDB->THETA3, SRHSs[0], I_B);     
     
      //for given c^{old},
      MatAdd(SqMat_IFace[1], SqMat_IFace[0], -tau*TDatabase::TimeDB->THETA2);
      memset(I_defect, 0, N_IsurfactDOF*SizeOfDouble);
      MatVectActive(SqMat_IFace[1], Isol_old, I_defect);
      Daxpy(N_IActive, 1, I_defect, I_B);  
      
      // save it for further iterative steps
      memset(IRhs_old, 0, N_IsurfactDOF*SizeOfDouble);
      memcpy(IRhs_old, I_B, N_IActive*SizeOfDouble);
      
      //reset the mat and rhs for assembling
      SQMATRICES_IFace[0]->Reset();
      SQMATRICES_IFace[1]->Reset();
      memset(SRHSs[0], 0, N_IsurfactDOF*SizeOfDouble);
     }
    else
     {
      // copy the rhs from j==0
      memcpy(I_B, IRhs_old, N_IActive*SizeOfDouble);
     }  // if( j==0 

    //assemble with the new iterative values
    switch(surf_couple_var)
     {
      case 1: // all source terms on rhs, using old solution
       AssembleSurf1D_SolubleSurfact(N_FESpaces, fesp, fefct, N_FESpaces_low,
                   IFacefesp, IFaceFeFunct[0], N_SquareMatrices, SQMATRICES_IFace, N_Rhs, SRHSs, 
                   IFaceferhs, N_List[0], N_List[1], Csol_old);
      break;

      case 2: // all source terms on rhs, but the solution (c_{i-1}) from previous iteration
       AssembleSurf1D_SolubleSurfact(N_FESpaces, fesp, fefct, N_FESpaces_low,
                   IFacefesp, IFaceFeFunct[0], N_SquareMatrices, SQMATRICES_IFace, N_Rhs, SRHSs, 
                   IFaceferhs, N_List[0], N_List[1], Csol_nonlinearstep);

      break;

      case 3: //  all source terms on rhs, solution (c_i) from current iteration
       AssembleSurf1D_SolubleSurfact(N_FESpaces, fesp, fefct, N_FESpaces_low,
                   IFacefesp, IFaceFeFunct[0], N_SquareMatrices, SQMATRICES_IFace, N_Rhs, SRHSs, 
                   IFaceferhs, N_List[0], N_List[1], Sol_All[3]);

      break;

      case 4: // c^\Gamma terms in source are in LHS, but the bulk solution  (c_{i-1}) from previous iteration
       AssembleSurf1D_SolubleSurfact_Implicit(N_FESpaces, fesp, fefct, N_FESpaces_low,
                   IFacefesp, Isol_old, N_SquareMatrices, SQMATRICES_IFace, N_Rhs, SRHSs, 
                   IFaceferhs, N_List[0], N_List[1], Csol_nonlinearstep, tau);
      break;

      case 5: // c^\Gamma terms in source are in LHS, bulk solution (c_i) from current iteration
       AssembleSurf1D_SolubleSurfact_Implicit(N_FESpaces, fesp, fefct, N_FESpaces_low,
                     IFacefesp, Isol_old, N_SquareMatrices, SQMATRICES_IFace, N_Rhs, SRHSs, 
                     IFaceferhs, N_List[0], N_List[1], Sol_All[3], tau); 
      break;
      default:
       OutPut("error in selecting linerizatin type for coupled surfaact eqns " << endl);
       exit(1);
     }

    // copy the currtent iterative rhs, for given c^{n+1}
    Daxpy(N_IActive, tau*TDatabase::TimeDB->THETA4, SRHSs[0], I_B);
    
    // set Dirichlet values
    if( (N_IsurfactDOF - N_IActive)>0)
     {
      memcpy(I_B+N_IActive, SRHSs[0]+N_IActive, (N_IsurfactDOF - N_IActive)*SizeOfDouble);
      memcpy(Sol_All[4]+N_IActive, SRHSs[0]+N_IActive, (N_IsurfactDOF - N_IActive)*SizeOfDouble);
     }
   //=====================================================================
   // assembling of system matrix
   //======================================================================== 
    MatAdd(SqMat_IFace[1], SqMat_IFace[0], tau*TDatabase::TimeDB->THETA1);      

    //solve the system
    DirectSolver(SqMat_IFace[1], I_B, Sol_All[4]);
 
   //======================================================================
   // surfactant on interphase phase - end
   //======================================================================  
  } //  for(j=0;j<Max_It_scalar;j++)  
#endif      

    //  ========================================================================================== 
//     //    test
//       MapSurfToDomain(IFaceFeFunct[0], FEFunctions_All[7], N_List[0], N_List[1]);
//       GetSurfactMass(FEFunctions_All[7], IFaceFeFunct[0], N_List[0], N_List[1], Surf_Mass);
//       FEFunctions_All[6]->GetMassAndMean(Params);      
//    
//       OutPut( "Time, TransportedMass " <<TDatabase::TimeDB->CURRENTTIME<<
//               " " << ( Surf_Mass[1]  - old_c_mass) << 
//               " " <<TDatabase::ParamDB->REACTOR_P14*(Surf_Mass[0]-old_c_Gamma_mass)<< endl);
// //       
// //    OutPut("solution B" << sqrt(Ddot(N_IsurfactDOF,Sol_All[4],Sol_All[4])) << endl);       
// 
//       
//         OutPut( "Time, Surfactant_Mass, Da*InterfaceSurfactant_Mass, InterfaceSurfactant_Conc " <<TDatabase::TimeDB->CURRENTTIME<<
//                " " <<Params[0]<< " "<< TDatabase::ParamDB->REACTOR_P14*Surf_Mass[0]<< " "<<Surf_Mass[0]/Surf_Mass[1]<< " "<<endl);     
//       OutPut( "Time, Surfactant_Mass_dif, InterfaceSurfactant_Mass_diff " <<TDatabase::TimeDB->CURRENTTIME<<
//               " " <<Params[0] - Initial_SurfactMass<< " "<<TDatabase::ParamDB->REACTOR_P14*Surf_Mass[0]-Initial_IFaceSurfactMass<< endl);
//        OutPut( "Time, Total Mass, Total Mass diff, RelativeTotal Mass diff, sphericity, KE " 
//               <<TDatabase::TimeDB->CURRENTTIME
//               <<" " <<Params[0]+TDatabase::ParamDB->REACTOR_P14*Surf_Mass[0]<< " " << ((Params[0]+TDatabase::ParamDB->REACTOR_P14*Surf_Mass[0]) - 
//             (Initial_IFaceSurfactMass + Initial_SurfactMass)) << " " << ((Params[0]+TDatabase::ParamDB->REACTOR_P14*Surf_Mass[0]) - 
//             (Initial_IFaceSurfactMass + Initial_SurfactMass))
//             /(Initial_IFaceSurfactMass + Initial_SurfactMass) <<" " << sphericity << " " << KE<<endl);       
//      
//           
//       
      
      
   //  ========================================================================================== 
//       exit(0);
      
      
// #ifdef __SURFACTDIFFTEST__
//      OutPut(endl << "CURRENT TIME: ");
//      OutPut(TDatabase::TimeDB->CURRENTTIME << endl);
//      
//      
//      continue; // cont to next substep l
// #endif  

// #ifdef __SURFACTCONVTEST__
//      OutPut(endl << "CURRENT TIME: ");
//      OutPut(TDatabase::TimeDB->CURRENTTIME << endl);
//      
//      continue; // cont to next substep l, no mesh movement
// #endif  

// cout<< "test Main " << endl;
// exit(0);

 
  //======================================================================   
  // check freesurf point on the solid surface
  // if so, change bounddary description
  //======================================================================                
  MovBoundVert[2][1]->GetCoords(SLPX, SLPY);
    
  if(SLPY <= 1e-8  )
   {              
    SLPY = 0.;
    MovBoundVert[0][0]->GetCoords(x, y);    
    MovBoundVert[2][1]->SetCoords(SLPX, SLPY);
    UpdateSlipBound->SetParams(x, y, SLPX-x, SLPY);

    Me = Free_Cells[0];
    Joint = Me->GetJoint(IsoCellEdgeNos[1][0]);
      
    IsoJoint = (TIsoBoundEdge *)Joint;
    IsoJoint->DeleteVertices();
    IsoJoint->ChangeEdgeID(BoundaryEdge);      

    Solid_Joint = (TBoundEdge *)Joint;
    Solid_Joint->ChangeBoundComp((Domain->GetBdPart(0))->GetBdComp(0));
    Solid_Joint->UpdateParameters(MovBoundVert[2][0], MovBoundVert[2][1]);
    Joint->ChangeType(BoundaryEdge);
    BoundComp = Solid_Joint->GetBoundComp();
    BoundComp->ChangeType(Line); 
    OutPut("Free surface first vert changed as solid vert"<< " x : "<< SLPX<< "  y : "<<SLPY<< endl); 
 
     cout<<"GetID "<< BoundComp->GetID() << endl;
//     exit(0);
    
    
    // impose no-penetration condition
    FeId =  FESpaces_All[0]->GetFE2D(IsoCellEdgeNos[0][0], Me);
    FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FeId);
    JointDOF = FeDesc->GetJointDOF(IsoCellEdgeNos[1][0]); 
    N_DOF = FeDesc->GetN_JointDOF();
    DOF = GlobalNumbers+BeginIndex[IsoCellEdgeNos[0][0]];
     
    for(k=0;k<N_DOF;k++)
      Sol_All[0][N_U + DOF[JointDOF[k]] ] = 0.;       
    
    delete [] Bound_Joint[0];
    delete [] Bound_Joint[1];                      

    delete [] MovBoundVert[0];
    delete [] MovBoundVert[1];   
    delete [] MovBoundVert[2];  

    delete [] Free_Joint;    
    delete [] Free_Cells;  
    delete [] IsoCellEdgeNos[0];  
    delete [] IsoCellEdgeNos[1];  
    
    GetMovingBoundData(coll, N_MovVert, Bound_Joint, MovBoundVert, Free_Joint,
                       Free_Cells, IsoCellEdgeNos, SLPX, y);
    
   //=====================================================================
   // Updating the Quard points on the solid surfaces
   //======================================================================         
   GridPos->GridToData();
   RefGridPos->GridToData();
  
   // Updating solid boundary 
   MovBoundVert[0][0]->GetCoords(x, y);
   MovBoundVert[2][0]->GetCoords(SLPX, SLPY);  
   UpdateSlipBound->SetParams(x, y, SLPX-x, SLPY-y);
   for(k=0;k<N_MovVert[0];k++)
    if(k==N_MovVert[0]-1)
     { Bound_Joint[0][k]->UpdateParameters(MovBoundVert[0][k], MovBoundVert[2][0]); }
    else
     { Bound_Joint[0][k]->UpdateParameters(MovBoundVert[0][k], MovBoundVert[0][k+1]); }
         
    // Updating axis boundary
    MovBoundVert[1][0]->GetCoords(x, y);
    MovBoundVert[0][0]->GetCoords(SLPX, SLPY);
    UpdateAxialBound->SetParams(x, y, SLPX-x, SLPY-y);
    for(k=0;k<N_MovVert[1];k++)
     if(k==N_MovVert[1]-1)
      { Bound_Joint[1][k]->UpdateParameters(MovBoundVert[1][k], MovBoundVert[0][0]); }
     else
      { Bound_Joint[1][k]->UpdateParameters(MovBoundVert[1][k], MovBoundVert[1][k+1]); }

   //======================================================================          
   // end Updating the Quard points on the solid surface
   //======================================================================       
   }  // if( SLPY <= 1e-8  )
       
  //======================================================================          
  // end change boundary description     
  // Remeshing Begin 
  //======================================================================      
//    if((l==0) && ((m % 1) == 0)) // if reparam or BD change ver happen the angle may become 0
    {
     Getcellangle(FESpaces_All[2], Angle);

    }
    
   if( Angle[0]<10.0 || Angle[1]>165.0) // || fabs(TDatabase::TimeDB->CURRENTTIME-0.005)<1e-5
    {
     OutPut( "MinAngle : "<< Angle[0]<< "  MaxAngle : "<<Angle[1]<< endl);    
#ifdef __SURFACT__  
     // interface surfactant
     MapSurfToDomain(IFaceFeFunct[0], FEFunctions_All[7], N_List[0], N_List[1]);    
     
#ifdef __SOLUBLE__  
//       Mass correction 
     GetSurfactMass(FEFunctions_All[7], IFaceFeFunct[0], N_List[0], N_List[1], Surf_Mass);
     FEFunctions_All[6]->GetMassAndMean(Params);      
     //save the old surfactant mass for correction
     old_c_mass = Params[0];
     old_c_Gamma_mass = Surf_Mass[0];
      
       OutPut( "Time, Surfactant_Mass, Da*InterfaceSurfactant_Mass, InterfaceSurfactant_Conc " <<TDatabase::TimeDB->CURRENTTIME<<
               " " <<Params[0]<< " "<< TDatabase::ParamDB->REACTOR_P14*Surf_Mass[0]<< " "<<Surf_Mass[0]/Surf_Mass[1]<< " "<<endl);     
      OutPut( "Time, Surfactant_Mass_dif, InterfaceSurfactant_Mass_diff " <<TDatabase::TimeDB->CURRENTTIME<<
              " " <<Params[0] - Initial_SurfactMass<< " "<<TDatabase::ParamDB->REACTOR_P14*Surf_Mass[0]-Initial_IFaceSurfactMass<< endl);
       OutPut( "Time, Total Mass, Total Mass diff, RelativeTotal Mass diff, sphericity, KE " <<TDatabase::TimeDB->CURRENTTIME<<
               " " <<Params[0]+TDatabase::ParamDB->REACTOR_P14*Surf_Mass[0]<< " " << ((Params[0]+TDatabase::ParamDB->REACTOR_P14*Surf_Mass[0]) - 
            (Initial_IFaceSurfactMass + Initial_SurfactMass)) << " " << ((Params[0]+TDatabase::ParamDB->REACTOR_P14*Surf_Mass[0]) - 
            (Initial_IFaceSurfactMass + Initial_SurfactMass))
            /(Initial_IFaceSurfactMass + Initial_SurfactMass) <<" " << sphericity << " " << KE<<endl);       
     
     
     
#endif       
     
     RemeshAxial3D_ImpDrop(Domain, FESpaces_All, FEVectFuncts_All, FEFunctions_All,
                           N_MovVert, Bound_Joint, MovBoundVert, Free_Joint, Free_Cells, IsoCellEdgeNos,
                           Sol_All, Rhs_All, SquareStructure_All, Structure_All, SqMat_All, Mat_All, 
                           IFaceDomain, IFaceFeSpaces, N_List, IFaceFeFunct, SqMat_IFace,
                           IFaceStruct, FE1D_List);      
#else

     RemeshAxial3D_ImpDrop(Domain, FESpaces_All, FEVectFuncts_All, FEFunctions_All,
                           N_MovVert, Bound_Joint, MovBoundVert, Free_Joint, Free_Cells, IsoCellEdgeNos,
                           Sol_All, Rhs_All, SquareStructure_All, Structure_All, SqMat_All, Mat_All, 
                           NULL, NULL, NULL, NULL, NULL,
                           NULL, NULL);
#endif 
//       t1 = GetTime();

   
     
#ifdef __SURFACT__     
     // copy to interface after bulk interpolation
     MapDomainToSurf(FEFunctions_All[7], IFaceFeFunct[0], N_List[0], N_List[1]);  
     N_S =  FESpaces_All[5]->GetN_DegreesOfFreedom();        
     memset(Sol_All[5], 0, N_S*SizeOfDouble);    
#endif 
     
     coll=FESpaces_All[0]->GetCollection();
     N_Cells = coll->GetN_Cells();
  
     Getcellangle(FESpaces_All[0], Angle);
      OutPut( "MinAngle : "<< Angle[0]<< "  MaxAngle : "<<Angle[1]<< endl);
   
     N_Remesh ++;
     remeshed=TRUE;

     GlobalNumbers = FESpaces_All[0]->GetGlobalNumbers();
     BeginIndex = FESpaces_All[0]->GetBeginIndex();
     N_Active =  FESpaces_All[0]->GetActiveBound();
     N_U = FESpaces_All[0]->GetN_DegreesOfFreedom();
     N_P = FESpaces_All[1]->GetN_DegreesOfFreedom();
     N_Unknowns = 2*N_U + N_P;
  
     N_G = FESpaces_All[2]->GetN_DegreesOfFreedom();
     N_GActive = FESpaces_All[2]->GetActiveBound();
     N_GBoundaryNodes = N_G - N_GActive;     
 
     delete [] B; delete [] defect; 
     B = new double[N_Unknowns];
     defect = new double[N_Unknowns];

#ifdef __ENERGY__    
     N_thermalDOF = FESpaces_All[3]->GetN_DegreesOfFreedom();
     N_thermalActive = FESpaces_All[3]->GetActiveBound();
     N_thermalNonActive = N_thermalDOF - N_thermalActive;  
     
     delete [] oldsol_T;
     oldsol_T = new double[N_thermalDOF];       
#endif       
  
#ifdef __SURFACT__     
   N_surfactDOF = FESpaces_All[4]->GetN_DegreesOfFreedom();
   N_surfactActive = FESpaces_All[4]->GetActiveBound();
   N_surfactNonActive = N_surfactDOF - N_surfactActive;     
   N_IsurfactDOF = IFaceFeSpaces[0]->GetN_DegreesOfFreedom();  
   N_IActive = IFaceFeSpaces[0]->GetActiveBound();    
   
   delete [] Csol_old;
   delete [] CRhs_old;
   delete [] C_B;
   delete [] C_defect;
   delete [] Csol_nonlinearstep;
   delete [] I_defect;
   delete [] Isol_old;
   delete [] I_B;
   delete [] IRhs_old;
    
   Csol_old = new double[N_surfactDOF];
   CRhs_old = new double[N_surfactDOF];
   C_B = new double[N_surfactDOF];  
   C_defect = new double[N_surfactDOF];  
   Csol_nonlinearstep = new double[N_surfactDOF];
   I_defect = new double[N_IsurfactDOF];
   Isol_old = new double[N_IsurfactDOF];
   I_B = new double[N_IsurfactDOF];
   IRhs_old = new double[N_IsurfactDOF];  
#endif
     
     GridKCol = SquareStructure_All[1]->GetKCol();
     GridRowPtr = SquareStructure_All[1]->GetRowPtr();
          
     delete [] refpos; delete []  auxpos;  delete [] pos;   delete [] ReparamMeshVelo;
     delete [] ReparamDisp; delete [] tmp_Gsol; delete [] tmp_Gd;

     tmp_Gd = new double[2*N_G];     
     tmp_Gsol = new double[2*N_G];
     refpos = new double[2*N_G];
     auxpos = new double[2*N_G];
     pos = new double[2*N_G];
     ReparamMeshVelo = new double[2*N_G];
     ReparamDisp = new double[2*N_G]; 
     
     delete RefGridPos; delete AuxGridPos; delete GridPos; delete ReparamPos;  

     ReparamPos = new TFEVectFunct2D(FESpaces_All[2], WString, WString, ReparamDisp, N_G, 2);       
     RefGridPos = new TFEVectFunct2D(FESpaces_All[2], WString, WString, refpos, N_G, 2);
     AuxGridPos = new TFEVectFunct2D(FESpaces_All[2], WString, WString, auxpos, N_G, 2);
     GridPos = new TFEVectFunct2D(FESpaces_All[2], WString, WString, pos, N_G, 2);

     delete Output;
     // prepare output (maxn_fespaces,  maxn_scalar,  maxn_vect, maxn_parameters, domain)
     Output = new TOutput2D(1, 3, 1, 2, Domain);
     Output->AddFEVectFunct(FEVectFuncts_All[0]);
     Output->AddFEFunction(FEFunctions_All[2]);      
#ifdef __SURFACT__      
  Output->AddFEFunction(FEFunctions_All[6]);
  Output->AddFEFunction(FEFunctions_All[7]);
#endif            
     
      if(TDatabase::ParamDB->WRITE_VTK)
       { 
#ifdef __SURFACT__   
        // interface surfactant
        MapSurfToDomain(IFaceFeFunct[0], FEFunctions_All[7], N_List[0], N_List[1]);
#endif
        os.seekp(std::ios::beg);
        if(img<10) os << "VTK/"<<VtkBaseName<<"_remesh.0000"<<RemeshImg<<".vtk" << ends;
         else if(img<100) os <<"VTK/"<< VtkBaseName<<"_remesh.000"<<RemeshImg<<".vtk" << ends;
          else if(img<1000) os <<"VTK/"<< VtkBaseName<<"_remesh.00"<<RemeshImg<<".vtk" << ends;
           else if(img<10000) os <<"VTK/"<< VtkBaseName<<"_remesh.0"<<RemeshImg<<".vtk" << ends;
            else  os << "VTK/"<<VtkBaseName<<"_remesh."<<RemeshImg<<".vtk" << ends;
        Output->WriteVtk(os.str().c_str());
        RemeshImg++;
       }    
       
#ifdef __SOLUBLE__           
      //======================================================================       
      // bulk surfactant mass correction
      FEFunctions_All[6]->CorrectMass(old_c_mass);
      
//       interface surfactant mass correction
      SurfSurfactCorrectMass(FEFunctions_All[7], IFaceFeFunct[0], N_List[0], N_List[1], old_c_Gamma_mass);
      
      MapSurfToDomain(IFaceFeFunct[0], FEFunctions_All[7], N_List[0], N_List[1]);     
      GetSurfactMass(FEFunctions_All[7], IFaceFeFunct[0], N_List[0], N_List[1], Surf_Mass);
      FEFunctions_All[6]->GetMassAndMean(Params);       
       
              OutPut( "Time, Surfactant_Mass, Da*InterfaceSurfactant_Mass, InterfaceSurfactant_Conc " <<TDatabase::TimeDB->CURRENTTIME<<
               " " <<Params[0]<< " "<< TDatabase::ParamDB->REACTOR_P14*Surf_Mass[0]<< " "<<Surf_Mass[0]/Surf_Mass[1]<< " "<<endl);     
      OutPut( "Time, Surfactant_Mass_dif, InterfaceSurfactant_Mass_diff " <<TDatabase::TimeDB->CURRENTTIME<<
              " " <<Params[0] - Initial_SurfactMass<< " "<<TDatabase::ParamDB->REACTOR_P14*Surf_Mass[0]-Initial_IFaceSurfactMass<< endl);
       OutPut( "Time, Total Mass, Total Mass diff, RelativeTotal Mass diff, sphericity, KE " <<TDatabase::TimeDB->CURRENTTIME<<
               " " <<Params[0]+TDatabase::ParamDB->REACTOR_P14*Surf_Mass[0]<< " " << ((Params[0]+TDatabase::ParamDB->REACTOR_P14*Surf_Mass[0]) - 
            (Initial_IFaceSurfactMass + Initial_SurfactMass)) << " " << ((Params[0]+TDatabase::ParamDB->REACTOR_P14*Surf_Mass[0]) - 
            (Initial_IFaceSurfactMass + Initial_SurfactMass))
            /(Initial_IFaceSurfactMass + Initial_SurfactMass) <<" " << sphericity << " " << KE<<endl);  
       
       
//        exit(0);
       //======================================================================  
#endif 

    } // if((Angle[0]<10.0) ||(Angle[1]>165.0))
  //======================================================================  
  // end Remeshing Begin 
  // Assembeling the grid matrix - Begin
  //======================================================================  

    fesp[0] = FESpaces_All[2];
    SQMATRICES_GRID[0] = SqMat_All[10];
    SQMATRICES_GRID[0]->Reset();
    SQMATRICES_GRID[1] = SqMat_All[11];
    SQMATRICES_GRID[1]->Reset();
    SQMATRICES_GRID[2] = SqMat_All[12];
    SQMATRICES_GRID[2]->Reset();
    SQMATRICES_GRID[3] = SqMat_All[13];
    SQMATRICES_GRID[3]->Reset();
    aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);  
       
    Assemble2D(1, fesp,
             4, SQMATRICES_GRID,
             0, NULL,
             0, NULL, NULL,
             DiscreteFormGrid,
             GridBoundaryConditions,
             GridBoundValues,
             aux);
    delete aux;   
   
     Entries[0] = SqMat_All[10]->GetEntries();
     Entries[1] = SqMat_All[11]->GetEntries();
     Entries[2] = SqMat_All[12]->GetEntries();
     Entries[3] = SqMat_All[13]->GetEntries();

   // for Dirichlet rows in off-diagonal matrices
   memset(Entries[1] + GridRowPtr[N_GActive], 0, (GridRowPtr[N_G] - GridRowPtr[N_GActive])*SizeOfDouble);
   memset(Entries[2] + GridRowPtr[N_GActive], 0, (GridRowPtr[N_G] - GridRowPtr[N_GActive])*SizeOfDouble);  
   
  //======================================================================  
  // end Assembeling the grid matrix
  //  Reparametrization of free surface - Begin
  //======================================================================  

  
//                 os.seekp(std::ios::beg);
//         if(img<10) os << "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
//          else if(img<100) os <<"VTK/"<< VtkBaseName<<".000"<<img<<".vtk" << ends;
//           else if(img<1000) os <<"VTK/"<< VtkBaseName<<".00"<<img<<".vtk" << ends;
//            else if(img<10000) os << "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
//             else  os << "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
//         Output->WriteVtk(os.str().c_str());
// 	img++;
    //======================================================================
    // nonlinear loop due to remeshing
    //first update the convection matrix with the new mesh velocity        
    //======================================================================
   if(remeshed)
    {           
     // working array for rhs is B, initialize B
     memset(B, 0, N_Unknowns*SizeOfDouble);

     if(remeshed)
       memset(Sol_All[1], 0, 2*N_G*SizeOfDouble);  
     
        DiscreteForm = DiscreteFormGalerkin;

        SQMATRICES[0] = SqMat_All[4];
        SQMATRICES[1] = SqMat_All[5];
        SQMATRICES[2] = SqMat_All[6];
        SQMATRICES[3] = SqMat_All[7];
        SQMATRICES[4] = SqMat_All[0];
        SQMATRICES[5] = SqMat_All[3];

        SQMATRICES[6] = SqMat_All[1];
        SQMATRICES[7] = SqMat_All[2];

        MATRICES[0] = Mat_All[0];
        MATRICES[1] = Mat_All[1];
        MATRICES[2] = Mat_All[2];
        MATRICES[3] = Mat_All[3];

        SQMATRICES[0]->Reset();
        SQMATRICES[1]->Reset();
        SQMATRICES[2]->Reset();
        SQMATRICES[3]->Reset();
        SQMATRICES[4]->Reset();
        SQMATRICES[5]->Reset();
        SQMATRICES[6]->Reset();
        SQMATRICES[7]->Reset();

        MATRICES[0]->Reset();
        MATRICES[1]->Reset();
        MATRICES[2]->Reset();
        MATRICES[3]->Reset();

        N_SquareMatrices = 6;
        N_RectMatrices = 4;

       // parameters which are the same for all NSTYPEs
        N_Rhs = 2;
        N_FESpaces = 3;

        fesp[0] = FESpaces_All[0];
        fesp[1] = FESpaces_All[1];
        fesp[2] = FESpaces_All[2];

        fefct[0] = FEFunctions_All[0];
        fefct[1] = FEFunctions_All[1];
        fefct[2] = FEFunctions_All[3];
        fefct[3] = FEFunctions_All[4];
 
        ferhs[0] = FESpaces_All[0];
        ferhs[1] = FESpaces_All[0];

        RHSs[0] = Rhs_All[0];
        RHSs[1] = Rhs_All[0] + N_U;
        RHSs[2] = Rhs_All[0] + 2*N_U;

        memset(Rhs_All[0], 0, N_Unknowns*SizeOfDouble);

       // 4 parameters are needed for assembling (u1_old, u2_old)
        aux =  new TAuxParam2D(MovingTNSN_FESpaces_Axial3D, MovingTNSN_Fct_Axial3D,
                               MovingTNSN_ParamFct_Axial3D,
                               MovingTNSN_FEValues_Axial3D,
                               fesp, fefct,
                               MovingTNSFct_Axial3D,
                               MovingTNSFEFctIndex_Axial3D,
                               MovingTNSFEMultiIndex_Axial3D,
                               MovingTNSN_Params_Axial3D, MovingTNSBeginParam_Axial3D);

       //======================================================================
       // assembling of matrices for each level
       // A_11 , (A_12), (A_21), (A_22), M_11, (M_22)
       //======================================================================
       Assemble2D(N_FESpaces, fesp,
                  N_SquareMatrices, SQMATRICES,
                  N_RectMatrices, MATRICES,
                  N_Rhs, RHSs, ferhs,
                  DiscreteForm,
                  BoundaryConditions,
                  BoundValues,
                  aux);

      delete aux;

     SqMat_All[8]->Reset(); // Matrix entries for freesurf int;
     SqMat_All[9]->Reset(); // no need to calculate in nonlinear steps			       

#ifdef __SURFACT__  
     // interface surfactant
     MapSurfToDomain(IFaceFeFunct[0], FEFunctions_All[7], N_List[0], N_List[1]);  
     
     FreeSurf_axial3D_new(SqMat_All[8], SqMat_All[9],  RHSs[0], RHSs[1], BoundCondition, tau,
                          FEFunctions_All[0]->GetValues(), FEFunctions_All[7], Params);     
#else     
     FreeSurf_axial3D_new(SqMat_All[8], SqMat_All[9],  RHSs[0], RHSs[1], BoundCondition, tau,
                          FEFunctions_All[0]->GetValues(), NULL, Params);
#endif  
    
     // Adding freesurf entries to A11 and A22
     MatAdd(SqMat_All[4], SqMat_All[8], 1);
     MatAdd(SqMat_All[7], SqMat_All[9], 1);

     // set rows of Dirichlet dof in off diagonal matrix blocks
     // to zero    
     // N_Active =  FESpaces_All[0]->GetActiveBound();
     // get row in off diagonal matrix where the Dirichlet nodes start
     RowPtr = SqMat_All[5]->GetRowPtr();
     // compute number of entries starting from this row to the end
     // of the matrix
     j = RowPtr[N_Active];
     k = RowPtr[N_U]-j;
     // get number of active dof
     // set these entries to zero
     memset(SqMat_All[5]->GetEntries()+j, 0, SizeOfDouble*k);
     memset(SqMat_All[6]->GetEntries()+j, 0, SizeOfDouble*k);

       if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
        {
          if (TDatabase::ParamDB->NSTYPE <4)
          {
            OutPut("For slip with friction bc NSTYPE 4 is ");
            OutPut("necessary !!!!! " << endl);
            exit(4711);
          }

          // prepare everything for the assembling of slip with friction bc
          // on all levels
          N_FESpaces = 1;
          N_SquareMatrices = 8;
          N_RectMatrices = 2;
          N_Rhs = 2;
          DiscreteForm = NULL;

          SQMATRICES[0] = SqMat_All[4];
          SQMATRICES[1] = SqMat_All[7];
          SQMATRICES[2] = SqMat_All[5];
          SQMATRICES[3] = SqMat_All[6];
          SQMATRICES[4] = SqMat_All[0];
          SQMATRICES[5] = SqMat_All[3];
          SQMATRICES[6] = SqMat_All[1];
          SQMATRICES[7] = SqMat_All[2];

          MATRICES[0] = Mat_All[2];
          MATRICES[1] = Mat_All[3];

          fesp[0] = FESpaces_All[0];
          ferhs[0] = FESpaces_All[0];
          ferhs[1] = FESpaces_All[0];

          RHSs[0] = Rhs_All[0];
          RHSs[1] = Rhs_All[0]+N_U;

          aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);

          Assemble2DSlipBC(N_FESpaces, fesp,
                           N_SquareMatrices, SQMATRICES,
                           N_RectMatrices, MATRICES,
                           N_Rhs, RHSs, ferhs,
                           DiscreteForm,
                           BoundaryConditions,
                           BoundValues,
                           aux, FEFunctions_All[0], FEFunctions_All[1]);

       delete aux;

      } // if (TDatabase::ParamDB->INTERNA

     //    scale the pressure matrices
     Dscal(Mat_All[2]->GetN_Entries(), tau, Mat_All[2]->GetEntries());
     Dscal(Mat_All[3]->GetN_Entries(), tau, Mat_All[3]->GetEntries());
     Dscal(Mat_All[0]->GetN_Entries(), tau, Mat_All[0]->GetEntries());
     Dscal(Mat_All[1]->GetN_Entries(), tau, Mat_All[1]->GetEntries());
     
     // update rhs
     Daxpy(N_Active, tau, Rhs_All[0], B);
     Daxpy(N_Active, tau, Rhs_All[0]+N_U, B+N_U);

     // update rhs by Laplacian and convective term initialy by current time step
     // scaled by current sub time step length and theta2
     // currently : M := M + gamma A
     // M = M + (-gamma - tau*TDatabase::TimeDB->THETA2) A
     MatAdd(SqMat_All[0], SqMat_All[4], - tau*TDatabase::TimeDB->THETA2);
     MatAdd(SqMat_All[1], SqMat_All[5], - tau*TDatabase::TimeDB->THETA2);
     MatAdd(SqMat_All[2], SqMat_All[6], - tau*TDatabase::TimeDB->THETA2);
     MatAdd(SqMat_All[3], SqMat_All[7], - tau*TDatabase::TimeDB->THETA2);

     // set current factor of steady state matrix
     gamma = -tau*TDatabase::TimeDB->THETA2;		     

     // defect = M * Sol
     // B:= B + defect (rhs)     
     MatVectActive(SqMat_All[0], Sol_All[0], defect);
     Daxpy(N_Active, 1, defect, B);
     MatVectActive(SqMat_All[1], Sol_All[0]+N_U, defect);
     Daxpy(N_Active, 1, defect, B);
     MatVectActive(SqMat_All[2], Sol_All[0], defect+N_U);
     Daxpy(N_Active, 1, defect+N_U, B+N_U);
     MatVectActive(SqMat_All[3], Sol_All[0]+N_U, defect+N_U);
     Daxpy(N_Active, 1, defect+N_U, B+N_U);     
     
    // set Dirichlet values
    // RHSs[0] still available from assembling
    memcpy(B+N_Active, RHSs[0]+N_Active, (N_U-N_Active)*SizeOfDouble);
    memcpy(B+N_Active+N_U, RHSs[1]+N_Active,(N_U-N_Active)*SizeOfDouble);

    // copy Dirichlet values from rhs into Sol[0][mg_level-1]
    memcpy(Sol_All[0]+N_Active, RHSs[0]+N_Active, (N_U-N_Active)*SizeOfDouble);
    memcpy(Sol_All[0]+N_Active+N_U, RHSs[1]+N_Active, (N_U-N_Active)*SizeOfDouble);

    //=====================================================================
    // the stiffness matrix is stored on M11, (M12, M21, M22)
    // assembling of system matrix
    //========================================================================
    // M = M + (-gamma + tau*TDatabase::TimeDB->THETA1) A
     MatAdd(SqMat_All[0], SqMat_All[4], -gamma + tau*TDatabase::TimeDB->THETA1);
     MatAdd(SqMat_All[1], SqMat_All[5], -gamma + tau*TDatabase::TimeDB->THETA1);
     MatAdd(SqMat_All[2], SqMat_All[6], -gamma + tau*TDatabase::TimeDB->THETA1);
     MatAdd(SqMat_All[3], SqMat_All[7], -gamma + tau*TDatabase::TimeDB->THETA1);

     // set current factor of steady state matrix
     gamma = tau*TDatabase::TimeDB->THETA1;     
  
     //======================================================================
     // nonlinear loop
     //======================================================================   
     N_LinIterCurr = 0;
     solver_time_curr = 0;      

     for(j=0;j<Max_It;j++)
      {
       memset(defect, 0, N_Unknowns*SizeOfDouble);

       SQMATRICES[0] = SqMat_All[0];
       SQMATRICES[1] = SqMat_All[1];
       SQMATRICES[2] = SqMat_All[2];
       SQMATRICES[3] = SqMat_All[3];
       MATRICES[0] = Mat_All[0];
       MATRICES[1] = Mat_All[1];
       MATRICES[2] = Mat_All[2];
       MATRICES[3] = Mat_All[3];      
       
      // compute defect
      Defect(sqmatrices, matrices, Sol_All[0], B, defect);

      residual =  Ddot(N_Unknowns, defect, defect);
      impuls_residual = Ddot(2*N_U, defect, defect);
      OutPut("nonlinear step " << setw(3) << j);
      OutPut(setw(14) << impuls_residual);
      OutPut(setw(14) << Ddot(N_P,defect+2*N_U,defect+2*N_U));
      OutPut(setw(14) << sqrt(residual));
      
      if(j>0)
       {
        OutPut(setw(14) << sqrt(residual)/oldresidual << endl);
       }
      else
       {
        OutPut(endl);
       }
       
      oldresidual = sqrt(residual);

      if ((((sqrt(residual)<=limit)||(j==Max_It-1)))  && (j>=TDatabase::ParamDB->SC_MINIT))
       {
        if (j==Max_It-1)
        j++;
        OutPut("ITE : " << setw(3) << j);
        OutPut(" (" << setw(3) << N_LinIterCurr << "/");
        OutPut(setw(3) << N_LinIter << " LINITE)");
        OutPut("  TIME FOR SOLVER : " << solver_time_curr << "/" << solver_time << "s");
        OutPut("  RES : " <<  sqrt(residual) << endl);
        // count total running time
        t4 =  GetTime();
        total_time += t4 - t3;
        t3 = t4;
        OutPut(TDatabase::TimeDB->CURRENTTIME << " total elapsed running time "<< total_time << endl);
        break;
       }

       //======================================================================
       // solve linear system
       //======================================================================
        t1 = GetTime();
        DirectSolver(SQMATRICES[0], SQMATRICES[1], SQMATRICES[2], SQMATRICES[3],
                     MATRICES[2], MATRICES[3], MATRICES[0], MATRICES[1],
                     B, Sol_All[0]);
        t2 = GetTime();
        solver_time_curr = t2-t1;
        solver_time += solver_time_curr;
 
       //======================================================================
       // end solve linear system
       //======================================================================
       // restore mass matrices by subtracting the A-matrices
       MatAdd(SqMat_All[0], SqMat_All[4], -gamma);
       MatAdd(SqMat_All[3], SqMat_All[7], -gamma);

       //======================================================================
       // assemble new matrix due to nonlinearity
       //======================================================================
       DiscreteForm = DiscreteFormNLGalerkin;	 
       N_RectMatrices = 0;
       N_Rhs = 0;
       N_FESpaces = 3;

       SQMATRICES[0] = SqMat_All[4];
       SQMATRICES[1] = SqMat_All[7];
       SQMATRICES[0]->Reset();
       SQMATRICES[1]->Reset();

       N_SquareMatrices = 2;
       last_sq = 1;
       
       fesp[0] = FESpaces_All[0];
       fesp[1] = FESpaces_All[1];
       fesp[2] = FESpaces_All[2];

       fefct[0] = FEFunctions_All[0];
       fefct[1] = FEFunctions_All[1];
       fefct[2] = FEFunctions_All[3];
       fefct[3] = FEFunctions_All[4];
 
       //======================================================================
       // assembling of matrices for each level due to nonlinearity
       // A_11, (A_22), no assembling of rhs
       //======================================================================
        aux =  new TAuxParam2D(MovingTNSN_FESpaces_Axial3D, MovingTNSN_Fct_Axial3D,
                               MovingTNSN_ParamFct_Axial3D,
                               MovingTNSN_FEValues_Axial3D,
                               fesp, fefct,
                               MovingTNSFct_Axial3D,
                               MovingTNSFEFctIndex_Axial3D,
                               MovingTNSFEMultiIndex_Axial3D,
                               MovingTNSN_Params_Axial3D, MovingTNSBeginParam_Axial3D);

         Assemble2D(N_FESpaces, fesp,
                    N_SquareMatrices, SQMATRICES,
                    N_RectMatrices, MATRICES,
                    N_Rhs, RHSs, ferhs,
                    DiscreteForm,
                    BoundaryConditions,
                    BoundValues,
                    aux); 
    
       // Adding freesurf entries to A11 and A22
       MatAdd(SqMat_All[4], SqMat_All[8], 1.);
       MatAdd(SqMat_All[7], SqMat_All[9], 1.); 
       

       if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
        {
          // prepare everything for the assembling of slip with friction bc
          // on all levels
          N_FESpaces = 1;
          N_SquareMatrices = 2;
          N_RectMatrices = 0;
          N_Rhs = 2;
          DiscreteForm = NULL;

          SQMATRICES[0] = SqMat_All[4];
          SQMATRICES[1] = SqMat_All[7];

          fesp[0] = FESpaces_All[0];
          ferhs[0] = FESpaces_All[0];
          ferhs[1] = FESpaces_All[0];

          RHSs[0] = Rhs_All[0];
          RHSs[1] = Rhs_All[0]+N_U;

          aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);

          Assemble2DSlipBC(N_FESpaces, fesp,
                           N_SquareMatrices, SQMATRICES,
                           N_RectMatrices, MATRICES,
                           N_Rhs, RHSs, ferhs,
                           DiscreteForm,
                           BoundaryConditions,
                           BoundValues,
                           aux, FEFunctions_All[0], FEFunctions_All[1]);

       delete aux;

      } // if (TDatabase::ParamDB->INTERN    
       
     MatAdd(SqMat_All[0], SqMat_All[4], tau*TDatabase::TimeDB->THETA1);
     MatAdd(SqMat_All[3], SqMat_All[7], tau*TDatabase::TimeDB->THETA1);       
    } //   for(j=0;j<Max_It;j++)   
            
    remeshed = FALSE;
         
   } //if (remeshed)  
  }// for(l=0;l<N_SubSteps;l++) 

//                 os.seekp(std::ios::beg);
//         if(img<10) os << "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
//          else if(img<100) os <<"VTK/"<< VtkBaseName<<".000"<<img<<".vtk" << ends;
//           else if(img<1000) os <<"VTK/"<< VtkBaseName<<".00"<<img<<".vtk" << ends;
//            else if(img<10000) os << "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
//             else  os << "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
//         Output->WriteVtk(os.str().c_str());
// 	img++;
// 	
// 	exit(0);


#ifndef __SURFACTDIFFTEST__
   MovBoundVert[0][0]->GetCoords(Lx, Ly);
   MovBoundVert[2][0]->GetCoords(Rx, Ry);
   
   if(MaxWetD<Rx-Lx)
   {
    MaxWetD = Rx-Lx;
    T_MaxWetD = TDatabase::TimeDB->CURRENTTIME;
   }
   

   //find the apex height of the droplet
   MovBoundVert[1][0]->GetCoords(x1, ApexHeight);
   for(k=0;k<N_MovVert[2];k++) // no need to set end vertices again
    {
     MovBoundVert[2][k]->GetCoords(x1, y1);
      if(y1>ApexHeight)
      { ApexHeight=y1;}
     }   
   
   if(MaxApexHeight<ApexHeight)
    {MaxApexHeight=ApexHeight;}
   
 if((m % 10 ) == 0   || m==1 )
  {
   MovBoundVert[0][N_MovVert[0]-1]->GetCoords(x1, y1);
   MovBoundVert[2][1]->GetCoords(x2, y2);
   MovBoundVert[2][2]->GetCoords(x3, y3);   
   MovBoundVert[2][3]->GetCoords(x4, y4);
   
   tx = x1-Rx;
   sx = x2-Rx;
   ty = y1-Ry;
   sy = y2-Ry;
   R_Theta[0] = acos( (tx*sx+ty*sy)/(sqrt(tx*tx+ty*ty)* sqrt(sx*sx+sy*sy)) )*(180/3.141592654);

   sx = x3-Rx;
   sy = y3-Ry;
   R_Theta[1] = acos( (tx*sx+ty*sy)/(sqrt(tx*tx+ty*ty)* sqrt(sx*sx+sy*sy)) )*(180./3.141592654);

   sx = ((x4))-Rx;
   sy = ((y4))-Ry;
   R_Theta[2] = acos( (tx*sx+ty*sy)/(sqrt(tx*tx+ty*ty)* sqrt(sx*sx+sy*sy)) )*(180./3.141592654); 
   
   MovBoundVert[1][0]->GetCoords(x1, y1);
   
   Get_KE(FEVectFuncts_All[0], Params);
   
   OutPut(setw(25)<<"T_MaxWetD, MaxWetD: " << T_MaxWetD << "   "<< MaxWetD << endl);
   
   if(!remeshed)
    OutPut(setw(25)<<"T, wd, MaxApexHeight, ApexHeight AxialZ,Ucl,RAng 1,2,3,   Volume, Diff, Rel. Diff : " 
                   << TDatabase::TimeDB->CURRENTTIME<<"   "<< Rx-Lx << "  "<< MaxApexHeight<< " " << ApexHeight
                   <<"   "<< y1<<"   "<<Params[2]<<"   "<<R_Theta[0]<<"   "<<R_Theta[1]<<"   "<<R_Theta[2]<<"   "<< Params[0]
                  <<"   "<< Params[0] - InitVolume<<"   "<< (Params[0] - InitVolume)/InitVolume << endl);
// #endif 
     CurrVolume = Params[0];
     KE = Params[1]; 
//    increase the TIMESTEPLENGTH  
   
     if( (Rx-Lx) >1 &&  AdaptTimeStep)
      {
//        TDatabase::TimeDB->TIMESTEPLENGTH = TDatabase::TimeDB->TIMESTEPLENGTH/10.;
       AdaptTimeStep = FALSE;
      }
      
#ifdef __SURFACT__  
      // interface surfactant
      MapSurfToDomain(IFaceFeFunct[0], FEFunctions_All[7], N_List[0], N_List[1]);
     
      // surfactant output 
      GetSurfactMass(FEFunctions_All[7], IFaceFeFunct[0], N_List[0], N_List[1], Surf_Mass);
      radius = pow((3.*CurrVolume), (1./3.)); // 2Pi get acanelled to to volume calculation in Get_KE
//       Output("CurrVolume : " << CurrVolume << " radius : " << radius <<endl);
      sphericity =  2.*Pi*radius*radius/(Surf_Mass[1]);  // hemisphire, so 2. * only   
      
//        Get_FeFunction2DMass(FEFunctions_All[6], Params);  
      FEFunctions_All[6]->GetMassAndMean(Params);      
// #ifndef __SURFACTCONVTEST__       
// #ifndef __SURFACTDIFFTEST__
#ifdef __SOLUBLE__     
       OutPut( "Time, GammaMaxR, GammaMax " <<TDatabase::TimeDB->CURRENTTIME<<
               " " <<GammaXmaxVal[0]<< " "<<GammaXmaxVal[1]<< " "<<endl);  
#endif       
// #endif          
       OutPut( "Time, Surfactant_Mass, Da*InterfaceSurfactant_Mass, InterfaceSurfactant_Conc " <<TDatabase::TimeDB->CURRENTTIME<<
               " " <<Params[0]<< " "<< TDatabase::ParamDB->REACTOR_P14*Surf_Mass[0]<< " "<<Surf_Mass[0]/Surf_Mass[1]<< " "<<endl);     
      OutPut( "Time, Surfactant_Mass_dif, InterfaceSurfactant_Mass_diff " <<TDatabase::TimeDB->CURRENTTIME<<
              " " <<Params[0] - Initial_SurfactMass<< " "<<TDatabase::ParamDB->REACTOR_P14*Surf_Mass[0]-Initial_IFaceSurfactMass<< endl);
       OutPut( "Time, Total Mass, Total Mass diff, RelativeTotal Mass diff, sphericity, KE " <<TDatabase::TimeDB->CURRENTTIME<<
               " " <<Params[0]+TDatabase::ParamDB->REACTOR_P14*Surf_Mass[0]<< " " << ((Params[0]+TDatabase::ParamDB->REACTOR_P14*Surf_Mass[0]) - 
            (Initial_IFaceSurfactMass + Initial_SurfactMass)) << " " << ((Params[0]+TDatabase::ParamDB->REACTOR_P14*Surf_Mass[0]) - 
            (Initial_IFaceSurfactMass + Initial_SurfactMass))
            /(Initial_IFaceSurfactMass + Initial_SurfactMass) <<" " << sphericity << " " << KE<<endl);     
// exit(0); 
          
/*#ifdef __SURFACTCONVTEST__         
      l_2_l_2u =  2.*Pi - Surf_Mass[0]; // not really L2,
      if(fabs(T_inftyL2)<fabs(l_2_l_2u))
        {
         T_inftyL2=fabs(l_2_l_2u);
         T_inftyL2_time=TDatabase::TimeDB->CURRENTTIME;
        }    
        
        Terrors[0] += (l_2_l_2u*l_2_l_2u +olderror * olderror)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
        olderror = l_2_l_2u;
      
       OutPut(TDatabase::TimeDB->CURRENTTIME <<  " L2(0,T;L2) " << sqrt(Terrors[0]) << " ");
       OutPut(" Mass Error: " <<  l_2_l_2u << " Rel Mass Error: " <<  l_2_l_2u/(2.*Pi) << " ");
       OutPut(T_inftyL2_time <<  " L2error_Max " << T_inftyL2 << endl);   
#endif    */   
       
      if(((m % (10*TDatabase::TimeDB->STEPS_PER_IMAGE)) == 0) || m==1)       
       {       
        PrintSurfSurfactant(N_MovVert[2], MovBoundVert[2], FEFunctions_All[7], N_BData_F, Surf_F);   
        PrintSurfSurfactant(N_MovVert[0], MovBoundVert[0], FEFunctions_All[7], N_BData_S, Surf_S);   
       } 
#endif     
   
  }
#endif


// #ifdef __SURFACTDIFFTEST__     
//       FEFunctions_All[7]->Interpolate(ExactS); 
//       MapDomainToSurf(FEFunctions_All[7], IFaceFeFunct[0], N_List[0], N_List[1]);   
//       memset(Sol_All[5], 0, N_S*SizeOfDouble); 
// #endif
      
      
 if((m % (int)(TDatabase::TimeDB->STEPS_PER_IMAGE) ) == 0   || m==1 )
  {
#ifdef __SURFACT__      
     // interface surfactant
     MapSurfToDomain(IFaceFeFunct[0], FEFunctions_All[7], N_List[0], N_List[1]);
#endif
     if(TDatabase::ParamDB->WRITE_VTK)
      { 
        os.seekp(std::ios::beg);
        if(img<10) os << "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
         else if(img<100) os <<"VTK/"<< VtkBaseName<<".000"<<img<<".vtk" << ends;
          else if(img<1000) os <<"VTK/"<< VtkBaseName<<".00"<<img<<".vtk" << ends;
           else if(img<10000) os << "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
            else  os << "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
        Output->WriteVtk(os.str().c_str());

//        os.seekp(std::ios::beg);
//        if(img<10) os << "VTK/"<<VtkBaseName<<".0000"<<img<<".plt" << ends;
//        else if(img<100) os << "VTK/"<<VtkBaseName<<".000"<<img<<".plt" << ends;
//        else if(img<1000) os << "VTK/"<<VtkBaseName<<".00"<<img<<".plt" << ends;
//        else if(img<10000) os << "VTK/"<<VtkBaseName<<".0"<<img<<".plt" << ends;
//        else  os << "VTK/"<<VtkBaseName<<"."<<img<<".plt" << ends; 
//        Output->WriteBinaryPlt(os.str().c_str());
 
        img++;
      }      
   } // if((m % (int)(TData

//    exit(0);
#ifdef __SURFACTDIFFTEST__    
// #ifndef __SURFACTCONVTEST__    
      MapSurfToDomain(IFaceFeFunct[0], FEFunctions_All[7], N_List[0], N_List[1]);
 
      GetSurfErrors(FEFunctions_All[7], IFaceFeFunct[0], errors, N_List[0], N_List[1]);


      Terrors[0] += (errors[0]*errors[0] + olderror * olderror)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
      olderror = errors[0];
//       Terrors[1] += (errors[1]*errors[1]);

      if(fabs(T_inftyL2)<fabs(errors[0]))
        {
         T_inftyL2=fabs(errors[0]);
         T_inftyL2_time=TDatabase::TimeDB->CURRENTTIME;
        }    

//       if(T_inftyH1<errors[1]) T_inftyH1=errors[1];      
      
      OutPut(TDatabase::TimeDB->CURRENTTIME <<  " L2(0,T;L2) " << sqrt(Terrors[0]) << " ");
      OutPut(" L2: " <<  errors[0] << " ");
      OutPut(T_inftyL2_time <<  " L2_Max " << T_inftyL2 << endl);   
      
//       OutPut( "L2: " << errors[0] << endl);
//       OutPut( "H1-semi: " << errors[1] << endl);
#endif      
// #else
      
   if( m % 200== 0 )
    {
     OutPut(setw(25)<<TDatabase::TimeDB->CURRENTTIME<<" No. ReParam : " << N_ReParam <<endl);  
     OutPut(setw(25)<<TDatabase::TimeDB->CURRENTTIME<<" No. Remeshed : " << N_Remesh <<endl);  
    }
  
// #endif   
   
// cout << " tau " << tau << endl;
// exit(0); 
   
//     if( ((m % 200)== 0) || (m==1) )
//      {
//       os.seekp(std::ios::beg);
//       if(N_BData<10) os << "BDData/Boundary.0000"<<N_BData<<".data" << ends;
//       else if(N_BData<100) os << "BDData/Boundary.000"<<N_BData<<".data" << ends;
//       else if(N_BData<1000) os << "BDData/Boundary.00"<<N_BData<<".data" << ends;
//       else if(N_BData<10000) os << "BDData/Boundary.0"<<N_BData<<".data" << ends;
//       else  os << "BDData/Boundary."<<N_BData<<".data" << ends;
// 
//       std::ofstream dat(os.str().c_str());
//       if (!dat)
//        {
//         cerr << "cannot open file for output" << endl;
//         return -1;
//        }
//       dat << "# Boundary data created for droplet by MooNMD" << endl;
//       dat << "# Current Time :" << TDatabase::TimeDB->CURRENTTIME << endl;
//       for(k=0;k<N_MovVert[0];k++) // no need to set end vertices again
//        {
//         MovBoundVert[0][k]->GetCoords(x1, y1);
// 
//         dat << x1 << " " <<  y1<< endl;
//        }
//       for(k=0;k<N_MovVert[2];k++) // no need to set end vertices again
//        {
//         MovBoundVert[2][k]->GetCoords(x1, y1);
//         dat << x1 << " " <<  y1<< endl;
//        }
//         MovBoundVert[1][0]->GetCoords(x1, y1);
//         dat << x1 << " " <<  y1<< endl;
//         
//         MovBoundVert[0][0]->GetCoords(x1, y1);
//         dat << x1 << " " <<  y1<< endl;
//         
//       dat.close();
//       cout << endl;
//       cout << "Boundary wrote output into file " << endl;
//       N_BData++;  
//       
//    } //if( ((m % 20)== 0) || (m==1) )
  
//     exit(0);
  } // while(TDatabase::TimeDB->CURRENTTIME< end_time)

  //======================================================================
  // end of time cycle
  //======================================================================   
     if(TDatabase::ParamDB->WRITE_VTK)
       { 
#ifdef __SURFACT__  	 
       // interface surfactant
       MapSurfToDomain(IFaceFeFunct[0], FEFunctions_All[7], N_List[0], N_List[1]);
#endif      
       os.seekp(std::ios::beg);
        if(img<10) os << "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
         else if(img<100) os << "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
          else if(img<1000) os << "VTK/"<<VtkBaseName<<".00"<<img<<".vtk" << ends;
           else if(img<10000) os << "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
            else  os << "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
        Output->WriteVtk(os.str().c_str());

//        os.seekp(std::ios::beg);
//        if(img<10) os << "VTK/"<<VtkBaseName<<".0000"<<img<<".plt" << ends;
//        else if(img<100) os << "VTK/"<<VtkBaseName<<".000"<<img<<".plt" << ends;
//        else if(img<1000) os << "VTK/"<<VtkBaseName<<".00"<<img<<".plt" << ends;
//        else if(img<10000) os << "VTK/"<<VtkBaseName<<".0"<<img<<".plt" << ends;
//        else  os << "VTK/"<<VtkBaseName<<"."<<img<<".plt" << ends; 
//        Output->WriteBinaryPlt(os.str().c_str());
// 	
        img++;
       }    
       
#ifdef __SURFACTDIFFTEST__    
#ifndef __SURFACTCONVTEST__         
      OutPut(TDatabase::TimeDB->CURRENTTIME <<  " L2(0,T;L2) " << sqrt(Terrors[0]) << " ");
      OutPut(" L2:: " <<  errors[0] << " ");
      OutPut(T_inftyL2_time <<  " L2_Max " << T_inftyL2 << endl);  
#endif
#endif      
      
    // count total running time
    OutPut(TDatabase::TimeDB->CURRENTTIME << " total elapsed running time "<< total_time << endl);

   CloseFiles();  
   return 0;  
}

