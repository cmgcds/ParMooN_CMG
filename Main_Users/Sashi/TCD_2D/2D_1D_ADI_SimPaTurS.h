#include <TimeConvDiff2D.h>
#include <MacroCell.h>

void ExampleFile()
{
  OutPut("Example: 2D_1D_ADI_SimPaTurS.h" << endl) ;

  TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD=1;
  #define __SIMPATURS__
}

void InitialU1(double x, double y, double *values)
{

  values[0] = 0.;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}


void InitialU2(double x, double y, double *values)
{
  values[0] = 0.;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}


// exact solution
void Exact(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}


void InitialCondition_Heat(double x, double y, double *values)
{
  double val;
  double T_D = TDatabase::ParamDB->REACTOR_P23;
  double T_W = TDatabase::ParamDB->REACTOR_P24;

  val = T_W/T_D;

  values[0] = val;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

// void Exact_Psd_Intl(double x, double y, double L, double *values)
void Exact_Psd_Intl(int N_Coord, double *X, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}
// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition_Heat(int i, double t, BoundCond &cond)
{

  switch(i)
  {
   case 0:
   case 2:
   case 4:
     cond = DIRICHLET;
   break;

   case 1:
   case 3: 
   case 5:
     cond = NEUMANN;
   break;

   default: 
      cout << "wrong boundary part number "  << i << endl;
      exit(0);
   break;
  }

}

void BoundValue_Heat(int BdComp, double Param, double &value)
{
  double val, dt, t,  y;
  double T_D = TDatabase::ParamDB->REACTOR_P23;
  double T_W = TDatabase::ParamDB->REACTOR_P24;

  val = T_W/T_D;

  switch(BdComp)
  {
   case 1:
   case 3:
   case 5:
     value=0.;
   break;

  case 0:
  case 2:
        value=val;
  break;
  

  case 4:
     y = 1 -  Param;
     dt = 1. - val;
//  cout << y << " t " << Param << endl;
     if( y>0 && y<(1./8.))
      {
	t = -6. + 12.*y*8.*(1. + 1./8. - y);
	t = val + dt*(1./(1. + exp(-t)));
//              cout << y << " t " << t << endl;
	value = t;  
      }
     else if(y<1. &&   y>(1. - 1./8.))
      {
	t = 6. + 12.*(1. - 1./8. - y)*8.*(1. + 1. -y);
	t = 1./(1. + exp(-t));
	t = val + dt*(1./(1. + exp(-t)));
//          cout << y << " t " << t << endl;
	value = t;      
      }
      else if(y==0. || y==1.)
      {
       value = val;  
      }   
      else
      {
	value = 1.;  
      }
     
  break;

  default:
     cout << "wrong boundary part number "  << BdComp << endl;
    break;
  }
}

// ========================================================================
// for scalar equations
// ========================================================================
void BilinearCoeffs_Heat(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  double eps;
  int i;
  double *coeff, *param;
  double x, y;

  if(TDatabase::ParamDB->REACTOR_P0)
    eps = 1.0/TDatabase::ParamDB->REACTOR_P0;
  else
    eps = 0.;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    x = X[i];
    y = Y[i];

    coeff[0] = eps;
    if(TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD)
     {
      coeff[1] = param[0];  // u1
      coeff[2] = param[1];  // u2
//       cout<< "coeff[0] eps " << eps << " u1 " << param[0]  << "u2 " << param[1] << endl;
     }
    else
     {
    coeff[1] = 0;  // u1
    coeff[2] = 0;  // u2
     }
    coeff[3] = 0.;

    coeff[4] = 0.; // f
  }
}



void InitialCondition_Conc(double x, double y, double *values)
{

  values[0] = 1.;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}


// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition_Conc(int i, double t, BoundCond &cond)
{
  switch(i)
  {
   case 0:
   case 1:
   case 2:
   case 3:
   case 5:
     cond = NEUMANN;
   break;

   case 4:
     cond = DIRICHLET;
   break;

   default: 
      cout << "wrong boundary part number "  << i << endl;
      exit(0);
   break;
  }
}

void BoundValue_Conc(int BdComp, double Param, double &value)
{
   double t,  y;
   
 switch(BdComp)
  {
  case 0:
  case 1:
  case 2:
  case 3:
  case 5:
     value=0.;
  break;

  case 4:
//      y = 1 -  Param;
// 
//      if( y>0 && y<(1./8.))
//       {
// 	t = -6. + 12.*y*8.*(1. + 1./8. - y);
// 	t = 1./(1. + exp(-t));
// //              cout << y << " t " << t << endl;
// 	value = t;  
//       }
//       else if(y<1. &&   y>(1. - 1./8.))
//       {
// 	t = 6. + 12.*(1. - 1./8. - y)*8.*(1. + 1. -y);
// 	t = 1./(1. + exp(-t));
// //          cout << y << " t " << t << endl;
// 	value = t;      
//       }
//       else if(y==0. || y==1.)
//       {
// 	value = 0.;  
//       }       
//       else
      {
	value = 1.;  
      }
   
      break;

      default: cout << "wrong boundary part number " << BdComp << endl;
    break;
  }
}


// ========================================================================
// for scalar equations
// ========================================================================
void BilinearCoeffs_Conc(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  double eps;
//   double a=0.2, b=0., c=0.;
  int i;
  double *coeff, *param;
  double x, y;

  if(TDatabase::ParamDB->REACTOR_P1)
    eps = 1.0/TDatabase::ParamDB->REACTOR_P1;
  else
    eps = 0.;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    x = X[i];
    y = Y[i];

    coeff[0] = eps;
    if(TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD)
     {
      coeff[1] = param[0];  // u1
      coeff[2] = param[1];  // u2
     }
    else
     {
      coeff[1] = 0;  // u1
      coeff[2] = 0;  // u2
     }
    coeff[3] = 0.;

    coeff[4] = 0.; // f
  }
}

void InitialCondition_Psd(double x, double y, double *values)
{
  values[0] = 0.;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition_Psd(int i, double t, BoundCond &cond)
{
  switch(i)
  {
   case 0:
   case 1:
   case 2:
   case 3:
   case 5:
     cond = NEUMANN;
   break;

   case 4:
     cond = DIRICHLET;
   break;

   default: 
      cout << "wrong boundary part number " << i << endl;
      exit(0);
   break;
  }
}

// void BoundValue_Psd(int BdComp, double Param, double &value)
// {
//   
//  double t, L = TDatabase::ParamDB->REACTOR_P29;
//   
//   
// //   if(L < 0.025) // condition at L_Min
//   if(L == 0 ) // condition at L_Min   
//    {
//     switch(BdComp)
//      {
//       case 0:
//       case 1:
//       case 2:
//       case 3:
//       case 5:
// 	value = 0.;
//       break;
//       case 4:
// 	value = 1.;
// //        value = exp(-L/0.0025);
//        
//       break;
//       
//       default: cout << "wrong boundary part number "  << BdComp  << endl;
// 	break;
//       }     
//    }
//   else if(L == TDatabase::ParamDB->REACTOR_P13) // condition at L_Max
//    {
//     switch(BdComp)
//      {
//       case 0:
//       case 1:
//       case 2:
//       case 3:
//       case 5:
// 	value = 0.;
//       break;
//       case 4:
// 	value = 0.;
//       break;
//       
//       default: cout << "wrong boundary part number "  << BdComp  << endl;
// 	break;
//       }     
//      }
//   else
//    {     
//     switch(BdComp)
//      {
//       case 0:
//       case 1:
//       case 2:
//       case 3:
//       case 5:
// 	value = 0.;
//       break;
//       case 4:
// 	value = 0.;
//       break;
//       
//       default: cout << "wrong boundary part number "  << BdComp  << endl;
// 	break;
//       }     
//      }
// 
//  }
//     

void BoundValue_Psd(int BdComp, double Param, double &value)
{  
 double t, L = TDatabase::ParamDB->REACTOR_P29;
 
    switch(BdComp)
     {
      case 0:
      case 1:
      case 2:
      case 3:
      case 5:
	value = 0.;
      break;
      case 4:
	if(L <0.025 && L>0.005)
	 {
	  t = (L - 0.005)/0.02;
	  value = 4.*t*(1.-t);
	 }
	else
	 value = 0.;
      break;
      
      default: cout << "wrong boundary part number "  << BdComp  << endl;
	break;
      }    

}

void BoundCondition_LminLMax(BoundCond &cond_Lmin, BoundCond &cond_Lmax)
{
  cond_Lmin = DIRICHLET;
//   cond_Lmin = NEUMANN;
  cond_Lmax = NEUMANN;
}


void BoundValue_LMin(double x, double y,  double *values)
 {
   double t;
  
   // no nucleation, so f = 0 at L=0
   
//  if((x==0) && y>(1./3) &&  y<(2./3))
//   {
//   if((x==0) && y>(1./3) &&  y<(1./3. + 1./24.))
//    {
//     t = -6. + 12.*(y - 1./3.)*24.*(1. + 1./3. + 1./24. - y);
//     t = 1./(1. + exp(-t));
// //              cout << y << " t " << t << endl;
//     values[0] = t;
//     values[1] = 0.;
//     values[2] = 0.;    
//    }
//   else if((x==0) && y>(2./3 - 1./24.)  &&  y<(2./3))
//    {
//      t = 6. + 12.*(2./3. - 1./24. - y)*24.*(1. + 2./3. -y);
//      t = 1./(1. + exp(-t));
// //          cout << y << " t " << t << endl;
//     values[0] = t;
//     values[1] = 0.;
//     values[2] = 0.;    
//    }  
//   else
//    {
//     values[0] = 1.;
//     values[1] = 0.;
//     values[2] = 0.;
//    }
//   }
//  else
  {
   values[0] = 0.;
   values[1] = 0.;
   values[2] = 0.;
  } 

 }


void BoundValue_LMax(double x, double y,  double *values)
 {
  values[0] = 0.;
  values[1] = 0.;
  values[2] = 0.;  
 }


// ========================================================================
// for scalar equations
// ========================================================================
void BilinearCoeffs_Psd(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  double eps;
  int i;
  double *coeff, *param;
  double x, y;

  if(TDatabase::ParamDB->REACTOR_P2)
    eps = 1.0/TDatabase::ParamDB->REACTOR_P2;
  else
    eps = 0.;


  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    x = X[i];
    y = Y[i];

    coeff[0] = eps;
    if(TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD)
     {
      coeff[1] = param[0];  // u1
      coeff[2] = param[1];  // u2
      //cout<< "coeff[0] eps " << eps << " u1 " << param[0]  << "u2 " << param[1] << endl;

     }
    else
     {
    coeff[1] = 0;  // u1
    coeff[2] = 0;  // u2
     }
    coeff[3] = 0.;

    coeff[4] = 0.; // f
  }
}


// // initial conditon
// void PsdAtLmin(double x, double y, double z, double &value)
// {
// //   //testing the implementation
// //   double t = TDatabase::TimeDB->CURRENTTIME;
// //   value = x*t/5.;
// 
// //   value = value:
// }


// initial conditon
// void InitialCondition_Psd_Intl(double x, double y, double z, double *values)
void InitialCondition_Psd_Intl(int N_Coord, double *X, double *values)
{

  //testing the implementation
//   double t = TDatabase::TimeDB->CURRENTTIME;
//    if(z < 0.025) 
//     values[0] = exp(-z/0.0025);  
//    else
     values[0] = 0;

//   values[0] = 0;
//    cout<< " test example "<<  values[0] <<endl;
}

// void BilinearCoeffs_Psd_Intl(int n_points, double *X, double *Y, double *Z,
//         double **parameters, double **coeffs)
void BilinearCoeffs_Psd_Intl(int n_points, int N_Dim, double **Coords, double **parameters, double **coeffs)
{
  int i;
  double eps, *coeff;                                  // *param;
  double x, y, z, c, a[3], b[3], s[3], h;
  double t = TDatabase::TimeDB->CURRENTTIME;

  double *X, *Y, *Z;
  
   X = coeffs[0];
   Y = coeffs[1];
   Z = coeffs[2];   
  
  b[0] = -1e8;// negative, so that C will be taken from the PBS growth term
  c = 0; 

  if(TDatabase::ParamDB->REACTOR_P3)
    eps = 1.0/TDatabase::ParamDB->REACTOR_P3;
  else
    eps = 0.;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    // param = parameters[i];

    x = X[i];
    y = Y[i];
    z = Z[i];

    // diffusion
    coeff[0] = eps;
    // convection in z direction
    coeff[1] = b[0];
    // reaction term
    coeff[2] = c;
    coeff[3] = 0.;
  }
}

void GetExampleFileData(BoundCondFunct2D **BoundaryConditions, BoundValueFunct2D **BoundValues, 
                        DoubleFunct2D **InitiaValues, CoeffFct2D **BilinearCoeffs, 
                        int &N_PBEqns, int &N_IndepntScalarEqns, int *Disctypes)
{

//   TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD=1;


   N_IndepntScalarEqns = 2;
   N_PBEqns = 1;
   #define __PBS__

   BilinearCoeffs[0] = BilinearCoeffs_Heat;
   BilinearCoeffs[1] = BilinearCoeffs_Conc;
   BilinearCoeffs[2] = BilinearCoeffs_Psd;

   BoundaryConditions[0] = BoundCondition_Heat;
   BoundaryConditions[1] = BoundCondition_Conc;
   BoundaryConditions[2] = BoundCondition_Psd;

   BoundValues[0] = BoundValue_Heat;
   BoundValues[1] = BoundValue_Conc;
   BoundValues[2] = BoundValue_Psd;

   InitiaValues[0] = InitialCondition_Heat;
   InitiaValues[1] = InitialCondition_Conc;
   InitiaValues[2] = InitialCondition_Psd;

   Disctypes[0] = GALERKIN;
   Disctypes[1] = GALERKIN;
//    Disctypes[2] = SDFEM;
   Disctypes[2] = GALERKIN;
}


void Generate1DMesh(TDomain *Domain, double Start, double End, int N_Cells)
{
  int i, j, N_V;
  int *Lines;
  double len, h, x, y, *X;
  TVertex **Vetrex;
  TJoint *Joint;
  TBaseCell  **CellTree;

  N_V = N_Cells+1;
  X = new double[N_V];

  if(Start!=0. || End!=1.)
   {
    cout << "Domain should be (0,1) for this example " << endl;
    exit(0);
   }

  for(i=0; i<N_V; i++)
   X[i] = 1. + (1. - Start)*(tanh(2.75*(double(i)/double(N_Cells) - 1.)))/tanh(2.75);

//   h = (End - Start)/N_Cells;
//   for(i=1; i<N_V; i++)
//    X[i] =  h*(double)i;

  X[0] = Start;
  X[N_V-1] = End;

//   for(i=0; i<N_V; i++)
//    cout<< " X[i] " << X[i] <<endl;

  Lines = new int[2*N_Cells];
  Vetrex = new TVertex*[N_V]; 

  y=0.;

  for(i=0; i<N_Cells; i++)
   {
    Lines[2*i]=i;
    Lines[2*i+1]=i+1;
    Vetrex[i] = new TVertex(X[i], y);
   }

  Vetrex[N_Cells] = new TVertex(X[N_V-1], y);

  CellTree = new TBaseCell*[N_Cells];

   for (i=0;i<N_Cells;i++)
   {
//     Vetrex[ i ]->GetCoords(x, y);
//     cout<< " x " << x<< " y " << y<<endl;
    CellTree[i] = new TMacroCell(TDatabase::RefDescDB[S_Line], 0);
    CellTree[i]->SetVertex(0, Vetrex[ Lines[ 2*i       ]]);
    CellTree[i]->SetVertex(1, Vetrex[ Lines[ 2*i + 1]]);
    ((TMacroCell *) CellTree[i])->SetSubGridID(0);

//  cout<< " x " <<CellTree[i]->GetN_Edges()<<endl;;
//  cout<< " x " <<TDatabase::RefDescDB[S_Line]->GetN_OrigEdges()<<endl;;

   }


//     Vetrex[ i ]->GetCoords(x, y);
//     cout<< " x " << x<< " y " << y<<endl;
//     exit(0);

   Domain->SetTreeInfo(CellTree, N_Cells);

   TDatabase::IteratorDB[It_EQ]->SetParam(Domain);
   TDatabase::IteratorDB[It_LE]->SetParam(Domain);
   TDatabase::IteratorDB[It_Finest]->SetParam(Domain);
   TDatabase::IteratorDB[It_Between]->SetParam(Domain);
   TDatabase::IteratorDB[It_OCAF]->SetParam(Domain);

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


