#include <TCD3D_Surface.h>


void TCD3D_Surf::MatricesAssemble (double Mult, double *coeff, double *param, double hK,
			 double **OrigValues, int *N_BaseFuncts,
			 double ***LocMatrices, double **LocRhs)
{
  int N_;
  double val, c0;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  double test100, test010, test001, ansatz100, ansatz010, ansatz001;
  double test000, ansatz000, u0, u1, u2, u0x, u1y, u2z;
  double **MatrixA, **MatrixM, *MatrixARow, *MatrixMRow;
  
  N_ = N_BaseFuncts[0];
  
  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];
  Orig3 = OrigValues[3];
  
  c0 = coeff[0];  // eps
  
  u0x = param[0]; // u1_x
  u1y = param[1]; // u2_y
  u2z = param[2]; // u3_z
  u0  = param[3]; // u1
  u1  = param[4]; // u2
  u2  = param[5]; // u3
  
  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  
  for (int i=0;i<N_;++i)
  {
    MatrixARow = MatrixA[i];
    MatrixMRow = MatrixM[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];
    
    for (int j=0;j<N_;++j)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val = c0 * (test100*ansatz100 + test010*ansatz010 + test001*ansatz001); // laplace
      //val +=  (u0*ansatz100 + u1*ansatz010 + u2*ansatz001)*test000	      // U.grad c
      //       /*+(u0x+u1y+u2z)*ansatz000*test000*/ ;                         // c divu;
      val *= Mult;
      
      MatrixARow[j] += val;
      
      val = test000*ansatz000;
      val *= Mult;
      
      MatrixMRow[j] += val;
    }
  }  
}
