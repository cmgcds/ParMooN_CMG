/** ************************************************************************ 
* @brief     source file for FE3D_ALE
             
* @brief  Description :

            Contains all the functions that are needed to calculate the Normals of the
            Surface in the GPU ( Device ).

            Functions and formulas are same as implemented in the sequential Implementation of ParMooN.
            They are just simplified ( in Datastructure wise ) for implementing in GPU

* @author  Thivin Anandh D
* @date      4-Sep-2019
* @History   Jan -8   - Primary Implementation
 ************************************************************************  */
#include <cuda.h>
#include "helper_cuda.h"
#include <FE3D_ALE.h>

// void (*hxTri_normal_func1[6]) (double xi , double eta, double zeta,double& n1,double& n2, double& n3 ,
//                             int ThreadId_new, double* X_cord, double* Y_cord, double* Z_cord, int* m_d_FreeSurfCellMapping,
//                             int N_VertexPerCell);




__global__ void C_calculate_normals_freeSurf(int N_FreeSurfvertex, int N_VertexPerCell,
                                             int* m_d_FreeSurfCellMapping,double* X_cord, double* Y_cord, double* Z_cord,
                                             double* m_d_nodalFunctionalRefValues,
                                             int* m_d_FreeSurfaceCells,int* m_d_freeSurfaceVertexLocal,int* m_d_freeSurfaceJoints,
                                             int* m_d_FreeSurfaceVertex,double* m_d_freeSurfaceNormal_1, double* m_d_freeSurfaceNormal_2, double* m_d_freeSurfaceNormal_3
                                             )
{
   double x[8];
   double y[8];
   double z[8];
   double a,b,c;
   double xc[8];
   double yc[8];
   double zc[8];


   // printf("%d\n", threadIdx.x);
   // Each Thread Will Copy its set of XYZ Co-Ordinates to the GPU
   int ThreadId = blockIdx.x*blockDim.x + threadIdx.x; 

   // printf("%d \t", ThreadId);

 
   if(ThreadId < N_FreeSurfvertex)
   {       
      int cellMapIndex    =  m_d_FreeSurfCellMapping[ThreadId];
      int shiftPosition   =  cellMapIndex * N_VertexPerCell;
         
      for ( int i = 0 ; i < N_VertexPerCell ; i++)
      {
         x[i] =   X_cord[shiftPosition + i];
         y[i] =   Y_cord[shiftPosition + i];
         z[i] =   Z_cord[shiftPosition + i];
         
      }
      

      // if(ThreadId == 80 )
      // {
      //    printf("\n Y Values : %f    %f   %f   %lf  %lf  %lf  %f   %f \n\n",  y[0], y[1],y[2], y[3],y[4],y[5],y[6],y[7]);
      //    printf("\n X Values : %f    %f   %f   %lf  %lf  %lf  %f   %f \n\n",  x[0], x[1],x[2], x[3],x[4],x[5],x[6],x[7]);
      //    printf("\n Z Values : %f    %f   %f   %lf  %lf  %lf  %f   %f \n\n",  z[0], z[1],z[2], z[3],z[4],z[5],z[6],z[7]);

      // }
      // Setting up Parameters for FE Cell -- HEXAtrilinear Mapping // 
      
      xc[0] = (x[0] + x[1] + x[3] + x[2] + x[4] + x[5] + x[7] + x[6]) * 0.125;
      xc[1] = (-x[0] + x[1] - x[3] + x[2] - x[4] + x[5] - x[7] + x[6]) * 0.125;
      xc[2] = (-x[0] - x[1] + x[3] + x[2] - x[4] - x[5] + x[7] + x[6]) * 0.125;
      xc[3] = (-x[0] - x[1] - x[3] - x[2] + x[4] + x[5] + x[7] + x[6]) * 0.125;
      xc[4] = (x[0] - x[1] - x[3] + x[2] + x[4] - x[5] - x[7] + x[6]) * 0.125;
      xc[5] = (x[0] - x[1] + x[3] - x[2] - x[4] + x[5] - x[7] + x[6]) * 0.125;
      xc[6] = (x[0] + x[1] - x[3] - x[2] - x[4] - x[5] + x[7] + x[6]) * 0.125;
      xc[7] = (-x[0] + x[1] + x[3] - x[2] + x[4] - x[5] - x[7] + x[6]) * 0.125;

      yc[0] = (y[0] + y[1] + y[3] + y[2] + y[4] + y[5] + y[7] + y[6]) * 0.125;
      yc[1] = (-y[0] + y[1] - y[3] + y[2] - y[4] + y[5] - y[7] + y[6]) * 0.125;
      yc[2] = (-y[0] - y[1] + y[3] + y[2] - y[4] - y[5] + y[7] + y[6]) * 0.125;
      yc[3] = (-y[0] - y[1] - y[3] - y[2] + y[4] + y[5] + y[7] + y[6]) * 0.125;
      yc[4] = (y[0] - y[1] - y[3] + y[2] + y[4] - y[5] - y[7] + y[6]) * 0.125;
      yc[5] = (y[0] - y[1] + y[3] - y[2] - y[4] + y[5] - y[7] + y[6]) * 0.125;
      yc[6] = (y[0] + y[1] - y[3] - y[2] - y[4] - y[5] + y[7] + y[6]) * 0.125;
      yc[7] = (-y[0] + y[1] + y[3] - y[2] + y[4] - y[5] - y[7] + y[6]) * 0.125;

      zc[0] = (z[0] + z[1] + z[3] + z[2] + z[4] + z[5] + z[7] + z[6]) * 0.125;
      zc[1] = (-z[0] + z[1] - z[3] + z[2] - z[4] + z[5] - z[7] + z[6]) * 0.125;
      zc[2] = (-z[0] - z[1] + z[3] + z[2] - z[4] - z[5] + z[7] + z[6]) * 0.125;
      zc[3] = (-z[0] - z[1] - z[3] - z[2] + z[4] + z[5] + z[7] + z[6]) * 0.125;
      zc[4] = (z[0] - z[1] - z[3] + z[2] + z[4] - z[5] - z[7] + z[6]) * 0.125;
      zc[5] = (z[0] - z[1] + z[3] - z[2] - z[4] + z[5] - z[7] + z[6]) * 0.125;
      zc[6] = (z[0] + z[1] - z[3] - z[2] - z[4] - z[5] + z[7] + z[6]) * 0.125;
      zc[7] = (-z[0] + z[1] + z[3] - z[2] + z[4] - z[5] - z[7] + z[6]) * 0.125;

      // hello();

      // Get the local Reference co-ordinates of the Vertex (xi eta ,zeta)
      double xi, eta , zeta;
      int localNodeNumber = m_d_freeSurfaceVertexLocal[ThreadId];

      if(localNodeNumber > N_VertexPerCell || localNodeNumber < 0)
      {
         printf(" Error : LocalVertex Number not valid");
         printf(" INFO   :  Function : calculate_normals_freeSurf , Class : FE3D_ALE, File : FE3D_ALE.cu");
         //Should exit the program -- CODE later                // THIVIN - TODO
      }

      xi    = m_d_nodalFunctionalRefValues[localNodeNumber];
      eta   = m_d_nodalFunctionalRefValues[localNodeNumber + N_VertexPerCell];
      zeta  = m_d_nodalFunctionalRefValues[localNodeNumber + 2*N_VertexPerCell];

   
      // Obtain the xi eta and zeta values required for 

      // Calculate Surface Normals  
      // Each Cuda Thread Will calculate Normal based on the value from its shared Array
      // Get the Tangent Vectors , and use their Cross product to optain the normal vector for the surface

      double xi_1, xi_2;
      
      int Jointno = m_d_freeSurfaceJoints[ThreadId];
      switch (Jointno)
      {
         case 0:
            xi_1 = xi; xi_2 = eta; break;
         case 1:
            xi_1 = zeta; xi_2 = xi; break;
         case 2:
            xi_1 = zeta; xi_2 = eta; break;
         case 3:
            xi_1 = zeta; xi_2 = -xi; break;
         case 4:
            xi_1 = eta; xi_2 = zeta;break;
         case 5:
            xi_1 = eta; xi_2 = xi;break;
         default:
         {
            printf(" Error in Class : Hexatrilinear.C - Function : GetRefValuesfromJointid \n");
            printf(" The Given Joint Id for the HEXAHEADRAL cell %d Does not exist ",Jointno);
         }
      }
         // if(ThreadId == 80)
         // {  
         //    printf(" Cell Mapping Number :  %d", cellMapIndex);
         //    printf(" xi_1 and xi_ 2 (dev):   %lf      %lf\n", xi_1 , xi_2);
         //    printf(" y1 and y4 , yc5 , yc7 (dev):   %lf    %lf  %lf    %lf \n", yc[1] ,yc[4],yc[5],yc[7] );
         //    printf(" Joint:   %d    \n", Jointno);

         // }
      // Find the Tangent Vectors for the Node 
      double p1 = xi_1; double p2 = xi_2;
      double t11,t12,t13,t21,t22,t23;
      double Xi,Eta,Zeta;
      switch(Jointno)
      {
         case 0:
            Xi = p1; Eta = p2; Zeta = -1;
            t11 = xc[2] + xc[4]*Xi + xc[6]*Zeta + xc[7]*Xi*Zeta;
            t12 = yc[2] + yc[4]*Xi + yc[6]*Zeta + yc[7]*Xi*Zeta;
            t13 = zc[2] + zc[4]*Xi + zc[6]*Zeta + zc[7]*Xi*Zeta;
      
            t21 = xc[1] + xc[4]*Eta + xc[5]*Zeta + xc[7]*Eta*Zeta;
            t22 = yc[1] + yc[4]*Eta + yc[5]*Zeta + yc[7]*Eta*Zeta;
            t23 = zc[1] + zc[4]*Eta + zc[5]*Zeta + zc[7]*Eta*Zeta;
         break;
      
         case 1:
            Xi = p2; Eta = -1; Zeta = p1;
            t11 = xc[1] + xc[4]*Eta + xc[5]*Zeta + xc[7]*Eta*Zeta;
            t12 = yc[1] + yc[4]*Eta + yc[5]*Zeta + yc[7]*Eta*Zeta;
            t13 = zc[1] + zc[4]*Eta + zc[5]*Zeta + zc[7]*Eta*Zeta;
      
            t21 = xc[3] + xc[5]*Xi + xc[6]*Eta + xc[7]*Xi*Eta;
            t22 = yc[3] + yc[5]*Xi + yc[6]*Eta + yc[7]*Xi*Eta;
            t23 = zc[3] + zc[5]*Xi + zc[6]*Eta + zc[7]*Xi*Eta;
         break;
      
         case 2:
            Xi = 1; Eta = p2; Zeta = p1;
            t11 = xc[2] + xc[4]*Xi + xc[6]*Zeta + xc[7]*Xi*Zeta;
            t12 = yc[2] + yc[4]*Xi + yc[6]*Zeta + yc[7]*Xi*Zeta;
            t13 = zc[2] + zc[4]*Xi + zc[6]*Zeta + zc[7]*Xi*Zeta;
      
            t21 = xc[3] + xc[5]*Xi + xc[6]*Eta + xc[7]*Xi*Eta;
            t22 = yc[3] + yc[5]*Xi + yc[6]*Eta + yc[7]*Xi*Eta;
            t23 = zc[3] + zc[5]*Xi + zc[6]*Eta + zc[7]*Xi*Eta;
         break;
      
         case 3:
            Xi = -p2; Eta = 1; Zeta = p1;
            t11 = xc[3] + xc[5]*Xi + xc[6]*Eta + xc[7]*Xi*Eta;
            t12 = yc[3] + yc[5]*Xi + yc[6]*Eta + yc[7]*Xi*Eta;
            t13 = zc[3] + zc[5]*Xi + zc[6]*Eta + zc[7]*Xi*Eta;
      
            t21 = xc[1] + xc[4]*Eta + xc[5]*Zeta + xc[7]*Eta*Zeta;
            t22 = yc[1] + yc[4]*Eta + yc[5]*Zeta + yc[7]*Eta*Zeta;
            t23 = zc[1] + zc[4]*Eta + zc[5]*Zeta + zc[7]*Eta*Zeta;
         break;
      
         case 4:
            Xi = -1; Eta = p1; Zeta = p2;
            t11 = xc[3] + xc[5]*Xi + xc[6]*Eta + xc[7]*Xi*Eta;
            t12 = yc[3] + yc[5]*Xi + yc[6]*Eta + yc[7]*Xi*Eta;
            t13 = zc[3] + zc[5]*Xi + zc[6]*Eta + zc[7]*Xi*Eta;
            
            t21 = xc[2] + xc[4]*Xi + xc[6]*Zeta + xc[7]*Xi*Zeta;
            t22 = yc[2] + yc[4]*Xi + yc[6]*Zeta + yc[7]*Xi*Zeta;
            t23 = zc[2] + zc[4]*Xi + zc[6]*Zeta + zc[7]*Xi*Zeta;
         break;
      
         case 5:
            Xi = p2; Eta = p1; Zeta = 1;
            t11 = xc[1] + xc[4]*Eta + xc[5]*Zeta + xc[7]*Eta*Zeta;
            t12 = yc[1] + yc[4]*Eta + yc[5]*Zeta + yc[7]*Eta*Zeta;
            t13 = zc[1] + zc[4]*Eta + zc[5]*Zeta + zc[7]*Eta*Zeta;
      
            t21 = xc[2] + xc[4]*Xi + xc[6]*Zeta + xc[7]*Xi*Zeta;
            t22 = yc[2] + yc[4]*Xi + yc[6]*Zeta + yc[7]*Xi*Zeta;
            t23 = zc[2] + zc[4]*Xi + zc[6]*Zeta + zc[7]*Xi*Zeta;
         break;
   
         default:
            printf("Wrong local joint number \n" );
            return;
      }
      

      double n1,  n2,n3;
      n1 = t12*t23 - t13*t22;
      n2 = t13*t21 - t11*t23;
      n3 = t11*t22 - t12*t21;

      double normLen = sqrt(n1*n1 + n2*n2 + n3*n3);

      n1 /= normLen;
      n2 /= normLen;
      n3 /= normLen;

      m_d_freeSurfaceNormal_1[ThreadId] = n1;
      m_d_freeSurfaceNormal_2[ThreadId] = n2;
      m_d_freeSurfaceNormal_3[ThreadId] = n3;

   }



}



__global__ void C_calculate_normals_freeSlip(int N_FreeSurfvertex, int N_VertexPerCell,
                                             int* m_d_FreeSurfCellMapping,double* X_cord, double* Y_cord, double* Z_cord,
                                             double* m_d_nodalFunctionalRefValues,
                                             int* m_d_FreeSurfaceCells,int* m_d_freeSurfaceVertexLocal,int* m_d_freeSurfaceJoints,
                                             int* m_d_FreeSurfaceVertex,double* m_d_freeSurfaceNormal_1, double* m_d_freeSurfaceNormal_2, double* m_d_freeSurfaceNormal_3
                                             )
{
   double x[8];
   double y[8];
   double z[8];
   double xc[8];
   double yc[8];
   double zc[8];


   // printf("%d\n", threadIdx.x);
   // Each Thread Will Copy its set of XYZ Co-Ordinates to the GPU
   int ThreadId = blockIdx.x*blockDim.x + threadIdx.x; 

 
   if(ThreadId < N_FreeSurfvertex)
   {       
      int cellMapIndex    =  m_d_FreeSurfCellMapping[ThreadId];
      int shiftPosition   =  cellMapIndex * N_VertexPerCell;
         
      for ( int i = 0 ; i < N_VertexPerCell ; i++)
      {
         x[i] =   X_cord[shiftPosition + i];
         y[i] =   Y_cord[shiftPosition + i];
         z[i] =   Z_cord[shiftPosition + i];
         
      }
      

      // if(ThreadId == 80 )
      // {
      //    printf("\n Y Values : %f    %f   %f   %lf  %lf  %lf  %f   %f \n\n",  y[0], y[1],y[2], y[3],y[4],y[5],y[6],y[7]);
      //    printf("\n X Values : %f    %f   %f   %lf  %lf  %lf  %f   %f \n\n",  x[0], x[1],x[2], x[3],x[4],x[5],x[6],x[7]);
      //    printf("\n Z Values : %f    %f   %f   %lf  %lf  %lf  %f   %f \n\n",  z[0], z[1],z[2], z[3],z[4],z[5],z[6],z[7]);

      // }
      // Setting up Parameters for FE Cell -- HEXAtrilinear Mapping // 
      
      xc[0] = (x[0] + x[1] + x[3] + x[2] + x[4] + x[5] + x[7] + x[6]) * 0.125;
      xc[1] = (-x[0] + x[1] - x[3] + x[2] - x[4] + x[5] - x[7] + x[6]) * 0.125;
      xc[2] = (-x[0] - x[1] + x[3] + x[2] - x[4] - x[5] + x[7] + x[6]) * 0.125;
      xc[3] = (-x[0] - x[1] - x[3] - x[2] + x[4] + x[5] + x[7] + x[6]) * 0.125;
      xc[4] = (x[0] - x[1] - x[3] + x[2] + x[4] - x[5] - x[7] + x[6]) * 0.125;
      xc[5] = (x[0] - x[1] + x[3] - x[2] - x[4] + x[5] - x[7] + x[6]) * 0.125;
      xc[6] = (x[0] + x[1] - x[3] - x[2] - x[4] - x[5] + x[7] + x[6]) * 0.125;
      xc[7] = (-x[0] + x[1] + x[3] - x[2] + x[4] - x[5] - x[7] + x[6]) * 0.125;

      yc[0] = (y[0] + y[1] + y[3] + y[2] + y[4] + y[5] + y[7] + y[6]) * 0.125;
      yc[1] = (-y[0] + y[1] - y[3] + y[2] - y[4] + y[5] - y[7] + y[6]) * 0.125;
      yc[2] = (-y[0] - y[1] + y[3] + y[2] - y[4] - y[5] + y[7] + y[6]) * 0.125;
      yc[3] = (-y[0] - y[1] - y[3] - y[2] + y[4] + y[5] + y[7] + y[6]) * 0.125;
      yc[4] = (y[0] - y[1] - y[3] + y[2] + y[4] - y[5] - y[7] + y[6]) * 0.125;
      yc[5] = (y[0] - y[1] + y[3] - y[2] - y[4] + y[5] - y[7] + y[6]) * 0.125;
      yc[6] = (y[0] + y[1] - y[3] - y[2] - y[4] - y[5] + y[7] + y[6]) * 0.125;
      yc[7] = (-y[0] + y[1] + y[3] - y[2] + y[4] - y[5] - y[7] + y[6]) * 0.125;

      zc[0] = (z[0] + z[1] + z[3] + z[2] + z[4] + z[5] + z[7] + z[6]) * 0.125;
      zc[1] = (-z[0] + z[1] - z[3] + z[2] - z[4] + z[5] - z[7] + z[6]) * 0.125;
      zc[2] = (-z[0] - z[1] + z[3] + z[2] - z[4] - z[5] + z[7] + z[6]) * 0.125;
      zc[3] = (-z[0] - z[1] - z[3] - z[2] + z[4] + z[5] + z[7] + z[6]) * 0.125;
      zc[4] = (z[0] - z[1] - z[3] + z[2] + z[4] - z[5] - z[7] + z[6]) * 0.125;
      zc[5] = (z[0] - z[1] + z[3] - z[2] - z[4] + z[5] - z[7] + z[6]) * 0.125;
      zc[6] = (z[0] + z[1] - z[3] - z[2] - z[4] - z[5] + z[7] + z[6]) * 0.125;
      zc[7] = (-z[0] + z[1] + z[3] - z[2] + z[4] - z[5] - z[7] + z[6]) * 0.125;

      // hello();
      
      // Get the local Reference co-ordinates of the Vertex (xi eta ,zeta)
      double xi, eta , zeta;
      int localNodeNumber = m_d_freeSurfaceVertexLocal[ThreadId];

      if(localNodeNumber > N_VertexPerCell || localNodeNumber < 0)
      {
         printf(" Error : LocalVertex Number not valid");
         printf(" INFO   :  Function : calculate_normals_freeSurf , Class : FE3D_ALE, File : FE3D_ALE.cu");
         //Should exit the program -- CODE later                // THIVIN - TODO
      }

      xi    = m_d_nodalFunctionalRefValues[localNodeNumber];
      eta   = m_d_nodalFunctionalRefValues[localNodeNumber + N_VertexPerCell];
      zeta  = m_d_nodalFunctionalRefValues[localNodeNumber + 2*N_VertexPerCell];

   
      // Obtain the xi eta and zeta values required for 

      // Calculate Surface Normals  
      // Each Cuda Thread Will calculate Normal based on the value from its shared Array
      // Get the Tangent Vectors , and use their Cross product to optain the normal vector for the surface

      double xi_1, xi_2;
      
      int Jointno = m_d_freeSurfaceJoints[ThreadId];
      switch (Jointno)
      {
         case 0:
            xi_1 = xi; xi_2 = eta; break;
         case 1:
            xi_1 = zeta; xi_2 = xi; break;
         case 2:
            xi_1 = zeta; xi_2 = eta; break;
         case 3:
            xi_1 = zeta; xi_2 = -xi; break;
         case 4:
            xi_1 = eta; xi_2 = zeta;break;
         case 5:
            xi_1 = eta; xi_2 = xi;break;
         default:
         {
            printf(" Error in Class : Hexatrilinear.C - Function : GetRefValuesfromJointid \n");
            printf(" The Given Joint Id for the HEXAHEADRAL cell %d Does not exist ",Jointno);
         }
      }
         // if(ThreadId == 80)
         // {  
         //    printf(" Cell Mapping Number :  %d", cellMapIndex);
         //    printf(" xi_1 and xi_ 2 (dev):   %lf      %lf\n", xi_1 , xi_2);
         //    printf(" y1 and y4 , yc5 , yc7 (dev):   %lf    %lf  %lf    %lf \n", yc[1] ,yc[4],yc[5],yc[7] );
         //    printf(" Joint:   %d    \n", Jointno);

         // }
      // Find the Tangent Vectors for the Node 
      double p1 = xi_1; double p2 = xi_2;
      double t11,t12,t13,t21,t22,t23;
      double Xi,Eta,Zeta;
      switch(Jointno)
      {
         case 0:
            Xi = p1; Eta = p2; Zeta = -1;
            t11 = xc[2] + xc[4]*Xi + xc[6]*Zeta + xc[7]*Xi*Zeta;
            t12 = yc[2] + yc[4]*Xi + yc[6]*Zeta + yc[7]*Xi*Zeta;
            t13 = zc[2] + zc[4]*Xi + zc[6]*Zeta + zc[7]*Xi*Zeta;
      
            t21 = xc[1] + xc[4]*Eta + xc[5]*Zeta + xc[7]*Eta*Zeta;
            t22 = yc[1] + yc[4]*Eta + yc[5]*Zeta + yc[7]*Eta*Zeta;
            t23 = zc[1] + zc[4]*Eta + zc[5]*Zeta + zc[7]*Eta*Zeta;
         break;
      
         case 1:
            Xi = p2; Eta = -1; Zeta = p1;
            t11 = xc[1] + xc[4]*Eta + xc[5]*Zeta + xc[7]*Eta*Zeta;
            t12 = yc[1] + yc[4]*Eta + yc[5]*Zeta + yc[7]*Eta*Zeta;
            t13 = zc[1] + zc[4]*Eta + zc[5]*Zeta + zc[7]*Eta*Zeta;
      
            t21 = xc[3] + xc[5]*Xi + xc[6]*Eta + xc[7]*Xi*Eta;
            t22 = yc[3] + yc[5]*Xi + yc[6]*Eta + yc[7]*Xi*Eta;
            t23 = zc[3] + zc[5]*Xi + zc[6]*Eta + zc[7]*Xi*Eta;
         break;
      
         case 2:
            Xi = 1; Eta = p2; Zeta = p1;
            t11 = xc[2] + xc[4]*Xi + xc[6]*Zeta + xc[7]*Xi*Zeta;
            t12 = yc[2] + yc[4]*Xi + yc[6]*Zeta + yc[7]*Xi*Zeta;
            t13 = zc[2] + zc[4]*Xi + zc[6]*Zeta + zc[7]*Xi*Zeta;
      
            t21 = xc[3] + xc[5]*Xi + xc[6]*Eta + xc[7]*Xi*Eta;
            t22 = yc[3] + yc[5]*Xi + yc[6]*Eta + yc[7]*Xi*Eta;
            t23 = zc[3] + zc[5]*Xi + zc[6]*Eta + zc[7]*Xi*Eta;
         break;
      
         case 3:
            Xi = -p2; Eta = 1; Zeta = p1;
            t11 = xc[3] + xc[5]*Xi + xc[6]*Eta + xc[7]*Xi*Eta;
            t12 = yc[3] + yc[5]*Xi + yc[6]*Eta + yc[7]*Xi*Eta;
            t13 = zc[3] + zc[5]*Xi + zc[6]*Eta + zc[7]*Xi*Eta;
      
            t21 = xc[1] + xc[4]*Eta + xc[5]*Zeta + xc[7]*Eta*Zeta;
            t22 = yc[1] + yc[4]*Eta + yc[5]*Zeta + yc[7]*Eta*Zeta;
            t23 = zc[1] + zc[4]*Eta + zc[5]*Zeta + zc[7]*Eta*Zeta;
         break;
      
         case 4:
            Xi = -1; Eta = p1; Zeta = p2;
            t11 = xc[3] + xc[5]*Xi + xc[6]*Eta + xc[7]*Xi*Eta;
            t12 = yc[3] + yc[5]*Xi + yc[6]*Eta + yc[7]*Xi*Eta;
            t13 = zc[3] + zc[5]*Xi + zc[6]*Eta + zc[7]*Xi*Eta;
            
            t21 = xc[2] + xc[4]*Xi + xc[6]*Zeta + xc[7]*Xi*Zeta;
            t22 = yc[2] + yc[4]*Xi + yc[6]*Zeta + yc[7]*Xi*Zeta;
            t23 = zc[2] + zc[4]*Xi + zc[6]*Zeta + zc[7]*Xi*Zeta;
         break;
      
         case 5:
            Xi = p2; Eta = p1; Zeta = 1;
            t11 = xc[1] + xc[4]*Eta + xc[5]*Zeta + xc[7]*Eta*Zeta;
            t12 = yc[1] + yc[4]*Eta + yc[5]*Zeta + yc[7]*Eta*Zeta;
            t13 = zc[1] + zc[4]*Eta + zc[5]*Zeta + zc[7]*Eta*Zeta;
      
            t21 = xc[2] + xc[4]*Xi + xc[6]*Zeta + xc[7]*Xi*Zeta;
            t22 = yc[2] + yc[4]*Xi + yc[6]*Zeta + yc[7]*Xi*Zeta;
            t23 = zc[2] + zc[4]*Xi + zc[6]*Zeta + zc[7]*Xi*Zeta;
         break;
   
         default:
            printf("Wrong local joint number \n" );
            return;
      }
      

      double n1,  n2,n3;
      n1 = t12*t23 - t13*t22;
      n2 = t13*t21 - t11*t23;
      n3 = t11*t22 - t12*t21;

      double normLen = sqrt(n1*n1 + n2*n2 + n3*n3);

      n1 /= normLen;
      n2 /= normLen;
      n3 /= normLen;

      m_d_freeSurfaceNormal_1[ThreadId] = n1;
      m_d_freeSurfaceNormal_2[ThreadId] = n2;
      m_d_freeSurfaceNormal_3[ThreadId] = n3;
   }



}



// ------------------ HOST FUNCTIONS -------------------------- //

void FE3D_ALE::C_calculateNormal_HostWrapper()
{
   int MAX_THREAD_PER_BLOCK = 64;
   int N_threads; 
   if(N_freeSurfaceVertex >= MAX_THREAD_PER_BLOCK) 
      N_threads = MAX_THREAD_PER_BLOCK;
   else
      N_threads = N_freeSurfaceVertex;
   int C_NUM_BLOCKS = std::ceil(double(N_freeSurfaceVertex)/MAX_THREAD_PER_BLOCK);


   dim3 dimGrid(C_NUM_BLOCKS);
   dim3 dimBlock(N_threads);

   // cout << " Entered Wrapperrrrrr for CUDA Function - GRID : " << C_NUM_BLOCKS << " Threads : " << N_threads << " numvert : " <<N_bd_FreeSlip_Vertex  <<endl;


   C_calculate_normals_freeSurf<<<dimGrid,N_threads, 0, FE3d_stream>>>(N_freeSurfaceVertex,m_h_N_VertexPerCell,m_d_FreeSurfCellMapping,
                                                                   X_cord, Y_cord,Z_cord,m_d_nodalFunctionalRefValues,
                                                                   m_d_FreeSurfaceCells,m_d_freeSurfaceVertexLocal, m_d_freeSurfaceJoints,
                                                                   m_d_FreeSurfaceVertex,m_d_freeSurfaceNormal_1,m_d_freeSurfaceNormal_2,m_d_freeSurfaceNormal_3
                                                                  );
   // cudaDeviceSynchronize();   
   // cout << " COmpleed kerner call"<<endl;
   checkCudaErrors(cudaMemcpyAsync(&Bd_FreeSurf_normal_1[0], m_d_freeSurfaceNormal_1 ,N_freeSurfaceVertex*sizeof(double),cudaMemcpyDeviceToHost,FE3d_stream )  );
   checkCudaErrors(cudaMemcpyAsync(&Bd_FreeSurf_normal_2[0], m_d_freeSurfaceNormal_2, N_freeSurfaceVertex*sizeof(double),cudaMemcpyDeviceToHost,FE3d_stream )  );
   checkCudaErrors(cudaMemcpyAsync(&Bd_FreeSurf_normal_3[0], m_d_freeSurfaceNormal_3 ,N_freeSurfaceVertex*sizeof(double),cudaMemcpyDeviceToHost,FE3d_stream )  );

}


void FE3D_ALE::C_calculateNormal_freeSlip_HostWrapper()
{

   int MAX_THREAD_PER_BLOCK = 64;
   int N_threads; 
   if(N_bd_FreeSlip_Vertex >= MAX_THREAD_PER_BLOCK) 
      N_threads = MAX_THREAD_PER_BLOCK;
   else
      N_threads = N_bd_FreeSlip_Vertex;
   int C_NUM_BLOCKS = std::ceil(double(N_bd_FreeSlip_Vertex)/MAX_THREAD_PER_BLOCK);


   dim3 dimGrid(C_NUM_BLOCKS);
   dim3 dimBlock(N_threads);


   // cout << " Entered Wrapperrrrrr for CUDA Function free Slip  - GRID : " << C_NUM_BLOCKS << " Threads : " << N_threads << " numvert : " <<N_bd_FreeSlip_Vertex  <<endl;


   C_calculate_normals_freeSlip <<< C_NUM_BLOCKS, N_threads,0, FE3d_stream_slip>>>(N_bd_FreeSlip_Vertex,m_h_N_VertexPerCell,m_d_Free_Slip_CellMapping,
                                                                   X_cord_slip, Y_cord_slip,Z_cord_slip,m_d_nodalFunctionalRefValues_slip,
                                                                   m_d_Free_Slip_Cells,m_d_free_Slip_VertexLocal, m_d_free_Slip_Joints,
                                                                   m_d_Free_Slip_Vertex,m_d_free_Slip_Normal_1,m_d_free_Slip_Normal_2,m_d_free_Slip_Normal_3
                                                                   );

   checkCudaErrors(cudaMemcpyAsync(Bd_normal_1.data(), m_d_free_Slip_Normal_1 ,N_bd_FreeSlip_Vertex*sizeof(double),cudaMemcpyDeviceToHost,FE3d_stream_slip )  );
   checkCudaErrors(cudaMemcpyAsync(Bd_normal_2.data(), m_d_free_Slip_Normal_2, N_bd_FreeSlip_Vertex*sizeof(double),cudaMemcpyDeviceToHost,FE3d_stream_slip )  );
   checkCudaErrors(cudaMemcpyAsync(Bd_normal_3.data(), m_d_free_Slip_Normal_3 ,N_bd_FreeSlip_Vertex*sizeof(double),cudaMemcpyDeviceToHost,FE3d_stream_slip )  );

   // cout << " NOrm test : " << Bd_normal_1[0] << "  " << Bd_normal_2[0] << "   " << Bd_normal_3[0] <<endl;
}



void FE3D_ALE::clearCudaVariables()
{
   checkCudaErrors(cudaStreamDestroy(FE3d_stream));
   if(m_d_nodalFunctionalRefValues) checkCudaErrors(cudaFree( m_d_nodalFunctionalRefValues));
   if(m_d_FreeSurfCellMapping) checkCudaErrors(cudaFree( m_d_FreeSurfCellMapping));
   if(m_d_FreeSurfaceVertex) checkCudaErrors(cudaFree( m_d_FreeSurfaceVertex));
   if(m_d_FreeSurfaceCells) checkCudaErrors(cudaFree( m_d_FreeSurfaceCells));
   if(m_d_freeSurfaceJoints) checkCudaErrors(cudaFree( m_d_freeSurfaceJoints));
   if(m_d_freeSurfaceVertexLocal) checkCudaErrors(cudaFree( m_d_freeSurfaceVertexLocal));
   if(m_d_freeSurfaceNormal_1) checkCudaErrors(cudaFree( m_d_freeSurfaceNormal_1));
   if(m_d_freeSurfaceNormal_2) checkCudaErrors(cudaFree( m_d_freeSurfaceNormal_2));
   if(m_d_freeSurfaceNormal_3) checkCudaErrors(cudaFree( m_d_freeSurfaceNormal_3));
   if(X_cord) checkCudaErrors(cudaFree( X_cord));
   if(Y_cord) checkCudaErrors(cudaFree( Y_cord));
   if(Z_cord) checkCudaErrors(cudaFree( Z_cord));

   checkCudaErrors(cudaStreamDestroy(FE3d_stream_slip));
   
   if(m_d_nodalFunctionalRefValues_slip) checkCudaErrors(cudaFree( m_d_nodalFunctionalRefValues_slip));
   if(m_d_Free_Slip_CellMapping) checkCudaErrors(cudaFree( m_d_Free_Slip_CellMapping));
   if(m_d_Free_Slip_Vertex) checkCudaErrors(cudaFree( m_d_Free_Slip_Vertex));
   if(m_d_Free_Slip_Cells) checkCudaErrors(cudaFree( m_d_Free_Slip_Cells));
   if(m_d_free_Slip_Joints) checkCudaErrors(cudaFree( m_d_free_Slip_Joints));
   if(m_d_free_Slip_VertexLocal) checkCudaErrors(cudaFree( m_d_free_Slip_VertexLocal));
   if(m_d_free_Slip_Normal_1) checkCudaErrors(cudaFree( m_d_free_Slip_Normal_1));
   if(m_d_free_Slip_Normal_2) checkCudaErrors(cudaFree( m_d_free_Slip_Normal_2));
   if(m_d_free_Slip_Normal_3) checkCudaErrors(cudaFree( m_d_free_Slip_Normal_3));
   if(X_cord_slip) checkCudaErrors(cudaFree( X_cord_slip));
   if(Y_cord_slip) checkCudaErrors(cudaFree( Y_cord_slip));
   if(Z_cord_slip) checkCudaErrors(cudaFree( Z_cord_slip));


}