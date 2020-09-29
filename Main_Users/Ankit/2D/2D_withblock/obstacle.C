#include <iostream>
#include <math.h>

using namespace std;

#define B 1 // { 0 1 1.35 }
#define sigma 3.46e6
#define pho 6360
#define mu 3.4e-7
#define a_duct .04

int main(int argc, char*argv[])
{
   double d = .01;
   double w =4;      // w = w/d
   double a = 4;      // a = a/d
   double Lu = 12;   // Lu = 12d/d
   double Ld = 42;   // Lu = 42d/d

   double xostart = -1.0/2;        // xostart = (-d/2)/d 
   double xoend = 1.0/2;           // xoend = (d/2)/d 
   double yostart = -1.0/2;        // yostart = (-d/2)/d 
   double yoend = 1.0/2 ;          // yoend = (d/2)/d 

   int division = 30;
   double x[division + 1] , y[division + 1] ; 
   int i = 0;
   
   for(i=0;i<=division;i++)
   {
     y[i] = yostart + i*(yoend-yostart)/division ;
     cout<<xostart<<"       "<<y[i]<<endl;
   }   
   
   for(i=1;i<=division;i++)
   {
     x[i] = xostart + i*(xoend-xostart)/division  ;
     cout<<x[i]<<"       "<<yoend<<endl;
   }
   
   for(i=1;i<=division;i++)
   {
     y[i] = yoend - i*(yoend-yostart)/division;
     cout<<xoend<<"       "<<y[i]<<endl;
   }   
   
   for(i=1;i<=division;i++)
   {
     x[i] = xoend - i*(xoend-xostart)/division ;
     cout<<x[i]<<"       "<<yostart<<endl;
   }
   

   



   
 /*  for(i=0;i<=division;i++)
   {
     x[i] = i*(xoend-xostart)/division + xostart ;
     cout<<x[i]<<"       "<<yostart<<endl;
   }
   
   for(i=0;i<=division;i++)
   {
     y[i] = i*(yoend-yostart)/division + yostart ;
     cout<<xoend<<"       "<<y[i]<<endl;
   }
   
   for(i=0;i<=division;i++)
   {
     x[i] = xoend - i*(xoend-xostart)/division  ;
     cout<<x[i]<<"       "<<yoend<<endl;
   }

   for(i=0;i<=division;i++)
   {
     y[i] = yoend - i*(yoend-yostart)/division ;
     cout<<xostart<<"       "<<y[i]<<endl;
   }
   */
}
