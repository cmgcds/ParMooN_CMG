double inner(double* a, double* b, int size)
{
    double val = 0;
    for( int k = 0 ; k < size ; k++)
        val += a[k]*b[k];
    return val;
}


double* addVec(double* a, double* b , int size,double* c)
{
    for ( int k = 0 ;k < size ; k++)
    {
        c[k] = a[k] + b[k];
    }
    return c;
}

double* addVec4(double* a, double* b ,double* c,double* d,int size, double* ans)
{
    for ( int k = 0 ;k < size ; k++)
    {
        ans[k] = a[k] + b[k] + c[k] + d[k];
    }
    return ans;
}


double* multVec (double* a, double* b , int size , double* c)
{
    for (int k = 0 ; k < size ; k++)
        c[k] = a[k]*b[k];
       
}

void Function(double* ans,double h, double* CoeffVector, double* t0, double* u_tilde,double* v_tilde, 
				  double* Dx_u_tilde,double* Dy_u_tilde,double* DDx_u_tilde,double* DDy_u_tilde,
				  double* DDx_v_tilde,double* DDy_v_tilde,
				  double* Dx_v_tilde, double* Dy_v_tilde,double* u_bar , double* v_bar, double* Dx_u_bar, double* Dy_u_bar,
				  double* Dx_v_bar, double* Dy_v_bar, double* p_tilde, double* Dx_p_tilde, double* Dy_p_tilde ,double* covariance)
{

    int N_R = 2;
    int N_S = 2;
    int R   = 1.0;
    int N_U = 4;
    double* temp = new double[N_U]();
    double* temp1 = new double[N_U]();
    double* temp2 = new double[N_U]();
    double* temp3 = new double[N_U]();
    double* temp4 = new double[N_U]();

    for ( int i=0 ; i < N_S ; i++)
    {
        // Get the ith component of Coeff vector
        double* localCoeff = ans + N_R*i;     // Answer vector
        double* u_tilde_i  = u_tilde + N_R*i;
        double* v_tilde_i  = v_tilde + N_R*i;

        for (int g=0 ; g<N_R ; g++)
        {
            for ( int a=0 ; a < N_S ; a++)
            {
                double* dx_p_tilde_a    = Dx_p_tilde + N_U*a;
                double* dy_p_tilde_a    = Dx_p_tilde + N_U*a;
                double* ddx_u_tilde_a   = DDx_u_tilde + N_U*a;
                double* ddx_v_tilde_a   = DDx_v_tilde + N_U*a;
                double* ddy_u_tilde_a   = DDy_u_tilde + N_U*a;
                double* ddy_v_tilde_a   = DDy_v_tilde + N_U*a;
                double* dx_u_tilde_a    = Dx_u_tilde + N_U*a;
                double* dx_v_tilde_a    = Dx_v_tilde + N_U*a;
                double* dy_u_tilde_a    = Dy_u_tilde + N_U*a;
                double* dy_v_tilde_a    = Dy_v_tilde + N_U*a;
                double* phi_a           = CoeffVector + N_R*a;
                double* u_tilde_a       = u_tilde + N_U*a;
                double* v_tilde_a       = v_tilde + N_U*a;

                localCoeff[g] +=  - phi_a[g] * inner(dx_p_tilde_a,u_tilde_i,N_U)  
                                  - phi_a[g] * inner(dy_p_tilde_a,v_tilde_i,N_U) ;
                
                    temp   = addVec(ddx_u_tilde_a,ddy_u_tilde_a,N_U,temp);

                localCoeff[g] +=    R * inner(temp,u_tilde_i,N_U);

                    temp =   addVec(ddx_v_tilde_a,ddy_v_tilde_a,N_U,temp);

                localCoeff[g] +=    phi_a[g] * R * inner(temp,v_tilde_i,N_U);
                
                    temp1 = multVec(u_bar,dx_u_tilde_a,N_U,temp1);
                    temp2 = multVec(u_tilde_a,Dx_u_bar,N_U,temp2);
                    temp3 = multVec(v_bar,dy_u_tilde_a,N_U,temp3);
                    temp4 = multVec(v_tilde_a,Dy_u_bar,N_U,temp4);

                    temp = addVec4(temp1,temp2,temp3,temp4,N_U,temp);

                localCoeff[g] += phi_a[g]*inner(temp,u_tilde_i,N_U);

                    temp1 = multVec(v_bar,dy_v_tilde_a,N_U,temp1);
                    temp2 = multVec(v_tilde_a,Dy_v_bar,N_U,temp2);
                    temp3 = multVec(u_bar,dx_v_tilde_a,N_U,temp3);
                    temp4 = multVec(u_tilde_a,Dx_v_bar,N_U,temp4);

                    temp = addVec4(temp1,temp2,temp3,temp4,N_U,temp);

                localCoeff[g] += phi_a[g]*inner(temp,v_tilde_i,N_U);

                for ( int b = 0 ; b < N_S ; b++)
                {
                    double* phi_b           = CoeffVector + N_R*b;
                    double* dx_u_tilde_b    = Dx_u_tilde + N_U*b;
                    double* dx_v_tilde_b    = Dx_v_tilde + N_U*b;
                    double* dy_u_tilde_b    = Dy_u_tilde + N_U*b;
                    double* dy_v_tilde_b    = Dy_v_tilde + N_U*b;
                    
                    
                    temp1 = multVec(u_tilde_a,dx_u_tilde_b,N_U,temp1);
                    temp2 = multVec(v_tilde_a,dy_u_tilde_b,N_U,temp2);
                    
                    temp = addVec(temp1,temp2,N_U,temp);


                    localCoeff[g] += - ( phi_b[g]*phi_a[g] - covariance[a*N_S + b])*inner(temp,u_tilde_i,N_U);

                    temp1 = multVec(v_tilde_a,dy_v_tilde_b,N_U,temp1);
                    temp2 = multVec(u_tilde_a,dx_v_tilde_b,N_U,temp2);
                    
                    temp = addVec(temp1,temp2,N_U,temp);

                    localCoeff[g] += - ( phi_b[g]*phi_a[g] - covariance[a*N_S + b])*inner(temp,v_tilde_i,N_U);
                }
                    

            }

        }
    }
  
    
}


void Daxpy_local(double* y, double a, double* x ,double* b, int size ,int columns )
{
    for ( int col = 0 ; col < columns; col++)
    {
        double* yy = y + size*col;
        double* xx = x + size*col;
        double* bb = b + size*col;

        for ( int i = 0 ; i < size; i++)
        {
            b[i] = a*xx[i] + yy[i];
        }
    }

}


void Solve_RK4(double h, double* CoeffVector, double* t0, double* u_tilde,double* v_tilde, 
				  double* Dx_u_tilde,double* Dy_u_tilde,double* DDx_u_tilde,double* DDy_u_tilde,
				  double* DDx_v_tilde,double* DDy_v_tilde,
				  double* Dx_v_tilde, double* Dy_v_tilde,double* u_bar , double* v_bar, double* Dx_u_bar, double* Dy_u_bar,
				  double* Dx_v_bar, double* Dy_v_bar, double* p_tilde, double* Dx_p_tilde, double* Dy_p_tilde ,double* covariance)
{

    int N_R = 2;
    int N_S = 2;
    int R   = 1.0;
    int N_U = 4;

    double* originalx0 = new double[N_R*N_S];
    memcpy(originalx0,CoeffVector,N_R*N_S*SizeOfDouble);
    double* k1 = new double[N_R*N_S];
    double* k2 = new double[N_R*N_S];
    double* k3 = new double[N_R*N_S];
    double* k4 = new double[N_R*N_S];

    double x0 = 0;
    
    Function(k1, h, CoeffVector, t0, u_tilde, v_tilde, 
				 Dx_u_tilde, Dy_u_tilde, DDx_u_tilde, DDy_u_tilde,
				 DDx_v_tilde, DDy_v_tilde,
				 Dx_v_tilde, Dy_v_tilde ,u_bar , v_bar, Dx_u_bar,Dy_u_bar,
				 Dx_v_bar,Dy_v_bar, p_tilde, Dx_p_tilde, Dy_p_tilde , covariance);


    // Perform y1 = y0 + k1/2   // N_S Columns available 
    Daxpy_local(CoeffVector,0.5,k1,originalx0,N_R, N_S);


    Function(k2, h, CoeffVector, t0, u_tilde, v_tilde, 
				 Dx_u_tilde, Dy_u_tilde, DDx_u_tilde, DDy_u_tilde,
				 DDx_v_tilde, DDy_v_tilde,
				 Dx_v_tilde, Dy_v_tilde ,u_bar , v_bar, Dx_u_bar,Dy_u_bar,
				 Dx_v_bar,Dy_v_bar, p_tilde, Dx_p_tilde, Dy_p_tilde , covariance);


        // Perform y1 = y0 + k2/2   // N_S Columns available 
    Daxpy_local(CoeffVector,0.5,k2,originalx0,N_R, N_S);

    Function(k3, h, CoeffVector, t0, u_tilde, v_tilde, 
                Dx_u_tilde, Dy_u_tilde, DDx_u_tilde, DDy_u_tilde,
                DDx_v_tilde, DDy_v_tilde,
                Dx_v_tilde, Dy_v_tilde ,u_bar , v_bar, Dx_u_bar,Dy_u_bar,
                Dx_v_bar,Dy_v_bar, p_tilde, Dx_p_tilde, Dy_p_tilde , covariance);

    // Perform y1 = y0 + k2/2   // N_S Columns available 
    Daxpy_local(CoeffVector,1,k2,originalx0,N_R, N_S);
    Function(k4, h, CoeffVector, t0, u_tilde, v_tilde, 
                Dx_u_tilde, Dy_u_tilde, DDx_u_tilde, DDy_u_tilde,
                DDx_v_tilde, DDy_v_tilde,
                Dx_v_tilde, Dy_v_tilde ,u_bar , v_bar, Dx_u_bar,Dy_u_bar,
                Dx_v_bar,Dy_v_bar, p_tilde, Dx_p_tilde, Dy_p_tilde , covariance);
    

    // Add the Final Solution to the original x0
    
    for ( int i = 0 ; i < N_R*N_S ; i++)
        CoeffVector[i] = originalx0[i] + h*k4[i];





    delete[] k1;
    delete[] k2;
    delete[] k3;
    delete[] k4;
    delete[] originalx0;


}




