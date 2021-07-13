
# %%
import numpy as np
import math

from numpy.lib.npyio import save


def RK4(t0,phi0,timeStep,N_S):
    x0 = t0
    y0 = phi0
    h  = timeStep

    ### Use the global Variables
    global Dx_p_titlde
    global Dy_p_titlde

    global u_titlde  
    global Dx_u_titlde
    global DDx_u_titlde
    global Dy_u_titlde
    global DDy_u_titlde

    global v_titlde  
    global Dx_v_titlde
    global DDx_v_titlde
    global Dy_v_titlde
    global DDy_v_titlde

    global u_bar
    global v_bar
    global Dx_u_bar
    global Dy_u_bar
    global Dx_v_bar
    global Dy_v_bar

    global N_U
    global R_E

    R_Bar = 1/R_E

    ### -- Get the ith Component from the Flat Array --- ####
    def comp(vec,i):
        return vec[i]

    # define the Inner product function
    def inner(vec_a,vec_b):
        val = 0.0
        size = vec_a.shape[0]
        if(vec_a.shape[0] == vec_b.shape[0]):
            for i in range(size):
                val += vec_a[i]*vec_b[i]
        else:
            print ("Error in dimensions of the vectors")
        return val 
        

    ## -- Define the derivative function -- #
    def function(x0):
        #Size of Phi0 - N_S - Declare a array
        val = np.zeros((N_R,N_S))


        for i in range(N_S):
            for a in range(N_S):
                val[:,i] +=  - phi0[:,a]*inner(Dx_p_titlde[:,a],u_titlde[:,i])  \
                           - phi0[:,a]*inner(Dy_p_titlde[:,a],v_titlde[:,i])  \
                           + R_Bar*phi0[:,a]*inner( (DDx_u_titlde[:,a] + DDy_u_titlde[:,a]),u_titlde[:,i] ) \
                           + R_Bar*phi0[:,a]*inner( (DDx_v_titlde[:,a] + DDy_v_titlde[:,a]),v_titlde[:,i] )  \
                           + phi0[:,a] * inner ( ( (u_bar*Dx_u_titlde[:,a][i]) + (u_titlde[:,a]*Dx_u_bar[i]) + (v_bar*Dy_u_titlde[:,a][i]) + (v_titlde[:,a],Dy_u_bar[i]) ) ,u_titlde[:,i] )  \
                           + phi0[:,a] * inner ( ( (v_bar*Dy_v_titlde[:,a][i]) + (v_titlde[:,a]*Dy_v_bar[i]) + (u_bar*Dx_v_titlde[:,a][i]) + (u_titlde[:,a],Dx_v_bar[i]) ) ,v_titlde[:,i] )  \
                
                for b in range(N_S):
                    


    
    k1 = function(x0,y0)
    k2 = function(x0+(h/2), y0+k1/2)
    k3 = function(x0+(h/2), y0+k2/2)
    k4 = function(x0+h, y0+k3)

    k = h*(k1 + 2*k2 + 2*k3 + k4)/ 6
    y1 = y0 + k

# Declare the Global Array Sizes
N_U = 2
N_S = 2
R_E = 1
N_R = 4


# Declare the Global Arrays
Dx_p_titlde   = np.ones((N_U , N_S))
Dy_p_titlde   = np.ones((N_U , N_S))

u_titlde     = np.ones((N_U , N_S))
Dx_u_titlde   = np.ones((N_U , N_S))
DDx_u_titlde  = np.ones((N_U , N_S))
Dy_u_titlde   = np.ones((N_U , N_S))
DDy_u_titlde  = np.ones((N_U , N_S))

v_titlde     = np.ones((N_U , N_S))
Dx_v_titlde   = np.ones((N_U , N_S))
DDx_v_titlde  = np.ones((N_U , N_S))
Dy_v_titlde   = np.ones((N_U , N_S))
DDy_v_titlde  = np.ones((N_U , N_S))

phi0          = np.ones((N_R,N_S))

u_bar        = np.ones(N_U )
v_bar        = np.ones(N_U )
Dx_u_bar     = np.ones(N_U )
Dy_u_bar     = np.ones(N_U)
Dx_v_bar     = np.ones(N_U )
Dy_v_bar     = np.ones(N_U)




# %%
