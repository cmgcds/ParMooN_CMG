

##################################################################################
############# @ Author : Thivin Anandh                          ###################
###########

import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
import sys
import plotly.express as px
import plotly.graph_objects as go


from SetParameters import *

## Returns two lists of the Quadrature points and Quadrature Values 
## Based on the given quad rule
def gaussianQuadrature(quadRule):

    if(quadRule == 2): 
        quadPoints = [0] * 2
        quadWeights = [0] * 2
        quadWeights[0] = 1.0000000000000000
        quadWeights[1] = 1.0000000000000000
        quadPoints[0]  = -0.5773502691896257
        quadPoints[1]  = 0.5773502691896257
    
    elif(quadRule == 3): 
        quadPoints = [0] * 3
        quadWeights = [0] * 3
        quadWeights[0] = 0.8888888888888888
        quadWeights[1] = 0.5555555555555556
        quadWeights[2] = 0.5555555555555556
        quadPoints[0]  = 0.0000000000000000
        quadPoints[1]  = -0.7745966692414834
        quadPoints[2] = 0.7745966692414834
    
    elif(quadRule == 4):
        quadPoints = [0] * 4
        quadWeights = [0] * 4
        quadWeights[0] = 0.6521451548625461
        quadWeights[1] = 0.6521451548625461
        quadWeights[2] = 0.3478548451374538
        quadWeights[3] = 0.3478548451374538
        quadPoints[0]  = -0.3399810435848563
        quadPoints[1]  = 0.3399810435848563
        quadPoints[2] = -0.8611363115940526
        quadPoints[3] = 0.8611363115940526

    else:
        print(" Choose a correct Quadrature Formula ") 
    

    return quadPoints,quadWeights




## Define the fuinction value at a quadrature point

def funcValQuad(q,val1,val2):
    shapeFunc2 =     ((q + 1) * 0.5)
    shapeFunc1 =    1.0 - shapeFunc2

    ## HARDCODED -- Assumed that the FE is always first order
    value  = val1*shapeFunc1 + val2*shapeFunc2
    return value


######  --- L2 Error -------- ########
def l2Error(AnalyticalSolution,eps):
    N_Cells = len(AnalyticalSolution) - 1
    x =  np.linspace(0,1,N_Cells+1)

    l2ErrorVal = 0.0
    l1ErrorVal = 0.0

    # get the Quadrature point and QuadratureW
    quadPoints,quadWeights = gaussianQuadrature(QUAD_RULE_NO)

    for i in range(N_Cells):
        start_x = x[i]
        end_x   = x[i+1]
        detJk   = (end_x - start_x)/2.0
        for q in range(len(quadPoints)):
            
            predSol   = funcValQuad(quadPoints[q],AnalyticalSolution[i],AnalyticalSolution[i+1])
            xtemp     = ((quadPoints[q] + 1) * 0.5)   # interpolates the q to [0,1] region
            xVal      = start_x + (end_x - start_x)*xtemp
            actualSol = ExactSolution(xVal,eps)
            l2ErrorVal += quadWeights[q] *detJk * (predSol-actualSol)**2
            l1ErrorVal += quadWeights[q] *detJk * abs(actualSol - predSol)
        # print( l2ErrorVal)
    l2ErrorVal =  np.sqrt(l2ErrorVal)                 

    print(" THE TOTAL L2 ERROR : ",  l2ErrorVal)
    print(" THE TOTAL L1 ERROR : ",  l1ErrorVal)
    return l2ErrorVal,l1ErrorVal


##### ----- RMSE ERROR -----  #####
def RMSE_Error(AnalyticalSolution,eps):
    N_Cells = len(AnalyticalSolution) - 1
    x =  np.linspace(0,1,N_Cells+1)

    MaxError  = -9999999999;
    RMSE_ErrorVal = 0.0
    smallL2   = 0.0
    for i in range(len(AnalyticalSolution)):
        actualSol = ExactSolution(x[i],eps)
        error = actualSol - AnalyticalSolution[i]
        RMSE_ErrorVal += (error)**2 /  len(AnalyticalSolution)
        smallL2    += (error)**2
        if( abs(actualSol - AnalyticalSolution[i]) > MaxError):
            MaxError = abs(actualSol - AnalyticalSolution[i])

    RMSE_ErrorVal = np.sqrt(RMSE_ErrorVal)
    smallL2       = np.sqrt(smallL2)
    print(" SMALL L2 ERROR : " , smallL2 )
    print(" RMSE_ERROR     : " , RMSE_ErrorVal )
    print(" THE TOTAL L-inf ERROR : " , MaxError )

    return smallL2, RMSE_ErrorVal, MaxError





########################### _________ MAIN CODE _____________ #################################
if __name__ == "__main__":
    ###GRAPH PLOTTING LISTS
    # epsnumbersList = [100000,1000000,1000000,10000000,100000000,1000000000,10000000000,100000000000]

    folderName = str(sys.argv[1])

    if(len(sys.argv) > 2):
        iterName  = str(sys.argv[2])

    ann_l2Error_list = [] 
    ann_l1Error_list = [] 
    ann_linfError_list = [] 
    ann_l2SmallError_list = [] 
    ann_RMSE_list = [] 


    analytical_l2Error_list = [] 
    analytical_l1Error_list = [] 
    analytical_linfError_list = [] 
    analytical_l2SmallError_list = [] 
    analytical_RMSE_list = [] 


    ann_l2Error_list_worst = [] 
    ann_l1Error_list_worst = [] 
    ann_linfError_list_worst = [] 
    ann_l2SmallError_list_worst = [] 
    ann_RMSE_list_worst = [] 


    analytical_l2Error_list_worst = [] 
    analytical_l1Error_list_worst = [] 
    analytical_linfError_list_worst = [] 
    analytical_l2SmallError_list_worst = [] 
    analytical_RMSE_list_worst = [] 

    epsValueList = []

    for list in range(len(epsnumbersList)):
        print (" ===========================================================")
        print(  " FOR PE NR : " , str( epsnumbersList[list])  )
        print (" ===========================================================")
        pe_nr = int(epsnumbersList[list])
        eps   = float(1.0/pe_nr)
        no = int(10**list)
        epsValueList.append(str(eps))
        

        filename  = folderName +  "/supgSol_" + str(epsnumbersList[list]) + ".csv"

        if(len(sys.argv) > 2):
            iterName  = str(sys.argv[2])
            filename  = folderName +  "/supgSol_" + str(epsnumbersList[list]) + "_" + str(iterName) + ".csv"


        # rawData_worst = genfromtxt(filenameWorst,delimiter=',')
        rawData = genfromtxt(filename, delimiter=',')

        # Remove the Last columns
        rawData = np.delete(rawData,rawData.shape[1]-1,1)
        # rawData_worst = np.delete(rawData_worst,rawData_worst.shape[1]-1,1)

        AnalyticalTauSolution = rawData[0]
        ComputedTauSolution   = rawData[1]
        GalerkingSolution     = rawData[2]

        print(" ERROR -- Computed ANN vs the Exact Solution - Best Model")
        a,b = l2Error(rawData[1],eps)
        c,d,e = RMSE_Error(rawData[1],eps)

        ann_l2Error_list.append(a)
        ann_l1Error_list.append(b)
        ann_l2SmallError_list.append(c)
        ann_RMSE_list.append(d)
        ann_linfError_list.append(e)

        print(" ERROR -- Analytical Tau vs the Exact Solution - Best Model")
        a,b = l2Error(rawData[0],eps)
        c,d,e = RMSE_Error(rawData[0],eps)

        analytical_l2Error_list.append(a)
        analytical_l1Error_list.append(b)
        analytical_l2SmallError_list.append(c)
        analytical_RMSE_list.append(d)
        analytical_linfError_list.append(e)

        ExactSolutionList = []

        X = np.linspace(0, 1,rawData.shape[1] , endpoint=True)
        
        for k in range(len(ComputedTauSolution)):
            ExactSolutionList.append(ExactSolution(X[k],eps))
        
        ############### GRAPH PLOTTERS #########################
    
        savefilename = folderName + "/" + str(pe_nr) + ".pdf"
        fig = plt.figure()
        plt.title("SUPG Solution using ANN Predicted Tau vs Analytical Tau ")
        plt.xlabel("eps")
        plt.ylabel("solution")
        plt.xlim(0,1)
        # plt.ylim(0,1.02)
        plt.plot(X,rawData[0] , color="blue", linewidth=1.5 ,  linestyle="-", label="Analytical_Tau")
        plt.plot(X,rawData[1] , color="red", linewidth=1.5 , linestyle=":",  label="ANN_Predicted")
        plt.plot(X,rawData[2] , color="purple", linewidth=1.5 , linestyle=":",  label="Galerkin")
        # plt.plot(X,rawData_worst[1] , color="black", linewidth=1.5 , linestyle="-",  label="ANN_Predicted-worst")
        plt.plot(X,ExactSolutionList , color="green", linewidth=1.5 , linestyle="-.",  label="Exact Solution")
        plt.legend(loc='upper left', frameon=False,bbox_to_anchor=(1.05, 1))
        plt.grid()
        
        plt.savefig(savefilename,dpi=400,bbox_inches='tight')
        plt.close(fig)

        savefilename = folderName + "/" + str(pe_nr) + ".html"
        fig = go.Figure()
        
        fig.add_trace(go.Scatter(x=X, y=rawData[0], name='Analytical_Tau',mode='lines+markers',marker_symbol="circle",
                         line=dict(color='firebrick', width=2)))
        fig.add_trace(go.Scatter(x=X, y=rawData[1], name='ANN_Predicted',mode='lines+markers',marker_symbol="square",
                         line=dict(color='blue', width=2)))
        fig.add_trace(go.Scatter(x=X, y=rawData[2], name='Galerkin',mode='lines+markers',marker_symbol="star",
                         line=dict(color='purple', width=2)))
        fig.add_trace(go.Scatter(x=X, y=ExactSolutionList, name='Exact_soln',mode='lines+markers',marker_symbol="diamond-open",
                         line=dict(color='green', width=2)))
        fig.update_traces( showlegend = True)
        fig.update_xaxes(title_text='x',ticks="outside",linecolor='black',mirror=True,showgrid=True, gridwidth=0.7, gridcolor='grey')
        fig.update_yaxes(title_text='Solution',ticks="inside",linecolor='black',mirror=True,showgrid=True, gridwidth=0.7, gridcolor='grey')
        titlestring = "SUPG Solution - ANN Predicted Tau vs Analytical Tau - " +  str("PE_NR ") + str(epsnumbersList[list])
        fig.update_layout(title=titlestring)
        
        fig.write_html(savefilename,auto_open=False)


    ########### -- Write to a File -- #########
    ListName            = folderName
    listExtension       = ".list" 
    savename            = ListName + listExtension
    MainArray           = [ann_l2Error_list,ann_l1Error_list,ann_l2SmallError_list,ann_RMSE_list,ann_linfError_list]
    MainArrayAnalytical = [analytical_l2Error_list,analytical_l1Error_list,analytical_l2SmallError_list,analytical_RMSE_list,analytical_linfError_list]

    np.savetxt(savename, 
            MainArray,
            delimiter =",", 
            fmt ='%10.9f')


    np.savetxt('AnalyticalTau', 
            MainArrayAnalytical,
            delimiter =",", 
            fmt ='%10.9f')