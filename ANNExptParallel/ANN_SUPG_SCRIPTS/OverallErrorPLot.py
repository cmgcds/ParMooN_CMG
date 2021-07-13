
import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
from itertools import product

from SetParameters import *

##### --- Plotly Scripts for Adding the minor Grids ----

def matplotter_zoom(analyticalError,Best100Data,Best200Data,Best400Data,Worst100Data,Worst200Data,Worst400Data,epsValueList,index,XX,title):
    ## Declare global Variables


    plottitle = title + " Error Convergence - Analytical vs ANN Models"
    savefilename = "plots/" + title + "_error_zoom" + ".pdf"
    fig = plt.figure()
    plt.title(plottitle)
    plt.xlabel("eps")
    plt.ylabel("Error")
    plt.plot(XX,analyticalError, color="red", linewidth=1.5 , marker="o", linestyle="-", label="Analytical_Tau")
    plt.plot(XX,Best100Data[index] , color="yellow", linewidth=1.5 , marker="v", linestyle=":",  label="Best100")
    plt.plot(XX,Best200Data[index] , color="brown", linewidth=1.5 , marker="^", linestyle=":",  label="Best200")
    plt.plot(XX,Best400Data[index] , color="black", linewidth=1.5 , marker="x", linestyle=":",  label="Best400")
    plt.legend(loc='upper left', frameon=False,bbox_to_anchor=(1.05, 1))
    plt.xticks(XX, epsValueList)
    plt.grid()
    plt.savefig(savefilename,dpi=400,bbox_inches='tight')
    plt.close(fig)


### Define the matplot lib ploting function
def matplotter(analyticalError,Best100Data,Best200Data,Best400Data,Worst100Data,Worst200Data,Worst400Data,epsValueList,index,XX,title):
   
    plottitle = title + " Error Convergence - Analytical vs ANN Models"
    savefilename = "plots/" + title + "_error" + ".pdf"
    fig = plt.figure()
    plt.title(plottitle)
    plt.xlabel("eps")
    plt.ylabel("Error")
    plt.plot(XX,analyticalError, color="red", linewidth=1.5 , marker="o", linestyle="-", label="Analytical_Tau")
    plt.plot(XX,Best100Data[index] , color="yellow", linewidth=1.5 , marker="v", linestyle=":",  label="Best100")
    plt.plot(XX,Best200Data[index] , color="brown", linewidth=1.5 , marker="^", linestyle=":",  label="Best200")
    plt.plot(XX,Best400Data[index] , color="black", linewidth=1.5 , marker="x", linestyle=":",  label="Best400")
    plt.plot(XX,Worst400Data[index] , color="purple", linewidth=1.5 , marker="d", linestyle=":",  label="Worst400")
    plt.plot(XX,Worst200Data[index] , color="cyan", linewidth=1.5 , marker="D", linestyle=":",  label="Worst200")
    plt.plot(XX,Worst100Data[index] , color="pink", linewidth=1.5 , marker="x", linestyle=":",  label="Worst100")
    plt.legend(loc='upper left', frameon=False,bbox_to_anchor=(1.05, 1))
    plt.xticks(XX, epsValueList)
    plt.grid()
    plt.savefig(savefilename,dpi=400,bbox_inches='tight')
    plt.close(fig)


def plotly_plot(analyticalError,Best100Data,Best200Data,Best400Data,Worst100Data,Worst200Data,Worst400Data,epsValueList,index,XX,title):
    ## Declare global Variables


    savefilename = savefilename = "plots/" + title  + ".html"
    plottitle = title + " Error Convergence - Analytical vs ANN Models"
    fig = go.Figure()

    fig.add_trace(go.Scatter(x=epsValueList, y=analyticalError, name='Analytical_tau_Error',mode='lines+markers',marker_symbol="circle",
                            line=dict(color='firebrick', width=3)))
    fig.add_trace(go.Scatter(x=epsValueList, y=Best100Data[index], name='Best100',mode='lines+markers',marker_symbol="star",
                            line=dict(color='magenta', width=3)))
    fig.add_trace(go.Scatter(x=epsValueList, y=Best200Data[index], name='Best200',mode='lines+markers',marker_symbol="triangle-up",
                            line=dict(color='goldenrod', width=3)))
    fig.add_trace(go.Scatter(x=epsValueList, y=Best400Data[index], name='Best400',mode='lines+markers',marker_symbol="triangle-down",
                            line=dict(color='cyan', width=3)))
    fig.add_trace(go.Scatter(x=epsValueList, y=Worst100Data[index], name='Worst100',mode='lines+markers',marker_symbol="octagon",
                            line=dict(color='yellow', width=3)))
    fig.add_trace(go.Scatter(x=epsValueList, y=Worst200Data[index], name='Worst200',mode='lines+markers',marker_symbol="hexagram",
                            line=dict(color='violet', width=3)))
    fig.add_trace(go.Scatter(x=epsValueList, y=Worst400Data[index], name='Worst400',mode='lines+markers',marker_symbol="pentagon",
                            line=dict(color='purple', width=3)))
    fig.update_traces( showlegend = True , marker=dict(size=14))
    fig.update_xaxes(title_text='Epsilon Values',ticks="outside",ticktext=epsValueList,linecolor='black',mirror=True,showgrid=True, gridwidth=0.7, gridcolor='grey')
    fig.update_yaxes(title_text='Error',ticks="inside",linecolor='black',mirror=True,showgrid=True, gridwidth=0.7, gridcolor='grey')
    titlestring = plottitle
    fig.update_layout(title=titlestring,plot_bgcolor='rgb(250,250,250)')
    fig.write_html(savefilename,auto_open=False)
    

if __name__ == "__main__":
    range_end = 12
    range_start = 5

    


    ##########  - READ FROM CSV FILES -- 
    Best100Data = genfromtxt('BestTDS100.list',delimiter=',')
    Best200Data = genfromtxt('BestTDS200.list',delimiter=',')
    Best400Data = genfromtxt('BestTDS400.list',delimiter=',')
    Worst100Data = genfromtxt('WorstTDS100.list',delimiter=',')
    Worst200Data = genfromtxt('WorstTDS200.list',delimiter=',')
    Worst400Data = genfromtxt('WorstTDS400.list',delimiter=',')
    analytical_Tau = genfromtxt('AnalyticalTau',delimiter=',')

    

    analytical_l2Error_list =   analytical_Tau[0]
    analytical_l1Error_list  =  analytical_Tau[1]
    analytical_linfError_list = analytical_Tau[4]
    analytical_rmseError_list = analytical_Tau[3]
    analytical_smallError_list = analytical_Tau[2]

    epsValueList = []

    for list in range(len(epsnumbersList)):
        pe_nr = int(epsnumbersList[list])
        eps   = 1.0/pe_nr
        epsValueList.append(str(eps))


    XX = np.linspace(0,1,len(epsnumbersList)) 
    

    ## --- Call the plotting Scripts for PDF Generation --- ##
    matplotter(analytical_l2Error_list,Best100Data,Best200Data,Best400Data,Worst100Data,Worst200Data,Worst400Data,epsValueList,0,XX,"L2")
    matplotter(analytical_l1Error_list ,Best100Data,Best200Data,Best400Data,Worst100Data,Worst200Data,Worst400Data,epsValueList,1,XX,"L1")
    matplotter(analytical_linfError_list ,Best100Data,Best200Data,Best400Data,Worst100Data,Worst200Data,Worst400Data,epsValueList,4,XX,"Linf")
    matplotter(analytical_rmseError_list ,Best100Data,Best200Data,Best400Data,Worst100Data,Worst200Data,Worst400Data,epsValueList,3,XX,"RMSE")
    matplotter(analytical_smallError_list ,Best100Data,Best200Data,Best400Data,Worst100Data,Worst200Data,Worst400Data,epsValueList,2,XX,"l2small")


    ## --- Call the Plotting function for zoomed in PDF generation --#
    matplotter_zoom(analytical_l2Error_list ,Best100Data,Best200Data,Best400Data,Worst100Data,Worst200Data,Worst400Data,epsValueList,0,XX,"L2")
    matplotter_zoom(analytical_l1Error_list ,Best100Data,Best200Data,Best400Data,Worst100Data,Worst200Data,Worst400Data,epsValueList,1,XX,"L1")
    matplotter_zoom(analytical_linfError_list ,Best100Data,Best200Data,Best400Data,Worst100Data,Worst200Data,Worst400Data,epsValueList,4,XX,"Linf")
    matplotter_zoom(analytical_rmseError_list ,Best100Data,Best200Data,Best400Data,Worst100Data,Worst200Data,Worst400Data,epsValueList,3,XX,"RMSE")
    matplotter_zoom(analytical_smallError_list ,Best100Data,Best200Data,Best400Data,Worst100Data,Worst200Data,Worst400Data,epsValueList,2,XX,"l2small")


    ###--- Call the PLOTLY FUNCTION to generate HTML Plots ---- #####
    plotly_plot(analytical_l2Error_list ,Best100Data,Best200Data,Best400Data,Worst100Data,Worst200Data,Worst400Data,epsValueList,0,XX,"L2")
    plotly_plot(analytical_l1Error_list ,Best100Data,Best200Data,Best400Data,Worst100Data,Worst200Data,Worst400Data,epsValueList,1,XX,"L1")
    plotly_plot(analytical_linfError_list ,Best100Data,Best200Data,Best400Data,Worst100Data,Worst200Data,Worst400Data,epsValueList,4,XX,"Linf")
    plotly_plot(analytical_rmseError_list ,Best100Data,Best200Data,Best400Data,Worst100Data,Worst200Data,Worst400Data,epsValueList,3,XX,"RMSE")
    plotly_plot(analytical_smallError_list ,Best100Data,Best200Data,Best400Data,Worst100Data,Worst200Data,Worst400Data,epsValueList,2,XX,"l2small")

    # plt.show()




   