# %%

############### PLOTS FOR VARIOUS EPSILON VALUES ######################################

import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go

## import epsnumber list
from SetParameters import *

########################### _________ MAIN CODE _____________ #################################
if __name__ == "__main__":
    # N_iterations = 5
    # epsnumbersList = [50,100,200,500,550,600,650,700]

    for index in range(len(epsnumbersList)):
        
        # LISTFOLDERNAMES = ['BestTDS100','BestTDS200','BestTDS400','WorstTDS100','WorstTDS200','WorstTDS400']
        # print (" ===========================================================")
        # print(  " FOR PE NR : " , str(int(epsnumbersList[index]))   )
        # print (" ===========================================================")
        pe_nr = int(epsnumbersList[index])
        eps   = float(1.0/pe_nr)
        
        tauXticks = np.linspace(1,N_iterations,N_iterations)
        

        TauValueList = [[] for n in range(len(LISTFOLDERNAMES) + 1)]
        LastIndex = len(LISTFOLDERNAMES) 

        for list in range(len(LISTFOLDERNAMES)):

            for iter in range(N_iterations):
                filename  = str(LISTFOLDERNAMES[list]) + "/supgtau_" + str(epsnumbersList[index]) + "_" + str(iter) + ".csv"
                k = np.genfromtxt(filename,delimiter=',')
                predTau = k[0]
                analyticalTau = k[1]
                TauValueList[list].append(predTau)
                if(list==0):
                    TauValueList[LastIndex].append(analyticalTau)
            
        
        X = np.linspace(0, 1,N_iterations , endpoint=True)
        ## PLotting Scripts
        # Graph Plotting 
        
        title = "Tau - Various ANN Models - Eps : " + str(eps) + str(index)
        filename1 = "plots/tau" + str(pe_nr) +".pdf"
        fig = plt.figure()
        plt.title(title)
        plt.xlabel("Iteration")
        plt.ylabel("Tau")
        # plt.xlim(0,5)
        # plt.ylim(0,1.02)
        plt.plot(X,TauValueList[0] , color="blue", linewidth=1.0 ,  linestyle="-",marker="." ,markersize="0.75", label=LISTFOLDERNAMES[0])
        plt.plot(X,TauValueList[1] , color="red", linewidth=1.0 , linestyle="-",  marker="," ,markersize="0.75", label=LISTFOLDERNAMES[1])
        plt.plot(X,TauValueList[2] , color="green", linewidth=1.0 , linestyle="-", marker="o" ,markersize="0.75",  label=LISTFOLDERNAMES[2])
        plt.plot(X,TauValueList[3] , color="purple", linewidth=1.0 , linestyle="-",  marker="v" ,markersize="0.75", label=LISTFOLDERNAMES[3])
        plt.plot(X,TauValueList[4] , color="yellow", linewidth=1.0 , linestyle="-",marker="^" ,markersize="0.75",   label=LISTFOLDERNAMES[4])
        plt.plot(X,TauValueList[5] , color="grey", linewidth=1.0 , linestyle="-", marker="<" ,markersize="0.75",  label=LISTFOLDERNAMES[5])
        plt.plot(X, TauValueList[6], color="black", linewidth=1.0 , linestyle="-", marker="<" ,markersize="0.75",  label="Analytical Tau")
        plt.legend(loc='upper left', frameon=False,bbox_to_anchor=(1.05, 1))
        plt.grid()
        plt.savefig(filename1,dpi=400,bbox_inches='tight')
        # plt.show()
        plt.close(fig)


        savefilename = savefilename = "plots/tau" + str(pe_nr)   + ".html"
        plottitle = "Sol - Various ANN Models - Eps : " + str(eps) 
        fig = go.Figure()

        fig.add_trace(go.Scatter(x=X, y=TauValueList[0], name=LISTFOLDERNAMES[0],mode='lines+markers',marker_symbol="circle",
                                line=dict(color='firebrick', width=3)))
        fig.add_trace(go.Scatter(x=X, y=TauValueList[1], name=LISTFOLDERNAMES[1],mode='lines+markers',marker_symbol="triangle-right",
                                line=dict(color='blue', width=3)))
        fig.add_trace(go.Scatter(x=X, y=TauValueList[2], name=LISTFOLDERNAMES[2],mode='lines+markers',marker_symbol="triangle-left",
                                line=dict(color='green', width=3)))
        fig.add_trace(go.Scatter(x=X, y=TauValueList[3], name=LISTFOLDERNAMES[3],mode='lines+markers',marker_symbol="star",
                                line=dict(color='magenta', width=3)))
        fig.add_trace(go.Scatter(x=X, y=TauValueList[4], name=LISTFOLDERNAMES[4],mode='lines+markers',marker_symbol="triangle-up",
                                line=dict(color='goldenrod', width=3)))
        fig.add_trace(go.Scatter(x=X, y=TauValueList[5], name=LISTFOLDERNAMES[5],mode='lines+markers',marker_symbol="triangle-down",
                                line=dict(color='cyan', width=3)))
        fig.add_trace(go.Scatter(x=X, y=TauValueList[6], name='AnalyticalTau',mode='lines+markers',marker_symbol="star-diamond",
                                line=dict(color='purple', width=3)))
        fig.update_traces( showlegend = True , marker=dict(size=14))
        fig.update_xaxes(title_text='Iteration',ticks="outside",ticktext=tauXticks,linecolor='black',mirror=True,showgrid=True, gridwidth=0.7, gridcolor='grey')
        fig.update_yaxes(title_text='Tau',ticks="inside",linecolor='black',mirror=True,showgrid=True, gridwidth=0.7, gridcolor='grey')
        titlestring = plottitle
        fig.update_layout(title=titlestring,plot_bgcolor='rgb(250,250,250)')
        fig.write_html(savefilename,auto_open=False)
# %%
