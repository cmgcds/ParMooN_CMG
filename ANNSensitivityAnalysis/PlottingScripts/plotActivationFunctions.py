import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

def lrelu(x):
    y = np.zeros_like(x);
    for i in range(np.size(x)):
        y[i] = max(0,x[i]) - 0.1*max(0,-x[i]);
        pass;
    return y

def sigmoid(x):
    return 1./(1. + np.exp(-10.*x));

def tanH(x):
    return (np.exp(x) - np.exp(-x))/ (np.exp(x) + np.exp(-x));

def plotAF():
    x = np.linspace(-1,1,200);

    fig, axs = plt.subplots(1,3, figsize=(6.4,4.8), dpi=300, constrained_layout=True);
    
    # Sigmoid function
    location = (0);
    y = sigmoid(x);
    axs[location].set_title("(a)Sigmoid Function", fontsize=10,y=0, pad=-40);
    axs[location].plot(x,y,color='red');
    axs[location].set_xlim(-1,1);
    axs[location].set_xticks([-1,-0.5,0,0.5,1]);
    axs[location].set_xticklabels([r"$-1$",r"$-0.5$",r"$0$",r"$0.5$",r"$1$"]);
    axs[location].set_ylim(-1,1);
    axs[location].set_yticks([-1,-0.5,0,0.5,1]);
    axs[location].set_yticklabels([r"$-1$",r"$-0.5$",r"$0$",r"$0.5$",r"$1$"]);
    axs[location].set_xlabel(r"$x$");
    axs[location].set_ylabel(r"$f(x)$", labelpad=-5);
    axs[location].axhline(0, color='black', linewidth = 1);
    axs[location].axvline(0, color='black', linewidth = 1);
    axs[location].grid();

    # Leaky-ReLU function
    y = lrelu(x);
    location = (1);
    axs[location].set_title("(b)Leaky-ReLU Function", fontsize=10, y=0, pad=-40);
    axs[location].plot(x,y,color='red');
    axs[location].set_xlim(-1,1);
    axs[location].set_xticks([-1,-0.5,0,0.5,1]);
    axs[location].set_xticklabels([r"$-1$",r"$-0.5$",r"$0$",r"$0.5$",r"$1$"]);
    axs[location].set_ylim(-1,1);
    axs[location].set_yticks([-1,-0.5,0,0.5,1]);
    axs[location].set_yticklabels([r"$-1$",r"$-0.5$",r"$0$",r"$0.5$",r"$1$"]);
    axs[location].set_xlabel(r"$x$");
    axs[location].set_ylabel(r"$f(x)$", labelpad=-5);
    axs[location].axhline(0, color='black', linewidth = 1);
    axs[location].axvline(0, color='black', linewidth = 1);
    axs[location].grid();

    # Hyperbolic Tangent function
    location = (2);
    y = tanH(5.*x);
    axs[location].set_title("(c)TanH Function", fontsize=10, y=0, pad=-40);
    axs[location].plot(x,y,color='red');
    axs[location].set_xlim(-1,1);
    axs[location].set_xticks([-1,-0.5,0,0.5,1]);
    axs[location].set_xticklabels([r"$-1$",r"$-0.5$",r"$0$",r"$0.5$",r"$1$"]);
    axs[location].set_ylim(-1,1);
    axs[location].set_yticks([-1,-0.5,0,0.5,1]);
    axs[location].set_yticklabels([r"$-1$",r"$-0.5$",r"$0$",r"$0.5$",r"$1$"]);
    axs[location].set_xlabel(r"$x$");
    axs[location].set_ylabel(r"$f(x)$", labelpad=-5);
    axs[location].axhline(0, color='black', linewidth = 1);
    axs[location].axvline(0, color='black', linewidth = 1);
    axs[location].grid();

    fig.tight_layout(rect=[0.02,0,0.98,0.45])
    plt.savefig("ActivationFunctions.pdf");

    os.system("pdfcrop ActivationFunctions.pdf ActivationFunctions.pdf");
    pass;


if __name__=="__main__":
    plt.rcParams["text.usetex"] = True
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = ["Computer Modern"]
    plotAF();
    pass;
