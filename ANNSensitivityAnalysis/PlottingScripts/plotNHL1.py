import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Name of the project
projectName = "ANN";

# Enter the run number i.e. folder output/ANN/runNUmber
runNumber = 1;

#_______________________________________________________
# Set paths and variables
#_______________________________________________________
# Output path for sensitivity analysis (./output)
outputPath = os.getcwd()+'/output';
# Location for storing test case (name is specified in inputData.py)
projectOutputDir = outputPath+'/'+projectName; 
runDir = projectOutputDir+'/'+str(runNumber); 

# Change into the run directory inside output/caseName
os.chdir(runDir);


#_______________________________________________________
# Prepare Sample data
#_______________________________________________________

# Load input space for sensitivity analysis
inputData = np.loadtxt("inputSpace.dat");
# Total number of samples
numberOfSamples = 72;
inputData = inputData[:numberOfSamples];
# Total number of input parameters i.e. dimension of Input space for sensitivity analysis
numberOfInputParam = 4;
# Total number of output results (4) 
# i.e. min error, max error, L1 error, L2 error
numberOfOutputParam = 4;

# Extract useful information about arrays etc. Refer Run.py to find out the sequence of data
sampleNumber = inputData[:,0]; # Will NOT be used
NHL = inputData[:,1]; # Will NOT be used
OPLTYPE = inputData[:,2];
HL_0_DIM = inputData[:,3];
HL_0_TYPE = inputData[:,4];
HL_1_DIM = inputData[:,5]; # Will NOT be used
HL_1_TYPE = inputData[:,6]; # Will NOT be used
HL_2_DIM = inputData[:,7]; # Will NOT be used
HL_2_TYPE = inputData[:,8]; # Will NOT be used
EPOCHS = inputData[:,9];

# Load the output space for sensitivity analysis
outputData = np.loadtxt("outputSpace_NHL1.dat");
L1Error = outputData[:,0];
L2Error = outputData[:,1];
MinError = outputData[:,2];
MaxError = outputData[:,3];

# Find out the NHL1, NHL2 and NHL3 data
flag1 = 0;
flag2 = 0;
flag3 = 0;
for i in range(numberOfSamples):
    if (inputData[i,2] == 0):
        flag1 = i+1;
    elif (inputData[i,2] == 1):
        flag2 = i+1;
    elif (inputData[i,2] == 3):
        flag3 = i+1;

#_______________________________________________________
# Process the data 
#_______________________________________________________
p5 = np.percentile(outputData,5, axis=0);
p5={'L1Error':p5[0], 'L2Error':p5[1], 'MinError':p5[2], 'MaxError':p5[3]};

p95 = np.percentile(outputData,95, axis=0);
p95={'L1Error':p95[0], 'L2Error':p95[1], 'MinError':p95[2], 'MaxError':p95[3]};

#_______________________________________________________
# Set color codes
#_______________________________________________________


shadeColor = 'beige';
s1color = 'red';
s2color = 'green';
s3color = 'darkblue';
s4color = 'darkorange';

s1size = '2';
s2size = '1';
s3size = '1.5';
s4size = '1.5';

#_______________________________________________________
# Plot the errors for all the networks (design space)
#_______________________________________________________

fig, axs = plt.subplots(2,2, figsize=(6.4,4.8), dpi=300, constrained_layout=True);

#_______________________________________________________
# L1 Error plot, location 0,0
#_______________________________________________________
location = (0,0)
axs[location].plot(np.arange(0,flag1), L1Error[0:flag1],"^", color=s1color, markersize=s1size, label="Sigmoid");
axs[location].plot(np.arange(flag1,flag2), L1Error[flag1:flag2],"s",color=s2color, markersize=s2size, label="Identity");
axs[location].plot(np.arange(flag2,flag3), L1Error[flag2:flag3],"v",color=s3color, markersize=s3size, label="Leaky-ReLU");
axs[location].plot(np.arange(flag3, numberOfSamples), L1Error[flag3:numberOfSamples],"o", color=s4color, markersize=s4size, label="TanH");

axs[location].legend(loc=0,ncol=3, fontsize=6.6);

axs[location].set_ylabel(r"$L_1$ Error");
axs[location].set_ylim(0.01,0.04);
axs[location].set_yticks([0.01,0.02, 0.03, 0.04]);

# highlight 5 and 95 percentile zones with demarkation line
axs[location].axhspan(0.01,p5['L1Error'],facecolor=shadeColor)
axs[location].axhspan(p95['L1Error'],0.04,facecolor=shadeColor)
axs[location].axhline(p5['L1Error'],color=s4color,linewidth='0.3');
axs[location].axhline(p95['L1Error'],color=s4color,linewidth='0.3');

#_______________________________________________________
# L2 Error plot
#_______________________________________________________
location = (1,0)
axs[location].plot(np.arange(0,flag1), L2Error[0:flag1],"^", color=s1color, markersize=s1size, label="Sigmoid");
axs[location].plot(np.arange(flag1,flag2), L2Error[flag1:flag2],"s",color=s2color, markersize=s2size, label="Identity");
axs[location].plot(np.arange(flag2,flag3), L2Error[flag2:flag3],"v",color=s3color, markersize=s3size, label="Leaky-ReLU");
axs[location].plot(np.arange(flag3, numberOfSamples), L2Error[flag3:numberOfSamples],"o", color=s4color, markersize=s4size, label="TanH");


axs[location].legend(loc=0,ncol=3, fontsize=6.6);

axs[location].set_xlabel(r"ANN Samples");

axs[location].set_ylabel(r"$L_2$ Error");
axs[location].set_ylim([0.3,1]);
axs[location].set_yticks([0.4,0.6,0.8,1]);

# highlight 5 and 95 percentile zones with demarkation line
axs[location].axhspan(0.3,p5['L2Error'],facecolor=shadeColor)
axs[location].axhspan(p95['L2Error'],1,facecolor=shadeColor)
axs[location].axhline(p5['L2Error'],color=s4color,linewidth='0.3');
axs[location].axhline(p95['L2Error'],color=s4color,linewidth='0.3');


#_______________________________________________________
# L_Min Error plot
#_______________________________________________________
location = (0,1);
axs[location].semilogy(np.arange(0,flag1), MinError[0:flag1],"^", color=s1color, markersize=s1size, label="Sigmoid");
axs[location].semilogy(np.arange(flag1,flag2), MinError[flag1:flag2],"s",color=s2color, markersize=s2size, label="Identity");
axs[location].semilogy(np.arange(flag2,flag3), MinError[flag2:flag3],"v",color=s3color, markersize=s3size, label="Leaky-ReLU");
axs[location].semilogy(np.arange(flag3, numberOfSamples), MinError[flag3:numberOfSamples],"o", color=s4color, markersize=s4size, label="TanH");


axs[location].legend(loc=1,ncol=3, fontsize=6.6);


axs[location].set_ylabel(r"Min Error");
axs[location].set_ylim([1e-8,1e-1]);
axs[location].set_yticks([1e-8,1e-6,1e-4,1e-2]);

# highlight 5 and 95 percentile zones with demarkation line
axs[location].axhspan(1e-8,p5['MinError'],facecolor=shadeColor)
axs[location].axhspan(p95['MinError'],1e-1,facecolor=shadeColor)
axs[location].axhline(p5['MinError'],color=s4color,linewidth='0.3');
axs[location].axhline(p95['MinError'],color=s4color,linewidth='0.3');


#_______________________________________________________
# L_Max Error plot
#_______________________________________________________
location = (1,1);
axs[location].plot(np.arange(0,flag1), MaxError[0:flag1],"^", color=s1color, markersize=s1size, label="Sigmoid");
axs[location].plot(np.arange(flag1,flag2), MaxError[flag1:flag2],"s",color=s2color, markersize=s2size, label="Identity");
axs[location].plot(np.arange(flag2,flag3), MaxError[flag2:flag3],"v",color=s3color, markersize=s3size, label="Leaky-ReLU");
axs[location].plot(np.arange(flag3, numberOfSamples), MaxError[flag3:numberOfSamples],"o", color=s4color, markersize=s4size, label="TanH");


axs[location].legend(loc=1,ncol=3, fontsize=6.6);

axs[location].set_xlabel(r"ANN Samples");

axs[location].set_ylabel(r"Max Error");
axs[location].set_ylim([1e-1,0.8]);
axs[location].set_yticks([2e-1,0.4,0.6]);

# highlight 5 and 95 percentile zones with demarkation line
axs[location].axhspan(1e-1,p5['MaxError'],facecolor=shadeColor)
axs[location].axhspan(p95['MaxError'],1e0,facecolor=shadeColor)
axs[location].axhline(p5['MaxError'],color=s4color,linewidth='0.3');
axs[location].axhline(p95['MaxError'],color=s4color,linewidth='0.3');

plt.savefig("ErrorD2NH1.pdf");




