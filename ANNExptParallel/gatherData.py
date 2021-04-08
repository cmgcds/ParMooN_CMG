from dataFunctions import *


if __name__=="__main__":
    # Create the list of the directories

    numberOfExperiments = 44;

    avgDirName = "Avg";

    dirList = [];

    for i in range(numberOfExperiments):
        dirList.append("Expt"+str(i));
        pass;


    # Delete the existing Avg folder if present
    curDir = os.getcwd();

    avgDirPath = curDir + '/output/' + avgDirName;

    # If average directory exists inside /output, it will be deleted and a new one will be created.
    if os.path.exists(avgDirPath):
        print("Deleting existing Avg directory...");
        os.system('rm -rf '+avgDirPath);
        pass;

    # Following function is defined in dataFunctions.py
    getAveragePerformance(dirList, avgDirName);

    pass;




    

