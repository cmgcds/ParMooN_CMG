import os
from multiprocessing import Process
from multiprocessing import Pool
from multiprocessing import cpu_count
import dataFunctions as DF




def runExperiment(ExptName):
    curDir = os.getcwd();
    os.chdir(ExptName);
    import RunExptAll as Rn
    exeName = ExptName+".exe";
    Rn.run(ExptName, exeName);
    os.chdir(curDir);
    pass;


def runAllExperiments(totalNumberOfExpts, numberOfProcesses):

    '''
    # Process approach
    for i in range(totalNumberOfExpts):
        exptName = "Expt"+str(i);
        p = Process(target=runExperiment, args=(exptName,))
        p.start();
        pass;
    '''


    # Pool approach
    exptList = []
    for i in range(totalNumberOfExpts):
        exptList.append("Expt"+str(int(i)));
        pass;

    p = Pool(numberOfProcesses);
    p.map(runExperiment, exptList);
    p.close();
    p.join();

    pass;
    





def preprocessingStep(totalNumberOfExpts):
    # Preprocessing step:
    # Create the directories and copy the scripts and dataset
    for i in range(totalNumberOfExpts):
        curDir = os.getcwd();

        name = "Expt"+str(i);
        exeName = name + ".exe";

        # Create a new directory if not existing
        if not os.path.exists(curDir+'/'+name):
            os.makedirs(name);
            pass;

        # Copy the executable file and rename it 
        os.system("cp parmoon_2D_SEQUENTIAL.exe "+name+"/"+exeName);

        # Copy the library
        os.system("cp -r lib "+name+"/.");

        # Copy the dataset
        os.system("cp allData.csv "+name+"/.");

        # Copy the python scripts
        os.system("cp *.py "+name+"/.");

        pass;
    pass;


def postProcessingStep(totalNumberOfExpts):

    # Create a list of the experiments
    exptList = []
    for i in range(totalNumberOfExpts):
        exptList.append("Expt"+str(int(i)));
        pass;

    # Create the output directory
    curDir = os.getcwd();
    outputDir = curDir+"/output";
    if not os.path.exists(outputDir):
        os.makedirs(outputDir);
        pass;


    # Move all the data files
    for i in range(totalNumberOfExpts):
        exptName = exptList[i];
        os.system("mv "+exptName+"/output/"+exptName+" output/.");
        pass;

    # Get average directory
    DF.getAveragePerformance(exptList, "Avg");

    pass;



if __name__=="__main__":

    numberOfExperiments = 2;

    # Find out the CPU count
    numberOfProcesses = cpu_count();

    preprocessingStep(numberOfExperiments);

    runAllExperiments(numberOfExperiments, numberOfProcesses);

    postProcessingStep(numberOfExperiments);
    pass;
