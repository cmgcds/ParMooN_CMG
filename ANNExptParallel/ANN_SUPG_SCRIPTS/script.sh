#!/bin/bash


#############################################################################
#######  Author  : Thivin Anandh D                              ############
#######  Purpose : To run each ANN model 'n' times for each     ############
#######            Epsilon Values                               ############
##############################################################################

####### -- USER INPUTS ------- ###


N_iterations=30
# Declare the ArrayNames for the folders
declare -a foldername=("BestTDS100" "BestTDS200" "BestTDS400" "WorstTDS100" "WorstTDS200" "WorstTDS400")
declare -a datfilename=( "BestMSEScriptTDS100.dat" "BestMSEScriptTDS200.dat" "BestMSEScriptTDS400.dat" "WorstMSEScriptTDS100.dat" "WorstMSEScriptTDS200.dat" "WorstMSEScriptTDS400.dat")
epsnumbersList=(100000 1000000 10000000 100000000 1000000000 10000000000 100000000000 1000000000000)

WorkingDir="cd1dANN_1"


#### export userinputs
export N_iterations





# get length of an array for dat file and other parameters
arraylength=${#foldername[@]}
echo ${arraylength}

#get length of array for epsvalueslist
epsLength=${#epsnumbersList[@]}
echo $epsLength



## get he Confirmation of all the Workign Directory
echo "[Question : ] Working Directory : ${WorkingDir}"
read -p "Is the Working Directory Correct ?" yn
    case $yn in
        [Yy]* ) echo .;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no.";;
    esac


## get he Confirmation of all the parameters in the file
echo "${epsnumbersList[@]}"
read -p "[Question : ]  Is the Peclet number values Correct ?" yn
    case $yn in
        [Yy]* ) echo .;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no.";;
    esac


## ---------------- get confirmation for N_Elements --------------------- ##
for i in `ls *.dat`
do
    echo " File : $i" 
    cat $i | grep "N_ELEM" 
done

read -p "[Question : ]  Is the H value Correct in all files ?" yn
    case $yn in
        [Yy]* ) echo .;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no.";;
    esac

## ---------------- get confirmation for Disctype --------------------- ##
for i in `ls *.dat`
do
    echo " File : $i" 
    cat $i | grep  "DISCTYPE" 
done

read -p "[Question : ] Is the DiscType value Correct in all files ?" yn
    case $yn in
        [Yy]* ) echo .;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no.";;
    esac


## ---------------- get confirmation for Exact Solution --------------------- ##
python3 -c 'from SetParameters import printExactSolution; printExactSolution()'

read -p "[Question : ] Is the Exact Solution value Correct in 'SetParameters.py' file?" yn
    case $yn in
        [Yy]* ) echo .;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no.";;
    esac



# FileName

for ((iter=0; iter<${N_iterations};iter++))
do
    for((i=0; i<${arraylength} ; i++));
    do
        modelName="model${foldername[$i]}.txt"
        datfilename=${datfilename[$i]}
        echo "Index : $i , FolderName : ${foldername[$i]}, fileName : ${filename}"

        mkdir -p ${foldername[$i]}

        ## Enter the Epsilon loop here
        for((j=0; j<${epsLength} ; j++));
        do
            eps=$((${epsnumbersList[j]}))
            echo " ----------------     Eps value : $eps  --------------------------------- "
            sed -i "s/PE_NR: .*/PE_NR: ${epsnumbersList[j]}/g" ${datfilename}        ## Change the pe_nr in the .dat file
            echo "j val ${epsnumbersList[j]} "
            ./parmoon_2D_SEQUENTIAL.exe ${datfilename}  ${modelName}

            ## transfer and rename the files accordingly
            csvname="supgSol_${epsnumbersList[j]}.csv"
            tauname="supgtau_${epsnumbersList[j]}.csv"
            newcsvname="supgSol_${epsnumbersList[j]}_${iter}.csv"
            newtauname="supgtau_${epsnumbersList[j]}_${iter}.csv"

            mv "${csvname}" "${newcsvname}" 
            mv "${tauname}" "${newtauname}" 
            mv "${newcsvname}"  "${foldername[$i]}"
            mv "${newtauname}"  "${foldername[$i]}"
        done

        ### -- rename the .list file with the suffix ----
        python3 mainPlot.py ${foldername[$i]} ${iter}
        mv "${foldername[$i]}.list" "${foldername[$i]}_${iter}.list"

    done
        # python3 mainPlot.py ${foldername[$i]}

    ########## ---- Change Training Data ---- ##################
    cd ../cd1dANN_ann
    python3 datagen.py
    cp *.txt ../${WorkingDir}/
    cd ../${WorkingDir}
done



echo " ----- RUN TIME COMPLETED -----------------"


##mkdir plots
mkdir -p plots

## Run the Average Error Script
python3 averageError.py
echo "[MESSAGE] :  Average Error Computed "

## Run the OverallError Plots
python3 OverallErrorPLot.py
echo "[MESSAGE] :  Overall Error pLotted "

## Run the EpsValue Plots
python3 epsPlots.py
echo "[MESSAGE] : EpsPlots Error plotted"


## Run the tau Value plots
python3 tauPlots.py
echo "[MESSAGE] : tau Error plotted"