# -*- coding: utf-8 -*-
"""
Project: MCFR_HyS_BOP IES NEUP
Author: Nicholas Dunkle
Advisor: Dr. Sandra Bogetic
PI: Dr. Nicholas Brown

The script performs following for a given numberOfRuns
1. create a folder 
2. create .mos file to load model and simulate 
3. run the simulation with timing using bash time and omc 
4. extract times using regex and split str list to make a number list using map
5. print time results to a .m file for easy plotting using matlab 

Steps:
- "python3 runFreq.py"
- control-Z
- bg
- disown -h
- Then wait for it to complete
- Then "python3 findPhase.py"

"""

# import packages
import re
import numpy as np
import os
import subprocess
import random
from shutil import copyfile
from scipy import interpolate

simBasePath = os.getcwd()

numberOfRuns = 150
rndRangeLow = 1000
rndRangeHigh = 1400

freq_space = np.logspace(
    -4, 2, num=numberOfRuns
)  # must be the same in both runFreq.py and findPhase.py

run = []
FreqInput = []
FreqInputVal = []
# HTCcoefCoreVal = []
# HTCcoefPHXVal = []

for runNumber in range(1, numberOfRuns + 1):

    NewPower = random.randint(rndRangeLow, rndRangeHigh)

    RePower = str(NewPower)
    Omega = freq_space[runNumber - 1]  # " + "{:.5f}".format(i) + "

    # Make dir
    workPath = "../run" + "{:.0f}".format(runNumber)  # "../f%f" % Omega #
    os.mkdir(workPath)

    with open("IES.mo", "r") as file:
        filedata = file.read()
        filedata = filedata.replace("%%freq%%", format(Omega))

    with open(workPath + "/IES.mo", "w") as batch:
        batch.write(filedata)
        batch.close()

    runtime = float(1400 + 4 / (Omega / (2 * np.pi)))
    numintervals = 100000  # float(runtime*2*(Omega/(2*np.pi)))
    # Write depl script
    runSim = (
        '//MCFR_HyS_BOP \n\
    setModelicaPath("/opt/dymola-2021x-x86_64/Modelica/Library/"); \n\
    getErrorString(); \n\
    loadModel(Modelica); getErrorString(); \n\
    loadFile("MCFR_HyS_BOP.mo"); \n\
    loadFile("IES.mo"); \n\
    simulate(IES,startTime=0,stopTime='
        + "{:.0f}".format(runtime)
        + ",numberOfIntervals="
        + "{:.0f}".format(numintervals)
        + ',tolerance=1E-8,method=dassl,outputFormat = "csv",fileNamePrefix ="run'
        + "{:.0f}".format(runNumber)
        + '",variableFilter="time|Core.NomPower|RBOP.Q_turbine|HyS.mdot_H2_out"); \n\
    getErrorString(); \n'
    )

    # Save depl script in work dir
    runSimFileName = "{}/runMCFR.mos".format(workPath)
    deplFile = open(os.path.join(workPath, runSimFileName), mode="w+")
    deplFile.write(runSim)
    deplFile.close()

    # Copy MCFR_HyS_BOP Library
    modelicaLibrary = "{}/MCFR_HyS_BOP.mo".format(workPath)
    copyfile("MCFR_HyS_BOP.mo", modelicaLibrary)

    # Give permission to create an executable
    bashCommand1 = "chmod +x {}".format(workPath)
    os.system(bashCommand1)

    # Execute omc with bash time as a python subprocess at workPath
    simOutput, timeOutput = subprocess.Popen(
        ["time", "-p", "omc", "runMCFR.mos"],
        cwd=workPath,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    ).communicate()

    run.append(runNumber)
    FreqInputVal.append(Omega)
    # HTCcoefCoreVal.append(HTCcoefCore)
    # HTCcoefPHXVal.append(HTCcoefPHX)

# Write timeResults as a .m file in proper format to plot using matlab
with open("senResultsX.m", "w") as text_file:
    print("Run = {}".format(run) + ";", file=text_file)
    print("Frequency = {}".format(FreqInputVal) + ";", file=text_file)
    # print("HTCcoefCore = {}".format(HTCcoefCoreVal)+";", file=text_file)
    # print("HTCcoefPHX = {}".format(HTCcoefPHXVal)+";", file=text_file)
