~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~ File Directory ~~~~~
MCFR_HyS_BOP.mo
- Explanation: The IES dynamic model

IES.mo
- Explanation: The Simbase class from the full model. It is necessary to have these separate for the 
modelica shell script (.mos) to work. If told to run MCFR_HyS_BOP.mo on its own, the computer will be 
confused and not know to run the simbase.
- NOTE: It is useful to have the parameters you wish to modify in this file. Have the parameter 
clearly marked (example: %%freq%%) so the search only finds what you wish to change.

runFreq.py
- Explanation: This script makes a desired number of directories for each run, takes the model files, 
copies them to each directory with the desired modification, writes runMCFR.mos, and runs that script 
to run the model. The desired outputs are recorded in a .csv file in that run's folder.

findPhase.py
- Explanation: This script goes into each of the directories, takes the relevant data, and fits it to 
a sinusoid function for the frequency analysis.

findThermal.py
- Explanation: Same general idea as to the findPhase.py but for the thermal power output.

findHydro.py
- Explanation: Same general idea as to the findPhase.py but for the hydrogen production output.

findElectro.py
- Explanation: Same general idea as to the findPhase.py but for the electrical power output.

~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~      Misc	    ~~~~~
- IMPORTANT! Have these files one directory below where you wish for the run directories to be. 
Those directories will be made one above.
- Limit the number of variables recorded to the .csv files. If nondefined, it will record all of them 
and significantly bog down processing time and bloat storage sizes.
- IMPORTANT! Wait for the runs to finish before using any of the findX.py scripts. They cannot find 
data that hasn't been written yet.

~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~    Commands  ~~~~~
If not using a queue submission system, multiple sets of runs can be processed concurrently via this method--
Steps:
- "python3 runFreq.py"
- control-Z
- "bg"
- "disown -h"
- Then wait for it to complete
- Then "python3 findPhase.py"

Explanation: When you run the python script, it holds the command terminal. Control-Z stops the script, 
allowing access to the command terminal again. Next, "bg" puts the run into the background, where it 
continues to run. Finally, "disown -h" will allow the script to continue running even when the SSH 
terminal is closed.

~~~~~~~~~~~~~~~~~~~~~~~~~