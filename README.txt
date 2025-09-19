The dynamic modeling work in this repository is part of a Nuclear Energy University Program (NEUP) project funded by the U.S. DOE. The objective of the project is to explore the application of advanced reactors within an Integrated Energy System (IES). In this context, an IES is the coupling of nuclear power to various heat intensive applications (e.g., electricity generation and hydrogen production).

~~~~ File and Folder Directory Descriptions ~~~~
"individualModels(BETA)"
- Explanation: The three systems (MCFR, HyS, and RBOP) as their own individual models. 
These files were made when developing the models. The complete IES model with the most refinement, 
features, and comments is MCFR_HyS_BOP.mo 

"results.zip"
- Explanation: The output data and plots for different transients and accident scenarios.

"scripts"
- Explanation: The files necessary for bulk-running the model numerous times with a slight modification. 
These were used for the frequency response analysis, which required 450 runs with a variable modified.

"MCFR_HyS_BOP.mo"
- Explanation: This is the completed IES dynamic model of the MCFR powering a Hybrid Sulfur cycle for 
hydrogen production and a Regenerative Rankine cycle for electricity production. Requires the OMEdit 
OpenModelica Connection Editor to run or similar.

"MCFR_IES_heat_diverter.mo"
- Explanation: This is an updated version of the completed IES model. This differs from the version in
"MCFR_HyS_BOP.mo" by the addition of steam tables in the Regenerative Rankine cycle’s steam generator
and by the separation of the tertiary loop into two sub loops. This change allows the heat to go to the
decomposer in the Hybrid Sulfur cycle and to the steam generator in the Regenerative Rankine cycle to
be altered independently of one another.
~~~~~~~~~~~~~~~~~~~~~~~~~


~~~~ Relevant Papers ~~~~
"Dynamic Modeling of a Fast Spectrum Molten Salt Reactor Integrated Energy System" journal article by N. Dunkle et Al (2025)
"Expanded Accident Scenarios and Frequency Characteristics for a Fast Spectrum Molten Salt Reactor in an Integrated Energy System" journal article by N. Dunkle et Al (2025)
"Safety and Safeguards for Molten Salt Reactors in Integrated Energy Systems" dissertation by N. Dunkle (2025)
"A reduced order sulfuric acid decomposition model for a nuclear-powered hybrid sulfur cycle" journal article by CT Callaway (2024)
~~~~~~~~~~~~~~~~~~~~~~~~~

