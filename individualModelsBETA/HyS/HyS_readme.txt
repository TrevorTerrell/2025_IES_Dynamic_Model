===============================================
Hybrid Sulfur Dynamic Model in OpenModelica ===
===============================================
OUTLINE:
- File overview
- Inputs and outputs
- Parameters
- Variables
- Misc

=====================
=== FILE OVERVIEW ===
=====================
For the standalone HyS model, the only file is HyS_Model.mo

Within HyS_Model.mo, there are three packages: Hybrid_Sulfur, Units, and Connectors.

The Hybrid_Sulfur package has the HyS model within it, which has all of the equations and parameters. The Units package provides the units for each parameter and variable. The Connectors package has all of the possible connectors that are used in the full IES model. With the standalone HyS model, no connectors are used.

==========================
=== INPUTS AND OUTPUTS ===
==========================
In the standalone HyS model, there are no actual inputs or outputs via connectors. Instead, the primary relevant inputs are: the flow rate on the hot side of the decomposer (mdot_H_nom), the temperature of the inlet on the hotside of the decomposer (T_H_in), the pressure on the secondary side of the decomposer (Pressure), and if the sulfur trioxide conversion is being benchmarked then the decomposer secondary side outlet temperature (T_vap) can be switched from a variable to a parameter.

Likewise, the relevant outputs are: the total heat energy used by the HyS system (Heat_used), the makeup water necessary to fully complete the electrolyis reaction (mdot_H2O_makeup), the hydrogen produced (mdot_H2_out)*, and the total electrical demand of the electrolysis reaction (P_demand). 

*[NOTE] - The mass flow rate of hydrogen produced may seem small, but is the lightest element after all. Assuming atmospheric pressure, a conversion from mass to volume shows a suitibly sized product of tens of cubic meters per second. 

==================
=== PARAMETERS ===
==================
- HyS_Model.Units.MassFlowRate mdot_H_nom = 5400;// Mass flow rate of hot fluid on primary side
- HyS_Model.Units.MolFlowRate Mdot_SA_nom = 62.309;// Mol flow rate of sulfur entering the decomposer
- HyS_Model.Units.MolFlowRate Mdot_H2O_nom = 37.691;// Mol flow rate of excess water in the HyS cycle
- HyS_Model.Units.SpecificHeat Cp_H = 0.00154912;// Specific heat of the hitec salt on the hot side of the decomposer
- HyS_Model.Units.SpecificHeat Cp_W = 7.5e-4;// Specific heat of the tubing material (4.35e-4 for Inconel)
- HyS_Model.Units.Mass M_H_node = 5000;// Mass of fluid in one of the hot fluid nodes [kg]
- HyS_Model.Units.Mass M_W_node = 0.0118;// Mass of one of the wall nodes [kg]
- HyS_Model.Units.Area A1_outer = 10;// Contact area between the tubes and the hot salt at node 1
- HyS_Model.Units.Area A2_outer = 3.668;// Contact area between the tubes and the hot salt at node 2
- HyS_Model.Units.Area A3_outer = 0.402;// Contact area between the tubes and the hot salt at node 3
- HyS_Model.Units.HeatTransferCoeff h_htc = 0.036;//1.2e-4;// Heat transfer coefficient between the tube and the two fluids
- HyS_Model.Units.Temperature T_boil = 227;// Saturated boiling temperature of secondary side [C]
- HyS_Model.Units.Temperature T_H_in = 800;// Decomposer hot side temperature [C]
- HyS_Model.Units.EnthalpyOfVaporization Hvap_SA = 0.073331;// Energy required to fully boil one mol of sulfuric acid (SA)
- HyS_Model.Units.EnthalpyOfVaporization Hvap_H2O = 0.035970;// Energy required to fully boil one mol of water
- HyS_Model.Units.UniversalGasConstant R = 0.00831447/1e3;// Universal Gas Constant [MJ/(mol*C)]
- HyS_Model.Units.Temperature T_ref_rxn = 298;//[K]
- HyS_Model.Units.Temperature T_ref_Kc = 973;//[K]
- HyS_Model.Units.Number epsilon = 0.25;
- HyS_Model.Units.Number Kc = 0.041;// Sulfur trioxide conversion constant
- HyS_Model.Units.Pressure Pressure = 100;// 1500// Secondary side pressure [kPa]
- HyS_Model.Units.HeatOfReaction H_rxn1 = 0.09753;//[MJ/mol]
- HyS_Model.Units.HeatOfReaction H_rxn2 = 0.09889;//0.09893;// [MJ/mol]
- HyS_Model.Units.MolarMass MM_SA = 0.0980785;// [kg/mol]
- HyS_Model.Units.MolarMass MM_H2O = 0.01801528;// [kg/mol]
- HyS_Model.Units.MolarMass MM_SO3 = 0.0800632;// [kg/mol]
- HyS_Model.Units.MolarMass MM_SO2 = 0.0640638;// [kg/mol]
- HyS_Model.Units.MolarMass MM_O2 = 0.0319988;// [kg/mol]
- HyS_Model.Units.FlowFraction FF_Decomposition = 1;// Flow fraction
- HyS_Model.Units.FlowFraction FF_HyS_Hot = 1;// Flow fraction
- HyS_Model.Units.FlowFraction FF_SA_recycle = 1;// Flow fraction
- HyS_Model.Units.Length pipe_thickness = 0.025;// Thickness of the tubes
- HyS_Model.Units.Time TauSArecycle_nom = 30;// Transit time for sulfuric acid exitting the decomposer then reentering it

=================
=== VARIABLES ===
=================
- HyS_Model.Units.Temperature T_vap;// Decomposer outlet temperature [C]
- HyS_Model.Units.Temperature T_H1;// Decomposer hot side node 1 temperature [C]
- HyS_Model.Units.Temperature T_H2;// Decomposer hot side node 2 temperature [C]
- HyS_Model.Units.Temperature T_H3;// Decomposer hot side node 3 temperature [C]
- HyS_Model.Units.MassFlowRate mdot_H;// Decomposer hot side mass flow rate 
- HyS_Model.Units.MassFlowRate mdot_S;// Decomposer secondary side mass flow rate 
- HyS_Model.Units.Temperature T_W1;// Tube node 1 temp [C]
- HyS_Model.Units.Temperature T_W2;// Tube node 2 temp [C]
- HyS_Model.Units.Temperature T_W3;// Tube node 3 temp [C]
- HyS_Model.Units.SpecificHeat Cp_SA;// Specific heat of sulfuric acid in gas phase
- HyS_Model.Units.SpecificHeat Cp_SA_liq;// Specific heat of sulfuric acid in liq phase
- HyS_Model.Units.SpecificHeat Cp_SO3;// Specific heat of sulfur trioxide in gas phase
- HyS_Model.Units.SpecificHeat Cp_SO2;// Specific heat of sulfur dioxide in gas phase
- HyS_Model.Units.SpecificHeat Cp_H2O;// Specific heat of water in gas phase
- HyS_Model.Units.SpecificHeat Cp_H2O_liq;// Specific heat of water in liq phase
- HyS_Model.Units.SpecificHeat Cp_O2;// Specific heat of oxygen gas
- HyS_Model.Units.SpecificHeat Cp_liq;// Specific heat of the liquid mix
- HyS_Model.Units.SpecificHeat Cp_vap;// Specific heat of the gas mix
- HyS_Model.Units.Power E_dot_boil;// Energy needed to fully boil the inflow
- HyS_Model.Units.Conversion X;// Ratio of SO3 that converts to SO2
- HyS_Model.Units.HeatOfReaction Q_rxn1;// Energy needed per mol of SA decomposition
- HyS_Model.Units.HeatOfReaction Q_rxn2;// Energy needed per mol of SO3 conversion
- HyS_Model.Units.HeatOfReaction Q_tot;// Heat used per mol of SO2 liberated. Used for efficiency comparisons
- HyS_Model.Units.HeatOfReaction Q_tot_CT;// Heat used per mol of SO2 liberated. This one uses CT's specific heat method
- HyS_Model.Units.HeatOfReaction dH_Rx;// Q_rxn2
- HyS_Model.Units.Temperature T_liq_in;// Temperature of fluid entering the decomposer secondary side
- HyS_Model.Units.Time vTauSArecycle;// Residency time
- HyS_Model.Units.SpecificHeat Cp_SO3_CT;// SO3 specific heat using CT's method
- HyS_Model.Units.SpecificHeat Cp_SO2_CT;// SO2 specific heat using CT's method
- HyS_Model.Units.SpecificHeat Cp_O2_CT;// O2 specific heat using CT's method
- HyS_Model.Units.HeatOfReaction Q_rxn2_CT;// SO3 conversion energy needed
- HyS_Model.Units.Conversion X_CT;// SO3 conversion ratio with CT's specific heat method
- HyS_Model.Units.Temperature delta_T_cond;// Change of temp via the condenser
- HyS_Model.Units.StandardEntropy S_SA;
- HyS_Model.Units.StandardEntropy S_H2;
- HyS_Model.Units.StandardEntropy S_SO2;
- HyS_Model.Units.StandardEntropy S_H2O;
- HyS_Model.Units.EnthalpyOfReaction dH0;
- HyS_Model.Units.Power P_demand;
- HyS_Model.Units.MolFlowRate Mdot_SA_rcy;
- HyS_Model.Units.MassFlowRate mdot_O2_out;
- HyS_Model.Units.MolFlowRate Mdot_SO2_out;
- HyS_Model.Units.MolFlowRate Mdot_H2O_Dout;
- HyS_Model.Units.MolFlowRate Mdot_rxn;
- HyS_Model.Units.MassFlowRate mdot_H2_out;
- HyS_Model.Units.MassFlowRate mdot_H2O_makeup;
- HyS_Model.Units.WeightPercent WtP_SA;
- HyS_Model.Units.WeightPercent WtP_SO2;
- HyS_Model.Units.Power Heat_used;

===========
== MISC ===
===========
~ Misc notes and tips will go here

























