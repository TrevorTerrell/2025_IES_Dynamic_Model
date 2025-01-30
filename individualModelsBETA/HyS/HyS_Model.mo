package HyS_Model
  package Hybrid_Sulfur
    model HyS
      // Parameter declarations
      // ===============
      // DECOMPOSISITION
      parameter HyS_Model.Units.MassFlowRate mdot_H_nom = 5400;// Mass flow rate of hot fluid on primary side
      parameter HyS_Model.Units.MolFlowRate Mdot_SA_nom = 62.309;
      parameter HyS_Model.Units.MolFlowRate Mdot_H2O_nom = 37.691;
      parameter HyS_Model.Units.SpecificHeat Cp_H = 0.00154912;
      parameter HyS_Model.Units.SpecificHeat Cp_W = 7.5e-4;// Specific heat of the tubing material (4.35e-4 for Inconel)
      parameter HyS_Model.Units.Mass M_H_node = 5000;// Mass of fluid in one of the hot fluid nodes [kg]
      parameter HyS_Model.Units.Mass M_W_node = 0.0118;// Mass of one of the wall nodes [kg]
    //  parameter HyS_Model.Units.Mass m_vap = 100;
      parameter HyS_Model.Units.Area A1_outer = 10;
      parameter HyS_Model.Units.Area A2_outer = 3.668;
      parameter HyS_Model.Units.Area A3_outer = 0.402;
    //  parameter HyS_Model.Units.Area A1_inner = 5;
      parameter HyS_Model.Units.HeatTransferCoeff h_htc = 0.036;//1.2e-4;// Heat transfer coefficient
      parameter HyS_Model.Units.Temperature T_boil = 227;// Saturated boiling temperature of secondary side [C]
      parameter HyS_Model.Units.Temperature T_H_in = 800;//600;
      parameter HyS_Model.Units.EnthalpyOfVaporization Hvap_SA = 0.073331;// Energy required to fully boil one mol of sulfuric acid (SA)
      parameter HyS_Model.Units.EnthalpyOfVaporization Hvap_H2O = 0.035970;// Energy required to fully boil one mol of water
      parameter HyS_Model.Units.UniversalGasConstant R = 0.00831447/1e3;// Universal Gas Constant [MJ/(mol*C)]
      parameter HyS_Model.Units.Temperature T_ref_rxn = 298;//[K]
      parameter HyS_Model.Units.Temperature T_ref_Kc = 973;//[K]
      parameter HyS_Model.Units.Number epsilon = 0.25;
      parameter HyS_Model.Units.Number Kc = 0.041;//0.041;
      parameter HyS_Model.Units.Pressure Pressure = 100;// 1500// Secondary side pressure [kPa]
      parameter HyS_Model.Units.HeatOfReaction H_rxn1 = 0.09753;//[MJ/mol]
      parameter HyS_Model.Units.HeatOfReaction H_rxn2 = 0.09889;//0.09893;// [MJ/mol]
      parameter HyS_Model.Units.MolarMass MM_SA = 0.0980785;// [kg/mol]
      parameter HyS_Model.Units.MolarMass MM_H2O = 0.01801528;
      parameter HyS_Model.Units.MolarMass MM_SO3 = 0.0800632;
      parameter HyS_Model.Units.MolarMass MM_SO2 = 0.0640638;
      parameter HyS_Model.Units.MolarMass MM_O2 = 0.0319988;
      parameter HyS_Model.Units.FlowFraction FF_Decomposition = 1;
      parameter HyS_Model.Units.FlowFraction FF_HyS_Hot = 1;
      parameter HyS_Model.Units.FlowFraction FF_SA_recycle = 1;
      parameter HyS_Model.Units.Length pipe_thickness = 0.025;
      parameter HyS_Model.Units.Time TauSArecycle_nom = 30;
      // ============
      // ELECTROLYSIS
      parameter HyS_Model.Units.Temperature T_elec = 127;
      parameter HyS_Model.Units.MolarMass MM_H2 = 0.002016;
      parameter HyS_Model.Units.Number z = 2;// Number of electrons
      parameter HyS_Model.Units.FaradayConstant F = 96485.3399;// Coulombs per mol of electrons
      parameter HyS_Model.Units.Volt E0_cell = 0.6;// From Gorensek
      parameter HyS_Model.Units.WeightPercent WtP_SA_nom = 0.9;
      parameter HyS_Model.Units.MolFlowRate Mdot_H2O_Eout = 200;//Nominal excess water flow rate
      parameter HyS_Model.Units.MolFlowRate Mdot_SA_in = 2000;//Nominal net sulfur flow rate
      // Variable declaration
      // ===============
      // DECOMPOSISITION
      HyS_Model.Units.Temperature T_vap;
      HyS_Model.Units.Temperature T_H1;
      HyS_Model.Units.Temperature T_H2;
      HyS_Model.Units.Temperature T_H3;
      HyS_Model.Units.MassFlowRate mdot_H;
      HyS_Model.Units.MassFlowRate mdot_S;
      HyS_Model.Units.Temperature T_W1;
      HyS_Model.Units.Temperature T_W2;
      HyS_Model.Units.Temperature T_W3;
      HyS_Model.Units.SpecificHeat Cp_SA;
      HyS_Model.Units.SpecificHeat Cp_SA_liq;
      HyS_Model.Units.SpecificHeat Cp_SO3;
      HyS_Model.Units.SpecificHeat Cp_SO2;
      HyS_Model.Units.SpecificHeat Cp_H2O;
      HyS_Model.Units.SpecificHeat Cp_H2O_liq;
      HyS_Model.Units.SpecificHeat Cp_O2;
      HyS_Model.Units.SpecificHeat Cp_liq;
      HyS_Model.Units.SpecificHeat Cp_vap;
      HyS_Model.Units.Power E_dot_boil;
      HyS_Model.Units.Conversion X;
      HyS_Model.Units.HeatOfReaction Q_rxn1;
      HyS_Model.Units.HeatOfReaction Q_rxn2;
      HyS_Model.Units.HeatOfReaction Q_tot;
      HyS_Model.Units.HeatOfReaction Q_tot_CT;
      HyS_Model.Units.HeatOfReaction dH_Rx;
      HyS_Model.Units.Temperature T_liq_in;
      HyS_Model.Units.Time vTauSArecycle;
      HyS_Model.Units.SpecificHeat Cp_SO3_CT;
      HyS_Model.Units.SpecificHeat Cp_SO2_CT;
      HyS_Model.Units.SpecificHeat Cp_O2_CT;
      HyS_Model.Units.HeatOfReaction Q_rxn2_CT;
      HyS_Model.Units.Conversion X_CT;
      // ============
      // ELECTROLYSIS
      HyS_Model.Units.Temperature delta_T_cond;
      HyS_Model.Units.StandardEntropy S_SA;
      HyS_Model.Units.StandardEntropy S_H2;
      HyS_Model.Units.StandardEntropy S_SO2;
      HyS_Model.Units.StandardEntropy S_H2O;
      HyS_Model.Units.EnthalpyOfReaction dH0;
      HyS_Model.Units.Power P_demand;
      // ============
      // MASS BALANCE
      HyS_Model.Units.MolFlowRate Mdot_SA_rcy;
      HyS_Model.Units.MassFlowRate mdot_O2_out;
      HyS_Model.Units.MolFlowRate Mdot_SO2_out;
      HyS_Model.Units.MolFlowRate Mdot_H2O_Dout;
      HyS_Model.Units.MolFlowRate Mdot_rxn;
      HyS_Model.Units.MassFlowRate mdot_H2_out;
      HyS_Model.Units.MassFlowRate mdot_H2O_makeup;
      HyS_Model.Units.WeightPercent WtP_SA;
      HyS_Model.Units.WeightPercent WtP_SO2;
      HyS_Model.Units.Power Heat_used;
    initial equation
      T_H1 = T_H_in;
      T_H2 = T_H_in - 25;
      T_H3 = T_H_in - 50;
      T_W1 = T_H_in - 10;
      T_W2 = T_H_in - 20;
      T_W3 = T_H_in - 30;
    equation
// ===============
// DECOMPOSISITION
      T_liq_in = T_elec;
      vTauSArecycle = TauSArecycle_nom/FF_SA_recycle;
    // =========
// Hot nodes (primary side)
      M_H_node*Cp_H*der(T_H1) = mdot_H*Cp_H*(T_H_in - T_H1) - h_htc*A1_outer/pipe_thickness*(((T_H_in + T_H1)/2) - T_W1);
      M_H_node*Cp_H*der(T_H2) = mdot_H*Cp_H*(T_H1 - T_H2) - h_htc*A2_outer/pipe_thickness*(((T_H1 + T_H2)/2) - T_W2);
      M_H_node*Cp_H*der(T_H3) = mdot_H*Cp_H*(T_H2 - T_H3) - h_htc*A3_outer/pipe_thickness*(((T_H2 + T_H3)/2) - T_W3);
    // ==========
// Wall nodes (tubing)
      M_W_node*Cp_W*der(T_W1) = h_htc*A1_outer/pipe_thickness*(((T_H_in + T_H1)/2) - T_W1) - mdot_S*Cp_vap*(T_vap - T_boil) - Mdot_rxn*(Q_rxn1 + Q_rxn2);
      M_W_node*Cp_W*der(T_W2) = h_htc*A2_outer/pipe_thickness*(((T_H1 + T_H2)/2) - T_W2) - E_dot_boil;
      M_W_node*Cp_W*der(T_W3) = h_htc*A3_outer/pipe_thickness*(((T_H2 + T_H3)/2) - T_W3) - mdot_S*Cp_liq*(T_boil - T_liq_in);
      E_dot_boil = Mdot_SA_in*Hvap_SA + Mdot_H2O_Eout*Hvap_H2O;
      Heat_used = mdot_S*Cp_liq*(T_boil - T_liq_in) + E_dot_boil + mdot_S*Cp_vap*(T_vap - T_boil) + Mdot_rxn*(Q_rxn1 + Q_rxn2);
    // ==========
// Cold nodes (secondary side)
    //  T_vap = 627;
  T_vap = T_H1 - 2*(T_H1 - T_W1);//2*T_W1-T_H1;
      Q_rxn1 = H_rxn1 + (Cp_H2O*MM_H2O + Cp_SO3*MM_SO3 - Cp_SA*MM_SA)*((T_vap + 273.15) - T_ref_rxn);
      Q_rxn2 = H_rxn2 + (-Cp_SO3*MM_SO3 + Cp_SO2*MM_SO2 + 0.5*Cp_O2*MM_O2)*((T_vap + 273.15) - T_ref_rxn);
      Q_tot = Q_rxn1 + Q_rxn2 + (mdot_S*Cp_liq*(T_boil - T_liq_in)/Mdot_SO2_out) + (E_dot_boil/Mdot_SO2_out) + (Cp_vap*mdot_S*(T_vap - T_boil)/Mdot_SO2_out);//HEAT USED PER MOL OF SO2 LIBERATED
    // ==========================
// Sulfur Trioxide Conversion
//  dH_Rx = (H_rxn2 + (-Cp_SO3*MM_SO3 + Cp_SO2*MM_SO2 + 0.5*Cp_O2*MM_O2)*((T_vap + 273.15) - T_ref_rxn))*(1e3);//[kJ/mol] CHECK UNITS
      dH_Rx  = H_rxn2 + (-Cp_SO3*MM_SO3 + Cp_SO2*MM_SO2 + 0.5*Cp_O2*MM_O2)*((T_vap + 273.15) - T_ref_rxn);//[MJ/mol]
//  Kc*exp((dH_Rx/R)*(((1/(T_ref_Kc))) - (1/(T_vap + 273.15)))) = ((4.86e-5*Pressure*(X/(1 + (epsilon*X)))*(1173/(T_vap + 273.15)))*(4.86e-5*Pressure*((0.5*X)/(1 + (epsilon*X)))*(1173/(T_vap + 273.15)))^(0.5))/(4.86e-5*Pressure*((1 - X)/(1 + (epsilon*X))));
//  Kc*exp((Q_rxn2/R)*(((1/(T_ref_Kc))) - (1/(T_vap + 273.15)))) = ((4.86e-5*Pressure*(X/(1 + (epsilon*X)))*(1173/(T_vap + 273.15)))*(4.86e-5*Pressure*((0.5*X)/(1 + (epsilon*X)))*(1173/(T_vap + 273.15)))^(0.5))/(4.86e-5*Pressure*((1 - X)/(1 + (epsilon*X))));
      Kc*exp((dH_Rx/R)*(((1/(T_ref_Kc))) - (1/(T_vap + 273.15)))) = ((4.86e-5*Pressure*(X/(1 + (epsilon*X)))*(1173/(T_vap + 273.15)))*(4.86e-5*Pressure*((0.5*X)/(1 + (epsilon*X)))*(1173/(T_vap + 273.15)))^(0.5))/(4.86e-5*Pressure*((1 - X)/(1 + (epsilon*X))));
    // ===========
// CT's method (just for comparison, not used elsewhere). The specific heats are per mol not per kg, like the Shomate ones
      Cp_SO3_CT = ((8.06) + ((1.056e-3)*(T_vap + 273.15)) - (2.028e5)/((T_vap + 273.15)^2))*R;
      Cp_SO2_CT = ((5.699) + ((8.01e-4)*(T_vap + 273.15)) - (1.015e5)/((T_vap + 273.15)^2))*R;
      Cp_O2_CT = ((3.693) + ((5.06e-4)*(T_vap + 273.15)) - (2.27e4)/((T_vap + 273.15)^2))*R;
      Q_rxn2_CT = H_rxn2 + (-Cp_SO3_CT + Cp_SO2_CT + 0.5*Cp_O2_CT)*((T_vap + 273.15) - T_ref_rxn);
      Kc*exp((Q_rxn2_CT/R)*(((1/(T_ref_Kc))) - (1/(T_vap + 273.15)))) = ((4.86e-5*Pressure*(X_CT/(1 + (epsilon*X_CT)))*(1173/(T_vap + 273.15)))*(4.86e-5*Pressure*((0.5*X_CT)/(1 + (epsilon*X_CT)))*(1173/(T_vap + 273.15)))^(0.5))/(4.86e-5*Pressure*((1 - X_CT)/(1 + (epsilon*X_CT))));
      Q_tot_CT = Q_rxn1 + Q_rxn2_CT + (mdot_S*Cp_liq*(T_boil - T_liq_in)/Mdot_SO2_out) + (E_dot_boil/Mdot_SO2_out) + (Cp_vap*mdot_S*(T_vap - T_boil)/Mdot_SO2_out);
    // =======================
// Specific Heat Constants [MJ/kg-K]. Uses the Shomate equations from NIST, then converts mols to kg
      Cp_SA = (47.28924 + 190.3314*((T_vap + 273.15)/1000) - 148.1299*(((T_vap + 273.15)/1000)^2) + 43.86631*(((T_vap + 273.15)/1000)^3) - 0.740016/(((T_vap + 273.15)/1000)^2))/(1e6*MM_SA);
      Cp_H2O = (30.092 + 6.832514*((T_vap + 273.15)/1000) + 6.793435*(((T_vap + 273.15)/1000)^2) - 2.534480*(((T_vap + 273.15)/1000)^3) + 0.082139/(((T_vap + 273.15)/1000)^2))/(1e6*MM_H2O);
      Cp_SO3 = (24.02503 + 119.4607*((T_vap + 273.15)/1000) - 94.38686*(((T_vap + 273.15)/1000)^2) + 26.96237*(((T_vap + 273.15)/1000)^3) - 0.117517/(((T_vap + 273.15)/1000)^2))/(1e6*MM_SO3);
      Cp_SO2 = (21.43049 + 74.35094*((T_vap + 273.15)/1000) - 57.75217*(((T_vap + 273.15)/1000)^2) + 16.35534*(((T_vap + 273.15)/1000)^3) + 0.086731/(((T_vap + 273.15)/1000)^2))/(1e6*MM_SO2);
      Cp_O2 = (30.03235 + 8.772972*((T_vap + 273.15)/1000) - 3.988133*(((T_vap + 273.15)/1000)^2) + 0.788313*(((T_vap + 273.15)/1000)^3) - 0.741599/(((T_vap + 273.15)/1000)^2))/(1e6*MM_O2);
      Cp_SA_liq = (47.28924 + 190.3314*((((T_boil - T_liq_in)/2) + 273.15)/1000) - 148.1299*(((((T_boil - T_liq_in)/2) + 273.15)/1000)^2) + 43.86631*(((((T_boil - T_liq_in)/2) + 273.15)/1000)^3) - 0.740016/(((((T_boil - T_liq_in)/2) + 273.15)/1000)^2))/(1e6*MM_SA);
      Cp_H2O_liq = (30.092 + 6.832514*((((T_boil - T_liq_in)/2) + 273.15)/1000) + 6.793435*(((((T_boil - T_liq_in)/2) + 273.15)/1000)^2) - 2.534480*(((((T_boil - T_liq_in)/2) + 273.15)/1000)^3) + 0.082139/(((((T_boil - T_liq_in)/2) + 273.15)/1000)^2))/(1e6*MM_H2O);
      Cp_liq = (WtP_SA*Cp_SA_liq) + ((1 - WtP_SA)*Cp_H2O_liq);
      Cp_vap = (WtP_SA*Cp_SA) + ((1 - WtP_SA)*Cp_H2O);
    // ====================
// Entropy calculations
      S_SA = (47.28924*log(((T_elec + 273.15)/1000)) + 190.3314*((T_elec + 273.15)/1000) - (148.1299*(((T_elec + 273.15)/1000)^2))/2 + (43.86631*(((T_elec + 273.15)/1000)^3))/3 + 0.740016/(2*(((T_elec + 273.15)/1000)^2)) + 301.2961)/(1e6);
      S_H2 = (33.066178*log(((T_elec + 273.15)/1000)) - 11.363417*((T_elec + 273.15)/1000) + (11.432816*(((T_elec + 273.15)/1000)^2))/2 - (2.772874*(((T_elec + 273.15)/1000)^3))/3 + 0.158558/(2*(((T_elec + 273.15)/1000)^2)) + 172.707974)/(1e6);
      S_SO2 = (21.43049*log(((T_elec + 273.15)/1000)) + 74.35094*((T_elec + 273.15)/1000) - (57.75217*(((T_elec + 273.15)/1000)^2))/2 + (16.35534*(((T_elec + 273.15)/1000)^3))/3 + 0.086731/(2*(((T_elec + 273.15)/1000)^2)) + 245.8872)/(1e6);
      S_H2O = (-203.6060*log(((T_elec + 273.15)/1000)) + 1523.290*((T_elec + 273.15)/1000) - (3196.413*(((T_elec + 273.15)/1000)^2))/2 + (2474.455*(((T_elec + 273.15)/1000)^3))/3 - 3.855326/(2*(((T_elec + 273.15)/1000)^2)) - 488.7163)/(1e6);
    // =================
// Electrolysis eqns
      dH0 = z*F*E0_cell/1e6 + (T_elec + 273.15)*((S_SA + S_H2) - (S_SO2 + 2*S_H2O));//Electric energy needed per mol of electrolysis rxn
      P_demand = dH0*Mdot_rxn;
      delta_T_cond = T_vap - T_elec;
    // ============
// Mass Balance
      mdot_H = mdot_H_nom*FF_HyS_Hot;//Hot salt flow rate
      mdot_S = Mdot_SA_in*MM_SA + Mdot_H2O_Eout*MM_H2O;
// ~~~
      Mdot_SA_rcy = Mdot_SA_in*(1 - X);
      Mdot_rxn = Mdot_SA_in*X;
      mdot_O2_out = 0.5*Mdot_rxn*MM_O2;
      Mdot_SO2_out = Mdot_rxn;
      Mdot_H2O_Dout = Mdot_H2O_Eout + Mdot_rxn;
// ~~~
      mdot_H2_out = Mdot_rxn*MM_H2;
      mdot_H2O_makeup = 2*Mdot_rxn*MM_H2O;//THIS TIMES A FLOW FRACTION
      WtP_SA = Mdot_SA_in*MM_SA/(Mdot_SA_in*MM_SA + Mdot_H2O_Eout*MM_H2O);
      WtP_SO2 = Mdot_SO2_out*MM_SO2/(Mdot_SO2_out*MM_SO2 + Mdot_H2O_Dout*MM_H2O);
      annotation(
        experiment(StartTime = 0, StopTime = 1000, Tolerance = 1e-06, Interval = 0.01));
    end HyS;
  end Hybrid_Sulfur;

  package Units
    type Temperature = Real(unit = "C", min = 0);
    type Length = Real(unit = "m", min = 0);
    type Area = Real(unit = "m^2", min = 0);
    type Volume = Real(unit = "m^3", min = 0);
    type Density = Real(unit = "kg/(m^3)", min = 0);
    type MolarMass = Real(unit = "kg/mol", min = 0);
    type Mass = Real(unit = "kg", min = 0);
    type Mol = Real(unit = "mol", min = 0);
    type MassFlowRate = Real(unit = "kg/s");
    type MolFlowRate = Real(unit = "mol/s");
    type Time = Real(unit = "s");
    type Pressure = Real(unit = "kPa", min = 0);
    // Reminder: 1Bar is equal to 100kPa
    type SpecificHeat = Real(unit = "MJ/(kg*C)");
    type HeatTransferCoeff = Real(unit = "MJ/(s*C*m)");
    type EnthalpyOfVaporization = Real(unit = "MJ/mol");
    type EnthalpyOfReaction = Real(unit = "MJ/mol");
    type EquilibriumConstant = Real(unit = "1");
    type Power = Real(unit = "MJ/s", min = 0);
    type UniversalGasConstant = Real(unit = "kJ/(mol*C)", min = 0);
    type Number = Real(unit = "1");
    type HeatOfReaction = Real(unit = "MJ/mol");
    type Hrxn_kJ = Real(unit = "kJ/mol");
    type WeightPercent = Real(unit = "1", min = 0, max = 1);
    type FlowFraction = Real(unit = "1", min = 0);
    type Energy = Real(unit = "MJ", min = 0);
    type NeutronDensity = Real(unit = "1/m^3", min = 0);
    type NominalNeutronPopulation = Real(unit = "1", min = 0);
    type Reactivity = Real(unit = "1");
    type Conversion = Real(unit = "1", min = 0, max = 1);
    type StandardEntropy = Real(unit = "MJ/(mol*K)");
    type FaradayConstant = Real(unit = "C/mol");
    type Volt = Real(unit = "J/C");
  end Units;

  package Connectors
    // Creates all inlets and outlets

    connector Temp_In
      HyS_Model.Units.Temperature T;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, fillPattern = FillPattern.Solid, extent = {{-51, 47}, {51, -47}}), Ellipse(origin = {-8, 36}, extent = {{0, -2}, {0, 2}})}),
        Icon(graphics = {Ellipse(origin = {1, -1}, fillPattern = FillPattern.Solid, extent = {{-51, 47}, {51, -47}})}));
    end Temp_In;

    connector Temp_Out
      HyS_Model.Units.Temperature T;
      annotation(
        Icon(graphics = {Ellipse(origin = {1, -1}, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}),
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}));
    end Temp_Out;

    connector DecayHeat_In
      HyS_Model.Units.SpecificHeat Q;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineColor = {0, 225, 255}, fillColor = {0, 225, 255}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}),
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {0, 225, 255}, fillColor = {0, 225, 255}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}));
    end DecayHeat_In;

    connector DecayHeat_Out
      HyS_Model.Units.SpecificHeat Q;
      annotation(
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {0, 225, 255}, fillColor = {0, 225, 255}, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}),
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineColor = {0, 225, 255}, fillColor = {0, 225, 255}, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}));
    end DecayHeat_Out;

    connector NominalNeutronPopulation
      HyS_Model.Units.NominalNeutronPopulation n;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}),
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}})}));
    end NominalNeutronPopulation;

    connector Reactivity
      HyS_Model.Units.Reactivity rho;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineColor = {78, 154, 6}, fillColor = {78, 154, 6}, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}})}),
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {78, 154, 6}, fillColor = {78, 154, 6}, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}})}));
    end Reactivity;

    connector ExternalReactivity
      HyS_Model.Units.Reactivity rho_Ex;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineColor = {78, 154, 6}, fillColor = {78, 154, 6}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}),
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {78, 154, 6}, fillColor = {78, 154, 6}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}));
    end ExternalReactivity;

    connector FlowFraction_In
      HyS_Model.Units.FlowFraction FF;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, fillColor = {245, 121, 0}, fillPattern = FillPattern.Solid, extent = {{-51, 47}, {51, -47}})}),
        Icon(graphics = {Ellipse(origin = {1, -1}, fillColor = {245, 121, 0}, fillPattern = FillPattern.Solid, extent = {{-51, 47}, {51, -47}})}));
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, fillColor = {245, 121, 0}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}),
        Icon(graphics = {Ellipse(origin = {1, -1}, fillColor = {245, 121, 0}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}));
    end FlowFraction_In;

    connector FlowFraction_Out
      HyS_Model.Units.FlowFraction FF;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineColor = {245, 121, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}),
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {245, 121, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}));
    end FlowFraction_Out;

    connector MassFlow_In
      HyS_Model.Units.MassFlowRate MF;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, fillColor = {85, 85, 255}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}),
        Icon(graphics = {Ellipse(origin = {1, -1}, fillColor = {85, 85, 255}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}));
    end MassFlow_In;

    connector MassFlow_Out
      HyS_Model.Units.MassFlowRate MF;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineColor = {85, 85, 255}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}),
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {85, 85, 255}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}));
    end MassFlow_Out;

    connector WeightPercent_In
      HyS_Model.Units.WeightPercent WP;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, fillColor = {255, 85, 255}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}),
        Icon(graphics = {Ellipse(origin = {1, -1}, fillColor = {255, 85, 255}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}));
    end WeightPercent_In;

    connector WeightPercent_Out
      HyS_Model.Units.WeightPercent WP;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineColor = {255, 85, 255}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}),
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {255, 85, 255}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}));
    end WeightPercent_Out;
  end Connectors;
end HyS_Model;
