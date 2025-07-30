package MCFR_HyS_BOP

  model Simbase
    //
    // Load parameters for Nuclear and HeatTransport packages
    // NOTE: These are the values to change for most transients. To turn a feature 'off', set the time for that feature (ex: tripPrimaryPump) to be much higher than the total simulation time. To turn it 'on', set the time to the desired value (usually 2000 or 5000 to allow settling after initialization).
    MCFR_HyS_BOP.HeatTransport.Core Core(P = 1200, TotalFuelMass = 638300, hAnom = 18, mdot_fuelNom = 45000) annotation(
      Placement(visible = true, transformation(origin = {-261, -89}, extent = {{-555, -555}, {555, 555}}, rotation = 0)));
    MCFR_HyS_BOP.HeatTransport.PrimaryHeatExchanger PHX(m_dot_pFluidNom_PHX = 45000, m_dot_sFluidNom_PHX = 27272.73, hApnNom_PHX = 20, hAsnNom_PHX = 20) annotation(
      Placement(visible = true, transformation(origin = {1237, 389}, extent = {{-427, -427}, {427, 427}}, rotation = 0)));
    MCFR_HyS_BOP.HeatTransport.SecondaryHeatExchanger SHX(m_dot_pFluidNom_SHX = 27272.73, m_dot_sFluidNom_SHX = 12910.56, hApnNom_SHX = 9.2014, hAsnNom_SHX = 9.2014) annotation(
      Placement(visible = true, transformation(origin = {1251, -302}, extent = {{-432, -432}, {432, 432}}, rotation = 0))); 
    MCFR_HyS_BOP.HeatTransport.OnceThroughSteamGenerator OTSG(m_dot_hitecNom = 6455.28, tauOTSGtoSHX = 10) annotation(Placement(visible = true, transformation(origin = {1386, -1133}, extent = {{-343, -343}, {343, 343}}, rotation = 0)));
    MCFR_HyS_BOP.HeatTransport.RBOP RBOP(m_dot_fdw_nom = 192.075, omega = 0, sinmag = 0.05, sinstart = 88888888, rampstart = 50000000, ramp_min = 0.2, m_dot_cool_nom = 20000, cooltriptime = 88888888) annotation(Placement(visible = true, transformation(origin = {279, -1258}, extent = {{-296, -296}, {296, 296}}, rotation = 0)));
    MCFR_HyS_BOP.HeatTransport.HyS HyS(mdot_H_nom = 6455.28, pressure_nom = 100, Mdot_SA_nom = 3453.3, Mdot_H2O_nom = 345.33, omega = 0, sinmag = 0.01, sinstart = 88888888, rampstart = 5000000, ramp_min = 0.2, sulfouttime = 88888888, prelieffail = 88888888) annotation(Placement(visible = true, transformation(origin = {2355, -608}, extent = {{-296, -296}, {296, 296}}, rotation = 0)));
    MCFR_HyS_BOP.Nuclear.mPKE mPKE(Lam = 8.77488E-06, lam = {1.24909E-02, 3.15594E-02, 1.10482E-01, 3.22830E-01, 1.34074, 9.03297}, beta = {2.12256E-04, 1.12099E-03, 1.11237E-03, 3.39814E-03, 1.19113E-03, 3.72466E-04}) annotation(
      Placement(visible = true, transformation(origin = {-1323, 315}, extent = {{-409, -409}, {409, 409}}, rotation = 0)));
    MCFR_HyS_BOP.Nuclear.DecayHeat decayheat(P = 1200, TotalFuelMass = 638300) annotation(
      Placement(visible = true, transformation(origin = {-1181, -471}, extent = {{-463, -463}, {463, 463}}, rotation = 0)));
    MCFR_HyS_BOP.Nuclear.Poisons Poisons(lam_bubble = 0.005) annotation(
      Placement(visible = true, transformation(origin = {-1265, -1015}, extent = {{-555, -555}, {555, 555}}, rotation = 0)));
    MCFR_HyS_BOP.Nuclear.ReactivityFeedback ReactivityFeedback(a_F = -4.4376e-05, a_R = 0.1164e-05, step_mag = -185*1e-5, omega = 0.1, sin_mag = 1*1e-5, ramp_rate = 11e-5, ramp_mag = 600e-5, stepInsertionTime = 500000000, sinInsertionTime = 88888888, rampInsertionTime = 500000000) annotation(
      Placement(visible = true, transformation(origin = {-519, -1063}, extent = {{545, -545}, {-545, 545}}, rotation = 0)));
    MCFR_HyS_BOP.Pumps.PrimaryLoopPump PLP(freeConvectionFF = 0.125, primaryPumpK = 0.092, tripPrimaryPump = 500000000) annotation(
      Placement(visible = true, transformation(origin = {364, 566}, extent = {{-262, -262}, {262, 262}}, rotation = 0)));
    MCFR_HyS_BOP.Pumps.SecondaryLoopPump SLP(freeConvectionFF = 0.05, secondaryPumpK = 0.092, tripSecondaryPump = 2000000000) annotation(
      Placement(visible = true, transformation(origin = {366, -96}, extent = {{-232, -232}, {232, 232}}, rotation = 0)));
    MCFR_HyS_BOP.Pumps.TertiaryLoopPump TLP(freeConvectionFF = 0.05, tertiaryPumpK = 0.092, tripTertiaryPump = 2000000000) annotation(
      Placement(visible = true, transformation(origin = {375, -583}, extent = {{-217, -217}, {217, 217}}, rotation = 0)));
    // Connections between modules. This is the information they send to each other.
  equation
    connect(PLP.primaryFF, Core.fuelFlowFraction) annotation(
      Line(points = {{364, 414}, {-483, 414}, {-483, -533}}, color = {245, 121, 0}));
    connect(PLP.primaryFF, mPKE.fuelFlowFrac) annotation(
      Line(points = {{364, 414}, {-1184, 414}, {-1184, 479.04}}, color = {245, 121, 0}));
    connect(PLP.primaryFF, PHX.primaryFF) annotation(
      Line(points = {{364, 414}, {938, 414}, {938, 602.5}}, color = {245, 121, 0}));
    connect(SLP.secondaryFF, SHX.secondaryFF) annotation(
      Line(points = {{366, -235}, {949, -235}, {949, -86}}, color = {245, 121, 0}));
    connect(SLP.secondaryFF, PHX.secondaryFF) annotation(
      Line(points = {{366, -235}, {1792, -235}, {1792, 176}}, color = {245, 121, 0}));
    connect(TLP.tertiaryFF, SHX.tertiaryFF) annotation(
      Line(points = {{375, -713}, {1813, -713}, {1813, -518}}, color = {245, 121, 0}));
    connect(Poisons.Poison_react, ReactivityFeedback.Poison_react) annotation(
      Line(points = {{-1154, -1237}, {-661, -1237}}, color = {78, 154, 6}));
    connect(Core.nPop, mPKE.n_population) annotation(
      Line(points = {{-705, 233}, {-1323, 233}, {-1323, 159.9}}, color = {20, 36, 248}));
    connect(mPKE.n_population, decayheat.nPop) annotation(
      Line(points = {{-1323, 159.58}, {-1329, 159.58}, {-1329, -369.42}}, color = {20, 36, 248}));
    connect(mPKE.n_population, Poisons.n_population) annotation(
      Line(points = {{-1323, 159.58}, {-1487, 159.58}, {-1487, -1237}}, color = {20, 36, 248}));
    connect(decayheat.decayHeat_Out, Core.P_decay) annotation(
      Line(points = {{-1014.32, -609.9}, {-802.82, -609.9}, {-802.82, -122}, {-705, -122}}, color = {0, 225, 255}));
    connect(decayheat.decayHeat_Out, PHX.P_decay) annotation(
      Line(points = {{-1014.32, -609.9}, {1357, -609.9}, {1357, 602.5}}, color = {0, 225, 255}));
    connect(ReactivityFeedback.feedback, mPKE.feedback) annotation(
      Line(points = {{-661, -954}, {-1322.7, -954}, {-1322.7, 479}}, color = {78, 154, 6}));
    connect(Core.fuelNode1, ReactivityFeedback.fuelNode1) annotation(
      Line(points = {{16.5, -489}, {-301, -489}, {-301, -878}}));
    connect(Core.fuelNode2, ReactivityFeedback.fuelNode2) annotation(
      Line(points = {{16.5, -255.5}, {-301, -255.5}, {-301, -1063}}));
    connect(Core.ReflcNode, ReactivityFeedback.ReflcNode) annotation(
      Line(points = {{16.5, 0}, {-301, 0}, {-301, -1226.5}}));
    connect(SHX.T_out_pFluid_SHX, PHX.T_in_sFluid_PHX) annotation(
      Line(points = {{1951, -172}, {1951, 40}, {1929, 40}, {1929, 252}}));
    connect(PHX.T_out_sFluid_PHX, SHX.T_in_pFluid_SHX) annotation(
      Line(points = {{810, 261}, {810, -172}}));
    connect(Core.temp_Out, PHX.T_in_pFluid_PHX) annotation(
      Line(points = {{16.5, 205}, {944, 205}, {944, 517}, {801, 517}}));
    connect(PHX.T_out_pFluid_PHX, Core.temp_In) annotation(
      Line(points = {{1929, 517}, {2446, 517}, {2446, 716}, {-878, 716}, {-878, -477.5}, {-705, -477.5}}));
  connect(TLP.tertiaryFF, OTSG.tertiaryFF) annotation(
      Line(points = {{375, -713}, {782, -713}, {782, -934}}, color = {245, 121, 0}));
  connect(SHX.T_out_sFluid_SHX, OTSG.T_in_hitec_OTSG) annotation(
      Line(points = {{819, -432}, {1082.5, -432}, {1082.5, -378}, {1082, -378}, {1082, -651}, {782, -651}, {782, -1030}}));
  connect(OTSG.T_out_hitec_OTSG, SHX.T_in_sFluid_OTSG) annotation(
      Line(points = {{1914, -1030}, {1951, -1030}, {1951, -440}}));
  connect(OTSG.m_dot_SH, RBOP.m_dot_OTSG) annotation(
      Line(points = {{1914, -1359}, {1540.5, -1359}, {1540.5, -1080}, {457, -1080}}, color = {170, 85, 255}));
  connect(RBOP.m_dot_fdw, OTSG.m_dot_fdw) annotation(
      Line(points = {{457, -1436}, {782, -1436}, {782, -1359}}, color = {170, 85, 255}));
  connect(OTSG.P_stm, RBOP.P_HP) annotation(
      Line(points = {{1359, -1359}, {1173.5, -1359}, {1173.5, -1080}, {279, -1080}}));
  connect(OTSG.T_SH1, RBOP.T_OTSG) annotation(
      Line(points = {{782, -1236}, {1467.5, -1236}, {1467.5, -1080}, {101, -1080}}));
  connect(RBOP.T_fdw, OTSG.T_fdw) annotation(
      Line(points = {{101, -1436}, {2278, -1436}, {2278, -1213}, {1914, -1213}, {1914, -1236}}));
  connect(SHX.T_out_sFluid_SHX, HyS.T_H_in) annotation(
      Line(points = {{820, -432}, {2177, -432}, {2177, -430}}));
  connect(TLP.tertiaryFF, HyS.FF_HyS_Hot) annotation(
      Line(points = {{376, -714}, {2355, -714}, {2355, -430}}, color = {245, 121, 0}));
  connect(HyS.T_hot_out, SHX.T_in_sFluid_HyS) annotation(
      Line(points = {{2177, -786}, {1950, -786}, {1950, -388}}));
    annotation(
      Diagram(coordinateSystem(extent = {{-1800, 1120}, {2740, -1740}})),
      version = "",
      uses,
      experiment(StartTime = 0, StopTime = 10000, Tolerance = 1e-06, Interval = 0.5));
  end Simbase;

  package Nuclear
    // Parameters and equations for reactivity effects (mPKE), decay heat, and poisons

    model mPKE
      // Modified Point Kinetic Equations. PKEs modified for circulation effects (advection of DNPs from the core)
      parameter MCFR_HyS_BOP.Units.NeutronGenerationTime Lam; // These neutronics parameters come from the MCFR-C neutronics model
      parameter MCFR_HyS_BOP.Units.PrecursorDecayConstant lam[6];
      parameter MCFR_HyS_BOP.Units.DelayedNeutronFrac beta[6];
      parameter MCFR_HyS_BOP.Units.ResidentTime tauCoreNom = 7.028; // Avg travel time within core
      parameter MCFR_HyS_BOP.Units.ResidentTime tauLoopNom = 7.028; // Avg travel time from core outlet to inlet
      //Variable declaration
      MCFR_HyS_BOP.Units.ResidentTime varTauCore;
      MCFR_HyS_BOP.Units.ResidentTime varTauLoop;
      MCFR_HyS_BOP.Units.PrecursorConc CG1;
      MCFR_HyS_BOP.Units.PrecursorConc CG2;
      MCFR_HyS_BOP.Units.PrecursorConc CG3;
      MCFR_HyS_BOP.Units.PrecursorConc CG4;
      MCFR_HyS_BOP.Units.PrecursorConc CG5;
      MCFR_HyS_BOP.Units.PrecursorConc CG6;
      MCFR_HyS_BOP.Units.PrecursorReturn CG1Return;
      MCFR_HyS_BOP.Units.PrecursorReturn CG2Return;
      MCFR_HyS_BOP.Units.PrecursorReturn CG3Return;
      MCFR_HyS_BOP.Units.PrecursorReturn CG4Return;
      MCFR_HyS_BOP.Units.PrecursorReturn CG5Return;
      MCFR_HyS_BOP.Units.PrecursorReturn CG6Return;
      MCFR_HyS_BOP.Units.Reactivity reactivity;
      MCFR_HyS_BOP.Units.DelayedNeutronFrac bterm[6];
      MCFR_HyS_BOP.Units.Reactivity rho0;
      input MCFR_HyS_BOP.Connectors.Reactivity feedback annotation(
        Placement(visible = true, transformation(origin = {0, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output MCFR_HyS_BOP.Connectors.NominalNeutronPopulation n_population annotation(
        Placement(visible = true, transformation(origin = {0, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input MCFR_HyS_BOP.Connectors.FlowFraction_In fuelFlowFrac annotation(
        Placement(visible = true, transformation(origin = {34, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {34, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    initial equation
      n_population.n = 1;
      CG1 = (beta[1]/Lam)*(1.0/(lam[1] - (exp(-lam[1]*tauLoopNom) - 1.0)/tauCoreNom));
      CG2 = (beta[2]/Lam)*(1.0/(lam[2] - (exp(-lam[2]*tauLoopNom) - 1.0)/tauCoreNom));
      CG3 = (beta[3]/Lam)*(1.0/(lam[3] - (exp(-lam[3]*tauLoopNom) - 1.0)/tauCoreNom));
      CG4 = (beta[4]/Lam)*(1.0/(lam[4] - (exp(-lam[4]*tauLoopNom) - 1.0)/tauCoreNom));
      CG5 = (beta[5]/Lam)*(1.0/(lam[5] - (exp(-lam[5]*tauLoopNom) - 1.0)/tauCoreNom));
      CG6 = (beta[6]/Lam)*(1.0/(lam[6] - (exp(-lam[6]*tauLoopNom) - 1.0)/tauCoreNom));
      bterm[1] = beta[1]/(1.0 + ((1.0 - exp(-lam[1]*tauLoopNom))/(lam[1]*tauCoreNom)));
      bterm[2] = beta[2]/(1.0 + ((1.0 - exp(-lam[2]*tauLoopNom))/(lam[2]*tauCoreNom)));
      bterm[3] = beta[3]/(1.0 + ((1.0 - exp(-lam[3]*tauLoopNom))/(lam[3]*tauCoreNom)));
      bterm[4] = beta[4]/(1.0 + ((1.0 - exp(-lam[4]*tauLoopNom))/(lam[4]*tauCoreNom)));
      bterm[5] = beta[5]/(1.0 + ((1.0 - exp(-lam[5]*tauLoopNom))/(lam[5]*tauCoreNom)));
      bterm[6] = beta[6]/(1.0 + ((1.0 - exp(-lam[6]*tauLoopNom))/(lam[6]*tauCoreNom)));
      rho0 = sum(beta) - sum(bterm);
    equation
      bterm[1] = beta[1]/(1.0 + ((1.0 - exp(-lam[1]*varTauLoop))/(lam[1]*varTauCore)));
      bterm[2] = beta[2]/(1.0 + ((1.0 - exp(-lam[2]*varTauLoop))/(lam[2]*varTauCore)));
      bterm[3] = beta[3]/(1.0 + ((1.0 - exp(-lam[3]*varTauLoop))/(lam[3]*varTauCore)));
      bterm[4] = beta[4]/(1.0 + ((1.0 - exp(-lam[4]*varTauLoop))/(lam[4]*varTauCore)));
      bterm[5] = beta[5]/(1.0 + ((1.0 - exp(-lam[5]*varTauLoop))/(lam[5]*varTauCore)));
      bterm[6] = beta[6]/(1.0 + ((1.0 - exp(-lam[6]*varTauLoop))/(lam[6]*varTauCore)));
      rho0 = sum(beta) - sum(bterm);
      reactivity = feedback.rho + rho0;
      der(n_population.n) = ((reactivity - sum(beta))/Lam*n_population.n) + (lam[1]*CG1) + (lam[2]*CG2) + (lam[3]*CG3) + (lam[4]*CG4) + (lam[5]*CG5) + (lam[6]*CG6);
      varTauCore = tauCoreNom/fuelFlowFrac.FF;
      varTauLoop = tauLoopNom/fuelFlowFrac.FF;
      CG1Return = delay(CG1, varTauLoop, 5000)*exp(-lam[1]*varTauLoop)/varTauCore;
      CG2Return = delay(CG2, varTauLoop, 5000)*exp(-lam[2]*varTauLoop)/varTauCore;
      CG3Return = delay(CG3, varTauLoop, 5000)*exp(-lam[3]*varTauLoop)/varTauCore;
      CG4Return = delay(CG4, varTauLoop, 5000)*exp(-lam[4]*varTauLoop)/varTauCore;
      CG5Return = delay(CG5, varTauLoop, 5000)*exp(-lam[5]*varTauLoop)/varTauCore;
      CG6Return = delay(CG6, varTauLoop, 5000)*exp(-lam[6]*varTauLoop)/varTauCore;
      der(CG1) = (beta[1]*n_population.n)/Lam - (lam[1]*CG1) - (CG1/varTauCore) + CG1Return;
      der(CG2) = (beta[2]*n_population.n)/Lam - (lam[2]*CG2) - (CG2/varTauCore) + CG2Return;
      der(CG3) = (beta[3]*n_population.n)/Lam - (lam[3]*CG3) - (CG3/varTauCore) + CG3Return;
      der(CG4) = (beta[4]*n_population.n)/Lam - (lam[4]*CG4) - (CG4/varTauCore) + CG4Return;
      der(CG5) = (beta[5]*n_population.n)/Lam - (lam[5]*CG5) - (CG5/varTauCore) + CG5Return;
      der(CG6) = (beta[6]*n_population.n)/Lam - (lam[6]*CG6) - (CG6/varTauCore) + CG6Return;
      annotation(
        Diagram(graphics = {Rectangle(origin = {2, 0}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, lineThickness = 0.5, extent = {{-50, 50}, {50, -50}}), Text(origin = {3, 1}, extent = {{37, -27}, {-37, 27}}, textString = "mPKE")}),
        Icon(graphics = {Rectangle(origin = {2, 0}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, lineThickness = 0.5, extent = {{-50, 50}, {50, -50}}), Text(origin = {3, 1}, extent = {{37, -27}, {-37, 27}}, textString = "mPKE")}));
    end mPKE;

    model DecayHeat
      // Decay Heat in fuel salt equations. Tracked through precursor groups similar to delayed neutrons. Three group method.
      parameter MCFR_HyS_BOP.Units.ReactorPower P;
      parameter MCFR_HyS_BOP.Units.DecayHeatYield DHYG1 = 9.635981959409105e-01; // Decay heat parameters from the NERTHUS dynamic model
      parameter MCFR_HyS_BOP.Units.DecayHeatYield DHYG2 = 3.560674858154914e-02;
      parameter MCFR_HyS_BOP.Units.DecayHeatYield DHYG3 = 7.950554775404400e-04;
      parameter MCFR_HyS_BOP.Units.Mass TotalFuelMass;
      parameter MCFR_HyS_BOP.Units.DecayHeatPrecursorDecayConstant DHlamG1 = 0.0945298;
      parameter MCFR_HyS_BOP.Units.DecayHeatPrecursorDecayConstant DHlamG2 = 0.00441957;
      parameter MCFR_HyS_BOP.Units.DecayHeatPrecursorDecayConstant DHlamG3 = 8.60979e-05;
      parameter MCFR_HyS_BOP.Units.FissFactor FissFactor = 2.464783802008740e-03; // Heat per fission relative to nominal power rating
      MCFR_HyS_BOP.Units.HeatTransferFraction DHG1;
      MCFR_HyS_BOP.Units.HeatTransferFraction DHG2;
      MCFR_HyS_BOP.Units.HeatTransferFraction DHG3;
      input MCFR_HyS_BOP.Connectors.NominalNeutronPopulation nPop annotation(
        Placement(visible = true, transformation(origin = {1, 33}, extent = {{-19, -19}, {19, 19}}, rotation = 0), iconTransformation(origin = {-32, 22}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output MCFR_HyS_BOP.Connectors.DecayHeat_Out decayHeat_Out annotation(
        Placement(visible = true, transformation(origin = {1, -31}, extent = {{-23, -23}, {23, 23}}, rotation = 0), iconTransformation(origin = {36, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    initial equation
      DHG1 = 0.025125;
      DHG2 = 0.0198578;
      DHG3 = 0.0227606;
    equation
      der(DHG1) = nPop.n*DHYG1*FissFactor - DHlamG1*DHG1;
      der(DHG2) = nPop.n*DHYG2*FissFactor - DHlamG2*DHG2;
      der(DHG3) = nPop.n*DHYG3*FissFactor - DHlamG3*DHG3;
      decayHeat_Out.Q = (DHG1 + DHG2 + DHG3)*P/TotalFuelMass;
      annotation(
        Diagram(graphics = {Rectangle(origin = {2, 0}, lineColor = {0, 225, 255}, lineThickness = 0.75, extent = {{-50, 50}, {50, -50}}), Text(origin = {3, 0}, extent = {{-45, 46}, {45, -46}}, textString = "Decay Heat")}),
        Icon(graphics = {Text(origin = {19, 41}, extent = {{-29, 33}, {29, -33}}, textString = "Decay Heat"), Rectangle(origin = {2, 0}, lineColor = {0, 225, 255}, lineThickness = 0.75, extent = {{-50, 50}, {50, -50}}), Text(origin = {1, 36}, extent = {{1, 2}, {-1, -2}}, textString = "text"), Text(origin = {-31, 9}, extent = {{-13, 5}, {13, -5}}, textString = "n(t)/n_0"), Text(origin = {19, -43}, extent = {{-31, 5}, {31, -5}}, textString = "Specific Decay Heat")}));
    end DecayHeat;

    model Poisons
      // Tracking of I-135, Xe-135, Pm-149, and Sm-149
      //Parameter declarations
      parameter MCFR_HyS_BOP.Units.IsotopicFissYield gamma_I = 0.06390; // Fraction of fissions that directly produce this isotope. Te135's fissyield is added to this one due to its short half live and that it decays into I135
      parameter MCFR_HyS_BOP.Units.IsotopicFissYield gamma_Xe = 0.0025761;
      parameter MCFR_HyS_BOP.Units.IsotopicFissYield gamma_Pm = 0.010816;
      parameter MCFR_HyS_BOP.Units.IsotopicFissYield gamma_Sm = 0.0;
      parameter MCFR_HyS_BOP.Units.IsotopicDecayConstant I135_lam = 2.926153244511758e-05; // Decay constants from https://periodictable.com/Isotopes
      parameter MCFR_HyS_BOP.Units.IsotopicDecayConstant Xe135_lam = 2.106574217602556e-05;
      parameter MCFR_HyS_BOP.Units.IsotopicDecayConstant Pm149_lam = 3.627371580423393e-06;
      parameter MCFR_HyS_BOP.Units.MicroAbsorptionCrossSection sig_Xe = 5.915e-24;//2.66449e-18; // Cross section for fast spectrum. Commented out is for thermal.
      parameter MCFR_HyS_BOP.Units.MicroAbsorptionCrossSection sig_Sm = 6.763e-24;//4.032e-20;
      parameter MCFR_HyS_BOP.Units.MacroAbsorptionCrossSection Sig_f = 0.00129788; // Fission cross section of fuel
      parameter MCFR_HyS_BOP.Units.MacroAbsorptionCrossSection Sig_a = 0.00302062; // Absorption cross section of fuel
      parameter MCFR_HyS_BOP.Units.NominalNeutronFlux phi_0 = 3.51513E+14; // Flux density in core
      parameter MCFR_HyS_BOP.Units.OffgasRemovalRate lam_bubble; // Removal rate of fission gases via bubble removal. See Dunkle_2023 "Effect of xenon removal rate on load following in high power thermal spectrum Molten-Salt Reactors (MSRs)" for more info.
      //Variable declarations
      MCFR_HyS_BOP.Units.IsotopicConc I135_conc;
      MCFR_HyS_BOP.Units.IsotopicConc Xe135_conc;
      MCFR_HyS_BOP.Units.IsotopicConc Pm149_conc;
      MCFR_HyS_BOP.Units.IsotopicConc Sm149_conc;
      MCFR_HyS_BOP.Units.IsotopicConc Xe135_0;
      MCFR_HyS_BOP.Units.IsotopicConc Sm149_0;
      MCFR_HyS_BOP.Units.Reactivity Xe_react;
      MCFR_HyS_BOP.Units.Reactivity Sm_react;
      input MCFR_HyS_BOP.Connectors.NominalNeutronPopulation n_population annotation(
        Placement(visible = true, transformation(origin = {-40, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-40, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output MCFR_HyS_BOP.Connectors.Reactivity Poison_react annotation(
        Placement(visible = true, transformation(origin = {20, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {20, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    initial equation
      I135_conc = gamma_I*Sig_f*phi_0/I135_lam;
      Xe135_conc = ((gamma_Xe + gamma_I)*Sig_f*phi_0)/(Xe135_lam + sig_Xe*phi_0 + lam_bubble);
      Pm149_conc = gamma_Pm*Sig_f*phi_0/Pm149_lam;
      Sm149_conc = ((gamma_Pm + gamma_Sm)*Sig_f)/(sig_Sm);
    equation
      Xe135_0 = ((gamma_Xe + gamma_I)*Sig_f*phi_0)/(Xe135_lam + sig_Xe*phi_0 + lam_bubble);
      Sm149_0 = ((gamma_Pm + gamma_Sm)*Sig_f)/(sig_Sm);
      der(I135_conc) = (n_population.n*gamma_I*Sig_f*phi_0) - (I135_lam*I135_conc);
      der(Xe135_conc) = (n_population.n*gamma_Xe*Sig_f*phi_0) + (I135_lam*I135_conc) - (sig_Xe*phi_0*n_population.n*Xe135_conc) - (Xe135_lam*Xe135_conc) - (lam_bubble*Xe135_conc);
      der(Pm149_conc) = (n_population.n*gamma_Pm*Sig_f*phi_0) - (Pm149_lam*Pm149_conc);
      der(Sm149_conc) = (n_population.n*gamma_Sm*Sig_f*phi_0) + (Pm149_lam*Pm149_conc) - (sig_Sm*phi_0*n_population.n*Sm149_conc);
      Xe_react = ((Xe135_0*sig_Xe)/Sig_a)*(1 - (Xe135_conc/Xe135_0));
      Sm_react = ((Sm149_0*sig_Sm)/Sig_a)*(1 - (Sm149_conc/Sm149_0));
      Poison_react.rho = Xe_react + Sm_react;
      annotation(
        Diagram(graphics = {Text(origin = {-35, -51}, extent = {{-17, 7}, {17, -7}}, textString = "Neutron Population"), Text(origin = {20, -51}, extent = {{-16, 5}, {16, -5}}, textString = "Poison Reactivity"), Text(origin = {-7, -7}, extent = {{-39, 17}, {39, -17}}, textString = "Poison Tracking"), Rectangle(origin = {-7, -28}, extent = {{-53, -32}, {53, 32}})}),
        Icon(graphics = {Text(origin = {-35, -51}, extent = {{-17, 7}, {17, -7}}, textString = "Neutron Population"), Text(origin = {20, -51}, extent = {{-16, 5}, {16, -5}}, textString = "Poison Reactivity"), Text(origin = {-7, -7}, extent = {{-39, 17}, {39, -17}}, textString = "Poison Tracking"), Rectangle(origin = {-7, -28}, extent = {{-53, -32}, {53, 32}})}));
    end Poisons;

    model ReactivityFeedback
      // Fuel temperature reactivity feedback, graphite temperature reactivity feedback, external reactivity, and total reactivity
      parameter MCFR_HyS_BOP.Units.Temperature FuelTempSetPointNode1 = 700; 
      parameter MCFR_HyS_BOP.Units.Temperature FuelTempSetPointNode2 = 720; // From Terrapower documents
      parameter MCFR_HyS_BOP.Units.Temperature ReflcTempSetPoint = 700;
      parameter MCFR_HyS_BOP.Units.TemperatureReactivityCoef a_F; // Fuel temp feedback coefficient. From MCFR-C neutronics modeling
      parameter MCFR_HyS_BOP.Units.TemperatureReactivityCoef a_R; // Reflector temp feedback coefficient. From MCFR-C neutronics modeling
      parameter MCFR_HyS_BOP.Units.Reactivity step_mag;
      parameter MCFR_HyS_BOP.Units.Frequency omega;
      parameter MCFR_HyS_BOP.Units.Reactivity sin_mag;
      parameter MCFR_HyS_BOP.Units.ReactivityRate ramp_rate;
      parameter MCFR_HyS_BOP.Units.Reactivity ramp_mag;
      parameter MCFR_HyS_BOP.Units.Reactivity dollar = 0.007407352; // Unused in current version. Can be used for reactivity in cents
      parameter MCFR_HyS_BOP.Units.InitiationTime stepInsertionTime;
      parameter MCFR_HyS_BOP.Units.InitiationTime sinInsertionTime;
      parameter MCFR_HyS_BOP.Units.InitiationTime rampInsertionTime;
      MCFR_HyS_BOP.Units.Reactivity FuelTempFeedbackNode1;
      MCFR_HyS_BOP.Units.Reactivity FuelTempFeedbackNode2;
      MCFR_HyS_BOP.Units.Reactivity ReflcTempFeedback;
      MCFR_HyS_BOP.Units.Reactivity TotalTempFeedback;
      MCFR_HyS_BOP.Units.Reactivity ReactExternalStep;
      MCFR_HyS_BOP.Units.Reactivity ReactExternalSin;
      MCFR_HyS_BOP.Units.Reactivity ReactExternalRamp;
      input MCFR_HyS_BOP.Connectors.Temp_In fuelNode1 annotation(
        Placement(visible = true, transformation(origin = {-40, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-40, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input MCFR_HyS_BOP.Connectors.Temp_In fuelNode2 annotation(
        Placement(visible = true, transformation(origin = {-40, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-40, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input MCFR_HyS_BOP.Connectors.Temp_In ReflcNode annotation(
        Placement(visible = true, transformation(origin = {-40, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-40, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output MCFR_HyS_BOP.Connectors.Reactivity feedback annotation(
        Placement(visible = true, transformation(origin = {26, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {26, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input MCFR_HyS_BOP.Connectors.Reactivity Poison_react annotation(
        Placement(visible = true, transformation(origin = {26, -32}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {26, -32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    initial equation
      ReactExternalStep = 0;
      ReactExternalSin = 0;
//  ReactExternalRamp = 0;
      FuelTempFeedbackNode1 = 0;
      FuelTempFeedbackNode2 = 0;
      ReflcTempFeedback = 0;
      TotalTempFeedback = 0;
    equation
      ReactExternalStep = 0 + (if time < stepInsertionTime then 0 else step_mag); // For step insertion
      ReactExternalSin = delay(sin_mag*sin(omega*time), sinInsertionTime); // For sinusoid insertion
      ReactExternalRamp = 0 + (if (time < rampInsertionTime) then 0 else if (time > rampInsertionTime + ramp_mag / ramp_rate) then ramp_mag else ramp_rate * (time - rampInsertionTime)); //For ramp insertion
      FuelTempFeedbackNode1 = (fuelNode1.T - FuelTempSetPointNode1)*(a_F/2); 
      FuelTempFeedbackNode2 = (fuelNode2.T - FuelTempSetPointNode2)*(a_F/2);
      ReflcTempFeedback = (ReflcNode.T - ReflcTempSetPoint)*a_R;
      TotalTempFeedback = FuelTempFeedbackNode1 + FuelTempFeedbackNode2 + ReflcTempFeedback + ReactExternalSin + ReactExternalStep + ReactExternalRamp + Poison_react.rho;
      feedback.rho = TotalTempFeedback;
      annotation(
        Diagram(graphics = {Rectangle(origin = {0, 2}, extent = {{-50, 50}, {50, -50}}), Text(origin = {13, 40}, extent = {{-31, 24}, {31, -24}}, textString = "Feedback"), Text(origin = {-32, 22}, extent = {{-16, 6}, {16, -6}}, textString = "Temp_FuelN1"), Text(origin = {-32, -12}, extent = {{-16, 6}, {16, -6}}, textString = "Temp_FuelN2"), Text(origin = {-32, -40}, extent = {{-16, 6}, {16, -6}}, textString = "Temp_GrapN1"), Text(origin = {26, 10}, extent = {{-16, 6}, {16, -6}}, textString = "Total_Feedback"), Text(origin = {26, -40}, extent = {{-12, 4}, {12, -4}}, textString = "Poison React")}),
        Icon(graphics = {Rectangle(origin = {0, 2}, extent = {{-50, 50}, {50, -50}}), Text(origin = {13, 40}, extent = {{-31, 24}, {31, -24}}, textString = "Feedback"), Text(origin = {-32, 22}, extent = {{-16, 6}, {16, -6}}, textString = "Temp_FuelN1"), Text(origin = {-32, -12}, extent = {{-16, 6}, {16, -6}}, textString = "Temp_FuelN2"), Text(origin = {-32, -40}, extent = {{-16, 6}, {16, -6}}, textString = "Temp_GrapN1"), Text(origin = {26, 10}, extent = {{-16, 6}, {16, -6}}, textString = "Total_Feedback"), Text(origin = {26, -40}, extent = {{-12, 4}, {12, -4}}, textString = "Poison React")}));
    end ReactivityFeedback;
  end Nuclear;
  package HeatTransport
    model Core
      // Fission power, decay power, nominal power, and core flow
      parameter MCFR_HyS_BOP.Units.ReactorPower P;
      parameter MCFR_HyS_BOP.Units.Mass TotalFuelMass;
      parameter MCFR_HyS_BOP.Units.Convection hAnom;
      parameter MCFR_HyS_BOP.Units.Mass m_FN1 = 159575; // From Terrapower documentation
      parameter MCFR_HyS_BOP.Units.Mass m_FN2 = 159575;
      parameter MCFR_HyS_BOP.Units.Mass m_RN = 969575;
      parameter MCFR_HyS_BOP.Units.SpecificHeatCapacity cP_fuel = 0.00066065; // From Terrapower documentation
      parameter MCFR_HyS_BOP.Units.SpecificHeatCapacity cP_reflector = 0.0013; 
      parameter MCFR_HyS_BOP.Units.MassFlowRate mdot_fuelNom;
      parameter MCFR_HyS_BOP.Units.VolumeImportance kFN1 = 0.465;
      parameter MCFR_HyS_BOP.Units.VolumeImportance kFN2 = 0.465;
      parameter MCFR_HyS_BOP.Units.VolumeImportance kR = 0.07;
      parameter MCFR_HyS_BOP.Units.HeatTransferFraction kHT_FN1 = 0.5;
      parameter MCFR_HyS_BOP.Units.HeatTransferFraction kHT_FN2 = 0.5;
      parameter MCFR_HyS_BOP.Units.ResidentTime tauCoreToPHX = 1; 
      parameter MCFR_HyS_BOP.Units.Temperature fuelNode1_initial = 700;
      parameter MCFR_HyS_BOP.Units.Temperature fuelNode2_initial = 720;
      parameter MCFR_HyS_BOP.Units.Temperature ReflcNode_initial = 700;
      parameter MCFR_HyS_BOP.Units.Mass m_DHRS = 31915; // Mass of fluid in the decay heat removal system (DHRS)
      parameter MCFR_HyS_BOP.Units.DHRStimeConstant DHRS_timeConstant = 10;
      parameter MCFR_HyS_BOP.Units.HeatFlowRate DHRS_MaxPowerRm = 0.1*1200;
      parameter MCFR_HyS_BOP.Units.HeatFlowRate DHRS_PowerBleed = 0.001*1200;
      parameter MCFR_HyS_BOP.Units.InitiationTime DHRS_time = 2000000000; // Time that the DHRS is engaged
      MCFR_HyS_BOP.Units.ResidentTime varTauCoreToPHX;
      MCFR_HyS_BOP.Units.NominalPower FissionPower;
      MCFR_HyS_BOP.Units.NominalPower DecayPower;
      MCFR_HyS_BOP.Units.NominalPower NomPower;
      MCFR_HyS_BOP.Units.ReactorPower PowerTotal;
      MCFR_HyS_BOP.Units.ReactorPower PowerFission;
      MCFR_HyS_BOP.Units.ReactorPower PowerDecay;
      MCFR_HyS_BOP.Units.MassFlowRate mdot_fuel;
      MCFR_HyS_BOP.Units.Convection hA;
      MCFR_HyS_BOP.Units.Temperature DHRS_Temp;
      MCFR_HyS_BOP.Units.HeatFlowRate DHRS_PowerRm;
      input MCFR_HyS_BOP.Connectors.Temp_In temp_In annotation(
        Placement(visible = true, transformation(origin = {-76, -70}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {-80, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input MCFR_HyS_BOP.Connectors.DecayHeat_In P_decay annotation(
        Placement(visible = true, transformation(origin = {-77, -7}, extent = {{-17, -17}, {17, 17}}, rotation = 0), iconTransformation(origin = {-80, -6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input MCFR_HyS_BOP.Connectors.NominalNeutronPopulation nPop annotation(
        Placement(visible = true, transformation(origin = {-80, 50}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {-80, 58}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output MCFR_HyS_BOP.Connectors.Temp_Out temp_Out annotation(
        Placement(visible = true, transformation(origin = {44, 60}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {50, 53}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
      output MCFR_HyS_BOP.Connectors.Temp_Out fuelNode1 annotation(
        Placement(visible = true, transformation(origin = {45, -78}, extent = {{-17, -17}, {17, 17}}, rotation = 0), iconTransformation(origin = {50, -72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output MCFR_HyS_BOP.Connectors.Temp_Out fuelNode2 annotation(
        Placement(visible = true, transformation(origin = {45, -40}, extent = {{-17, -17}, {17, 17}}, rotation = 0), iconTransformation(origin = {50, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output MCFR_HyS_BOP.Connectors.Temp_Out ReflcNode annotation(
        Placement(visible = true, transformation(origin = {44, 16}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {50, 16}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input MCFR_HyS_BOP.Connectors.FlowFraction_In fuelFlowFraction annotation(
        Placement(visible = true, transformation(origin = {-40, -80}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {-40, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    initial equation
      fuelNode2.T = fuelNode2_initial;
      fuelNode1.T = fuelNode1_initial;
      ReflcNode.T = ReflcNode_initial;
      DHRS_Temp = 718;
    equation
      FissionPower = nPop.n*(1 - 0.067);
      DecayPower = P_decay.Q/P*TotalFuelMass;
      NomPower = FissionPower + DecayPower - (DHRS_PowerRm/P);
      mdot_fuel = mdot_fuelNom*fuelFlowFraction.FF;
//  hA = hAnom*(0.8215*fuelFlowFraction.FF^6 - 4.108*fuelFlowFraction.FF^5 + 7.848*fuelFlowFraction.FF^4 - 7.165*fuelFlowFraction.FF^3 + 3.004*fuelFlowFraction.FF^2 + 0.5903*fuelFlowFraction.FF + 0.008537); // UNUSED method for determining hA as a function of flow fraction
      hA = hAnom;
      varTauCoreToPHX = tauCoreToPHX/fuelFlowFraction.FF; // Variable due to flow fraction. (Slow transport = high residency time)
      m_FN1*cP_fuel*der(fuelNode1.T) = mdot_fuel*cP_fuel*(temp_In.T - fuelNode1.T) + kFN1*FissionPower*P - hA*kHT_FN1*(fuelNode1.T - ReflcNode.T) + P_decay.Q*m_FN1;
      m_FN2*cP_fuel*der(fuelNode2.T) = mdot_fuel*cP_fuel*(fuelNode1.T - fuelNode2.T) + kFN2*FissionPower*P - hA*kHT_FN2*(fuelNode1.T - ReflcNode.T) + P_decay.Q*m_FN2;
      m_RN*cP_reflector*der(ReflcNode.T) = hA*(fuelNode1.T - ReflcNode.T) + kR*FissionPower*P;
      DHRS_PowerRm = (DHRS_MaxPowerRm - DHRS_PowerBleed)/(1 + exp(log(1/1E-3 - 1)*(1 - (time - DHRS_time)/DHRS_timeConstant))) + DHRS_PowerBleed;
      m_DHRS*cP_fuel*der(DHRS_Temp) = mdot_fuel*cP_fuel*(fuelNode2.T - DHRS_Temp) + P_decay.Q*m_DHRS - DHRS_PowerRm;
      temp_Out.T = delay(DHRS_Temp, varTauCoreToPHX, 5000); // Delayed to simulate transport time between components
      PowerTotal = NomPower * P;
      PowerFission = FissionPower * P;
      PowerDecay = DecayPower * P;
      annotation(
        Diagram(graphics = {Rectangle(origin = {-15.4297, -10.014}, extent = {{-75.8099, 90.0986}, {75.8099, -90.0986}}), Rectangle(origin = {-40, 20}, lineColor = {255, 63, 63}, fillColor = {255, 63, 63}, fillPattern = FillPattern.Solid, lineThickness = 0, extent = {{-20, 20}, {20, -20}}), Rectangle(origin = {-40, -40}, lineColor = {255, 63, 63}, fillColor = {255, 63, 63}, fillPattern = FillPattern.Solid, lineThickness = 0, extent = {{-20, 20}, {20, -20}}), Rectangle(origin = {15, -10}, lineColor = {94, 95, 92}, fillColor = {94, 95, 92}, fillPattern = FillPattern.Solid, lineThickness = 0, extent = {{-15, 50}, {15, -50}}), Text(origin = {-16, 70}, extent = {{-47, 12}, {47, -12}}, textString = "Core Heat Transfer"), Text(origin = {-76, 35}, extent = {{-13, 13}, {13, -13}}, textString = "n(t)/n_0"), Text(origin = {-76, -20}, extent = {{-13, 13}, {13, -13}}, textString = "Decay Heat"), Text(origin = {-75, -85}, extent = {{-13, 13}, {13, -13}}, textString = "Temp_In"), Text(origin = {46, 45}, extent = {{-13, 13}, {13, -13}}, textString = "Temp_Out"), Text(origin = {47, 2}, extent = {{-13, 13}, {13, -13}}, textString = "Temp_G1"), Text(origin = {45, -54}, extent = {{-13, 13}, {13, -13}}, textString = "Temp_FN2"), Text(origin = {45, -92}, extent = {{-13, 13}, {13, -13}}, textString = "Temp_FN1")}),
        Icon(graphics = {Rectangle(origin = {-15.4297, -10.014}, extent = {{-75.8099, 90.0986}, {75.8099, -90.0986}}), Rectangle(origin = {-40, 20}, lineColor = {255, 63, 63}, fillColor = {255, 63, 63}, fillPattern = FillPattern.Solid, lineThickness = 0, extent = {{-20, 20}, {20, -20}}), Rectangle(origin = {-40, -40}, lineColor = {255, 63, 63}, fillColor = {255, 63, 63}, fillPattern = FillPattern.Solid, lineThickness = 0, extent = {{-20, 20}, {20, -20}}), Rectangle(origin = {15, -10}, lineColor = {94, 95, 92}, fillColor = {94, 95, 92}, fillPattern = FillPattern.Solid, lineThickness = 0, extent = {{-15, 50}, {15, -50}}), Text(origin = {-16, 70}, extent = {{-47, 12}, {47, -12}}, textString = "Core Heat Transfer"), Text(origin = {-76, 35}, extent = {{-13, 13}, {13, -13}}, textString = "n(t)/n_0"), Text(origin = {-76, -20}, extent = {{-13, 13}, {13, -13}}, textString = "Decay Heat"), Text(origin = {-75, -85}, extent = {{-13, 13}, {13, -13}}, textString = "Temp_In"), Text(origin = {46, 45}, extent = {{-13, 13}, {13, -13}}, textString = "Temp_Out"), Text(origin = {47, 2}, extent = {{-13, 13}, {13, -13}}, textString = "Temp_G1"), Text(origin = {45, -54}, extent = {{-13, 13}, {13, -13}}, textString = "Temp_FN2"), Text(origin = {45, -92}, extent = {{-13, 13}, {13, -13}}, textString = "Temp_FN1")}),
        experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.002));
    end Core;

    // Parameters and equations for all heat exchangers and pumps

    model PrimaryHeatExchanger
      // PHX parameters and equations
      // Parameter declaration
      parameter MCFR_HyS_BOP.Units.Mass m_PN1_PHX = 63830;
      // Fuel Salt node 1 mass
      parameter MCFR_HyS_BOP.Units.Mass m_PN2_PHX = 63830;
      // Fuel Salt node 2 mass
      parameter MCFR_HyS_BOP.Units.Mass m_PN3_PHX = 63830;
      // Fuel Salt node 3 mass
      parameter MCFR_HyS_BOP.Units.Mass m_PN4_PHX = 63830;
      // Fuel Salt node 4 mass
      parameter MCFR_HyS_BOP.Units.Mass m_TN1_PHX = 2860;
      // PHX Tubing node 1 mass
      parameter MCFR_HyS_BOP.Units.Mass m_TN2_PHX = 2860;
      // PHX Tubing node 2 mass
      parameter MCFR_HyS_BOP.Units.Mass m_SN1_PHX = 54545;
      // Secondary Side node 1 mass
      parameter MCFR_HyS_BOP.Units.Mass m_SN2_PHX = 54545;
      // Secondary Side node 2 mass
      parameter MCFR_HyS_BOP.Units.Mass m_SN3_PHX = 54545;
      // Secondary Side node 3 mass
      parameter MCFR_HyS_BOP.Units.Mass m_SN4_PHX = 54545;
      // Secondary Side node 4 mass
      parameter MCFR_HyS_BOP.Units.Mass m_pipe = 31915;
      // Mass of fluid in pipe between PHX to Core. Necessary due to decay heat present even in transport
      parameter MCFR_HyS_BOP.Units.SpecificHeatCapacity cP_pFluid_PHX = 0.00066065; // From Terrapower documents
      // Fuel salt specific heat [MJ/(kg*C)]
      parameter MCFR_HyS_BOP.Units.SpecificHeatCapacity cP_Tube_PHX = 0.00057779; // specific heat capacity of tubes (MJ/(kg-C)) ORNL-TM-0728 p.20
      // PHX Tubing specific heat [MJ/(kg*C)]
      parameter MCFR_HyS_BOP.Units.SpecificHeatCapacity cP_sFluid_PHX = 0.0011; // Secondary salt. See Dunkle_2025 "Dynamic Modeling of a Fast Spectrum Molten Salt Reactor Integrated Energy System" Table 1
      // Coolant salt specific heat [MJ/(kg*C)]
      parameter MCFR_HyS_BOP.Units.MassFlowRate m_dot_pFluidNom_PHX;
      // Fuel salt nominal flow rate
      parameter MCFR_HyS_BOP.Units.MassFlowRate m_dot_sFluidNom_PHX;
      // Coolant salt nominal flow rate
      parameter MCFR_HyS_BOP.Units.Convection hApnNom_PHX;
      // PHX primary side heat transfer coefficient [MJ/(s*C)]
      parameter MCFR_HyS_BOP.Units.Convection hAsnNom_PHX;
      // PHX secondary side heat transfer coefficient [MJ/(s*C)]
      parameter MCFR_HyS_BOP.Units.ResidentTime tauPHXtoCore = 1;
      // Avg travel time from PHX to pipe
      parameter MCFR_HyS_BOP.Units.ResidentTime tauPHXtoSHX = 10;
      // Avg travel time from PHX to SHX
      //parameter MCFR_HyS_BOP.Units.ResidentTime tauPipeToCore = 0.6887;
      parameter MCFR_HyS_BOP.Units.Temperature T_PN1_PHX_initial = 708.6;
      parameter MCFR_HyS_BOP.Units.Temperature T_PN2_PHX_initial = 699.1;
      parameter MCFR_HyS_BOP.Units.Temperature T_PN3_PHX_initial = 689.3;
      parameter MCFR_HyS_BOP.Units.Temperature T_PN4_PHX_initial = 679.5;
      parameter MCFR_HyS_BOP.Units.Temperature T_TN1_PHX_initial = 693.8;
      parameter MCFR_HyS_BOP.Units.Temperature T_TN2_PHX_initial = 674.1;
      parameter MCFR_HyS_BOP.Units.Temperature T_SN1_PHX_initial = 659.0;
      parameter MCFR_HyS_BOP.Units.Temperature T_SN2_PHX_initial = 669.1;
      parameter MCFR_HyS_BOP.Units.Temperature T_SN3_PHX_initial = 679.0;
      parameter MCFR_HyS_BOP.Units.Temperature T_SN4_PHX_initial = 688.9;
      // Variable declaration
      MCFR_HyS_BOP.Units.ResidentTime varTauPHXtoCore;
      MCFR_HyS_BOP.Units.ResidentTime varTauPHXtoSHX;
      MCFR_HyS_BOP.Units.MassFlowRate m_dot_pFluid_PHX;
      MCFR_HyS_BOP.Units.MassFlowRate m_dot_sFluid_PHX;
      MCFR_HyS_BOP.Units.Convection hApn_PHX;
      MCFR_HyS_BOP.Units.Convection hAsn_PHX;
      MCFR_HyS_BOP.Units.Temperature T_PN1_PHX;
      MCFR_HyS_BOP.Units.Temperature T_PN2_PHX;
      MCFR_HyS_BOP.Units.Temperature T_PN3_PHX;
      MCFR_HyS_BOP.Units.Temperature T_PN4_PHX;
      MCFR_HyS_BOP.Units.Temperature T_TN1_PHX;
      MCFR_HyS_BOP.Units.Temperature T_TN2_PHX;
      MCFR_HyS_BOP.Units.Temperature T_SN1_PHX;
      MCFR_HyS_BOP.Units.Temperature T_SN2_PHX;
      MCFR_HyS_BOP.Units.Temperature T_SN3_PHX;
      MCFR_HyS_BOP.Units.Temperature T_SN4_PHX;
      MCFR_HyS_BOP.Units.ReactorPower P_primaryLoop;
      MCFR_HyS_BOP.Units.Temperature Pipe_Temp;
      // Connections to other modules
      input MCFR_HyS_BOP.Connectors.Temp_In T_in_pFluid_PHX annotation(
        Placement(visible = true, transformation(origin = {-100, 30}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-102, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input MCFR_HyS_BOP.Connectors.Temp_In T_in_sFluid_PHX annotation(
        Placement(visible = true, transformation(origin = {160, -30}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {162, -32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output MCFR_HyS_BOP.Connectors.Temp_Out T_out_sFluid_PHX annotation(
        Placement(visible = true, transformation(origin = {-100, -30}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-100, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output MCFR_HyS_BOP.Connectors.Temp_Out T_out_pFluid_PHX annotation(
        Placement(visible = true, transformation(origin = {160, 30}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {162, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input MCFR_HyS_BOP.Connectors.DecayHeat_In P_decay annotation(
        Placement(visible = true, transformation(origin = {30, 50}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {28, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input MCFR_HyS_BOP.Connectors.FlowFraction_In primaryFF annotation(
        Placement(visible = true, transformation(origin = {-70, 50}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {-70, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input MCFR_HyS_BOP.Connectors.FlowFraction_In secondaryFF annotation(
        Placement(visible = true, transformation(origin = {130, -50}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {130, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    initial equation
      T_PN1_PHX = T_PN1_PHX_initial;
      T_PN2_PHX = T_PN2_PHX_initial;
      T_PN3_PHX = T_PN3_PHX_initial;
      T_PN4_PHX = T_PN4_PHX_initial;
      T_TN1_PHX = T_TN1_PHX_initial;
      T_TN2_PHX = T_TN2_PHX_initial;
      T_SN1_PHX = T_SN1_PHX_initial;
      T_SN2_PHX = T_SN2_PHX_initial;
      T_SN3_PHX = T_SN3_PHX_initial;
      T_SN4_PHX = T_SN4_PHX_initial;
      Pipe_Temp = 680.8;
    equation
      m_dot_pFluid_PHX = m_dot_pFluidNom_PHX*primaryFF.FF;
      m_dot_sFluid_PHX = m_dot_sFluidNom_PHX*secondaryFF.FF;
      varTauPHXtoCore = tauPHXtoCore/primaryFF.FF;
      varTauPHXtoSHX = tauPHXtoSHX/secondaryFF.FF;
// These next two are the simpler version for constant htc
      hApn_PHX = hApnNom_PHX;
      hAsn_PHX = hAsnNom_PHX;
// These next two are optional and not currently used. They attempt to include some complexity of htc changing due to flow rate turbulence
//hApn_PHX = hApnNom_PHX*(0.8215*primaryFF.FF^6 - 4.108*primaryFF.FF^5 + 7.848*primaryFF.FF^4 - 7.165*primaryFF.FF^3 + 3.004*primaryFF.FF^2 + 0.5903*primaryFF.FF + 0.008537);
//hAsn_PHX = hAsnNom_PHX*(0.8215*secondaryFF.FF^6 - 4.108*secondaryFF.FF^5 + 7.848*secondaryFF.FF^4 - 7.165*secondaryFF.FF^3 + 3.004*secondaryFF.FF^2 + 0.5903*secondaryFF.FF + 0.008537);
      m_PN1_PHX*cP_pFluid_PHX*der(T_PN1_PHX) = m_dot_pFluid_PHX*cP_pFluid_PHX*(T_in_pFluid_PHX.T - T_PN1_PHX) - hApn_PHX*(T_PN1_PHX - T_TN1_PHX) + P_decay.Q*m_PN1_PHX;
      m_PN2_PHX*cP_pFluid_PHX*der(T_PN2_PHX) = m_dot_pFluid_PHX*cP_pFluid_PHX*(T_PN1_PHX - T_PN2_PHX) - hApn_PHX*(T_PN1_PHX - T_TN1_PHX) + P_decay.Q*m_PN2_PHX;
      m_PN3_PHX*cP_pFluid_PHX*der(T_PN3_PHX) = m_dot_pFluid_PHX*cP_pFluid_PHX*(T_PN2_PHX - T_PN3_PHX) - hApn_PHX*(T_PN3_PHX - T_TN2_PHX) + P_decay.Q*m_PN3_PHX;
      m_PN4_PHX*cP_pFluid_PHX*der(T_PN4_PHX) = m_dot_pFluid_PHX*cP_pFluid_PHX*(T_PN3_PHX - T_PN4_PHX) - hApn_PHX*(T_PN3_PHX - T_TN2_PHX) + P_decay.Q*m_PN4_PHX;
      
      m_TN1_PHX*cP_Tube_PHX*der(T_TN1_PHX) = 2*hApn_PHX*(T_PN1_PHX - T_TN1_PHX) - 2*hAsn_PHX*(T_TN1_PHX - T_SN3_PHX);
      m_TN2_PHX*cP_Tube_PHX*der(T_TN2_PHX) = 2*hApn_PHX*(T_PN3_PHX - T_TN2_PHX) - 2*hAsn_PHX*(T_TN2_PHX - T_SN1_PHX);
      
      m_SN1_PHX*cP_sFluid_PHX*der(T_SN1_PHX) = m_dot_sFluid_PHX*cP_sFluid_PHX*(T_in_sFluid_PHX.T - T_SN1_PHX) + hAsn_PHX*(T_TN2_PHX - T_SN1_PHX);
      m_SN2_PHX*cP_sFluid_PHX*der(T_SN2_PHX) = m_dot_sFluid_PHX*cP_sFluid_PHX*(T_SN1_PHX - T_SN2_PHX) + hAsn_PHX*(T_TN2_PHX - T_SN1_PHX);
      m_SN3_PHX*cP_sFluid_PHX*der(T_SN3_PHX) = m_dot_sFluid_PHX*cP_sFluid_PHX*(T_SN2_PHX - T_SN3_PHX) + hAsn_PHX*(T_TN1_PHX - T_SN3_PHX);
      m_SN4_PHX*cP_sFluid_PHX*der(T_SN4_PHX) = m_dot_sFluid_PHX*cP_sFluid_PHX*(T_SN3_PHX - T_SN4_PHX) + hAsn_PHX*(T_TN1_PHX - T_SN3_PHX);
      
      T_out_sFluid_PHX.T = delay(T_SN4_PHX, varTauPHXtoSHX, 5000); // Delay for transport time
      m_pipe*cP_pFluid_PHX*der(Pipe_Temp) = m_dot_pFluid_PHX*cP_pFluid_PHX*(T_PN4_PHX - Pipe_Temp) + P_decay.Q*m_pipe;
      T_out_pFluid_PHX.T = delay(Pipe_Temp, varTauPHXtoCore, 5000);
      P_primaryLoop = m_dot_pFluid_PHX*cP_pFluid_PHX*(T_in_pFluid_PHX.T - T_PN4_PHX) + P_decay.Q*4*m_PN1_PHX;
      annotation(
        Diagram(graphics = {Rectangle(origin = {30.2864, 0.0340212}, extent = {{-149.242, 59.993}, {149.242, -59.993}}), Rectangle(origin = {-60, 30}, lineColor = {170, 0, 0}, fillColor = {170, 0, 0}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {120, 30}, lineColor = {170, 0, 0}, fillColor = {170, 0, 0}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {0, 30}, lineColor = {170, 0, 0}, fillColor = {170, 0, 0}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {-59.95, -29.61}, lineColor = {0, 85, 255}, fillColor = {0, 85, 255}, fillPattern = FillPattern.Solid, extent = {{-19.63, -9.61}, {19.63, 9.61}}), Rectangle(origin = {-0.12, -29.82}, lineColor = {0, 85, 255}, fillColor = {0, 85, 255}, fillPattern = FillPattern.Solid, extent = {{-19.88, 9.82}, {19.88, -9.82}}), Rectangle(origin = {59.8, -29.85}, lineColor = {0, 85, 255}, fillColor = {0, 85, 255}, fillPattern = FillPattern.Solid, extent = {{-19.8, 9.85}, {19.8, -9.85}}), Rectangle(origin = {119.85, -29.89}, lineColor = {0, 85, 255}, fillColor = {0, 85, 255}, fillPattern = FillPattern.Solid, extent = {{-19.85, 9.89}, {19.85, -9.89}}), Rectangle(origin = {89.92, -1.05}, lineColor = {136, 138, 133}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-49.92, 8.83}, {49.92, -8.83}}), Rectangle(origin = {60, 30}, lineColor = {170, 0, 0}, fillColor = {170, 0, 0}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {-30.08, -1.05}, lineColor = {136, 138, 133}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-49.92, 8.83}, {49.92, -8.83}}), Text(origin = {-100, 8}, extent = {{-14, 6}, {14, -6}}, textString = "Tp_IN"), Text(origin = {-100, -50}, extent = {{-14, 6}, {14, -6}}, textString = "Ts_OUT"), Text(origin = {160, -50}, extent = {{-14, 6}, {14, -6}}, textString = "Ts_IN"), Text(origin = {162, 12}, extent = {{-14, 6}, {14, -6}}, textString = "Tp_OUT"), Text(origin = {147, 50}, extent = {{59, -14}, {-59, 14}}, textString = "PHX")}, coordinateSystem(extent = {{-120, 80}, {220, -80}})),
        Icon(graphics = {Rectangle(origin = {30.2864, 0.0340212}, extent = {{-149.242, 59.993}, {149.242, -59.993}}), Rectangle(origin = {-60, 30}, lineColor = {170, 0, 0}, fillColor = {170, 0, 0}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {120, 30}, lineColor = {170, 0, 0}, fillColor = {170, 0, 0}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {0, 30}, lineColor = {170, 0, 0}, fillColor = {170, 0, 0}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {-59.95, -29.61}, lineColor = {0, 85, 255}, fillColor = {0, 85, 255}, fillPattern = FillPattern.Solid, extent = {{-19.63, -9.61}, {19.63, 9.61}}), Rectangle(origin = {-0.12, -29.82}, lineColor = {0, 85, 255}, fillColor = {0, 85, 255}, fillPattern = FillPattern.Solid, extent = {{-19.88, 9.82}, {19.88, -9.82}}), Rectangle(origin = {59.8, -29.85}, lineColor = {0, 85, 255}, fillColor = {0, 85, 255}, fillPattern = FillPattern.Solid, extent = {{-19.8, 9.85}, {19.8, -9.85}}), Rectangle(origin = {119.85, -29.89}, lineColor = {0, 85, 255}, fillColor = {0, 85, 255}, fillPattern = FillPattern.Solid, extent = {{-19.85, 9.89}, {19.85, -9.89}}), Rectangle(origin = {89.92, -1.05}, lineColor = {136, 138, 133}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-49.92, 8.83}, {49.92, -8.83}}), Rectangle(origin = {60, 30}, lineColor = {170, 0, 0}, fillColor = {170, 0, 0}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {-30.08, -1.05}, lineColor = {136, 138, 133}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-49.92, 8.83}, {49.92, -8.83}}), Text(origin = {-100, 8}, extent = {{-14, 6}, {14, -6}}, textString = "Tp_IN"), Text(origin = {-100, -50}, extent = {{-14, 6}, {14, -6}}, textString = "Ts_OUT"), Text(origin = {160, -50}, extent = {{-14, 6}, {14, -6}}, textString = "Ts_IN"), Text(origin = {162, 12}, extent = {{-14, 6}, {14, -6}}, textString = "Tp_OUT"), Text(origin = {147, 50}, extent = {{59, -14}, {-59, 14}}, textString = "PHX")}));
    end PrimaryHeatExchanger;

    model SecondaryHeatExchanger
      // SHX parameters and equations
      // Parameter declaration
      parameter MCFR_HyS_BOP.Units.Mass m_PN1_SHX = 54545;
      // Hot side (coolant salt from PHX) node 1 mass
      parameter MCFR_HyS_BOP.Units.Mass m_PN2_SHX = 54545;
      // Hot side (coolant salt from PHX) node 2 mass
      parameter MCFR_HyS_BOP.Units.Mass m_PN3_SHX = 54545;
      // Hot side (coolant salt from PHX) node 3 mass
      parameter MCFR_HyS_BOP.Units.Mass m_PN4_SHX = 54545;
      // Hot side (coolant salt from PHX) node 4 mass
      parameter MCFR_HyS_BOP.Units.Mass m_TN1_SHX = 4268;
      // SHX Tubing node 1 mass
      parameter MCFR_HyS_BOP.Units.Mass m_TN2_SHX = 4268;
      // SHX Tubing node 2 mass
      parameter MCFR_HyS_BOP.Units.Mass m_SN1_SHX = 20460;
      // Cold Side (Hitec salt from OTSG) node 1 mass
      parameter MCFR_HyS_BOP.Units.Mass m_SN2_SHX = 20460;
      // Cold Side (Hitec salt from OTSG) node 2 mass
      parameter MCFR_HyS_BOP.Units.Mass m_SN3_SHX = 20460;
      // Cold Side (Hitec salt from OTSG) node 3 mass
      parameter MCFR_HyS_BOP.Units.Mass m_SN4_SHX = 20460;
      // Cold Side (Hitec salt from OTSG) node 4 mass
      parameter MCFR_HyS_BOP.Units.SpecificHeatCapacity cP_pFluid_SHX = 0.0011; // See IES paper Table 1 for details
      // Fuel salt specific heat [MJ/(kg*C)]
      parameter MCFR_HyS_BOP.Units.SpecificHeatCapacity cP_Tube_SHX = 0.00057779;
      // PHX Tubing specific heat [MJ/(kg*C)]
      parameter MCFR_HyS_BOP.Units.SpecificHeatCapacity cP_sFluid_SHX = 0.00154912; // Hitec salt specific heat capacity
      // Coolant salt specific heat [MJ/(kg*C)]
      parameter MCFR_HyS_BOP.Units.MassFlowRate m_dot_pFluidNom_SHX;
      // Fuel salt nominal flow rate
      parameter MCFR_HyS_BOP.Units.MassFlowRate m_dot_sFluidNom_SHX;
      // Coolant salt nominal flow rate
      parameter MCFR_HyS_BOP.Units.Convection hApnNom_SHX;
      // PHX primary side heat transfer coefficient [MJ/(s*C)]
      parameter MCFR_HyS_BOP.Units.Convection hAsnNom_SHX;
      // PHX secondary side heat transfer coefficient [MJ/(s*C)]
      parameter MCFR_HyS_BOP.Units.ResidentTime tauSHXtoPHX = 10;
      // Avg travel time from SHX to PHX
      parameter MCFR_HyS_BOP.Units.ResidentTime tauSHXtoOTSG = 10;
      // Avg travel time from SHX to OTSG
      parameter MCFR_HyS_BOP.Units.Temperature T_PN1_SHX_initial = 679.6;
      parameter MCFR_HyS_BOP.Units.Temperature T_PN2_SHX_initial = 670.4;
      parameter MCFR_HyS_BOP.Units.Temperature T_PN3_SHX_initial = 659.6;
      parameter MCFR_HyS_BOP.Units.Temperature T_PN4_SHX_initial = 648.9;
      parameter MCFR_HyS_BOP.Units.Temperature T_TN1_SHX_initial = 649.5;
      parameter MCFR_HyS_BOP.Units.Temperature T_TN2_SHX_initial = 624.5;
      parameter MCFR_HyS_BOP.Units.Temperature T_SN1_SHX_initial = 589.4;
      parameter MCFR_HyS_BOP.Units.Temperature T_SN2_SHX_initial = 605.6;
      parameter MCFR_HyS_BOP.Units.Temperature T_SN3_SHX_initial = 619.4;
      parameter MCFR_HyS_BOP.Units.Temperature T_SN4_SHX_initial = 633.3;
      // Variable declaration
      MCFR_HyS_BOP.Units.ResidentTime varTauSHXtoOTSG;
      MCFR_HyS_BOP.Units.ResidentTime varTauSHXtoPHX;
      MCFR_HyS_BOP.Units.MassFlowRate m_dot_pFluid_SHX;
      MCFR_HyS_BOP.Units.MassFlowRate m_dot_sFluid_SHX;
      MCFR_HyS_BOP.Units.Convection hApn_SHX;
      MCFR_HyS_BOP.Units.Convection hAsn_SHX;
      MCFR_HyS_BOP.Units.Temperature T_in_sFluid_SHX;
      MCFR_HyS_BOP.Units.Temperature T_PN1_SHX;
      MCFR_HyS_BOP.Units.Temperature T_PN2_SHX;
      MCFR_HyS_BOP.Units.Temperature T_PN3_SHX;
      MCFR_HyS_BOP.Units.Temperature T_PN4_SHX;
      MCFR_HyS_BOP.Units.Temperature T_TN1_SHX;
      MCFR_HyS_BOP.Units.Temperature T_TN2_SHX;
      MCFR_HyS_BOP.Units.Temperature T_SN1_SHX;
      MCFR_HyS_BOP.Units.Temperature T_SN2_SHX;
      MCFR_HyS_BOP.Units.Temperature T_SN3_SHX;
      MCFR_HyS_BOP.Units.Temperature T_SN4_SHX;
      MCFR_HyS_BOP.Units.ReactorPower P_TertiaryLoop;
      MCFR_HyS_BOP.Units.ReactorPower P_SecondaryLoop;
      // Connections to other modules
      input MCFR_HyS_BOP.Connectors.Temp_In T_in_pFluid_SHX annotation(
        Placement(visible = true, transformation(origin = {-100, 30}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-102, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input MCFR_HyS_BOP.Connectors.Temp_In T_in_sFluid_OTSG annotation(
        Placement(visible = true, transformation(origin = {160, -30}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {162, -32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input MCFR_HyS_BOP.Connectors.Temp_In T_in_sFluid_HyS annotation(
        Placement(visible = true, transformation(origin = {160, -10}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {162, -20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output MCFR_HyS_BOP.Connectors.Temp_Out T_out_sFluid_SHX annotation(
        Placement(visible = true, transformation(origin = {-100, -30}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-100, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output MCFR_HyS_BOP.Connectors.Temp_Out T_out_pFluid_SHX annotation(
        Placement(visible = true, transformation(origin = {160, 30}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {162, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input MCFR_HyS_BOP.Connectors.FlowFraction_In secondaryFF annotation(
        Placement(visible = true, transformation(origin = {-70, 50}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {-70, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input MCFR_HyS_BOP.Connectors.FlowFraction_In tertiaryFF annotation(
        Placement(visible = true, transformation(origin = {130, -50}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {130, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    initial equation
      T_PN1_SHX = T_PN1_SHX_initial;
      T_PN2_SHX = T_PN2_SHX_initial;
      T_PN3_SHX = T_PN3_SHX_initial;
      T_PN4_SHX = T_PN4_SHX_initial;
      T_TN1_SHX = T_TN1_SHX_initial;
      T_TN2_SHX = T_TN2_SHX_initial;
      T_SN1_SHX = T_SN1_SHX_initial;
      T_SN2_SHX = T_SN2_SHX_initial;
      T_SN3_SHX = T_SN3_SHX_initial;
      T_SN4_SHX = T_SN4_SHX_initial;
    equation
      m_dot_pFluid_SHX = m_dot_pFluidNom_SHX*secondaryFF.FF;
      m_dot_sFluid_SHX = m_dot_sFluidNom_SHX*tertiaryFF.FF;
      varTauSHXtoPHX = tauSHXtoPHX/secondaryFF.FF;
      varTauSHXtoOTSG = tauSHXtoOTSG/tertiaryFF.FF;
// These next two are optional. They attempt to include some complexity of htc changing due to flow rate turbulence
//hApn_SHX = hApnNom_SHX*(0.8215*secondaryFF.FF^6 - 4.108*secondaryFF.FF^5 + 7.848*secondaryFF.FF^4 - 7.165*secondaryFF.FF^3 + 3.004*secondaryFF.FF^2 + 0.5903*secondaryFF.FF + 0.008537);
//hAsn_SHX = hAsnNom_SHX*(0.8215*tertiaryFF.FF^6 - 4.108*tertiaryFF.FF^5 + 7.848*tertiaryFF.FF^4 - 7.165*tertiaryFF.FF^3 + 3.004*tertiaryFF.FF^2 + 0.5903*tertiaryFF.FF + 0.008537);
// These next two are the simpler version for constant htc
      hApn_SHX = hApnNom_SHX;
      hAsn_SHX = hAsnNom_SHX;
      T_in_sFluid_SHX = (T_in_sFluid_OTSG.T + T_in_sFluid_HyS.T)/2;
      
      m_PN1_SHX*cP_pFluid_SHX*der(T_PN1_SHX) = m_dot_pFluid_SHX*cP_pFluid_SHX*(T_in_pFluid_SHX.T - T_PN1_SHX) - hApn_SHX*(T_PN1_SHX - T_TN1_SHX);
      m_PN2_SHX*cP_pFluid_SHX*der(T_PN2_SHX) = m_dot_pFluid_SHX*cP_pFluid_SHX*(T_PN1_SHX - T_PN2_SHX) - hApn_SHX*(T_PN1_SHX - T_TN1_SHX);
      m_PN3_SHX*cP_pFluid_SHX*der(T_PN3_SHX) = m_dot_pFluid_SHX*cP_pFluid_SHX*(T_PN2_SHX - T_PN3_SHX) - hApn_SHX*(T_PN3_SHX - T_TN2_SHX);
      m_PN4_SHX*cP_pFluid_SHX*der(T_PN4_SHX) = m_dot_pFluid_SHX*cP_pFluid_SHX*(T_PN3_SHX - T_PN4_SHX) - hApn_SHX*(T_PN3_SHX - T_TN2_SHX);
      
      m_TN1_SHX*cP_Tube_SHX*der(T_TN1_SHX) = 2*hApn_SHX*(T_PN1_SHX - T_TN1_SHX) - 2*hAsn_SHX*(T_TN1_SHX - T_SN3_SHX);
      m_TN2_SHX*cP_Tube_SHX*der(T_TN2_SHX) = 2*hApn_SHX*(T_PN3_SHX - T_TN2_SHX) - 2*hAsn_SHX*(T_TN2_SHX - T_SN1_SHX);
      
      m_SN1_SHX*cP_sFluid_SHX*der(T_SN1_SHX) = m_dot_sFluid_SHX*cP_sFluid_SHX*(T_in_sFluid_SHX - T_SN1_SHX) + hAsn_SHX*(T_TN2_SHX - T_SN1_SHX);
      m_SN2_SHX*cP_sFluid_SHX*der(T_SN2_SHX) = m_dot_sFluid_SHX*cP_sFluid_SHX*(T_SN1_SHX - T_SN2_SHX) + hAsn_SHX*(T_TN2_SHX - T_SN1_SHX);
      m_SN3_SHX*cP_sFluid_SHX*der(T_SN3_SHX) = m_dot_sFluid_SHX*cP_sFluid_SHX*(T_SN2_SHX - T_SN3_SHX) + hAsn_SHX*(T_TN1_SHX - T_SN3_SHX);
      m_SN4_SHX*cP_sFluid_SHX*der(T_SN4_SHX) = m_dot_sFluid_SHX*cP_sFluid_SHX*(T_SN3_SHX - T_SN4_SHX) + hAsn_SHX*(T_TN1_SHX - T_SN3_SHX); 
           
      T_out_pFluid_SHX.T = delay(T_PN4_SHX, varTauSHXtoPHX, 5000);
      T_out_sFluid_SHX.T = delay(T_SN4_SHX, varTauSHXtoOTSG, 5000);
      P_SecondaryLoop = m_dot_pFluid_SHX*cP_pFluid_SHX*(T_in_pFluid_SHX.T - T_PN4_SHX);
      P_TertiaryLoop = m_dot_sFluid_SHX*cP_sFluid_SHX*(T_SN4_SHX - T_in_sFluid_SHX);
      annotation(
        Diagram(graphics = {Rectangle(origin = {30.2864, 0.0340212}, extent = {{-149.242, 59.993}, {149.242, -59.993}}), Rectangle(origin = {-60, 30}, lineColor = {255, 170, 0}, fillColor = {255, 170, 0}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {120, 30}, lineColor = {255, 170, 0}, fillColor = {255, 170, 0}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {0, 30}, lineColor = {255, 170, 0}, fillColor = {255, 170, 0}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {-59.95, -29.61}, lineColor = {170, 85, 255}, fillColor = {170, 85, 255}, fillPattern = FillPattern.Solid, extent = {{-19.63, -9.61}, {19.63, 9.61}}), Rectangle(origin = {-0.12, -29.82}, lineColor = {170, 85, 255}, fillColor = {170, 85, 255}, fillPattern = FillPattern.Solid, extent = {{-19.88, 9.82}, {19.88, -9.82}}), Rectangle(origin = {59.8, -29.85}, lineColor = {170, 85, 255}, fillColor = {170, 85, 255}, fillPattern = FillPattern.Solid, extent = {{-19.8, 9.85}, {19.8, -9.85}}), Rectangle(origin = {119.85, -29.89}, lineColor = {170, 85, 255}, fillColor = {170, 85, 255}, fillPattern = FillPattern.Solid, extent = {{-19.85, 9.89}, {19.85, -9.89}}), Rectangle(origin = {89.92, -1.05}, lineColor = {136, 138, 133}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-49.92, 8.83}, {49.92, -8.83}}), Rectangle(origin = {60, 30}, lineColor = {255, 170, 0}, fillColor = {255, 170, 0}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {-30.08, -1.05}, lineColor = {136, 138, 133}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-49.92, 8.83}, {49.92, -8.83}}), Text(origin = {-100, 8}, extent = {{-14, 6}, {14, -6}}, textString = "Tp_IN"), Text(origin = {-100, -50}, extent = {{-14, 6}, {14, -6}}, textString = "Ts_OUT"), Text(origin = {160, -50}, extent = {{-14, 6}, {14, -6}}, textString = "Ts_IN"), Text(origin = {162, 12}, extent = {{-14, 6}, {14, -6}}, textString = "Tp_OUT"), Text(origin = {147, 50}, extent = {{59, -14}, {-59, 14}}, textString = "SHX")}, coordinateSystem(extent = {{-120, 80}, {220, -80}})),
        Icon(graphics = {Rectangle(origin = {30.2864, 0.0340212}, extent = {{-149.242, 59.993}, {149.242, -59.993}}), Rectangle(origin = {-60, 30}, lineColor = {255, 170, 0}, fillColor = {255, 170, 0}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {120, 30}, lineColor = {255, 170, 0}, fillColor = {255, 170, 0}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {0, 30}, lineColor = {255, 170, 0}, fillColor = {255, 170, 0}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {-59.95, -29.61}, lineColor = {170, 85, 255}, fillColor = {170, 85, 255}, fillPattern = FillPattern.Solid, extent = {{-19.63, -9.61}, {19.63, 9.61}}), Rectangle(origin = {-0.12, -29.82}, lineColor = {170, 85, 255}, fillColor = {170, 85, 255}, fillPattern = FillPattern.Solid, extent = {{-19.88, 9.82}, {19.88, -9.82}}), Rectangle(origin = {59.8, -29.85}, lineColor = {170, 85, 255}, fillColor = {170, 85, 255}, fillPattern = FillPattern.Solid, extent = {{-19.8, 9.85}, {19.8, -9.85}}), Rectangle(origin = {119.85, -29.89}, lineColor = {170, 85, 255}, fillColor = {170, 85, 255}, fillPattern = FillPattern.Solid, extent = {{-19.85, 9.89}, {19.85, -9.89}}), Rectangle(origin = {89.92, -1.05}, lineColor = {136, 138, 133}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-49.92, 8.83}, {49.92, -8.83}}), Rectangle(origin = {60, 30}, lineColor = {255, 170, 0}, fillColor = {255, 170, 0}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {-30.08, -1.05}, lineColor = {136, 138, 133}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-49.92, 8.83}, {49.92, -8.83}}), Text(origin = {-100, 8}, extent = {{-14, 6}, {14, -6}}, textString = "Tp_IN"), Text(origin = {-100, -50}, extent = {{-14, 6}, {14, -6}}, textString = "Ts_OUT"), Text(origin = {160, -50}, extent = {{-14, 6}, {14, -6}}, textString = "Ts_IN"), Text(origin = {162, 12}, extent = {{-14, 6}, {14, -6}}, textString = "Tp_OUT"), Text(origin = {147, 50}, extent = {{59, -14}, {-59, 14}}, textString = "SHX")}),
        experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002));
    end SecondaryHeatExchanger;

    model OnceThroughSteamGenerator
      // OTSG parameters and equations
      // Parameter declaration
      parameter MCFR_HyS_BOP.Units.Density rho_hitec = 1680; // From http://www.coal2nuclear.com/MSR%20-%20HITEC%20Heat%20Transfer%20Salt.pdf
      parameter MCFR_HyS_BOP.Units.Density rho_tube = 8425.712;
      parameter MCFR_HyS_BOP.Units.Density rho_SC = 778.8648; // Density of subcooled fluid in neighbor node to saturated boiling
      parameter MCFR_HyS_BOP.Units.Density rho_sh = 45.05219; // Density of the superheated steam
      parameter MCFR_HyS_BOP.Units.SpecificHeatCapacity cP_hitec = 0.00154912; // Hitec salt specific heat [MJ/(kg*C)]
      parameter MCFR_HyS_BOP.Units.SpecificHeatCapacity cP_Tube_OTSG = 456.056e-6; // OTSG Tubing specific heat [MJ/(kg*C)]
      parameter MCFR_HyS_BOP.Units.SpecificHeatCapacity cP_SH1 = 0.002573333673413; // Specific heat approximation for the steam outlet
      parameter MCFR_HyS_BOP.Units.SpecificHeatCapacity cP_SH2 = 0.003687517297276; // Specific heat approximation for the steam between boiling and outlet
      parameter MCFR_HyS_BOP.Units.SpecificHeatCapacity cP_fdw = 0.004393935709417; // Specific heat approximation of subcooled fluid at OTSG inlet
      parameter MCFR_HyS_BOP.Units.SpecificHeatCapacity cP_SC = 0.005121564867136; // Specific heat of subcooled fluid in neighbor node to saturated boiling
      parameter MCFR_HyS_BOP.Units.MassFlowRate m_dot_hitecNom;// = 12910.56;//3379; // Hitec salt nominal flow rate
      parameter MCFR_HyS_BOP.Units.MassFlowRate m_dot_fdwNom = 1.925004320027648e+02;             // Feedwater nominal flow rate
      parameter MCFR_HyS_BOP.Units.Convection h_pw = 9101.343578935E-6; // OTSG primary side to tubes heat transfer coefficient [MJ/(s*C)]
      parameter MCFR_HyS_BOP.Units.Convection h_wsh = 5732.378387768E-6; // OTSG tubes to superheated steam heat transfer coefficient [MJ/(s*C)]
      parameter MCFR_HyS_BOP.Units.Convection h_wb = 13349.334646671E-6; // OTSG tubes to saturated boiling heat transfer coefficient [MJ/(s*C)]
      parameter MCFR_HyS_BOP.Units.Convection h_wsc = 8385.005556375E-6; // OTSG tubes to subcooled fluid heat transfer coefficient [MJ/(s*C)]
      parameter MCFR_HyS_BOP.Units.ResidentTime tauOTSGtoSHX;           // Avg travel time from PHX to SHX
      parameter MCFR_HyS_BOP.Units.Temperature T_PN1_OTSG_initial = 633;//605.5;
      parameter MCFR_HyS_BOP.Units.Temperature T_PN2_OTSG_initial = 614;//579.6;
      parameter MCFR_HyS_BOP.Units.Temperature T_PN3_OTSG_initial = 603;//557.5;
      parameter MCFR_HyS_BOP.Units.Temperature T_PN4_OTSG_initial = 592;//537.3;
      parameter MCFR_HyS_BOP.Units.Temperature T_PN5_OTSG_initial = 584;//523.9;
      parameter MCFR_HyS_BOP.Units.Temperature T_PN6_OTSG_initial = 575;//507.3;
      parameter MCFR_HyS_BOP.Units.Temperature T_TN1_OTSG_initial = 632;//583.4;
      parameter MCFR_HyS_BOP.Units.Temperature T_TN2_OTSG_initial = 605;//542.1;
      parameter MCFR_HyS_BOP.Units.Temperature T_TN3_OTSG_initial = 432;//408.7;
      parameter MCFR_HyS_BOP.Units.Temperature T_TN4_OTSG_initial = 428;//401.0;
      parameter MCFR_HyS_BOP.Units.Temperature T_TN5_OTSG_initial = 454;//416.0;
      parameter MCFR_HyS_BOP.Units.Temperature T_TN6_OTSG_initial = 410;//368.8;
      parameter MCFR_HyS_BOP.Units.Temperature T_SH1_OTSG_initial = 630;//594.7;
      parameter MCFR_HyS_BOP.Units.Temperature T_SH2_OTSG_initial = 592;//537.6;
      parameter MCFR_HyS_BOP.Units.Temperature T_sat_OTSG_initial = 327.8;
      parameter MCFR_HyS_BOP.Units.Temperature T_SC2_OTSG_initial = 251.2;
      parameter MCFR_HyS_BOP.Units.Temperature X_5 = 2.469804064618581e+02; // Saturated boiling temperature at zero pressure [C]
      parameter MCFR_HyS_BOP.Units.TsatConstant K_5 = 6.473415732352336; // Saturated boiling temperature change per pressure [C/Psat]
      parameter MCFR_HyS_BOP.Units.X_4 X_4 = 1.942396453493191; // Hfg = Hg - Hf (Specific enthalpy as a saturated gas/vapor minus specific enthalpy as a saturated fluid)
      parameter MCFR_HyS_BOP.Units.K_4 K_4 = -0.062356611803191; // Ptable polyfitted with Hfg_table
      parameter MCFR_HyS_BOP.Units.Number K_sc = 17.615718797133468;//
      parameter MCFR_HyS_BOP.Units.Number K_1 = 11.632248097704855;// Necessary for the pressure changes. See Vikram Singh's thesis for more information.
      parameter MCFR_HyS_BOP.Units.Number K_b = 12.834854076292217;
      parameter MCFR_HyS_BOP.Units.Length L_sc_initial = 0.3;//1.4469*(23)^(1/3); //Length of the subcooled domain
      parameter MCFR_HyS_BOP.Units.Length L_b_initial = 0.3;//2.3645*(23)^(1/3); //Length of the boiling domain
      parameter MCFR_HyS_BOP.Units.Length L_sh_initial = 9.4;//4.723*(23)^(1/3); //Length of the superheated domain
      parameter MCFR_HyS_BOP.Units.Length L_OTSG = 10;//8.5344*(23)^(1/3); //Total OTSG length. Doesn't change
      parameter MCFR_HyS_BOP.Units.Length D_inner = 0.0141478; // Inner diameter of the tubes
      parameter MCFR_HyS_BOP.Units.Length D_outer = 0.015875; // Outer diameter of the tubes
      parameter MCFR_HyS_BOP.Units.Area A_flow = 1.2245; // Cross-sectional area of flow on secondary side [m^2]
      parameter MCFR_HyS_BOP.Units.Number N_tubes = 10500;//6546*(23);
      parameter MCFR_HyS_BOP.Units.Number pi = 2*Modelica.Math.asin(1.0);
      parameter MCFR_HyS_BOP.Units.EnthalpyChange dHs_dPs = -0.042433622114357; // Superheated steam Specific enthalpy change per pressure change
      parameter MCFR_HyS_BOP.Units.CompressibilityFactor Z_ss = 0.76634; // Compressibility factor at 570K and 60 atm
      parameter MCFR_HyS_BOP.Units.UniversalGasConstant R_ugc = 8.3114462e-6; // Universal gas constant [MJ/(mol*C)]
      parameter MCFR_HyS_BOP.Units.MolarMass MM_stm = 0.018; // Molar mass of steam
      parameter MCFR_HyS_BOP.Units.Pressure P_setpoint = 12.5; // Nominal secondary side pressure [MPa]
      // Variable declaration
      MCFR_HyS_BOP.Units.ResidentTime varTauOTSGtoSHX;
      MCFR_HyS_BOP.Units.MassFlowRate m_dot_hitec;
      MCFR_HyS_BOP.Units.MassFlowRate m_dot_SC;
      MCFR_HyS_BOP.Units.MassFlowRate m_dot_B;
      MCFR_HyS_BOP.Units.Temperature T_PN1;
      MCFR_HyS_BOP.Units.Temperature T_PN2;
      MCFR_HyS_BOP.Units.Temperature T_PN3;
      MCFR_HyS_BOP.Units.Temperature T_PN4;
      MCFR_HyS_BOP.Units.Temperature T_PN5;
      MCFR_HyS_BOP.Units.Temperature T_PN6;
      MCFR_HyS_BOP.Units.Temperature T_TN1;
      MCFR_HyS_BOP.Units.Temperature T_TN2;
      MCFR_HyS_BOP.Units.Temperature T_TN3;
      MCFR_HyS_BOP.Units.Temperature T_TN4;
      MCFR_HyS_BOP.Units.Temperature T_TN5;
      MCFR_HyS_BOP.Units.Temperature T_TN6;
      MCFR_HyS_BOP.Units.Temperature T_SH2;
      MCFR_HyS_BOP.Units.Temperature T_sat;
      MCFR_HyS_BOP.Units.Temperature T_SC2;
      MCFR_HyS_BOP.Units.Density rho_b;
      MCFR_HyS_BOP.Units.Length L_sc;
      MCFR_HyS_BOP.Units.Length L_b;
      MCFR_HyS_BOP.Units.Length L_sh;
      MCFR_HyS_BOP.Units.Mass M_sh;
      MCFR_HyS_BOP.Units.ReactorPower P_OTSG_Primary;
      // Connections to other modules
      input MCFR_HyS_BOP.Connectors.Temp_In T_in_hitec_OTSG annotation(
        Placement(visible = true, transformation(origin = {-176, 30}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-176, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output MCFR_HyS_BOP.Connectors.Temp_Out T_out_hitec_OTSG annotation(
        Placement(visible = true, transformation(origin = {154, 30}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {154, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input MCFR_HyS_BOP.Connectors.FlowFraction_In tertiaryFF annotation(
        Placement(visible = true, transformation(origin = {-176, 58}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {-176, 58}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output MCFR_HyS_BOP.Connectors.Pressure P_stm annotation(
        Placement(visible = true, transformation(origin = {-8, -64}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-8, -66}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input MCFR_HyS_BOP.Connectors.FlowRate_In m_dot_fdw annotation(
        Placement(visible = true, transformation(origin = {-176, -62}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-176, -66}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output MCFR_HyS_BOP.Connectors.FlowRate_Out m_dot_SH annotation(
        Placement(visible = true, transformation(origin = {154, -62}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {154, -66}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input MCFR_HyS_BOP.Connectors.Temp_In T_fdw annotation(
        Placement(visible = true, transformation(origin = {154, -30}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {154, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output MCFR_HyS_BOP.Connectors.Temp_Out T_SH1 annotation(
        Placement(visible = true, transformation(origin = {-176, -30}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-176, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    initial equation
// Necessary: Initial lengths (L_sc, L_b, and L_sh), initial temps, initial pressures
      T_PN1 = T_PN1_OTSG_initial;
      T_PN2 = T_PN2_OTSG_initial;
      T_PN3 = T_PN3_OTSG_initial;
      T_PN4 = T_PN4_OTSG_initial;
      T_PN5 = T_PN5_OTSG_initial;
      T_PN6 = T_PN6_OTSG_initial;
      T_TN1 = T_TN1_OTSG_initial;
      T_TN2 = T_TN2_OTSG_initial;
      T_TN3 = T_TN3_OTSG_initial;
      T_TN4 = T_TN4_OTSG_initial;
      T_TN5 = T_TN5_OTSG_initial;
      T_TN6 = T_TN6_OTSG_initial;
      T_SH1.T = T_SH1_OTSG_initial;
      T_SH2 = T_SH2_OTSG_initial;
      T_SC2 = T_SC2_OTSG_initial;
      L_sc = L_sc_initial;
      L_b = L_b_initial;
      L_sh = L_sh_initial;
      P_stm.P = P_setpoint;
      M_sh = A_flow * rho_sh * L_sh_initial;
    equation
      m_dot_hitec = m_dot_hitecNom*tertiaryFF.FF;
      varTauOTSGtoSHX = tauOTSGtoSHX/tertiaryFF.FF;
// Primary Side Nodes
      der(T_PN1)*(N_tubes*pi*(D_inner^2)*0.25*L_sh*0.5*rho_hitec*cP_hitec) = ((m_dot_hitec*cP_hitec*(T_in_hitec_OTSG.T - T_PN1)) + (h_pw*N_tubes*pi*D_inner*L_sh*0.5*(T_TN1 - T_PN1))); // Primary side node #1 temperature eqn
      der(T_PN2)*(N_tubes*pi*(D_inner^2)*0.25*L_sh*0.5*rho_hitec*cP_hitec) = ((m_dot_hitec*cP_hitec*(T_PN1 - T_PN2)) + (h_pw*N_tubes*pi*D_inner*L_sh*0.5*(T_TN2 - T_PN2))); // Primary side node #2 temperature eqn
      der(T_PN3)*(N_tubes*pi*(D_inner^2)*0.25*L_b*0.5*rho_hitec*cP_hitec) = ((m_dot_hitec*cP_hitec*(T_PN2 - T_PN3)) + (h_pw*N_tubes*pi*D_inner*L_b*0.5*(T_TN3 - T_PN3))); // Primary side node #3 temperature eqn
      der(T_PN4)*(N_tubes*pi*(D_inner^2)*0.25*L_b*0.5*rho_hitec*cP_hitec) = ((m_dot_hitec*cP_hitec*(T_PN3 - T_PN4)) + (h_pw*N_tubes*pi*D_inner*L_b*0.5*(T_TN4 - T_PN4))); // Primary side node #4 temperature eqn
      der(T_PN5)*(N_tubes*pi*(D_inner^2)*0.25*L_sc*0.5*rho_hitec*cP_hitec) = ((m_dot_hitec*cP_hitec*(T_PN4 - T_PN5)) + (h_pw*N_tubes*pi*D_inner*L_sc*0.5*(T_TN5 - T_PN5))); // Primary side node #5 temperature eqn
      der(T_PN6)*(N_tubes*pi*(D_inner^2)*0.25*L_sc*0.5*rho_hitec*cP_hitec) = ((m_dot_hitec*cP_hitec*(T_PN5 - T_PN6)) + (h_pw*N_tubes*pi*D_inner*L_sc*0.5*(T_TN6 - T_PN6))); // Primary side node #6 temperature eqn
      T_out_hitec_OTSG.T = delay(T_PN6, varTauOTSGtoSHX, 5000);
//// Tube Nodes
      der(T_TN1)*(N_tubes*pi*(D_outer^2 - D_inner^2)*0.25*L_sh*0.5*rho_tube*cP_Tube_OTSG) = (h_pw*(N_tubes*pi*D_inner*L_sh*0.5)*(T_PN1 - T_TN1) + h_wsh*(N_tubes*pi*D_outer*L_sh*0.5)*(T_SH1.T - T_TN1)); // Tubing node #1 temperature eqn
      der(T_TN2)*(N_tubes*pi*(D_outer^2 - D_inner^2)*0.25*L_sh*0.5*rho_tube*cP_Tube_OTSG) = (h_pw*(N_tubes*pi*D_inner*L_sh*0.5)*(T_PN2 - T_TN2) + h_wsh*(N_tubes*pi*D_outer*L_sh*0.5)*(T_SH2 - T_TN2)); // Tubing node #2 temperature eqn
      der(T_TN3)*(N_tubes*pi*(D_outer^2 - D_inner^2)*0.25*L_b*0.5*rho_tube*cP_Tube_OTSG) = (h_pw*(N_tubes*pi*D_inner*L_b*0.5)*(T_PN3 - T_TN3) + h_wb*(N_tubes*pi*D_outer*L_b*0.5)*(T_sat - T_TN3)); // Tubing node #3 temperature eqn
      der(T_TN4)*(N_tubes*pi*(D_outer^2 - D_inner^2)*0.25*L_b*0.5*rho_tube*cP_Tube_OTSG) = (h_pw*(N_tubes*pi*D_inner*L_b*0.5)*(T_PN4 - T_TN4) + h_wb*(N_tubes*pi*D_outer*L_b*0.5)*(T_sat - T_TN4)); // Tubing node #4 temperature eqn
      der(T_TN5)*(N_tubes*pi*(D_outer^2 - D_inner^2)*0.25*L_sc*0.5*rho_tube*cP_Tube_OTSG) = (h_pw*(N_tubes*pi*D_inner*L_sc*0.5)*(T_PN5 - T_TN5) + h_wsc*(N_tubes*pi*D_outer*L_sc*0.5)*(T_sat - T_TN5)); // Tubing node #5 temperature eqn
      der(T_TN6)*(N_tubes*pi*(D_outer^2 - D_inner^2)*0.25*L_sc*0.5*rho_tube*cP_Tube_OTSG) = (h_pw*(N_tubes*pi*D_inner*L_sc*0.5)*(T_PN6 - T_TN6) + h_wsc*(N_tubes*pi*D_outer*L_sc*0.5)*(T_SC2 - T_TN6)); // Tubing node #6 temperature eqn
//
//// Secondary Side Nodes
      der(T_SH1.T)*(M_sh*0.5*cP_SH1) = (h_wsh*(N_tubes*pi*D_outer*L_sh*0.5)*(T_TN1 - T_SH1.T) + (m_dot_SH.mdot*cP_SH1*(T_SH2 - T_SH1.T)) + (A_flow*L_sh*0.5 + dHs_dPs)*der(P_stm.P)); // Superheated steam node #1 (outlet) temperature eqn
    der(T_SH2)*(M_sh*0.5*cP_SH2) = (h_wsh*(N_tubes*pi*D_outer*L_sh*0.5)*(T_TN2 - T_SH2) + (m_dot_B*cP_SH2*(T_sat - T_SH2)) + (A_flow*L_sh*0.5 - dHs_dPs)*der(P_stm.P));
// Superheated steam node #2 temperature eqn
      T_sat = X_5 + K_5*P_stm.P; // Saturated boiling temp eqn
        der(T_SC2) = ((h_wsc*(N_tubes*pi*D_outer*L_sc*0.5)*(T_TN6-T_SC2) + (cP_fdw*m_dot_fdw.mdot*T_fdw.T) - cP_SC*((m_dot_fdw.mdot+m_dot_SC)*0.5*T_SC2) + A_flow*L_sc*0.5*der(P_stm.P)) * (2/(A_flow * rho_SC *cP_SC)) - T_SC2*der(L_sc))*(1/(L_sc)); // Subcooled fluid node temperature eqn
// Secondary side pressure:
      der(P_stm.P)*A_flow*L_sh = ((Z_ss*R_ugc)*(der(M_sh)*(T_SH1.T + T_SH2) + M_sh*(der(T_SH1.T) + der(T_SH2)))/(2*MM_stm)) - (P_stm.P*A_flow*der(L_sh)); // Secondary side pressure. T_sh1 is OTSG outlet. Assumption: P_subcooled = P_saturatedboiling = P_superheated
//der(P_stm.P)*A_flow*L_sh = ((Z_ss*R_ugc)*(der(M_sh)*(T_SH2) + M_sh*(der(T_SH2)))/(MM_stm)) - (P_stm.P*A_flow*der(L_sh)); //Alternative method
//P_stm.P = 12.5; // Default pressure. Uncomment this and comment the der(P_stm.P) eqn to set pressure constant.
//  der(P_stm) = 0; //Alternative
//=======
// Density, length, and mass flow rates:
      rho_b = 25.884875705 + 12.835*P_stm.P; // Two-phase boiling density eqn. See Singh_2020 (https://www.sciencedirect.com/science/article/pii/S0029549319304881)
      der(L_sc) = (1/(A_flow * 0.5 * rho_SC * cP_SC * (T_SC2 - T_sat))) * (h_wsc * N_tubes * pi * D_outer * L_sc * 0.5 * (T_TN5 - T_sat) + cP_SC * (m_dot_fdw.mdot * T_SC2 - (m_dot_fdw.mdot + m_dot_SC) * 0.5 * T_sat) + A_flow * L_sc * 0.5 * (K_sc * cP_SC * (T_SC2 - 2 * T_sat) - 1) * der(P_stm.P) - (A_flow * rho_SC * cP_SC * K_1 * L_sc * 0.5 * der(P_stm.P))); // Subcooled domain length eqn
      der(L_b)*A_flow*rho_b = m_dot_SC - m_dot_B - (A_flow*L_b*K_b*der(P_stm.P)); // Boiling domain length eqn
      L_sh = L_OTSG - L_b - L_sc; // Total OTSG length eqn
//  der(L_sh) = -der(L_b) - der(L_sc); // Superheated domain length eqn
      m_dot_SH.mdot = m_dot_fdw.mdot * P_stm.P / P_setpoint; // Flow rate out of superheated domain (and OTSG)
      m_dot_SC = m_dot_fdw.mdot - (rho_SC*A_flow*der(L_sc)) - (A_flow*L_sc*K_sc*der(P_stm.P)); // Flow rate out of subcooled domain
      m_dot_B = (h_wb*(N_tubes*pi*D_outer*L_b)*(T_TN3*0.5 + T_TN4*0.5 - T_sat))/(X_4 + K_4*P_stm.P); // Flow rate out of boiling domain
//  m_dot_SC = 1.925004320027648e+02;
//  m_dot_B = 1.925004320027648e+02;
      der(M_sh) = m_dot_B - m_dot_SH.mdot; // Rate of steam mass change
      P_OTSG_Primary = m_dot_hitec*cP_hitec*(T_in_hitec_OTSG.T - T_PN6);
    annotation(
        Diagram(graphics = {Rectangle(origin = {-5, -2}, extent = {{-185, 74}, {185, -74}}), Rectangle(origin = {-128, 30}, lineColor = {236, 127, 211}, fillColor = {236, 127, 211}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {-128, -30}, lineColor = {127, 236, 152}, fillColor = {127, 236, 152}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {-80, 30}, lineColor = {236, 127, 211}, fillColor = {236, 127, 211}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {-32, 30}, lineColor = {236, 127, 211}, fillColor = {236, 127, 211}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {16, 30}, lineColor = {236, 127, 211}, fillColor = {236, 127, 211}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {64, 30}, lineColor = {236, 127, 211}, fillColor = {236, 127, 211}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {112, 30}, lineColor = {236, 127, 211}, fillColor = {236, 127, 211}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {-80, -30}, lineColor = {127, 236, 152}, fillColor = {127, 236, 152}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {16, -30}, lineColor = {127, 236, 152}, fillColor = {127, 236, 152}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {-32, -30}, lineColor = {127, 236, 152}, fillColor = {127, 236, 152}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {112, -30}, lineColor = {127, 236, 152}, fillColor = {127, 236, 152}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {64, -30}, lineColor = {127, 236, 152}, fillColor = {127, 236, 152}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {-128, 0}, lineColor = {136, 138, 133}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {-80, 0}, lineColor = {136, 138, 133}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {16, 0}, lineColor = {136, 138, 133}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {-32, 0}, lineColor = {136, 138, 133}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {112, 0}, lineColor = {136, 138, 133}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {64, 0}, lineColor = {136, 138, 133}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Text(origin = {-10, 58}, extent = {{-38, 16}, {38, -16}}, textString = "OTSG")}, coordinateSystem(extent = {{-200, 80}, {220, -80}})),
        Icon(graphics = {Rectangle(origin = {-5, -2}, extent = {{-185, 74}, {185, -74}}), Rectangle(origin = {-128, 30}, lineColor = {236, 127, 211}, fillColor = {236, 127, 211}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {-128, -30}, lineColor = {127, 236, 152}, fillColor = {127, 236, 152}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {-80, 30}, lineColor = {236, 127, 211}, fillColor = {236, 127, 211}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {-32, 30}, lineColor = {236, 127, 211}, fillColor = {236, 127, 211}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {16, 30}, lineColor = {236, 127, 211}, fillColor = {236, 127, 211}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {64, 30}, lineColor = {236, 127, 211}, fillColor = {236, 127, 211}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {112, 30}, lineColor = {236, 127, 211}, fillColor = {236, 127, 211}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {-80, -30}, lineColor = {127, 236, 152}, fillColor = {127, 236, 152}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {16, -30}, lineColor = {127, 236, 152}, fillColor = {127, 236, 152}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {-32, -30}, lineColor = {127, 236, 152}, fillColor = {127, 236, 152}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {112, -30}, lineColor = {127, 236, 152}, fillColor = {127, 236, 152}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {64, -30}, lineColor = {127, 236, 152}, fillColor = {127, 236, 152}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {-128, 0}, lineColor = {136, 138, 133}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {-80, 0}, lineColor = {136, 138, 133}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {16, 0}, lineColor = {136, 138, 133}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {-32, 0}, lineColor = {136, 138, 133}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {112, 0}, lineColor = {136, 138, 133}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {64, 0}, lineColor = {136, 138, 133}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Text(origin = {-10, 58}, extent = {{-38, 16}, {38, -16}}, textString = "OTSG")}));
    end OnceThroughSteamGenerator;

    model RBOP
    //======================
      //PARAMETERS:
    parameter MCFR_HyS_BOP.Units.Pressure P_LP = 1.15;
    parameter MCFR_HyS_BOP.Units.Pressure P_co = 0.008;
    parameter MCFR_HyS_BOP.Units.Pressure P_atm = 0.101325;
    parameter MCFR_HyS_BOP.Units.Number K_RH = 0.2;//0.16; //Ratio of steam taken from steam valve for reheat loop
    parameter MCFR_HyS_BOP.Units.Number K_HPB = 0.02;//0.1634; //Ratio of high pressure steam bled before turbine for feedwater heater
    parameter MCFR_HyS_BOP.Units.Number K_LPB_initial = 0.12;//0.2174; //Ratio of low pressure steam bled before turbine for feedwater heater. Controller changes this
    parameter MCFR_HyS_BOP.Units.Temperature T_setpoint = 180;//212; //Desired feedwater temp
    parameter MCFR_HyS_BOP.Units.Temperature T_cool_in = 22;//Temp of cooling water for condenser
    parameter MCFR_HyS_BOP.Units.MassFlowRate m_dot_cool_nom;
    parameter MCFR_HyS_BOP.Units.MassFlowRate m_dot_fdw_nom;
    parameter MCFR_HyS_BOP.Units.Time rampstart;
    parameter MCFR_HyS_BOP.Units.Number ramp_min;
    //May be up to 27700. Source: https://www.nuclear-power.com/nuclear-power-plant/turbine-generator-power-conversion-system/main-condenser-steam-condenser/
    parameter MCFR_HyS_BOP.Units.ResidentTime TauFDWtoOTSG = 10;
    parameter MCFR_HyS_BOP.Units.FlowFraction FF_fdw = 1;
    parameter MCFR_HyS_BOP.Units.Frequency omega;
    parameter MCFR_HyS_BOP.Units.InitiationTime sinstart;
    parameter MCFR_HyS_BOP.Units.Number sinmag;
    parameter MCFR_HyS_BOP.Units.Time cooltriptime;
    parameter MCFR_HyS_BOP.Units.Number turbine_efficiency = 0.9;
      //======================
      //VARIABLES:
    MCFR_HyS_BOP.Units.MassFlowRate m_dot_v;
    MCFR_HyS_BOP.Units.MassFlowRate m_dot_HPT_in;
    MCFR_HyS_BOP.Units.MassFlowRate m_dot_HPT_bleed;
    MCFR_HyS_BOP.Units.MassFlowRate m_dot_HPT_out;
    MCFR_HyS_BOP.Units.SpecificEnthalpy H_HPT_out;
    MCFR_HyS_BOP.Units.SpecificEnthalpy H_dHPT;
    MCFR_HyS_BOP.Units.SpecificEnthalpy H_HPT_bleed;
    MCFR_HyS_BOP.Units.HeatFlowRate Q_HPT;
    MCFR_HyS_BOP.Units.Temperature T_HPT_out;
    //===
    MCFR_HyS_BOP.Units.Temperature T_sat_MS;
    MCFR_HyS_BOP.Units.SpecificEnthalpy H_vap_MS;
    MCFR_HyS_BOP.Units.SpecificEnthalpy H_f_MS;
    MCFR_HyS_BOP.Units.Number X_stm;
    MCFR_HyS_BOP.Units.MassFlowRate m_dot_liq;
    MCFR_HyS_BOP.Units.MassFlowRate m_dot_MS_out;
    MCFR_HyS_BOP.Units.SpecificEnthalpy H_MS;
    MCFR_HyS_BOP.Units.Temperature T_MS;
    //===
    MCFR_HyS_BOP.Units.MassFlowRate m_dot_RH;
    MCFR_HyS_BOP.Units.MassFlowRate m_dot_vRH;
    MCFR_HyS_BOP.Units.HeatFlowRate Q_max_HP;
    MCFR_HyS_BOP.Units.HeatFlowRate Q_max_LP;
    MCFR_HyS_BOP.Units.HeatFlowRate Q_RH;
    MCFR_HyS_BOP.Units.SpecificEnthalpy H_RH;
    MCFR_HyS_BOP.Units.Temperature T_RH;
    MCFR_HyS_BOP.Units.SpecificEnthalpy H_vRH;
    MCFR_HyS_BOP.Units.Temperature T_vRH;
      //===
    MCFR_HyS_BOP.Units.MassFlowRate m_dot_LPT_in;
    MCFR_HyS_BOP.Units.MassFlowRate m_dot_LPT_bleed;
    MCFR_HyS_BOP.Units.MassFlowRate m_dot_LPT_out;
    MCFR_HyS_BOP.Units.SpecificEnthalpy H_LPT_out;
    MCFR_HyS_BOP.Units.SpecificEnthalpy H_dLPT;
    MCFR_HyS_BOP.Units.SpecificEnthalpy H_LPT_bleed;
    MCFR_HyS_BOP.Units.HeatFlowRate Q_LPT;
    MCFR_HyS_BOP.Units.Temperature T_LPT_out;
    //===
    MCFR_HyS_BOP.Units.Temperature T_sat_con;
    MCFR_HyS_BOP.Units.SpecificEnthalpy H_f_con;
    MCFR_HyS_BOP.Units.SpecificEnthalpy H_vap_con;
    MCFR_HyS_BOP.Units.Number X_inlet_con;
    MCFR_HyS_BOP.Units.HeatFlowRate Q_con;
    MCFR_HyS_BOP.Units.Temperature T_stm_in;
    MCFR_HyS_BOP.Units.SpecificEnthalpy H_cool_in;
    MCFR_HyS_BOP.Units.SpecificEnthalpy H_limit;
    MCFR_HyS_BOP.Units.HeatFlowRate Q_cool_max;
    MCFR_HyS_BOP.Units.SpecificEnthalpy H_cool_out;
    MCFR_HyS_BOP.Units.Temperature T_cool_out;
    MCFR_HyS_BOP.Units.SpecificEnthalpy H_con_out;
    MCFR_HyS_BOP.Units.HeatFlowRate Q_cool;
    MCFR_HyS_BOP.Units.MassFlowRate m_dot_con;
    MCFR_HyS_BOP.Units.MassFlowRate m_dot_cool;
    //===
    MCFR_HyS_BOP.Units.SpecificEnthalpy H_setpoint;
    MCFR_HyS_BOP.Units.SpecificEnthalpy H_fdw;
    MCFR_HyS_BOP.Units.Mass M_fdw;
    MCFR_HyS_BOP.Units.Number K_LPB_check;
    MCFR_HyS_BOP.Units.Number K_LPB;
    MCFR_HyS_BOP.Units.ResidentTime varTauFDWtoOTSG;
    MCFR_HyS_BOP.Units.Temperature T_feed;
    //===
    MCFR_HyS_BOP.Units.SpecificEnthalpy H_HPT_in;
    MCFR_HyS_BOP.Units.MassFlowRate m_dot_fdw2;
    MCFR_HyS_BOP.Units.HeatFlowRate Q_turbine;
  //======================
      input MCFR_HyS_BOP.Connectors.Pressure P_HP annotation(
        Placement(visible = true, transformation(origin = {0, 60}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input MCFR_HyS_BOP.Connectors.FlowRate_In m_dot_OTSG annotation(
        Placement(visible = true, transformation(origin = {60, 60}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {60, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output MCFR_HyS_BOP.Connectors.Temp_Out T_fdw annotation(
        Placement(visible = true, transformation(origin = {-60, -60}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-60, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input MCFR_HyS_BOP.Connectors.Temp_In T_OTSG annotation(
        Placement(visible = true, transformation(origin = {-60, 60}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-60, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output MCFR_HyS_BOP.Connectors.FlowRate_Out m_dot_fdw annotation(
        Placement(visible = true, transformation(origin = {60, -60}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {60, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    initial equation
    H_HPT_in = Water.specificEnthalpy_pT(P_HP.P*1e6, T_OTSG.T + 273.15)/1e6;
    P_HP.P =12.5;
    equation
       H_HPT_in = Water.specificEnthalpy_pT(P_HP.P*1e6, T_OTSG.T + 273.15)/1e6;
       varTauFDWtoOTSG = TauFDWtoOTSG/FF_fdw;
//======================
//VALVE AND NOZZLE CHEST EQNS:
       m_dot_v = (m_dot_OTSG.mdot+13.929976)*K_RH;
//======================
//HIGH PRESSURE TURBINE EQNS:
       m_dot_HPT_in = (m_dot_OTSG.mdot+13.929976) - m_dot_v;
       m_dot_HPT_bleed = K_HPB*m_dot_HPT_in;
       m_dot_HPT_out = m_dot_HPT_in - m_dot_HPT_bleed;
       H_HPT_out = Water.isentropicEnthalpy(P_LP*1e6, Water.setState_phX(P_HP.P*1e6, H_HPT_in*1e6))/1e6;
       H_dHPT = H_HPT_in - H_HPT_out;
       H_HPT_bleed = H_HPT_in;
       Q_HPT = H_dHPT*m_dot_HPT_out;
       T_HPT_out = Water.temperature_ph(P_LP*1e6, H_HPT_out*1e6) - 273.15;
//======================
//MOISTURE SEPARATOR EQNS:
       T_sat_MS = MoistAir.saturationTemperature(P_LP*1e6, 220, 600) - 273.15;
       H_vap_MS = MoistAir.enthalpyOfVaporization(T_sat_MS + 273.15)/1e6;
       H_f_MS = Water.specificEnthalpy_pT(P_LP*1e6, T_sat_MS + 273.15)/1e6;
//H_liq is H_f
       X_stm = (H_HPT_out - H_f_MS)/H_vap_MS;
    if X_stm >= 1 then
       m_dot_liq = 0;
       m_dot_MS_out = m_dot_HPT_out;
       H_MS = H_HPT_out;
    else
       m_dot_MS_out = X_stm*m_dot_HPT_out;
       m_dot_liq = m_dot_HPT_out - m_dot_MS_out;
       H_MS = H_f_MS + H_vap_MS;
//H_HPT_out/X_stm;
    end if;
       T_MS = Water.temperature_ph(P_LP*1e6, H_MS*1e6) - 273.15;
//======================
//REHEATER EQNS:
       m_dot_RH = m_dot_MS_out;
       m_dot_vRH = m_dot_v;
// ~~~~~~~~~~~~~
      if T_MS < 322 then
        Q_max_HP = m_dot_vRH*(Water.specificEnthalpy_pT(P_HP.P*1e6, T_OTSG.T + 273.15)/1e6 - Water.specificEnthalpy_pT(P_HP.P*1e6, (327.82) + 273.15)/1e6) + 1.4*((322 - T_MS)/(T_OTSG.T - T_MS))*(T_vRH - T_MS);
      else
        Q_max_HP = m_dot_vRH*(Water.specificEnthalpy_pT(P_HP.P*1e6, T_OTSG.T + 273.15)/1e6 - Water.specificEnthalpy_pT(P_HP.P*1e6, (T_MS + 5.82) + 273.15)/1e6);
    end if;
       Q_max_LP = m_dot_RH*(Water.specificEnthalpy_pT(P_LP*1e6, (T_OTSG.T - 10) + 273.15)/1e6 - Water.specificEnthalpy_pT(P_LP*1e6, T_MS + 273.15)/1e6);
    if Q_max_HP < Q_max_LP then
        Q_RH = Q_max_HP;
        H_RH = Water.specificEnthalpy_pT(P_LP*1e6, T_MS + 273.15)/1e6 + Q_RH/m_dot_RH;
        T_RH = Water.temperature_ph(P_LP*1e6, H_RH*1e6) - 273.15;
        H_vRH = Water.specificEnthalpy_pT(P_HP.P*1e6, T_OTSG.T + 273.15)/1e6 - Q_RH/m_dot_vRH;
        T_vRH = Water.temperature_ph(P_HP.P*1e6, H_vRH*1e6) - 273.15;
    else
        Q_RH = Q_max_LP;
        H_RH = Water.specificEnthalpy_pT(P_LP*1e6, T_MS + 273.15)/1e6 + Q_RH/m_dot_RH;
        T_RH = Water.temperature_ph(P_LP*1e6, H_RH*1e6) - 273.15;
        H_vRH = Water.specificEnthalpy_pT(P_HP.P*1e6, T_OTSG.T + 273.15)/1e6 - Q_RH/m_dot_vRH;
        T_vRH = Water.temperature_ph(P_HP.P*1e6, H_vRH*1e6) - 273.15;
    end if;
//Explanation: The limits of the reheater are that the reheated outlet can be heated nearly to the temp of the valve side inlet, or the valve side outlet can be cooled nearly down to the temp of the regeated side inlet. The reheated side cannot be heated beyond the valve side max temp nor can the valve side transfer heat when below the secondary side's minimum temp due to Newton's law of cooling.
//======================
//LOW PRESSURE TURBINE EQNS:
       m_dot_LPT_in = m_dot_RH;
       m_dot_LPT_bleed = K_LPB*m_dot_LPT_in;
       m_dot_LPT_out = m_dot_LPT_in - m_dot_LPT_bleed;
       H_LPT_out = Water.isentropicEnthalpy(P_co*1e6, Water.setState_phX(P_LP*1e6, H_RH*1e6))/1e6;
       H_dLPT = H_RH - H_LPT_out;
       H_LPT_bleed = H_RH;
       Q_LPT = H_dLPT*m_dot_LPT_out;
       T_LPT_out = Water.temperature_ph(P_co*1e6, H_LPT_out*1e6) - 273.15;
       Q_turbine = (Q_HPT + Q_LPT)*turbine_efficiency;
//======================
//CONDENSER EQNS:
       m_dot_cool = m_dot_cool_nom*((1 - 0.01)*exp(-0.092*delay(time, cooltriptime)) + 0.01);
       T_sat_con = MoistAir.saturationTemperature(P_co*1e6, 220, 600) - 273.15;
       H_f_con = Water.specificEnthalpy_pT(P_co*1e6, T_sat_con + 273.15)/1e6;
       H_vap_con = MoistAir.enthalpyOfVaporization(T_stm_in + 273.15)/1e6;
       X_inlet_con = (H_LPT_out - H_f_con)/H_vap_con;
       Q_con = m_dot_LPT_out*(H_LPT_out - H_f_con);
//Energy demand
       T_stm_in = Water.temperature_ph(P_co*1e6, H_LPT_out*1e6) - 273.15;
       H_cool_in = Water.specificEnthalpy_pT(P_atm*1e6, T_cool_in + 273.15)/1e6;
       H_limit = Water.specificEnthalpy_pT(P_atm*1e6, T_stm_in - 5 + 273.15)/1e6;
       Q_cool_max = m_dot_cool*(H_limit - H_cool_in);
    if Q_con < Q_cool_max then
        H_cool_out = (Q_con/m_dot_cool) + H_cool_in;
        T_cool_out = Water.temperature_ph(P_atm*1e6, H_cool_out*1e6) - 273.15;
        H_con_out = H_f_con;
    else
        H_cool_out = (Q_cool_max/m_dot_cool) + H_cool_in;
        T_cool_out = Water.temperature_ph(P_atm*1e6, H_cool_out*1e6) - 273.15;
        H_con_out = H_LPT_out - (Q_cool_max/m_dot_con);
    end if;
       Q_cool = (H_cool_out - H_cool_in)*m_dot_cool;
       m_dot_con = m_dot_LPT_out;
//Method of determining outlet flow based on enthalpy. Hot inlet and outlet enthalpy known, along with coolant inlet. Use coolant dH_co*mdot_in as demand, which is satisfied by the difference between coolant outlet (unknown) and inlet, which can meet it up to a certain fraction of the Hot inlet. This determines the flow rate as a ratio of how well the demand is satisfied times the inlet flow rate.
//======================
//FEEDWATER HEATER EQNS:
//   m_dot_fdw.mdot = m_dot_fdw_nom;//388.963;//(m_dot_HPT_bleed + m_dot_LPT_bleed + m_dot_liq + m_dot_vRH + m_dot_con);
       m_dot_fdw2 = (m_dot_HPT_bleed + m_dot_LPT_bleed + m_dot_liq + m_dot_vRH + m_dot_con);
       m_dot_fdw.mdot = m_dot_fdw_nom * (1 + delay(sinmag*sin(omega*time),sinstart)) * (1 + (if (time < rampstart) then 0 else if (time < rampstart + 500) then -(1-ramp_min)*(time - rampstart)/500 else if (time < rampstart + 1500) then -(1-ramp_min) else if (time < rampstart + 2000) then -(1-ramp_min)*(rampstart + 2000 - time)/500 else 0));
       H_setpoint = Water.specificEnthalpy_pT(P_HP.P*1e6, T_setpoint + 273.15)/1e6;
       K_LPB_check = (m_dot_fdw.mdot/(H_LPT_bleed*m_dot_LPT_in))*(H_setpoint - ((H_HPT_bleed*m_dot_HPT_bleed/m_dot_fdw.mdot) + (H_f_MS*m_dot_liq/m_dot_fdw.mdot) + (H_vRH*m_dot_vRH/m_dot_fdw.mdot) + (H_con_out*m_dot_con/m_dot_fdw.mdot)));
    if K_LPB_check > 0 then
        K_LPB = K_LPB_check;
    elseif K_LPB_check > 1 then
        K_LPB = 1;
    else
       K_LPB = 0;
    end if;
       T_feed = Water.temperature_ph(P_HP.P*1e6, H_fdw*1e6) - 273.15;
       T_fdw.T = delay(T_feed, varTauFDWtoOTSG, 5000);
       H_fdw = (H_HPT_bleed*m_dot_HPT_bleed/m_dot_fdw.mdot) + (H_LPT_bleed*m_dot_LPT_bleed/m_dot_fdw.mdot) + (H_f_MS*m_dot_liq/m_dot_fdw.mdot) + (H_vRH*m_dot_vRH/m_dot_fdw.mdot) + (H_con_out*m_dot_con/m_dot_fdw.mdot);
       der(M_fdw) = m_dot_HPT_bleed + m_dot_LPT_bleed + m_dot_liq + m_dot_vRH + m_dot_con - m_dot_fdw.mdot;
//======================
    annotation(
        experiment(StartTime = 0, StopTime = 5000, Tolerance = 1e-06, Interval = 0.1),
        Diagram(graphics = {Rectangle(extent = {{-80, 80}, {80, -80}}), Text(origin = {1, 29}, extent = {{-43, 15}, {43, -15}}, textString = "RBOP"), Polygon(lineColor = {49, 49, 49}, fillColor = {37, 116, 243}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, lineThickness = 1.25, points = {{0, 0}, {-60, -40}, {-60, 40}, {0, 0}, {60, 40}, {60, -40}, {60, -40}, {60, -40}, {0, 0}, {0, 0}, {0, 0}})}),
        Icon(graphics = {Polygon(lineColor = {49, 49, 49}, fillColor = {37, 116, 243}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, lineThickness = 1.25, points = {{0, 0}, {-60, -40}, {-60, 40}, {0, 0}, {60, 40}, {60, -40}, {60, -40}, {60, -40}, {0, 0}, {0, 0}, {0, 0}}), Text(origin = {1, 39}, extent = {{-43, 15}, {43, -15}}, textString = "RBOP"), Rectangle(extent = {{-80, 80}, {80, -80}})}));    
    end RBOP;

    model HyS
    // Parameter declarations
      // ===============
      // DECOMPOSISITION
      parameter MCFR_HyS_BOP.Units.MassFlowRate mdot_H_nom;// = 5400;// Mass flow rate of hot fluid on primary side
      parameter MCFR_HyS_BOP.Units.SpecificHeatCapacity Cp_H = 0.00154912; // Hitec salt specific heat [MJ/(kg*C)]
      parameter MCFR_HyS_BOP.Units.SpecificHeatCapacity Cp_W = 7.5e-4;// Specific heat of the tubing material (4.35e-4 for Inconel)
      parameter MCFR_HyS_BOP.Units.Mass M_H_node = 5000;// Mass of fluid in one of the hot fluid nodes [kg]
      parameter MCFR_HyS_BOP.Units.Mass M_W_node = 500;    //0.0118;// Mass of one of the wall nodes [kg]
      parameter MCFR_HyS_BOP.Units.Area A1_outer = 10;
      parameter MCFR_HyS_BOP.Units.Area A2_outer = 3.668;
      parameter MCFR_HyS_BOP.Units.Area A3_outer = 0.402;
      parameter MCFR_HyS_BOP.Units.HeatTransferCoeff h_htc = 0.036;//1.2e-4;// Heat transfer coefficient
      parameter MCFR_HyS_BOP.Units.Temperature T_boil = 227;    // Saturated boiling temperature of secondary side [C]
      parameter MCFR_HyS_BOP.Units.EnthalpyOfVaporization Hvap_SA = 0.073331;// Energy required to fully boil one mol of sulfuric acid (SA)
      parameter MCFR_HyS_BOP.Units.EnthalpyOfVaporization Hvap_H2O = 0.035970;// Energy required to fully boil one mol of water
      parameter MCFR_HyS_BOP.Units.UniversalGasConstant R = 0.00831447/1e3;// Universal Gas Constant [MJ/(mol*C)]
      parameter MCFR_HyS_BOP.Units.Temperature T_ref_rxn = 298;//[K]
      parameter MCFR_HyS_BOP.Units.Temperature T_ref_Kc = 973;//[K]
      parameter MCFR_HyS_BOP.Units.Number epsilon = 0.25;
      parameter MCFR_HyS_BOP.Units.Number Kc = 0.041;//0.041;
      parameter MCFR_HyS_BOP.Units.pressure_kPa pressure_nom;// = 100;// 1500// Secondary side pressure [kPa]
      parameter MCFR_HyS_BOP.Units.HeatOfReaction H_rxn1 = 0.09753;//[MJ/mol]
      parameter MCFR_HyS_BOP.Units.HeatOfReaction H_rxn2 = 0.09889;//0.09893;// [MJ/mol]
      parameter MCFR_HyS_BOP.Units.MolarMass MM_SA = 0.0980785;// [kg/mol]
      parameter MCFR_HyS_BOP.Units.MolarMass MM_H2O = 0.01801528;
      parameter MCFR_HyS_BOP.Units.MolarMass MM_SO3 = 0.0800632;
      parameter MCFR_HyS_BOP.Units.MolarMass MM_SO2 = 0.0640638;
      parameter MCFR_HyS_BOP.Units.MolarMass MM_O2 = 0.0319988;
      parameter MCFR_HyS_BOP.Units.FlowFraction FF_Decomposition = 1;
      parameter MCFR_HyS_BOP.Units.FlowFraction FF_SA_recycle = 1;
      parameter MCFR_HyS_BOP.Units.Length pipe_thickness = 0.025;
      parameter MCFR_HyS_BOP.Units.Time TauSArecycle_nom = 30;
      parameter MCFR_HyS_BOP.Units.Frequency omega;
      parameter MCFR_HyS_BOP.Units.InitiationTime sinstart;
      parameter MCFR_HyS_BOP.Units.Number sinmag;
    // ============
      // ELECTROLYSIS
      parameter MCFR_HyS_BOP.Units.Temperature T_elec = 127;
      parameter MCFR_HyS_BOP.Units.MolarMass MM_H2 = 0.002016;
      parameter MCFR_HyS_BOP.Units.Number z = 2;// Number of electrons
      parameter MCFR_HyS_BOP.Units.FaradayConstant F = 96485.3399;// Coulombs per mol of electrons
      parameter MCFR_HyS_BOP.Units.Volt E0_cell = 0.6;// From Gorensek
      parameter MCFR_HyS_BOP.Units.WeightPercent WtP_SA_nom = 0.9;
      parameter MCFR_HyS_BOP.Units.MolFlowRate Mdot_H2O_nom;//Nominal excess water flow rate
      parameter MCFR_HyS_BOP.Units.MolFlowRate Mdot_SA_nom;//Nominal net sulfur flow rate
      // ~~~
      parameter MCFR_HyS_BOP.Units.Time rampstart;
      parameter MCFR_HyS_BOP.Units.Number ramp_min;
      parameter MCFR_HyS_BOP.Units.ResidentTime TauHyStoSHX = 10;
      parameter MCFR_HyS_BOP.Units.FlowFraction FF_HyS = 1;
      parameter MCFR_HyS_BOP.Units.Time sulfouttime;
      parameter MCFR_HyS_BOP.Units.Time prelieffail;
    // Variable declaration
      // ===============
      // DECOMPOSISITION
      MCFR_HyS_BOP.Units.Temperature T_vap;
      MCFR_HyS_BOP.Units.Temperature T_H1;
      MCFR_HyS_BOP.Units.Temperature T_H2;
      MCFR_HyS_BOP.Units.Temperature T_H3;
      MCFR_HyS_BOP.Units.MassFlowRate mdot_H;
      MCFR_HyS_BOP.Units.MassFlowRate mdot_S;
      MCFR_HyS_BOP.Units.Temperature T_W1;
      MCFR_HyS_BOP.Units.Temperature T_W2;
      MCFR_HyS_BOP.Units.Temperature T_W3;
      MCFR_HyS_BOP.Units.SpecificHeatCapacity Cp_SA;
      MCFR_HyS_BOP.Units.SpecificHeatCapacity Cp_SA_liq;
      MCFR_HyS_BOP.Units.SpecificHeatCapacity Cp_SO3;
      MCFR_HyS_BOP.Units.SpecificHeatCapacity Cp_SO2;
      MCFR_HyS_BOP.Units.SpecificHeatCapacity Cp_H2O;
      MCFR_HyS_BOP.Units.SpecificHeatCapacity Cp_H2O_liq;
      MCFR_HyS_BOP.Units.SpecificHeatCapacity Cp_O2;
      MCFR_HyS_BOP.Units.SpecificHeatCapacity Cp_liq;
      MCFR_HyS_BOP.Units.SpecificHeatCapacity Cp_vap;
      MCFR_HyS_BOP.Units.Power E_dot_boil;
      MCFR_HyS_BOP.Units.Conversion X;
      MCFR_HyS_BOP.Units.HeatOfReaction Q_rxn1;
      MCFR_HyS_BOP.Units.HeatOfReaction Q_rxn2;
      MCFR_HyS_BOP.Units.HeatOfReaction Q_tot;
      MCFR_HyS_BOP.Units.HeatOfReaction Q_tot_CT;
      MCFR_HyS_BOP.Units.HeatOfReaction dH_Rx;
      MCFR_HyS_BOP.Units.Temperature T_liq_in;
      MCFR_HyS_BOP.Units.Time vTauSArecycle;
      MCFR_HyS_BOP.Units.SpecificHeatCapacity Cp_SO3_CT;
      MCFR_HyS_BOP.Units.SpecificHeatCapacity Cp_SO2_CT;
      MCFR_HyS_BOP.Units.SpecificHeatCapacity Cp_O2_CT;
      MCFR_HyS_BOP.Units.HeatOfReaction Q_rxn2_CT;
      MCFR_HyS_BOP.Units.Conversion X_CT;
      MCFR_HyS_BOP.Units.Power P_HyS_Primary;
      MCFR_HyS_BOP.Units.ResidentTime varTauHyStoSHX;
      // ============
      // ELECTROLYSIS
      MCFR_HyS_BOP.Units.Temperature delta_T_cond;
      MCFR_HyS_BOP.Units.StandardEntropy S_SA;
      MCFR_HyS_BOP.Units.StandardEntropy S_H2;
      MCFR_HyS_BOP.Units.StandardEntropy S_SO2;
      MCFR_HyS_BOP.Units.StandardEntropy S_H2O;
      MCFR_HyS_BOP.Units.EnthalpyOfReaction dH0;
      MCFR_HyS_BOP.Units.Power P_demand;
      // ============
      // MASS BALANCE
      MCFR_HyS_BOP.Units.MolFlowRate Mdot_SA_rcy;
      MCFR_HyS_BOP.Units.MassFlowRate mdot_O2_out;
      MCFR_HyS_BOP.Units.MolFlowRate Mdot_SO2_out;
      MCFR_HyS_BOP.Units.MolFlowRate Mdot_H2O_Dout;
      MCFR_HyS_BOP.Units.MolFlowRate Mdot_rxn;
      MCFR_HyS_BOP.Units.MassFlowRate mdot_H2_out;
      MCFR_HyS_BOP.Units.MassFlowRate mdot_H2O_makeup;
      MCFR_HyS_BOP.Units.WeightPercent WtP_SA;
      MCFR_HyS_BOP.Units.WeightPercent WtP_SO2;
      MCFR_HyS_BOP.Units.Power Heat_used;
      MCFR_HyS_BOP.Units.MolFlowRate Mdot_H2O_Eout;
      MCFR_HyS_BOP.Units.MolFlowRate Mdot_SA_in;
      MCFR_HyS_BOP.Units.pressure_kPa pressure;
    // ===============
      // CONNECTIONS
      input MCFR_HyS_BOP.Connectors.Temp_In T_H_in annotation(
        Placement(visible = true, transformation(origin = {-60, 60}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-60, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input MCFR_HyS_BOP.Connectors.FlowFraction_In FF_HyS_Hot annotation(
        Placement(visible = true, transformation(origin = {0, 60}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {0, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output MCFR_HyS_BOP.Connectors.Temp_Out T_hot_out annotation(
        Placement(visible = true, transformation(origin = {-60, -60}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-60, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    initial equation
      T_H1 = T_H_in.T;
      T_H2 = T_H_in.T - 25;
      T_H3 = T_H_in.T - 50;
      T_W1 = T_H_in.T - 10;
      T_W2 = T_H_in.T - 20;
      T_W3 = T_H_in.T - 30;
    equation
      pressure = pressure_nom /((1 - 0.11)*exp(-0.092*delay(time, prelieffail)) + 0.11);
// ===============
// DECOMPOSISITION
      T_liq_in = T_elec;
      vTauSArecycle = TauSArecycle_nom/FF_SA_recycle;
      varTauHyStoSHX = TauHyStoSHX/FF_HyS;
// =========
// Hot nodes (primary side)
      M_H_node*Cp_H*der(T_H1) = mdot_H*Cp_H*(T_H_in.T - T_H1) - h_htc*A1_outer/pipe_thickness*(((T_H_in.T + T_H1)/2) - T_W1);
      M_H_node*Cp_H*der(T_H2) = mdot_H*Cp_H*(T_H1 - T_H2) - h_htc*A2_outer/pipe_thickness*(((T_H1 + T_H2)/2) - T_W2);
      M_H_node*Cp_H*der(T_H3) = mdot_H*Cp_H*(T_H2 - T_H3) - h_htc*A3_outer/pipe_thickness*(((T_H2 + T_H3)/2) - T_W3);
      T_hot_out.T = delay(T_H3, varTauHyStoSHX, 5000);
// ==========
// Wall nodes (tubing)
      M_W_node*Cp_W*der(T_W1) = h_htc*A1_outer/pipe_thickness*(((T_H_in.T + T_H1)/2) - T_W1) - mdot_S*Cp_vap*(T_vap - T_boil) - Mdot_rxn*(Q_rxn1 + Q_rxn2);
      M_W_node*Cp_W*der(T_W2) = h_htc*A2_outer/pipe_thickness*(((T_H1 + T_H2)/2) - T_W2) - E_dot_boil;
      M_W_node*Cp_W*der(T_W3) = h_htc*A3_outer/pipe_thickness*(((T_H2 + T_H3)/2) - T_W3) - mdot_S*Cp_liq*(T_boil - T_liq_in);
      E_dot_boil = Mdot_SA_in*Hvap_SA + Mdot_H2O_Eout*Hvap_H2O;
      Heat_used = mdot_S*Cp_liq*(T_boil - T_liq_in) + E_dot_boil + mdot_S*Cp_vap*(T_vap - T_boil) + Mdot_rxn*(Q_rxn1 + Q_rxn2);
// ==========
// Cold nodes (secondary side)
      T_vap = T_H1 - 2*(T_H1 - T_W1);//Equation simplifies down to this. See the T_W1 eqn for the energy usage.
      Q_rxn1 = H_rxn1 + (Cp_H2O*MM_H2O + Cp_SO3*MM_SO3 - Cp_SA*MM_SA)*((T_vap + 273.15) - T_ref_rxn);
      Q_rxn2 = H_rxn2 + (-Cp_SO3*MM_SO3 + Cp_SO2*MM_SO2 + 0.5*Cp_O2*MM_O2)*((T_vap + 273.15) - T_ref_rxn);
      Q_tot = Q_rxn1 + Q_rxn2 + (mdot_S*Cp_liq*(T_boil - T_liq_in)/Mdot_SO2_out) + (E_dot_boil/Mdot_SO2_out) + (Cp_vap*mdot_S*(T_vap - T_boil)/Mdot_SO2_out);    //HEAT USED PER MOL OF SO2 LIBERATED
// ==========================
// Sulfur Trioxide Conversion
//  dH_Rx = (H_rxn2 + (-Cp_SO3*MM_SO3 + Cp_SO2*MM_SO2 + 0.5*Cp_O2*MM_O2)*((T_vap + 273.15) - T_ref_rxn))*(1e3);//[kJ/mol] CHECK UNITS
      dH_Rx  = H_rxn2 + (-Cp_SO3*MM_SO3 + Cp_SO2*MM_SO2 + 0.5*Cp_O2*MM_O2)*((T_vap + 273.15) - T_ref_rxn);    //[MJ/mol]
      Kc*exp((dH_Rx/R)*(((1/(T_ref_Kc))) - (1/(T_vap + 273.15)))) = ((4.86e-5*pressure*(X/(1 + (epsilon*X)))*(1173/(T_vap + 273.15)))*(4.86e-5*pressure*((0.5*X)/(1 + (epsilon*X)))*(1173/(T_vap + 273.15)))^(0.5))/(4.86e-5*pressure*((1 - X)/(1 + (epsilon*X))));
// ===========
// CT Calloways's method for calculating the specific heat capacities (just for comparison, not used elsewhere). The specific heats are per mol not per kg, like the Shomate ones.
      Cp_SO3_CT = ((8.06) + ((1.056e-3)*(T_vap + 273.15)) - (2.028e5)/((T_vap + 273.15)^2))*R;
      Cp_SO2_CT = ((5.699) + ((8.01e-4)*(T_vap + 273.15)) - (1.015e5)/((T_vap + 273.15)^2))*R;
      Cp_O2_CT = ((3.693) + ((5.06e-4)*(T_vap + 273.15)) - (2.27e4)/((T_vap + 273.15)^2))*R;
      Q_rxn2_CT = H_rxn2 + (-Cp_SO3_CT + Cp_SO2_CT + 0.5*Cp_O2_CT)*((T_vap + 273.15) - T_ref_rxn);
      Kc*exp((Q_rxn2_CT/R)*(((1/(T_ref_Kc))) - (1/(T_vap + 273.15)))) = ((4.86e-5*pressure*(X_CT/(1 + (epsilon*X_CT)))*(1173/(T_vap + 273.15)))*(4.86e-5*pressure*((0.5*X_CT)/(1 + (epsilon*X_CT)))*(1173/(T_vap + 273.15)))^(0.5))/(4.86e-5*pressure*((1 - X_CT)/(1 + (epsilon*X_CT))));
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
      Mdot_H2O_Eout = Mdot_H2O_nom * (1 + delay(sinmag*sin(omega*time),sinstart)) * (1 + (if (time < rampstart) then 0 else if (time < rampstart + 500) then -(1-ramp_min)*(time - rampstart)/500 else if (time < rampstart + 1500) then -(1-ramp_min) else if (time < rampstart + 2000) then -(1-ramp_min)*(rampstart + 2000 - time)/500 else 0));
      Mdot_SA_in = Mdot_SA_nom * ((1 - 0.01)*exp(-0.092*delay(time, sulfouttime)) + 0.01) * (1 + (if (time < rampstart) then 0 else if (time < rampstart + 500) then -(1-ramp_min)*(time - rampstart)/500 else if (time < rampstart + 1500) then -(1-ramp_min) else if (time < rampstart + 2000) then -(1-ramp_min)*(rampstart + 2000 - time)/500 else 0));
// ~~~
      mdot_H = mdot_H_nom*FF_HyS_Hot.FF;//Hot salt flow rate
      mdot_S = Mdot_SA_in*MM_SA + Mdot_H2O_Eout*MM_H2O;
// ~~~
      Mdot_SA_rcy = Mdot_SA_in*(1 - X);
      Mdot_rxn = Mdot_SA_in*X;
      mdot_O2_out = 0.5*Mdot_rxn*MM_O2;
      Mdot_SO2_out = Mdot_rxn;
      Mdot_H2O_Dout = Mdot_H2O_Eout + Mdot_rxn;
// ~~~
      mdot_H2_out = Mdot_rxn*MM_H2;
      mdot_H2O_makeup = 2*Mdot_rxn*MM_H2O;//Add a flow fraction here to modulate the makeup water feed rate
      WtP_SA = Mdot_SA_in*MM_SA/(Mdot_SA_in*MM_SA + Mdot_H2O_Eout*MM_H2O);
      WtP_SO2 = Mdot_SO2_out*MM_SO2/(Mdot_SO2_out*MM_SO2 + Mdot_H2O_Dout*MM_H2O);
      P_HyS_Primary = mdot_H_nom*Cp_H*(T_H_in.T - T_H3);
    annotation(
        Diagram(graphics = {Text(origin = {0, 23}, extent = {{-56, 21}, {56, -21}}, textString = "HyS"), Rectangle(extent = {{-80, 80}, {80, -80}}), Rectangle(origin = {0, -20}, fillColor = {255, 63, 63}, fillPattern = FillPattern.HorizontalCylinder, lineThickness = 1, extent = {{-60, 20}, {60, -20}}, radius = 10)}),
        Icon(graphics = {Text(origin = {0, 23}, extent = {{-56, 21}, {56, -21}}, textString = "HyS"), Rectangle(origin = {0, -20}, fillColor = {255, 63, 63}, fillPattern = FillPattern.HorizontalCylinder, lineThickness = 1, extent = {{-60, 20}, {60, -20}}, radius = 10), Rectangle(extent = {{-80, 80}, {80, -80}})}));
end HyS;

    model UHX
      parameter MCFR_HyS_BOP.Units.Mass m_PN1UHX;
      parameter MCFR_HyS_BOP.Units.ReactorPower UHXP;
      parameter MCFR_HyS_BOP.Units.SpecificHeatCapacity cP_pFluidUHX = 0.00154912;
      parameter MCFR_HyS_BOP.Units.MassFlowRate m_dot_pFluidUHXnom = 12910.56;
      parameter MCFR_HyS_BOP.Units.ResidentTime tauUHXtoSHX = 1;
      parameter MCFR_HyS_BOP.Units.InitiationTime tripUHX;
      parameter MCFR_HyS_BOP.Units.Demand SetDemand;
      MCFR_HyS_BOP.Units.Temperature UHXtemp;
      MCFR_HyS_BOP.Units.ResidentTime varTauUHXtoSHX;
      MCFR_HyS_BOP.Units.Demand PowerDemand;
      MCFR_HyS_BOP.Units.MassFlowRate m_dot_pFluidUHX;
      input MCFR_HyS_BOP.Connectors.Temp_In UHXtemp_In annotation(
        Placement(visible = true, transformation(origin = {-70, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-70, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output MCFR_HyS_BOP.Connectors.Temp_Out UHXtemp_Out annotation(
        Placement(visible = true, transformation(origin = {26, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {26, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input MCFR_HyS_BOP.Connectors.FlowFraction_In tertiaryFF annotation(
        Placement(visible = true, transformation(origin = {-70, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-70, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    initial equation
//
    equation
      PowerDemand = if time < tripUHX then 1 else SetDemand;
      varTauUHXtoSHX = tauUHXtoSHX / tertiaryFF.FF;
      m_dot_pFluidUHX = m_dot_pFluidUHXnom * tertiaryFF.FF;
//m_PN1UHX * cP_pFluidUHX * der(UHXtemp) = m_dot_pFluidUHX * cP_pFluidUHX * (UHXtemp_In.T - UHXtemp) - UHXP * PowerDemand;
      UHXtemp = UHXtemp_In.T - ((UHXP * PowerDemand)/(m_dot_pFluidUHX*cP_pFluidUHX));
      UHXtemp_Out.T = delay(UHXtemp, varTauUHXtoSHX, 5000);
      annotation(
        Diagram(graphics = {Rectangle(origin = {-20, 0}, extent = {{-60, 60}, {60, -60}}), Text(origin = {-64, -14}, extent = {{-14, 12}, {14, -12}}, textString = "Temp In"), Text(origin = {24, -14}, extent = {{-14, 12}, {14, -12}}, textString = "Temp Out"), Text(origin = {8, 44}, extent = {{-34, 18}, {34, -18}}, textString = "UHX")}),
        Icon(graphics = {Text(origin = {8, 44}, extent = {{-34, 18}, {34, -18}}, textString = "UHX"), Rectangle(origin = {-20, 0}, extent = {{-60, 60}, {60, -60}}), Text(origin = {24, -14}, extent = {{-14, 12}, {14, -12}}, textString = "Temp Out"), Text(origin = {-64, -14}, extent = {{-14, 12}, {14, -12}}, textString = "Temp In")}));
    end UHX;
    replaceable package Water = Modelica.Media.Water.StandardWater;
    replaceable package MoistAir = Modelica.Media.Air.MoistAir;
  end HeatTransport;

  package Pumps
    model PrimaryLoopPump
      // Fuel Salt Pump parameters and equations
      parameter MCFR_HyS_BOP.Units.FlowFraction freeConvectionFF;
      parameter MCFR_HyS_BOP.Units.PumpConstant primaryPumpK;
      parameter MCFR_HyS_BOP.Units.InitiationTime tripPrimaryPump;
      output MCFR_HyS_BOP.Connectors.FlowFraction_Out primaryFF annotation(
        Placement(visible = true, transformation(origin = {0, -40}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {0, -58}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    initial equation
      primaryFF.FF = 1;
    equation
      primaryFF.FF = (1 - freeConvectionFF)*exp(-primaryPumpK*delay(time, tripPrimaryPump)) + freeConvectionFF;
      annotation(
        Diagram(graphics = {Rectangle(lineColor = {245, 121, 0}, lineThickness = 1, extent = {{-60, 60}, {60, -60}}), Text(origin = {2, 6}, extent = {{-42, 16}, {42, -16}}, textString = "Primary Pump")}),
        Icon(graphics = {Rectangle(origin = {0, -20}, lineColor = {245, 121, 0}, lineThickness = 1, extent = {{-60, 60}, {60, -60}}), Text(origin = {4, -20}, extent = {{-42, 16}, {42, -16}}, textString = "Primary Pump")}));
    end PrimaryLoopPump;

    model SecondaryLoopPump
      // Coolant Salt Pump parameters and equations
      parameter MCFR_HyS_BOP.Units.FlowFraction freeConvectionFF;
      parameter MCFR_HyS_BOP.Units.PumpConstant secondaryPumpK;
      parameter MCFR_HyS_BOP.Units.InitiationTime tripSecondaryPump;
      output MCFR_HyS_BOP.Connectors.FlowFraction_Out secondaryFF annotation(
        Placement(visible = true, transformation(origin = {-1.77636e-15, -40}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {1.77636e-15, -60}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    initial equation
      secondaryFF.FF = 1;
    equation
      secondaryFF.FF = (1 - freeConvectionFF)*exp(-secondaryPumpK*delay(time, tripSecondaryPump)) + freeConvectionFF;
      annotation(
        Diagram(graphics = {Rectangle(lineColor = {0, 85, 0}, lineThickness = 1, extent = {{-60, 60}, {60, -60}}), Text(origin = {2, 6}, extent = {{-42, 16}, {42, -16}}, textString = "Secondary Pump")}),
        Icon(graphics = {Rectangle(origin = {0, -20}, lineColor = {245, 121, 0}, lineThickness = 1, extent = {{-60, 60}, {60, -60}}), Text(origin = {4, -20}, extent = {{-42, 16}, {42, -16}}, textString = "Secondary Pump")}));
    end SecondaryLoopPump;

    model TertiaryLoopPump
      // Hitec Salt Pump parameters and equations
      parameter MCFR_HyS_BOP.Units.FlowFraction freeConvectionFF;
      parameter MCFR_HyS_BOP.Units.PumpConstant tertiaryPumpK;
      parameter MCFR_HyS_BOP.Units.InitiationTime tripTertiaryPump;
      output MCFR_HyS_BOP.Connectors.FlowFraction_Out tertiaryFF annotation(
        Placement(visible = true, transformation(origin = {-1.77636e-15, -40}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {1.77636e-15, -60}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    initial equation
      tertiaryFF.FF = 1;
    equation
      tertiaryFF.FF = (1 - freeConvectionFF)*exp(-tertiaryPumpK*delay(time, tripTertiaryPump)) + freeConvectionFF;
      annotation(
        Diagram(graphics = {Rectangle(lineColor = {85, 85, 255}, lineThickness = 1, extent = {{-60, 60}, {60, -60}}), Text(origin = {2, 6}, extent = {{-42, 16}, {42, -16}}, textString = "Tertiary Pump")}),
        Icon(graphics = {Rectangle(origin = {0, -20}, lineColor = {245, 121, 0}, lineThickness = 1, extent = {{-60, 60}, {60, -60}}), Text(origin = {4, -20}, extent = {{-42, 16}, {42, -16}}, textString = "Tertiary Pump")}));
    end TertiaryLoopPump;
  end Pumps;

  package Units
    // Defines all units in the model
    type Frequency = Real(unit = "rad/s");
    // Used in: Nuclear > mPKE
    type NeutronDensity = Real(unit = "1/m^3", min = 0);
    // Neutron Density [#n/m^3]
    // Used in:
    type NominalNeutronPopulation = Real(unit = "1", min = 0);
    // Nominal Neutron Population [#n]
    // Used in:
    type NeutronGenerationTime = Real(unit = "s", min = 0);
    // Neutron Generation Time (Lambda) [s]
    // Used in:
    type Reactivity = Real(unit = "1");
    // Absolute Reactivity [unitless]
    // Used in:
    type PrecursorDecayConstant = Real(unit = "1/s", min = 0);
    // Delayed Precursor Decay Constant [s^-1]
    // Used in:
    type PrecursorConc = Real(unit = "1", min = 0);
    // Delayed Neutron Precursor
    // Used in:
    type PrecursorReturn = Real(unit = "1/s", min = 0);
    // Used in:
    type DelayedNeutronFrac = Real(unit = "1", min = 0);
    // Delayed Neutron Fraction [unitless]
    // Used in:
    type VolumeImportance = Real(unit = "1", min = 0);
    // Fraction of Fission Power Deposited in Core Nodes [m^3]
    // Used in: Nuclear > Temperature feedback
    type TemperatureReactivityCoef = Real(unit = "1/C");
    // Absolute Temperature Reactivity Coef [C^-1]
    // Used in: Nuclear > Poisons
    type IsotopicConc = Real(unit = "1/m^3", min = 0);
    // Isotopic Concentration [#atoms/m^3]
    // Used in:
    type IsotopicDecayConstant = Real(unit = "1/s", min = 0);
    // Isotopic Decay Constant [s^-1]
    // Used in:
    type IsotopicFissYield = Real(unit = "1", min = 0);
    // Isotopic Yield per fission
    // Used in: Nuclear > Poisons
    type NominalNeutronFlux = Real(unit = "1/(cm^2.s)");
    // OLD: "n/(cm^2.s)"
    // Used in:
    type MassFlowRate = Real(unit = "kg/s", min = 0);
    // Mass Flow Rate [kg/m^3]
    // Used in:
    type ResidentTime = Real(unit = "s", min = 0);
    // Fuel Salt Resident Time (tau) [s]
    // Used in:
    type FlowFraction = Real(unit = "1", min = 0);
    // Flow Fraction [unitless]
    // Used in:
    //type SpecificHeat = Real(unit = "MJ/(kg.s)");
    type DecayHeatGen = Real(unit = "MJ/(kg.s)");
    // Used in:
    type HeatFlowRate = Real(unit = "MJ/s");
    // Heat Flow Rate [MJ/kg.s]
    // Used in: Heat Transport
    type NominalPower = Real(unit = "1", min = 0);
    // Nominal reactor power as a fraction
    // Used in:
    type ReactorPower = Real(unit = "MW", min = 0);
    // Nominal Reactor Thermal Power
    // Used in:
    type Convection = Real(unit = "MJ/(s.C)", min = 0);
    // Rate of heat transfer by convection [J/(s.C) == W/C]
    // Used in:
    type HeatTransferFraction = Real(unit = "1", min = 0);
    // Fraction of Heat Trasfer from Node to an adjacent Node [unitless]
    // Used in: Nuclear > DecayHeat
    type DecayHeatPrecursorDecayConstant = Real(unit = "1/s", min = 0);
    // Decay heat precursor decay constant [s^-1]
    // Used in:
    type DecayHeatYield = Real(unit = "1/s", min = 0);
    // Used in:
    type DecayHeatFraction = Real(units = "1", min = 0);
    // Used in:
    type Mass = Real(unit = "kg", min = 0);
    // Used in:
    type Density = Real(unit = "kg/(m^3)", min = 0);
    type SpecificHeatCapacity = Real(unit = "MJ/(kg.C)", min = 0);
    // Used in:
    type Temperature = Real(unit = "C", min = 0);
    // Component Related
    // Used in:
    type DHRStimeConstant = Real(unit = "1/s", min = 0);
    // DHRS time constant [per seconds]
    // Used in:
    type InitiationTime = Real(unit = "s", min = 0);
    // Start component [seconds]
    // Used in:
    type Demand = Real(unit = "1", min = 0);
    // Component demand
    // Used in:
    type PumpConstant = Real(unit = "1/s", min = 0);
    // Used in:
    type Pressure = Real(unit = "MPa", min = 0);
    type Length = Real(unit = "m", min = 0);
    type Number = Real(unit = "1");
    type Area = Real(unit = "m^2", min = 0);
    type MicroAbsorptionCrossSection = Real(unit = "cm^2", min = 0);
    type MacroAbsorptionCrossSection = Real(unit = "1/cm", min = 0);
    type TsatConstant = Real(unit = "C/MPa", min = 0);
    type CompressibilityFactor = Real(unit = "1", min = 0);
    type UniversalGasConstant = Real(unit = "MJ/(mol.C)", min = 0);
    type MolarMass = Real(unit = "kg/mol", min = 0);
    type K_4 = Real(unit = "(m^6.s)/(kg^2)");
    type X_4 = Real(unit = "(m^5)/(kg.s)", min = 0);
    type OffgasRemovalRate = Real(unit = "1/s", min = 0);
    type FissFactor = Real(unit = "1", min = 0);
    type EnthalpyChange = Real(unit = "1");
    type SpecificEnthalpy = Real(unit = "MJ/kg");
    type MolFlowRate = Real(unit = "mol/s");
    type pressure_kPa = Real(unit = "kPa", min = 0);
    // Reminder: 1Bar is equal to 100kPa
    type HeatTransferCoeff = Real(unit = "MJ/(s*C*m)");
    type EnthalpyOfVaporization = Real(unit = "MJ/mol");
    type HeatOfReaction = Real(unit = "MJ/mol");
    type Time = Real(unit = "s");
    type FaradayConstant = Real(unit = "C/mol");
    type Volt = Real(unit = "J/C");
    type WeightPercent = Real(unit = "1", min = 0, max = 1);
    type Power = Real(unit = "MJ/s", min = 0);
    type Conversion = Real(unit = "1", min = 0, max = 1);
    type StandardEntropy = Real(unit = "MJ/(mol*K)");
    type EnthalpyOfReaction = Real(unit = "MJ/mol");

    type ReactivityRate = Real(unit = "1/s");
    
    
  end Units;

  package Connectors
    // Creates all inlets and outlets

    connector Temp_In
      MCFR_HyS_BOP.Units.Temperature T;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, fillPattern = FillPattern.Solid, extent = {{-51, 47}, {51, -47}}), Ellipse(origin = {-8, 36}, extent = {{0, -2}, {0, 2}})}),
        Icon(graphics = {Ellipse(origin = {1, -1}, fillPattern = FillPattern.Solid, extent = {{-51, 47}, {51, -47}})}));
    end Temp_In;

    connector Temp_Out
      MCFR_HyS_BOP.Units.Temperature T;
      annotation(
        Icon(graphics = {Ellipse(origin = {1, -1}, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}),
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}));
    end Temp_Out;

    connector DecayHeat_In
      MCFR_HyS_BOP.Units.DecayHeatGen Q;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineColor = {0, 225, 255}, fillColor = {0, 225, 255}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}),
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {0, 225, 255}, fillColor = {0, 225, 255}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}));
    end DecayHeat_In;

    connector DecayHeat_Out
      MCFR_HyS_BOP.Units.DecayHeatGen Q;
      annotation(
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {0, 225, 255}, fillColor = {0, 225, 255}, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}),
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineColor = {0, 225, 255}, fillColor = {0, 225, 255}, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}));
    end DecayHeat_Out;

    connector NominalNeutronPopulation
      MCFR_HyS_BOP.Units.NominalNeutronPopulation n;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}),
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}})}));
    end NominalNeutronPopulation;

    connector Reactivity
      MCFR_HyS_BOP.Units.Reactivity rho;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineColor = {78, 154, 6}, fillColor = {78, 154, 6}, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}})}),
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {78, 154, 6}, fillColor = {78, 154, 6}, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}})}));
    end Reactivity;

    connector Pressure
      MCFR_HyS_BOP.Units.Pressure P;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, fillColor = {255, 170, 255}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}})}),
        Icon(graphics = {Ellipse(origin = {1, -1}, fillColor = {255, 170, 255}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}})}));
    end Pressure;
  
    connector ExternalReactivity
      MCFR_HyS_BOP.Units.Reactivity rho_Ex;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineColor = {78, 154, 6}, fillColor = {78, 154, 6}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}),
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {78, 154, 6}, fillColor = {78, 154, 6}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}));
    end ExternalReactivity;

    connector FlowFraction_In
      MCFR_HyS_BOP.Units.FlowFraction FF;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, fillColor = {245, 121, 0}, fillPattern = FillPattern.Solid, extent = {{-51, 47}, {51, -47}})}),
        Icon(graphics = {Ellipse(origin = {1, -1}, fillColor = {245, 121, 0}, fillPattern = FillPattern.Solid, extent = {{-51, 47}, {51, -47}})}));
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, fillColor = {245, 121, 0}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}),
        Icon(graphics = {Ellipse(origin = {1, -1}, fillColor = {245, 121, 0}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}));
    end FlowFraction_In;

    connector FlowFraction_Out
      MCFR_HyS_BOP.Units.FlowFraction FF;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineColor = {245, 121, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}),
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {245, 121, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}));
    end FlowFraction_Out;
    
    connector FlowRate_In
      MCFR_HyS_BOP.Units.MassFlowRate mdot;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, fillColor = {170, 85, 255}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}),
        Icon(graphics = {Ellipse(origin = {1, -1}, fillColor = {170, 85, 255}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}));
    end FlowRate_In;
  
    connector FlowRate_Out
      MCFR_HyS_BOP.Units.MassFlowRate mdot;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineColor = {170, 85, 255}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}),
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {170, 85, 255}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}));
    end FlowRate_Out;
  end Connectors;

  model notes
  equation
// ===== Misc =====
// Errors/warnings at the start of a simulation are due to variables settling into steady state after initialization. Transients should thus begin a few thousand seconds in.
// Hard crashes of the model are often due to the OTSG. As is, the model cannot handle when the superheated steam domain length goes to zero. That happens when the temperatures in the OTSG go too low for too long. Improvement of the OTSG in this regard means that saturated steam is allowed through. Adjustments would possibly be needed in the RBOP too.
// ===== Organization =====
//// The "Simbase" class is what connects all the modules/components together. It has the values for the parameters that are most important and/or most likely to be changed to simulate different transients.
//// The "Nuclear" package has the modules necessary to model the reactor and give it improved capabilities/features. These include modified point kinetics (mPKE) for tracking neutron population, DecayHeat for tracking decay heat via precursor groups (similar to delayed neutrons), Poisons for tracking poison and poison parent concentrations, and ReactivityFeedback for determining all feedbacks and producing a net reactivity value for mPKE.
//// The "HeatTransport" package has all the components of the IES. These include the MCFR Core, the eight primary heat exchangers (which are combined as one component), the secondary heat exchanger (also called the Thermal Manifold), the Once-Through Steam Generator (OTSG), Rankine Balance of Plant (RBOP), and hybrid sulfur cycle (HyS). The RBOP and HyS are not purely heat exchangers, have interdepencies with their components, and thus are grouped into their own model class script in this organization. The 'Water' and 'MoistAir' packages are necessary for the RBOP to function. The Ultimate Heat Exchanger (UHX) is a simple heat removal component that is not used by default but was used for debugging.
//// The "Pumps" package has the modules for each of the main salt loops (primary, secondary, and tertiary). These are simple and only determine the flow fraction when the pump is tripped. By default, they output the flow fraction as 1.
//// The "Units" package lists and defines each of the units used throughout the model.
//// The "Connectors" package lists and defines each of the connectors used throughout the model.
  end notes;
end MCFR_HyS_BOP;
