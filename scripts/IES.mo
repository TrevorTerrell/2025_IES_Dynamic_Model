model IES
//
  // Load parameters for Nuclear and HeatTransport packages
  MCFR_HyS_BOP.HeatTransport.Core Core(P = 1200, TotalFuelMass = 638300, hAnom = 18, mdot_fuelNom = 45000) annotation(
    Placement(visible = true, transformation(origin = {-261, -89}, extent = {{-555, -555}, {555, 555}}, rotation = 0)));
  MCFR_HyS_BOP.HeatTransport.PrimaryHeatExchanger PHX(m_dot_pFluidNom_PHX = 45000, m_dot_sFluidNom_PHX = 27272.73, hApnNom_PHX = 20, hAsnNom_PHX = 20) annotation(
    Placement(visible = true, transformation(origin = {1237, 389}, extent = {{-427, -427}, {427, 427}}, rotation = 0)));
  MCFR_HyS_BOP.HeatTransport.SecondaryHeatExchanger SHX(m_dot_pFluidNom_SHX = 27272.73, m_dot_sFluidNom_SHX = 12910.56, hApnNom_SHX = 9.2014, hAsnNom_SHX = 9.2014) annotation(
    Placement(visible = true, transformation(origin = {1251, -302}, extent = {{-432, -432}, {432, 432}}, rotation = 0))); 
  MCFR_HyS_BOP.HeatTransport.OnceThroughSteamGenerator OTSG(m_dot_hitecNom = 6455.28, tauOTSGtoSHX = 10) annotation(Placement(visible = true, transformation(origin = {1386, -1133}, extent = {{-343, -343}, {343, 343}}, rotation = 0)));
  MCFR_HyS_BOP.HeatTransport.RBOP RBOP(m_dot_fdw_nom = 192.075, omega = 0, sinmag = 0.01, sinstart = 88888888, rampstart = 50000000, ramp_min = 0.2, m_dot_cool = 20000) annotation(Placement(visible = true, transformation(origin = {279, -1258}, extent = {{-296, -296}, {296, 296}}, rotation = 0)));
  MCFR_HyS_BOP.HeatTransport.HyS HyS(mdot_H_nom = 6455.28, pressure = 100, Mdot_SA_nom = 3453.3, Mdot_H2O_nom = 345.33, omega = 0, sinmag = 0.01, sinstart = 88888888, rampstart = 500000000, ramp_min = 0.2) annotation(Placement(visible = true, transformation(origin = {2355, -608}, extent = {{-296, -296}, {296, 296}}, rotation = 0)));
  MCFR_HyS_BOP.Nuclear.mPKE mPKE(Lam = 8.77488E-06, lam = {1.24909E-02, 3.15594E-02, 1.10482E-01, 3.22830E-01, 1.34074, 9.03297}, beta = {2.12256E-04, 1.12099E-03, 1.11237E-03, 3.39814E-03, 1.19113E-03, 3.72466E-04}) annotation(
    Placement(visible = true, transformation(origin = {-1323, 315}, extent = {{-409, -409}, {409, 409}}, rotation = 0)));
  MCFR_HyS_BOP.Nuclear.DecayHeat decayheat(P = 1200, TotalFuelMass = 638300) annotation(
    Placement(visible = true, transformation(origin = {-1181, -471}, extent = {{-463, -463}, {463, 463}}, rotation = 0)));
  MCFR_HyS_BOP.Nuclear.Poisons Poisons(lam_bubble = 0.005) annotation(
    Placement(visible = true, transformation(origin = {-1265, -1015}, extent = {{-555, -555}, {555, 555}}, rotation = 0)));
  MCFR_HyS_BOP.Nuclear.ReactivityFeedback ReactivityFeedback(a_F = -4.4376e-05, a_R = 0.1164e-05, step_mag = -185*1e-5, omega = %%freq%%, sin_mag = 1*1e-5, ramp_rate = -1e-5, ramp_low = 0.4, stepInsertionTime = 500000000, sinInsertionTime = 1000, rampInsertionTime = 500000000) annotation(
    Placement(visible = true, transformation(origin = {-519, -1063}, extent = {{545, -545}, {-545, 545}}, rotation = 0)));
  MCFR_HyS_BOP.Pumps.PrimaryLoopPump PLP(freeConvectionFF = 0.1, primaryPumpK = 0.092, tripPrimaryPump = 50000000) annotation(
    Placement(visible = true, transformation(origin = {364, 566}, extent = {{-262, -262}, {262, 262}}, rotation = 0)));
  MCFR_HyS_BOP.Pumps.SecondaryLoopPump SLP(freeConvectionFF = 0.05, secondaryPumpK = 0.092, tripSecondaryPump = 2000000000) annotation(
    Placement(visible = true, transformation(origin = {366, -96}, extent = {{-232, -232}, {232, 232}}, rotation = 0)));
  MCFR_HyS_BOP.Pumps.TertiaryLoopPump TLP(freeConvectionFF = 0.05, tertiaryPumpK = 0.092, tripTertiaryPump = 2000000000) annotation(
    Placement(visible = true, transformation(origin = {375, -583}, extent = {{-217, -217}, {217, 217}}, rotation = 0)));
  // Connections between modules
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
end IES;
