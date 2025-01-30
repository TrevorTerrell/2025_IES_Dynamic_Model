model MCFR_Simbase
  //
  // Load parameters for Nuclear and HeatTransport packages
  MCFR_Packages.HeatTransport.Core Core(TotalFuelMass = 638300, hAnom = 18, mdot_fuelNom = 45000) annotation(
    Placement(visible = true, transformation(origin = {-261, -89}, extent = {{-555, -555}, {555, 555}}, rotation = 0)));
  MCFR_Packages.HeatTransport.PrimaryHeatExchanger PHX(m_dot_pFluidNom_PHX = 45000, m_dot_sFluidNom_PHX = 27272.73, hApnNom_PHX = 20, hAsnNom_PHX = 20) annotation(
    Placement(visible = true, transformation(origin = {1581, 125}, extent = {{-427, -427}, {427, 427}}, rotation = 0)));
  MCFR_Packages.HeatTransport.SecondaryHeatExchanger SHX(m_dot_pFluidNom_SHX = 27272.73, m_dot_sFluidNom_SHX = 12910.56, hApnNom_SHX = 9.2014, hAsnNom_SHX = 9.2014) annotation(
    Placement(visible = true, transformation(origin = {1611, -610}, extent = {{-432, -432}, {432, 432}}, rotation = 0)));
  MCFR_Packages.HeatTransport.UHX UHX(UHXP = 1200, tripUHX = 2000000000, SetDemand = 0.5) annotation(
    Placement(visible = true, transformation(origin = {1829, -1383}, extent = {{479, 479}, {-479, -479}}, rotation = -180))); 
  MCFR_Packages.Nuclear.mPKE mPKE(Lam = 8.77488E-06, lam = {1.24909E-02, 3.15594E-02, 1.10482E-01, 3.22830E-01, 1.34074, 9.03297}, beta = {2.12256E-04, 1.12099E-03, 1.11237E-03, 3.39814E-03, 1.19113E-03, 3.72466E-04}) annotation(
    Placement(visible = true, transformation(origin = {-1323, 315}, extent = {{-409, -409}, {409, 409}}, rotation = 0)));
  MCFR_Packages.Nuclear.DecayHeat decayheat(P = 1200, TotalFuelMass = 638300) annotation(
    Placement(visible = true, transformation(origin = {-1181, -471}, extent = {{-463, -463}, {463, 463}}, rotation = 0)));
  MCFR_Packages.Nuclear.Poisons Poisons(lam_bubble = 0.005) annotation(
    Placement(visible = true, transformation(origin = {-1283, -1041}, extent = {{-555, -555}, {555, 555}}, rotation = 0)));
  MCFR_Packages.Nuclear.ReactivityFeedback ReactivityFeedback(a_F = -4.4376e-05, a_R = 0.1164e-05, step_mag = 185*1e-5, omega = 0.01, sin_mag = 50*1e-5, stepInsertionTime = 2000000000, sinInsertionTime = 5000) annotation(
    Placement(visible = true, transformation(origin = {-83, -1051}, extent = {{545, -545}, {-545, 545}}, rotation = 0)));
  MCFR_Packages.Pumps.PrimaryLoopPump PLP(freeConvectionFF = 0.05, primaryPumpK = 0.092, tripPrimaryPump = 2000000000) annotation(
    Placement(visible = true, transformation(origin = {570, 566}, extent = {{-262, -262}, {262, 262}}, rotation = 0)));
  MCFR_Packages.Pumps.SecondaryLoopPump SLP(freeConvectionFF = 0.05, secondaryPumpK = 0.092, tripSecondaryPump = 2000000000) annotation(
    Placement(visible = true, transformation(origin = {596, -92}, extent = {{-232, -232}, {232, 232}}, rotation = 0)));
  MCFR_Packages.Pumps.TertiaryLoopPump TLP(freeConvectionFF = 0.05, tertiaryPumpK = 0.092, tripTertiaryPump = 2000000000) annotation(
    Placement(visible = true, transformation(origin = {595, -815}, extent = {{-217, -217}, {217, 217}}, rotation = 0)));
  // Connections between modules
equation
  connect(PLP.primaryFF, Core.fuelFlowFraction) annotation(
    Line(points = {{570, 414}, {-483, 414}, {-483, -533}}, color = {245, 121, 0}));
  connect(PLP.primaryFF, mPKE.fuelFlowFrac) annotation(
    Line(points = {{570, 414}, {-1184, 414}, {-1184, 479.04}}, color = {245, 121, 0}));
  connect(PLP.primaryFF, PHX.primaryFF) annotation(
    Line(points = {{570, 414}, {1282, 414}, {1282, 338.04}}, color = {245, 121, 0}));
  connect(SLP.secondaryFF, SHX.secondaryFF) annotation(
    Line(points = {{596, -231.2}, {1309, -231.2}, {1309, -394.2}}, color = {245, 121, 0}));
  connect(SLP.secondaryFF, PHX.secondaryFF) annotation(
    Line(points = {{596, -231.2}, {2136, -231.2}, {2136, -88.2}}, color = {245, 121, 0}));
  connect(TLP.tertiaryFF, SHX.tertiaryFF) annotation(
    Line(points = {{595, -945}, {2173, -945}, {2173, -826.2}}, color = {245, 121, 0}));
  connect(Poisons.Poison_react, ReactivityFeedback.Poison_react) annotation(
    Line(points = {{-1172, -1263}, {-452.5, -1263}, {-452.5, -1225}, {-225, -1225}}, color = {78, 154, 6}));
  connect(Core.nPop, mPKE.n_population) annotation(
    Line(points = {{-705, 233}, {-1323, 233}, {-1323, 159.9}}, color = {20, 36, 248}));
  connect(mPKE.n_population, decayheat.nPop) annotation(
    Line(points = {{-1323, 159.58}, {-1329, 159.58}, {-1329, -369.42}}, color = {20, 36, 248}));
  connect(mPKE.n_population, Poisons.n_population) annotation(
    Line(points = {{-1323, 159.58}, {-1505, 159.58}, {-1505, -1263}}, color = {20, 36, 248}));
  connect(decayheat.decayHeat_Out, Core.P_decay) annotation(
    Line(points = {{-1014.32, -609.9}, {-802.82, -609.9}, {-802.82, -122}, {-705, -122}}, color = {0, 225, 255}));
  connect(decayheat.decayHeat_Out, PHX.P_decay) annotation(
    Line(points = {{-1014.32, -609.9}, {1700.68, -609.9}, {1700.68, 338.1}}, color = {0, 225, 255}));
  connect(ReactivityFeedback.feedback, mPKE.feedback) annotation(
    Line(points = {{-224.7, -942}, {-1322.7, -942}, {-1322.7, 479}}, color = {78, 154, 6}));
  connect(Core.fuelNode1, ReactivityFeedback.fuelNode1) annotation(
    Line(points = {{16.5, -489}, {135, -489}, {135, -865.6}}));
  connect(Core.fuelNode2, ReactivityFeedback.fuelNode2) annotation(
    Line(points = {{16.5, -255.5}, {135, -255.5}, {135, -1051}}));
  connect(Core.ReflcNode, ReactivityFeedback.ReflcNode) annotation(
    Line(points = {{16.5, 0}, {135, 0}, {135, -1215.2}}));
  connect(SHX.T_out_pFluid_SHX, PHX.T_in_sFluid_PHX) annotation(
    Line(points = {{2310.84, -480.4}, {2310.84, -206.4}, {2272.84, -206.4}, {2272.84, -12.4}}));
  connect(PHX.T_out_sFluid_PHX, SHX.T_in_pFluid_SHX) annotation(
    Line(points = {{1154, -3.1}, {1154, -199.6}, {1170, -199.6}, {1170, -480.1}}));
  connect(UHX.UHXtemp_Out, SHX.T_in_sFluid_SHX) annotation(
    Line(points = {{1954, -1393}, {1954, -748}, {2310, -748}}));
  connect(SHX.T_out_sFluid_SHX, UHX.UHXtemp_In) annotation(
    Line(points = {{1180, -740}, {1184, -740}, {1184, -1393}, {1494, -1393}}));
  connect(Core.temp_Out, PHX.T_in_pFluid_PHX) annotation(
    Line(points = {{16.5, 205}, {944, 205}, {944, 254}, {1146, 254}}));
  connect(PHX.T_out_pFluid_PHX, Core.temp_In) annotation(
    Line(points = {{2272, 254}, {2446, 254}, {2446, 716}, {-878, 716}, {-878, -477.5}, {-705, -477.5}}));
  connect(TLP.tertiaryFF, UHX.tertiaryFF) annotation(
    Line(points = {{595, -945}, {1494, -945}, {1494, -1143.5}}, color = {245, 121, 0}));
  annotation(
    Diagram(coordinateSystem(extent = {{-1800, 1120}, {2740, -1740}})),
    version = "",
    uses,
    experiment(StartTime = 0, StopTime = 10000, Tolerance = 1e-06, Interval = 0.1));
end MCFR_Simbase;
