package MCFR_Packages
  package HeatTransport
    model Core
      // Fission power, decay power, nominal power, and core flow
      parameter MCFR_Packages.Units.ReactorPower P = 1200;
      parameter MCFR_Packages.Units.Mass TotalFuelMass;
      parameter MCFR_Packages.Units.Convection hAnom;
      parameter MCFR_Packages.Units.Mass m_FN1 = 159575;
      parameter MCFR_Packages.Units.Mass m_FN2 = 159575;
      parameter MCFR_Packages.Units.Mass m_RN = 969575;
      parameter MCFR_Packages.Units.SpecificHeatCapacity cP_fuel = 0.00066065;
      parameter MCFR_Packages.Units.SpecificHeatCapacity cP_reflector = 0.0013;
      parameter MCFR_Packages.Units.MassFlowRate mdot_fuelNom;
      parameter MCFR_Packages.Units.VolumeImportance kFN1 = 0.465;
      parameter MCFR_Packages.Units.VolumeImportance kFN2 = 0.465;
      parameter MCFR_Packages.Units.VolumeImportance kR = 0.07;
      parameter MCFR_Packages.Units.HeatTransferFraction kHT_FN1 = 0.5;
      parameter MCFR_Packages.Units.HeatTransferFraction kHT_FN2 = 0.5;
      parameter MCFR_Packages.Units.ResidentTime tauCoreToPHX = 0.6887;
      parameter MCFR_Packages.Units.Temperature fuelNode1_initial = 700;
      parameter MCFR_Packages.Units.Temperature fuelNode2_initial = 720;
      parameter MCFR_Packages.Units.Temperature ReflcNode_initial = 700;
      parameter MCFR_Packages.Units.Mass m_DHRS = 31915;
      parameter MCFR_Packages.Units.DHRStimeConstant DHRS_timeConstant = 10;
      parameter MCFR_Packages.Units.HeatFlowRate DHRS_MaxPowerRm = 0.1*1200;
      parameter MCFR_Packages.Units.HeatFlowRate DHRS_PowerBleed = 0.001*1200;
      parameter MCFR_Packages.Units.InitiationTime DHRS_time = 2000000000;
      MCFR_Packages.Units.ResidentTime varTauCoreToPHX;
      MCFR_Packages.Units.NominalPower FissionPower;
      MCFR_Packages.Units.NominalPower DecayPower;
      MCFR_Packages.Units.NominalPower NomPower;
      MCFR_Packages.Units.MassFlowRate mdot_fuel;
      MCFR_Packages.Units.Convection hA;
      MCFR_Packages.Units.Temperature DHRS_Temp;
      MCFR_Packages.Units.HeatFlowRate DHRS_PowerRm;
      input MCFR_Packages.Connectors.Temp_In temp_In annotation(
        Placement(visible = true, transformation(origin = {-76, -70}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {-80, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input MCFR_Packages.Connectors.DecayHeat_In P_decay annotation(
        Placement(visible = true, transformation(origin = {-77, -7}, extent = {{-17, -17}, {17, 17}}, rotation = 0), iconTransformation(origin = {-80, -6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input MCFR_Packages.Connectors.NominalNeutronPopulation nPop annotation(
        Placement(visible = true, transformation(origin = {-80, 50}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {-80, 58}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output MCFR_Packages.Connectors.Temp_Out temp_Out annotation(
        Placement(visible = true, transformation(origin = {44, 60}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {50, 53}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
      output MCFR_Packages.Connectors.Temp_Out fuelNode1 annotation(
        Placement(visible = true, transformation(origin = {45, -78}, extent = {{-17, -17}, {17, 17}}, rotation = 0), iconTransformation(origin = {50, -72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output MCFR_Packages.Connectors.Temp_Out fuelNode2 annotation(
        Placement(visible = true, transformation(origin = {45, -40}, extent = {{-17, -17}, {17, 17}}, rotation = 0), iconTransformation(origin = {50, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output MCFR_Packages.Connectors.Temp_Out ReflcNode annotation(
        Placement(visible = true, transformation(origin = {44, 16}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {50, 16}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input MCFR_Packages.Connectors.FlowFraction_In fuelFlowFraction annotation(
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
    //  hA = hAnom*(0.8215*fuelFlowFraction.FF^6 - 4.108*fuelFlowFraction.FF^5 + 7.848*fuelFlowFraction.FF^4 - 7.165*fuelFlowFraction.FF^3 + 3.004*fuelFlowFraction.FF^2 + 0.5903*fuelFlowFraction.FF + 0.008537);
    hA = hAnom;
      varTauCoreToPHX = tauCoreToPHX/fuelFlowFraction.FF;
      m_FN1*cP_fuel*der(fuelNode1.T) = mdot_fuel*cP_fuel*(temp_In.T - fuelNode1.T) + kFN1*FissionPower*P - hA*kHT_FN1*(fuelNode1.T - ReflcNode.T) + P_decay.Q*m_FN1;
      m_FN2*cP_fuel*der(fuelNode2.T) = mdot_fuel*cP_fuel*(fuelNode1.T - fuelNode2.T) + kFN2*FissionPower*P - hA*kHT_FN2*(fuelNode1.T - ReflcNode.T) + P_decay.Q*m_FN2;
      m_RN*cP_reflector*der(ReflcNode.T) = hA*(fuelNode1.T - ReflcNode.T) + kR*FissionPower*P;
      DHRS_PowerRm = (DHRS_MaxPowerRm - DHRS_PowerBleed)/(1 + exp(log(1/1E-3 - 1)*(1 - (time - DHRS_time)/DHRS_timeConstant))) + DHRS_PowerBleed;
      m_DHRS*cP_fuel*der(DHRS_Temp) = mdot_fuel*cP_fuel*(fuelNode2.T - DHRS_Temp) + P_decay.Q*m_DHRS - DHRS_PowerRm;
      temp_Out.T = delay(DHRS_Temp, varTauCoreToPHX, 5000);
      annotation(
        Diagram(graphics = {Rectangle(origin = {-15.4297, -10.014}, extent = {{-75.8099, 90.0986}, {75.8099, -90.0986}}), Rectangle(origin = {-40, 20}, lineColor = {255, 63, 63}, fillColor = {255, 63, 63}, fillPattern = FillPattern.Solid, lineThickness = 0, extent = {{-20, 20}, {20, -20}}), Rectangle(origin = {-40, -40}, lineColor = {255, 63, 63}, fillColor = {255, 63, 63}, fillPattern = FillPattern.Solid, lineThickness = 0, extent = {{-20, 20}, {20, -20}}), Rectangle(origin = {15, -10}, lineColor = {94, 95, 92}, fillColor = {94, 95, 92}, fillPattern = FillPattern.Solid, lineThickness = 0, extent = {{-15, 50}, {15, -50}}), Text(origin = {-16, 70}, extent = {{-47, 12}, {47, -12}}, textString = "Core Heat Transfer"), Text(origin = {-76, 35}, extent = {{-13, 13}, {13, -13}}, textString = "n(t)/n_0"), Text(origin = {-76, -20}, extent = {{-13, 13}, {13, -13}}, textString = "Decay Heat"), Text(origin = {-75, -85}, extent = {{-13, 13}, {13, -13}}, textString = "Temp_In"), Text(origin = {46, 45}, extent = {{-13, 13}, {13, -13}}, textString = "Temp_Out"), Text(origin = {47, 2}, extent = {{-13, 13}, {13, -13}}, textString = "Temp_G1"), Text(origin = {45, -54}, extent = {{-13, 13}, {13, -13}}, textString = "Temp_FN2"), Text(origin = {45, -92}, extent = {{-13, 13}, {13, -13}}, textString = "Temp_FN1")}),
        Icon(graphics = {Rectangle(origin = {-15.4297, -10.014}, extent = {{-75.8099, 90.0986}, {75.8099, -90.0986}}), Rectangle(origin = {-40, 20}, lineColor = {255, 63, 63}, fillColor = {255, 63, 63}, fillPattern = FillPattern.Solid, lineThickness = 0, extent = {{-20, 20}, {20, -20}}), Rectangle(origin = {-40, -40}, lineColor = {255, 63, 63}, fillColor = {255, 63, 63}, fillPattern = FillPattern.Solid, lineThickness = 0, extent = {{-20, 20}, {20, -20}}), Rectangle(origin = {15, -10}, lineColor = {94, 95, 92}, fillColor = {94, 95, 92}, fillPattern = FillPattern.Solid, lineThickness = 0, extent = {{-15, 50}, {15, -50}}), Text(origin = {-16, 70}, extent = {{-47, 12}, {47, -12}}, textString = "Core Heat Transfer"), Text(origin = {-76, 35}, extent = {{-13, 13}, {13, -13}}, textString = "n(t)/n_0"), Text(origin = {-76, -20}, extent = {{-13, 13}, {13, -13}}, textString = "Decay Heat"), Text(origin = {-75, -85}, extent = {{-13, 13}, {13, -13}}, textString = "Temp_In"), Text(origin = {46, 45}, extent = {{-13, 13}, {13, -13}}, textString = "Temp_Out"), Text(origin = {47, 2}, extent = {{-13, 13}, {13, -13}}, textString = "Temp_G1"), Text(origin = {45, -54}, extent = {{-13, 13}, {13, -13}}, textString = "Temp_FN2"), Text(origin = {45, -92}, extent = {{-13, 13}, {13, -13}}, textString = "Temp_FN1")}),
        experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.002));
    end Core;

    // Parameters and equations for all heat exchangers and pumps

    model PrimaryHeatExchanger
      // PHX parameters and equations
      // Parameter declaration
      parameter MCFR_Packages.Units.Mass m_PN1_PHX = 63830;
      // Fuel Salt node 1 mass
      parameter MCFR_Packages.Units.Mass m_PN2_PHX = 63830;
      // Fuel Salt node 2 mass
      parameter MCFR_Packages.Units.Mass m_PN3_PHX = 63830;
      // Fuel Salt node 3 mass
      parameter MCFR_Packages.Units.Mass m_PN4_PHX = 63830;
      // Fuel Salt node 4 mass
      parameter MCFR_Packages.Units.Mass m_TN1_PHX = 2860;
      // PHX Tubing node 1 mass
      parameter MCFR_Packages.Units.Mass m_TN2_PHX = 2860;
      // PHX Tubing node 2 mass
      parameter MCFR_Packages.Units.Mass m_SN1_PHX = 54545;
      // Secondary Side node 1 mass
      parameter MCFR_Packages.Units.Mass m_SN2_PHX = 54545;
      // Secondary Side node 2 mass
      parameter MCFR_Packages.Units.Mass m_SN3_PHX = 54545;
      // Secondary Side node 3 mass
      parameter MCFR_Packages.Units.Mass m_SN4_PHX = 54545;
      // Secondary Side node 4 mass
      parameter MCFR_Packages.Units.Mass m_pipe = 31915;
      // Mass of fluid in pipe between PHX to Core
      parameter MCFR_Packages.Units.SpecificHeatCapacity cP_pFluid_PHX = 0.00066065;
      // Fuel salt specific heat [MJ/(kg*C)]
      parameter MCFR_Packages.Units.SpecificHeatCapacity cP_Tube_PHX = 0.00057779;
      // PHX Tubing specific heat [MJ/(kg*C)]
      parameter MCFR_Packages.Units.SpecificHeatCapacity cP_sFluid_PHX = 0.0011;
      // Coolant salt specific heat [MJ/(kg*C)]
      parameter MCFR_Packages.Units.MassFlowRate m_dot_pFluidNom_PHX;
      // Fuel salt nominal flow rate
      parameter MCFR_Packages.Units.MassFlowRate m_dot_sFluidNom_PHX;
      // Coolant salt nominal flow rate
      parameter MCFR_Packages.Units.Convection hApnNom_PHX;
      // PHX primary side heat transfer coefficient [MJ/(s*C)]
      parameter MCFR_Packages.Units.Convection hAsnNom_PHX;
      // PHX secondary side heat transfer coefficient [MJ/(s*C)]
      parameter MCFR_Packages.Units.ResidentTime tauPHXtoCore = 0.6887;
      // Avg travel time from PHX to pipe
      parameter MCFR_Packages.Units.ResidentTime tauPHXtoSHX = 2.5;
      // Avg travel time from PHX to SHX
      //parameter MCFR_Packages.Units.ResidentTime tauPipeToCore = 0.6887;
      parameter MCFR_Packages.Units.Temperature T_PN1_PHX_initial = 708.6;
      parameter MCFR_Packages.Units.Temperature T_PN2_PHX_initial = 699.1;
      parameter MCFR_Packages.Units.Temperature T_PN3_PHX_initial = 689.3;
      parameter MCFR_Packages.Units.Temperature T_PN4_PHX_initial = 679.5;
      parameter MCFR_Packages.Units.Temperature T_TN1_PHX_initial = 693.8;
      parameter MCFR_Packages.Units.Temperature T_TN2_PHX_initial = 674.1;
      parameter MCFR_Packages.Units.Temperature T_SN1_PHX_initial = 659.0;
      parameter MCFR_Packages.Units.Temperature T_SN2_PHX_initial = 669.1;
      parameter MCFR_Packages.Units.Temperature T_SN3_PHX_initial = 679.0;
      parameter MCFR_Packages.Units.Temperature T_SN4_PHX_initial = 688.9;
      // Variable declaration
      MCFR_Packages.Units.ResidentTime varTauPHXtoCore;
      MCFR_Packages.Units.ResidentTime varTauPHXtoSHX;
      MCFR_Packages.Units.MassFlowRate m_dot_pFluid_PHX;
      MCFR_Packages.Units.MassFlowRate m_dot_sFluid_PHX;
      MCFR_Packages.Units.Convection hApn_PHX;
      MCFR_Packages.Units.Convection hAsn_PHX;
      MCFR_Packages.Units.Temperature T_PN1_PHX;
      MCFR_Packages.Units.Temperature T_PN2_PHX;
      MCFR_Packages.Units.Temperature T_PN3_PHX;
      MCFR_Packages.Units.Temperature T_PN4_PHX;
      MCFR_Packages.Units.Temperature T_TN1_PHX;
      MCFR_Packages.Units.Temperature T_TN2_PHX;
      MCFR_Packages.Units.Temperature T_SN1_PHX;
      MCFR_Packages.Units.Temperature T_SN2_PHX;
      MCFR_Packages.Units.Temperature T_SN3_PHX;
      MCFR_Packages.Units.Temperature T_SN4_PHX;
      MCFR_Packages.Units.ReactorPower P_primaryLoop;
      MCFR_Packages.Units.Temperature Pipe_Temp;
      // Connections to other modules
      input MCFR_Packages.Connectors.Temp_In T_in_pFluid_PHX annotation(
        Placement(visible = true, transformation(origin = {-100, 30}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-102, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input MCFR_Packages.Connectors.Temp_In T_in_sFluid_PHX annotation(
        Placement(visible = true, transformation(origin = {160, -30}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {162, -32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output MCFR_Packages.Connectors.Temp_Out T_out_sFluid_PHX annotation(
        Placement(visible = true, transformation(origin = {-100, -30}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-100, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output MCFR_Packages.Connectors.Temp_Out T_out_pFluid_PHX annotation(
        Placement(visible = true, transformation(origin = {160, 30}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {162, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input MCFR_Packages.Connectors.DecayHeat_In P_decay annotation(
        Placement(visible = true, transformation(origin = {30, 50}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {28, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input MCFR_Packages.Connectors.FlowFraction_In primaryFF annotation(
        Placement(visible = true, transformation(origin = {-70, 50}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {-70, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input MCFR_Packages.Connectors.FlowFraction_In secondaryFF annotation(
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
    // These next two are optional. They attempt to include some complexity of htc changing due to flow rate turbulence
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
      T_out_sFluid_PHX.T = delay(T_SN4_PHX, varTauPHXtoSHX, 5000);
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
      parameter MCFR_Packages.Units.Mass m_PN1_SHX = 54545;
      // Hot side (coolant salt from PHX) node 1 mass
      parameter MCFR_Packages.Units.Mass m_PN2_SHX = 54545;
      // Hot side (coolant salt from PHX) node 2 mass
      parameter MCFR_Packages.Units.Mass m_PN3_SHX = 54545;
      // Hot side (coolant salt from PHX) node 3 mass
      parameter MCFR_Packages.Units.Mass m_PN4_SHX = 54545;
      // Hot side (coolant salt from PHX) node 4 mass
      parameter MCFR_Packages.Units.Mass m_TN1_SHX = 4268;
      // SHX Tubing node 1 mass
      parameter MCFR_Packages.Units.Mass m_TN2_SHX = 4268;
      // SHX Tubing node 2 mass
      parameter MCFR_Packages.Units.Mass m_SN1_SHX = 20460;
      // Cold Side (Hitec salt from OTSG) node 1 mass
      parameter MCFR_Packages.Units.Mass m_SN2_SHX = 20460;
      // Cold Side (Hitec salt from OTSG) node 2 mass
      parameter MCFR_Packages.Units.Mass m_SN3_SHX = 20460;
      // Cold Side (Hitec salt from OTSG) node 3 mass
      parameter MCFR_Packages.Units.Mass m_SN4_SHX = 20460;
      // Cold Side (Hitec salt from OTSG) node 4 mass
      parameter MCFR_Packages.Units.SpecificHeatCapacity cP_pFluid_SHX = 0.0011;
      // Fuel salt specific heat [MJ/(kg*C)]
      parameter MCFR_Packages.Units.SpecificHeatCapacity cP_Tube_SHX = 0.00057779;
      // PHX Tubing specific heat [MJ/(kg*C)]
      parameter MCFR_Packages.Units.SpecificHeatCapacity cP_sFluid_SHX = 0.00154912;
      // Coolant salt specific heat [MJ/(kg*C)]
      parameter MCFR_Packages.Units.MassFlowRate m_dot_pFluidNom_SHX;
      // Fuel salt nominal flow rate
      parameter MCFR_Packages.Units.MassFlowRate m_dot_sFluidNom_SHX;
      // Coolant salt nominal flow rate
      parameter MCFR_Packages.Units.Convection hApnNom_SHX;
      // PHX primary side heat transfer coefficient [MJ/(s*C)]
      parameter MCFR_Packages.Units.Convection hAsnNom_SHX;
      // PHX secondary side heat transfer coefficient [MJ/(s*C)]
      parameter MCFR_Packages.Units.ResidentTime tauSHXtoPHX = 1;
      // Avg travel time from SHX to pipe
      parameter MCFR_Packages.Units.ResidentTime tauSHXtoOTSG = 1;
      // Avg travel time from SHX to PHX
      parameter MCFR_Packages.Units.Temperature T_PN1_SHX_initial = 679.6;
      parameter MCFR_Packages.Units.Temperature T_PN2_SHX_initial = 670.4;
      parameter MCFR_Packages.Units.Temperature T_PN3_SHX_initial = 659.6;
      parameter MCFR_Packages.Units.Temperature T_PN4_SHX_initial = 648.9;
      parameter MCFR_Packages.Units.Temperature T_TN1_SHX_initial = 649.5;
      parameter MCFR_Packages.Units.Temperature T_TN2_SHX_initial = 624.5;
      parameter MCFR_Packages.Units.Temperature T_SN1_SHX_initial = 589.4;
      parameter MCFR_Packages.Units.Temperature T_SN2_SHX_initial = 605.6;
      parameter MCFR_Packages.Units.Temperature T_SN3_SHX_initial = 619.4;
      parameter MCFR_Packages.Units.Temperature T_SN4_SHX_initial = 633.3;
      // Variable declaration
      MCFR_Packages.Units.ResidentTime varTauSHXtoOTSG;
      MCFR_Packages.Units.ResidentTime varTauSHXtoPHX;
      MCFR_Packages.Units.MassFlowRate m_dot_pFluid_SHX;
      MCFR_Packages.Units.MassFlowRate m_dot_sFluid_SHX;
      MCFR_Packages.Units.Convection hApn_SHX;
      MCFR_Packages.Units.Convection hAsn_SHX;
      MCFR_Packages.Units.Temperature T_PN1_SHX;
      MCFR_Packages.Units.Temperature T_PN2_SHX;
      MCFR_Packages.Units.Temperature T_PN3_SHX;
      MCFR_Packages.Units.Temperature T_PN4_SHX;
      MCFR_Packages.Units.Temperature T_TN1_SHX;
      MCFR_Packages.Units.Temperature T_TN2_SHX;
      MCFR_Packages.Units.Temperature T_SN1_SHX;
      MCFR_Packages.Units.Temperature T_SN2_SHX;
      MCFR_Packages.Units.Temperature T_SN3_SHX;
      MCFR_Packages.Units.Temperature T_SN4_SHX;
      MCFR_Packages.Units.ReactorPower P_TertiaryLoop;
      MCFR_Packages.Units.ReactorPower P_SecondaryLoop;
      // Connections to other modules
      input MCFR_Packages.Connectors.Temp_In T_in_pFluid_SHX annotation(
        Placement(visible = true, transformation(origin = {-100, 30}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-102, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input MCFR_Packages.Connectors.Temp_In T_in_sFluid_SHX annotation(
        Placement(visible = true, transformation(origin = {160, -30}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {162, -32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output MCFR_Packages.Connectors.Temp_Out T_out_sFluid_SHX annotation(
        Placement(visible = true, transformation(origin = {-100, -30}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-100, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output MCFR_Packages.Connectors.Temp_Out T_out_pFluid_SHX annotation(
        Placement(visible = true, transformation(origin = {160, 30}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {162, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input MCFR_Packages.Connectors.FlowFraction_In secondaryFF annotation(
        Placement(visible = true, transformation(origin = {-70, 50}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {-70, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input MCFR_Packages.Connectors.FlowFraction_In tertiaryFF annotation(
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
      m_PN1_SHX*cP_pFluid_SHX*der(T_PN1_SHX) = m_dot_pFluid_SHX*cP_pFluid_SHX*(T_in_pFluid_SHX.T - T_PN1_SHX) - hApn_SHX*(T_PN1_SHX - T_TN1_SHX);
      m_PN2_SHX*cP_pFluid_SHX*der(T_PN2_SHX) = m_dot_pFluid_SHX*cP_pFluid_SHX*(T_PN1_SHX - T_PN2_SHX) - hApn_SHX*(T_PN1_SHX - T_TN1_SHX);
      m_PN3_SHX*cP_pFluid_SHX*der(T_PN3_SHX) = m_dot_pFluid_SHX*cP_pFluid_SHX*(T_PN2_SHX - T_PN3_SHX) - hApn_SHX*(T_PN3_SHX - T_TN2_SHX);
      m_PN4_SHX*cP_pFluid_SHX*der(T_PN4_SHX) = m_dot_pFluid_SHX*cP_pFluid_SHX*(T_PN3_SHX - T_PN4_SHX) - hApn_SHX*(T_PN3_SHX - T_TN2_SHX);
      m_TN1_SHX*cP_Tube_SHX*der(T_TN1_SHX) = 2*hApn_SHX*(T_PN1_SHX - T_TN1_SHX) - 2*hAsn_SHX*(T_TN1_SHX - T_SN3_SHX);
      m_TN2_SHX*cP_Tube_SHX*der(T_TN2_SHX) = 2*hApn_SHX*(T_PN3_SHX - T_TN2_SHX) - 2*hAsn_SHX*(T_TN2_SHX - T_SN1_SHX);
      m_SN1_SHX*cP_sFluid_SHX*der(T_SN1_SHX) = m_dot_sFluid_SHX*cP_sFluid_SHX*(T_in_sFluid_SHX.T - T_SN1_SHX) + hAsn_SHX*(T_TN2_SHX - T_SN1_SHX);
      m_SN2_SHX*cP_sFluid_SHX*der(T_SN2_SHX) = m_dot_sFluid_SHX*cP_sFluid_SHX*(T_SN1_SHX - T_SN2_SHX) + hAsn_SHX*(T_TN2_SHX - T_SN1_SHX);
      m_SN3_SHX*cP_sFluid_SHX*der(T_SN3_SHX) = m_dot_sFluid_SHX*cP_sFluid_SHX*(T_SN2_SHX - T_SN3_SHX) + hAsn_SHX*(T_TN1_SHX - T_SN3_SHX);
      m_SN4_SHX*cP_sFluid_SHX*der(T_SN4_SHX) = m_dot_sFluid_SHX*cP_sFluid_SHX*(T_SN3_SHX - T_SN4_SHX) + hAsn_SHX*(T_TN1_SHX - T_SN3_SHX);
      T_out_pFluid_SHX.T = delay(T_PN4_SHX, varTauSHXtoPHX, 5000);
      T_out_sFluid_SHX.T = delay(T_SN4_SHX, varTauSHXtoOTSG, 5000);
      P_SecondaryLoop = m_dot_pFluid_SHX*cP_pFluid_SHX*(T_in_pFluid_SHX.T - T_PN4_SHX);
      P_TertiaryLoop = m_dot_sFluid_SHX*cP_sFluid_SHX*(T_SN4_SHX - T_in_sFluid_SHX.T);
      annotation(
        Diagram(graphics = {Rectangle(origin = {30.2864, 0.0340212}, extent = {{-149.242, 59.993}, {149.242, -59.993}}), Rectangle(origin = {-60, 30}, lineColor = {255, 170, 0}, fillColor = {255, 170, 0}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {120, 30}, lineColor = {255, 170, 0}, fillColor = {255, 170, 0}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {0, 30}, lineColor = {255, 170, 0}, fillColor = {255, 170, 0}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {-59.95, -29.61}, lineColor = {170, 85, 255}, fillColor = {170, 85, 255}, fillPattern = FillPattern.Solid, extent = {{-19.63, -9.61}, {19.63, 9.61}}), Rectangle(origin = {-0.12, -29.82}, lineColor = {170, 85, 255}, fillColor = {170, 85, 255}, fillPattern = FillPattern.Solid, extent = {{-19.88, 9.82}, {19.88, -9.82}}), Rectangle(origin = {59.8, -29.85}, lineColor = {170, 85, 255}, fillColor = {170, 85, 255}, fillPattern = FillPattern.Solid, extent = {{-19.8, 9.85}, {19.8, -9.85}}), Rectangle(origin = {119.85, -29.89}, lineColor = {170, 85, 255}, fillColor = {170, 85, 255}, fillPattern = FillPattern.Solid, extent = {{-19.85, 9.89}, {19.85, -9.89}}), Rectangle(origin = {89.92, -1.05}, lineColor = {136, 138, 133}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-49.92, 8.83}, {49.92, -8.83}}), Rectangle(origin = {60, 30}, lineColor = {255, 170, 0}, fillColor = {255, 170, 0}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {-30.08, -1.05}, lineColor = {136, 138, 133}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-49.92, 8.83}, {49.92, -8.83}}), Text(origin = {-100, 8}, extent = {{-14, 6}, {14, -6}}, textString = "Tp_IN"), Text(origin = {-100, -50}, extent = {{-14, 6}, {14, -6}}, textString = "Ts_OUT"), Text(origin = {160, -50}, extent = {{-14, 6}, {14, -6}}, textString = "Ts_IN"), Text(origin = {162, 12}, extent = {{-14, 6}, {14, -6}}, textString = "Tp_OUT"), Text(origin = {147, 50}, extent = {{59, -14}, {-59, 14}}, textString = "SHX")}, coordinateSystem(extent = {{-120, 80}, {220, -80}})),
        Icon(graphics = {Rectangle(origin = {30.2864, 0.0340212}, extent = {{-149.242, 59.993}, {149.242, -59.993}}), Rectangle(origin = {-60, 30}, lineColor = {255, 170, 0}, fillColor = {255, 170, 0}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {120, 30}, lineColor = {255, 170, 0}, fillColor = {255, 170, 0}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {0, 30}, lineColor = {255, 170, 0}, fillColor = {255, 170, 0}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {-59.95, -29.61}, lineColor = {170, 85, 255}, fillColor = {170, 85, 255}, fillPattern = FillPattern.Solid, extent = {{-19.63, -9.61}, {19.63, 9.61}}), Rectangle(origin = {-0.12, -29.82}, lineColor = {170, 85, 255}, fillColor = {170, 85, 255}, fillPattern = FillPattern.Solid, extent = {{-19.88, 9.82}, {19.88, -9.82}}), Rectangle(origin = {59.8, -29.85}, lineColor = {170, 85, 255}, fillColor = {170, 85, 255}, fillPattern = FillPattern.Solid, extent = {{-19.8, 9.85}, {19.8, -9.85}}), Rectangle(origin = {119.85, -29.89}, lineColor = {170, 85, 255}, fillColor = {170, 85, 255}, fillPattern = FillPattern.Solid, extent = {{-19.85, 9.89}, {19.85, -9.89}}), Rectangle(origin = {89.92, -1.05}, lineColor = {136, 138, 133}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-49.92, 8.83}, {49.92, -8.83}}), Rectangle(origin = {60, 30}, lineColor = {255, 170, 0}, fillColor = {255, 170, 0}, fillPattern = FillPattern.Solid, extent = {{-20, 10}, {20, -10}}), Rectangle(origin = {-30.08, -1.05}, lineColor = {136, 138, 133}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-49.92, 8.83}, {49.92, -8.83}}), Text(origin = {-100, 8}, extent = {{-14, 6}, {14, -6}}, textString = "Tp_IN"), Text(origin = {-100, -50}, extent = {{-14, 6}, {14, -6}}, textString = "Ts_OUT"), Text(origin = {160, -50}, extent = {{-14, 6}, {14, -6}}, textString = "Ts_IN"), Text(origin = {162, 12}, extent = {{-14, 6}, {14, -6}}, textString = "Tp_OUT"), Text(origin = {147, 50}, extent = {{59, -14}, {-59, 14}}, textString = "SHX")}),
        experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002));
    end SecondaryHeatExchanger;

    model UHX
      parameter MCFR_Packages.Units.Mass m_PN1UHX;
      parameter MCFR_Packages.Units.ReactorPower UHXP;
      parameter MCFR_Packages.Units.SpecificHeatCapacity cP_pFluidUHX = 0.00154912;
      parameter MCFR_Packages.Units.MassFlowRate m_dot_pFluidUHXnom = 12910.56;
      parameter MCFR_Packages.Units.ResidentTime tauUHXtoSHX = 1;
      parameter MCFR_Packages.Units.InitiationTime tripUHX;
      parameter MCFR_Packages.Units.Demand SetDemand;
      MCFR_Packages.Units.Temperature UHXtemp;
      MCFR_Packages.Units.ResidentTime varTauUHXtoSHX;
      MCFR_Packages.Units.Demand PowerDemand;
      MCFR_Packages.Units.MassFlowRate m_dot_pFluidUHX;
      input MCFR_Packages.Connectors.Temp_In UHXtemp_In annotation(
        Placement(visible = true, transformation(origin = {-70, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-70, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output MCFR_Packages.Connectors.Temp_Out UHXtemp_Out annotation(
        Placement(visible = true, transformation(origin = {26, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {26, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input MCFR_Packages.Connectors.FlowFraction_In tertiaryFF annotation(
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
  end HeatTransport;

  package Nuclear
    // Parameters and equations for reactivity effects (mPKE), decay heat, and poisons

    model mPKE
      // MSR Point Kinetic Equations
      parameter MCFR_Packages.Units.NeutronGenerationTime Lam;
      parameter MCFR_Packages.Units.PrecursorDecayConstant lam[6];
      parameter MCFR_Packages.Units.DelayedNeutronFrac beta[6];
      parameter MCFR_Packages.Units.ResidentTime tauCoreNom = 7.028; // Avg travel time within core
      parameter MCFR_Packages.Units.ResidentTime tauLoopNom = 7.028; // Avg travel time from core outlet to inlet
      //Variable declaration
      MCFR_Packages.Units.ResidentTime varTauCore;
      MCFR_Packages.Units.ResidentTime varTauLoop;
      MCFR_Packages.Units.PrecursorConc CG1;
      MCFR_Packages.Units.PrecursorConc CG2;
      MCFR_Packages.Units.PrecursorConc CG3;
      MCFR_Packages.Units.PrecursorConc CG4;
      MCFR_Packages.Units.PrecursorConc CG5;
      MCFR_Packages.Units.PrecursorConc CG6;
      MCFR_Packages.Units.PrecursorReturn CG1Return;
      MCFR_Packages.Units.PrecursorReturn CG2Return;
      MCFR_Packages.Units.PrecursorReturn CG3Return;
      MCFR_Packages.Units.PrecursorReturn CG4Return;
      MCFR_Packages.Units.PrecursorReturn CG5Return;
      MCFR_Packages.Units.PrecursorReturn CG6Return;
      MCFR_Packages.Units.Reactivity reactivity;
      MCFR_Packages.Units.DelayedNeutronFrac bterm[6];
      MCFR_Packages.Units.Reactivity rho0;
      input MCFR_Packages.Connectors.Reactivity feedback annotation(
        Placement(visible = true, transformation(origin = {0, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output MCFR_Packages.Connectors.NominalNeutronPopulation n_population annotation(
        Placement(visible = true, transformation(origin = {0, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input MCFR_Packages.Connectors.FlowFraction_In fuelFlowFrac annotation(
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
//n_population.n = 1;
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
      // Decay Heat in fuel salt equations
      parameter MCFR_Packages.Units.ReactorPower P;
      parameter MCFR_Packages.Units.DecayHeatYield DHYG1 = 9.635981959409105e-01;
      parameter MCFR_Packages.Units.DecayHeatYield DHYG2 = 3.560674858154914e-02;
      parameter MCFR_Packages.Units.DecayHeatYield DHYG3 = 7.950554775404400e-04;
      parameter MCFR_Packages.Units.Mass TotalFuelMass;
      parameter MCFR_Packages.Units.DecayHeatPrecursorDecayConstant DHlamG1 = 0.0945298;
      parameter MCFR_Packages.Units.DecayHeatPrecursorDecayConstant DHlamG2 = 0.00441957;
      parameter MCFR_Packages.Units.DecayHeatPrecursorDecayConstant DHlamG3 = 8.60979e-05;
      parameter MCFR_Packages.Units.FissFactor FissFactor = 2.464783802008740e-03; // Heat per fission relative to nominal power rating
      MCFR_Packages.Units.HeatTransferFraction DHG1;
      MCFR_Packages.Units.HeatTransferFraction DHG2;
      MCFR_Packages.Units.HeatTransferFraction DHG3;
      input MCFR_Packages.Connectors.NominalNeutronPopulation nPop annotation(
        Placement(visible = true, transformation(origin = {1, 33}, extent = {{-19, -19}, {19, 19}}, rotation = 0), iconTransformation(origin = {-32, 22}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output MCFR_Packages.Connectors.DecayHeat_Out decayHeat_Out annotation(
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
      parameter MCFR_Packages.Units.IsotopicFissYield gamma_I = 0.06390; // Fraction of fissions that directly produce this isotope
      parameter MCFR_Packages.Units.IsotopicFissYield gamma_Xe = 0.0025761;
      parameter MCFR_Packages.Units.IsotopicFissYield gamma_Pm = 0.010816;
      parameter MCFR_Packages.Units.IsotopicFissYield gamma_Sm = 0.0;
      parameter MCFR_Packages.Units.IsotopicDecayConstant I135_lam = 2.926153244511758e-05;
      parameter MCFR_Packages.Units.IsotopicDecayConstant Xe135_lam = 2.106574217602556e-05;
      parameter MCFR_Packages.Units.IsotopicDecayConstant Pm149_lam = 3.627371580423393e-06;
      parameter MCFR_Packages.Units.MicroAbsorptionCrossSection sig_Xe = 2.66449e-18;
      parameter MCFR_Packages.Units.MicroAbsorptionCrossSection sig_Sm = 4.032e-20;
      parameter MCFR_Packages.Units.MacroAbsorptionCrossSection Sig_f = 0.00129788; // Fission cross section of fuel
      parameter MCFR_Packages.Units.MacroAbsorptionCrossSection Sig_a = 0.00302062; // Absorption cross section of fuel
      parameter MCFR_Packages.Units.NominalNeutronFlux phi_0 = 3.51513E+14; // Flux density in core
      parameter MCFR_Packages.Units.OffgasRemovalRate lam_bubble; // Removal rate of fission gases via bubble removal
      //Variable declarations
      MCFR_Packages.Units.IsotopicConc I135_conc;
      MCFR_Packages.Units.IsotopicConc Xe135_conc;
      MCFR_Packages.Units.IsotopicConc Pm149_conc;
      MCFR_Packages.Units.IsotopicConc Sm149_conc;
      MCFR_Packages.Units.IsotopicConc Xe135_0;
      MCFR_Packages.Units.IsotopicConc Sm149_0;
      MCFR_Packages.Units.Reactivity Xe_react;
      MCFR_Packages.Units.Reactivity Sm_react;
      input MCFR_Packages.Connectors.NominalNeutronPopulation n_population annotation(
        Placement(visible = true, transformation(origin = {-40, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-40, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output MCFR_Packages.Connectors.Reactivity Poison_react annotation(
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
      parameter MCFR_Packages.Units.Temperature FuelTempSetPointNode1 = 700;
      parameter MCFR_Packages.Units.Temperature FuelTempSetPointNode2 = 720;
      parameter MCFR_Packages.Units.Temperature ReflcTempSetPoint = 700;
      parameter MCFR_Packages.Units.TemperatureReactivityCoef a_F; // Fuel temp feedback coefficient
      parameter MCFR_Packages.Units.TemperatureReactivityCoef a_R; // Reflector temp feedback coefficient
      parameter MCFR_Packages.Units.Reactivity step_mag;
      parameter MCFR_Packages.Units.Frequency omega;
      parameter MCFR_Packages.Units.Reactivity sin_mag;
      parameter MCFR_Packages.Units.InitiationTime stepInsertionTime;
      parameter MCFR_Packages.Units.InitiationTime sinInsertionTime;
      MCFR_Packages.Units.Reactivity FuelTempFeedbackNode1;
      MCFR_Packages.Units.Reactivity FuelTempFeedbackNode2;
      MCFR_Packages.Units.Reactivity ReflcTempFeedback;
      MCFR_Packages.Units.Reactivity TotalTempFeedback;
      MCFR_Packages.Units.Reactivity ReactExternalStep;
      MCFR_Packages.Units.Reactivity ReactExternalSin;
      input MCFR_Packages.Connectors.Temp_In fuelNode1 annotation(
        Placement(visible = true, transformation(origin = {-40, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-40, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input MCFR_Packages.Connectors.Temp_In fuelNode2 annotation(
        Placement(visible = true, transformation(origin = {-40, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-40, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input MCFR_Packages.Connectors.Temp_In ReflcNode annotation(
        Placement(visible = true, transformation(origin = {-40, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-40, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      output MCFR_Packages.Connectors.Reactivity feedback annotation(
        Placement(visible = true, transformation(origin = {26, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {26, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      input MCFR_Packages.Connectors.Reactivity Poison_react annotation(
        Placement(visible = true, transformation(origin = {26, -32}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {26, -32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    initial equation
      ReactExternalStep = 0;
      ReactExternalSin = 0;
      FuelTempFeedbackNode1 = 0;
      FuelTempFeedbackNode2 = 0;
      ReflcTempFeedback = 0;
      TotalTempFeedback = 0;
    equation
      ReactExternalStep = 0 + (if time < stepInsertionTime then 0 else step_mag);
      ReactExternalSin = 0 + (if time < sinInsertionTime then 0 else sin_mag*sin(omega*time));
      FuelTempFeedbackNode1 = (fuelNode1.T - FuelTempSetPointNode1)*(a_F/2);
      FuelTempFeedbackNode2 = (fuelNode2.T - FuelTempSetPointNode2)*(a_F/2);
      ReflcTempFeedback = (ReflcNode.T - ReflcTempSetPoint)*a_R;
      TotalTempFeedback = FuelTempFeedbackNode1 + FuelTempFeedbackNode2 + ReflcTempFeedback + ReactExternalSin + ReactExternalStep + Poison_react.rho;
      feedback.rho = TotalTempFeedback;
      annotation(
        Diagram(graphics = {Rectangle(origin = {0, 2}, extent = {{-50, 50}, {50, -50}}), Text(origin = {13, 40}, extent = {{-31, 24}, {31, -24}}, textString = "Feedback"), Text(origin = {-32, 22}, extent = {{-16, 6}, {16, -6}}, textString = "Temp_FuelN1"), Text(origin = {-32, -12}, extent = {{-16, 6}, {16, -6}}, textString = "Temp_FuelN2"), Text(origin = {-32, -40}, extent = {{-16, 6}, {16, -6}}, textString = "Temp_GrapN1"), Text(origin = {26, 10}, extent = {{-16, 6}, {16, -6}}, textString = "Total_Feedback"), Text(origin = {26, -40}, extent = {{-12, 4}, {12, -4}}, textString = "Poison React")}),
        Icon(graphics = {Rectangle(origin = {0, 2}, extent = {{-50, 50}, {50, -50}}), Text(origin = {13, 40}, extent = {{-31, 24}, {31, -24}}, textString = "Feedback"), Text(origin = {-32, 22}, extent = {{-16, 6}, {16, -6}}, textString = "Temp_FuelN1"), Text(origin = {-32, -12}, extent = {{-16, 6}, {16, -6}}, textString = "Temp_FuelN2"), Text(origin = {-32, -40}, extent = {{-16, 6}, {16, -6}}, textString = "Temp_GrapN1"), Text(origin = {26, 10}, extent = {{-16, 6}, {16, -6}}, textString = "Total_Feedback"), Text(origin = {26, -40}, extent = {{-12, 4}, {12, -4}}, textString = "Poison React")}));
    end ReactivityFeedback;
  end Nuclear;

  package Pumps
    model PrimaryLoopPump
      // Fuel Salt Pump parameters and equations
      parameter MCFR_Packages.Units.FlowFraction freeConvectionFF;
      parameter MCFR_Packages.Units.PumpConstant primaryPumpK;
      parameter MCFR_Packages.Units.InitiationTime tripPrimaryPump;
      output MCFR_Packages.Connectors.FlowFraction_Out primaryFF annotation(
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
      parameter MCFR_Packages.Units.FlowFraction freeConvectionFF;
      parameter MCFR_Packages.Units.PumpConstant secondaryPumpK;
      parameter MCFR_Packages.Units.InitiationTime tripSecondaryPump;
      output MCFR_Packages.Connectors.FlowFraction_Out secondaryFF annotation(
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
      parameter MCFR_Packages.Units.FlowFraction freeConvectionFF;
      parameter MCFR_Packages.Units.PumpConstant tertiaryPumpK;
      parameter MCFR_Packages.Units.InitiationTime tripTertiaryPump;
      output MCFR_Packages.Connectors.FlowFraction_Out tertiaryFF annotation(
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
    type SpecificHeat = Real(unit = "MJ/(kg.s)");
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
    type DecayHeatFraction = Real(tinits = "1", min = 0);
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
    type Number = Real(unit = "1", min = 0);
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
    type K_sc = Real(unit = "1");
    type SpecificEnthalpy = Real(unit = "MJ/kg");
  end Units;

  package Connectors
    // Creates all inlets and outlets

    connector Temp_In
      MCFR_Packages.Units.Temperature T;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, fillPattern = FillPattern.Solid, extent = {{-51, 47}, {51, -47}}), Ellipse(origin = {-8, 36}, extent = {{0, -2}, {0, 2}})}),
        Icon(graphics = {Ellipse(origin = {1, -1}, fillPattern = FillPattern.Solid, extent = {{-51, 47}, {51, -47}})}));
    end Temp_In;

    connector Temp_Out
      MCFR_Packages.Units.Temperature T;
      annotation(
        Icon(graphics = {Ellipse(origin = {1, -1}, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}),
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}));
    end Temp_Out;

    connector DecayHeat_In
      MCFR_Packages.Units.SpecificHeat Q;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineColor = {0, 225, 255}, fillColor = {0, 225, 255}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}),
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {0, 225, 255}, fillColor = {0, 225, 255}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}));
    end DecayHeat_In;

    connector DecayHeat_Out
      MCFR_Packages.Units.SpecificHeat Q;
      annotation(
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {0, 225, 255}, fillColor = {0, 225, 255}, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}),
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineColor = {0, 225, 255}, fillColor = {0, 225, 255}, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}));
    end DecayHeat_Out;

    connector NominalNeutronPopulation
      MCFR_Packages.Units.NominalNeutronPopulation n;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}),
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}})}));
    end NominalNeutronPopulation;

    connector Reactivity
      MCFR_Packages.Units.Reactivity rho;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineColor = {78, 154, 6}, fillColor = {78, 154, 6}, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}})}),
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {78, 154, 6}, fillColor = {78, 154, 6}, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}})}));
    end Reactivity;

    connector ExternalReactivity
      MCFR_Packages.Units.Reactivity rho_Ex;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineColor = {78, 154, 6}, fillColor = {78, 154, 6}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}),
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {78, 154, 6}, fillColor = {78, 154, 6}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}));
    end ExternalReactivity;

    connector FlowFraction_In
      MCFR_Packages.Units.FlowFraction FF;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, fillColor = {245, 121, 0}, fillPattern = FillPattern.Solid, extent = {{-51, 47}, {51, -47}})}),
        Icon(graphics = {Ellipse(origin = {1, -1}, fillColor = {245, 121, 0}, fillPattern = FillPattern.Solid, extent = {{-51, 47}, {51, -47}})}));
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, fillColor = {245, 121, 0}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}),
        Icon(graphics = {Ellipse(origin = {1, -1}, fillColor = {245, 121, 0}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}));
    end FlowFraction_In;

    connector FlowFraction_Out
      MCFR_Packages.Units.FlowFraction FF;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineColor = {245, 121, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}),
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {245, 121, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}));
    end FlowFraction_Out;
  end Connectors;

  package ExampleSims
  // For test cases, performance issues, debugging, and etc

    model m_test
      parameter Modelica.Units.SI.Pressure Press = 12e6;
      //parameter Modelica.Units.SI.Temperature Temp = 596;
      //parameter SI.Time duration(min=0.0, start=2);
      Modelica.Units.SI.Temperature Temp;
      Modelica.Units.SI.SpecificEnergy h_test;
    replaceable package Medium = Modelica.Media.Water.StandardWater;
    equation
      //Temp = offset + (if time < startTime then 0 else if time < (startTime + duration) then (time - startTime)*height/duration else height);
      Temp = time + 400;
h_test = Medium.specificEnthalpy_pT(Press, Temp);
    end m_test;
  end ExampleSims;
end MCFR_Packages;
