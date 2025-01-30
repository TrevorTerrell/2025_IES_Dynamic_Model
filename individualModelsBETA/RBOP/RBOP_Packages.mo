package RBOP_Packages
  package Rankine_Cycle model
    BOP
      //======================
      //PARAMETERS:
      parameter RBOP_Packages.Units.Pressure P_HP = 12.5;
      //INPUT FROM OTSG
      parameter RBOP_Packages.Units.Pressure P_LP = 1.15;
      parameter RBOP_Packages.Units.Pressure P_co = 0.008;
      parameter RBOP_Packages.Units.Pressure P_atm = 0.101325;
      parameter RBOP_Packages.Units.MassFlowRate m_dot_OTSG = 192.5;
      //INPUT FROM OTSG
      parameter RBOP_Packages.Units.Temperature T_OTSG_initial = 600;
      //INPUT FROM OTSG
      parameter RBOP_Packages.Units.Number K_RH = 0.2;
      //0.16;
      parameter RBOP_Packages.Units.Number K_HPB = 0.02;
      //0.1634;
      parameter RBOP_Packages.Units.Number K_LPB_initial = 0.12;
      //0.2174;
      //parameter RBOP_Packages.Units.SpecificEnthalpy H_HPT_in = Water.specificEnthalpy_pT(P_HP*1e6, T_OTSG+273.15)/1e6;
      parameter RBOP_Packages.Units.SpecificEnthalpy H_HPT_in_initial = Water.specificEnthalpy_pT(P_HP*1e6, T_OTSG_initial + 273.15)/1e6;
      parameter RBOP_Packages.Units.Temperature T_setpoint = 212;
      parameter RBOP_Packages.Units.Temperature T_cool_in = 22;
      parameter RBOP_Packages.Units.MassFlowRate m_dot_cool = 20000;
      //May be up to 27700. Source: https://www.nuclear-power.com/nuclear-power-plant/turbine-generator-power-conversion-system/main-condenser-steam-condenser/
      //======================
      //VARIABLES:
      RBOP_Packages.Units.MassFlowRate m_dot_v;
      RBOP_Packages.Units.MassFlowRate m_dot_HPT_in;
      RBOP_Packages.Units.MassFlowRate m_dot_HPT_bleed;
      RBOP_Packages.Units.MassFlowRate m_dot_HPT_out;
      RBOP_Packages.Units.SpecificEnthalpy H_HPT_out;
      RBOP_Packages.Units.SpecificEnthalpy H_dHPT;
      RBOP_Packages.Units.SpecificEnthalpy H_HPT_bleed;
      RBOP_Packages.Units.HeatFlowRate Q_HPT;
      RBOP_Packages.Units.Temperature T_HPT_out;
      //===
      RBOP_Packages.Units.Temperature T_sat_MS;
      RBOP_Packages.Units.SpecificEnthalpy H_vap_MS;
      RBOP_Packages.Units.SpecificEnthalpy H_f_MS;
      RBOP_Packages.Units.Number X_stm;
      RBOP_Packages.Units.MassFlowRate m_dot_liq;
      RBOP_Packages.Units.MassFlowRate m_dot_MS_out;
      RBOP_Packages.Units.SpecificEnthalpy H_MS;
      RBOP_Packages.Units.Temperature T_MS;
      //===
      RBOP_Packages.Units.MassFlowRate m_dot_RH;
      RBOP_Packages.Units.MassFlowRate m_dot_vRH;
      RBOP_Packages.Units.HeatFlowRate Q_max_HP;
      RBOP_Packages.Units.HeatFlowRate Q_max_LP;
      RBOP_Packages.Units.HeatFlowRate Q_RH;
      RBOP_Packages.Units.SpecificEnthalpy H_RH;
      RBOP_Packages.Units.Temperature T_RH;
      RBOP_Packages.Units.SpecificEnthalpy H_vRH;
      RBOP_Packages.Units.Temperature T_vRH;
      //  RBOP_Packages.Units.SpecificEnthalpy H_test;
      //===
      RBOP_Packages.Units.MassFlowRate m_dot_LPT_in;
      RBOP_Packages.Units.MassFlowRate m_dot_LPT_bleed;
      RBOP_Packages.Units.MassFlowRate m_dot_LPT_out;
      RBOP_Packages.Units.SpecificEnthalpy H_LPT_out;
      RBOP_Packages.Units.SpecificEnthalpy H_dLPT;
      RBOP_Packages.Units.SpecificEnthalpy H_LPT_bleed;
      RBOP_Packages.Units.HeatFlowRate Q_LPT;
      RBOP_Packages.Units.Temperature T_LPT_out;
      //===
      RBOP_Packages.Units.Temperature T_sat_con;
      RBOP_Packages.Units.SpecificEnthalpy H_f_con;
      RBOP_Packages.Units.SpecificEnthalpy H_vap_con;
      RBOP_Packages.Units.Number X_inlet_con;
      RBOP_Packages.Units.HeatFlowRate Q_con;
      RBOP_Packages.Units.Temperature T_stm_in;
      RBOP_Packages.Units.SpecificEnthalpy H_cool_in;
      RBOP_Packages.Units.SpecificEnthalpy H_limit;
      RBOP_Packages.Units.HeatFlowRate Q_cool_max;
      RBOP_Packages.Units.SpecificEnthalpy H_cool_out;
      RBOP_Packages.Units.Temperature T_cool_out;
      RBOP_Packages.Units.SpecificEnthalpy H_con_out;
      RBOP_Packages.Units.HeatFlowRate Q_cool;
      RBOP_Packages.Units.MassFlowRate m_dot_con;
      //===
      RBOP_Packages.Units.MassFlowRate m_dot_fdw;
      RBOP_Packages.Units.SpecificEnthalpy H_setpoint;
      RBOP_Packages.Units.Temperature T_fdw;
      RBOP_Packages.Units.SpecificEnthalpy H_fdw;
      RBOP_Packages.Units.Mass M_fdw;
      RBOP_Packages.Units.Number K_LPB_check;
      RBOP_Packages.Units.Number K_LPB;
      //===
      RBOP_Packages.Units.Temperature T_OTSG;
      RBOP_Packages.Units.SpecificEnthalpy H_HPT_in;
      //======================
    initial equation
      T_OTSG = T_OTSG_initial;
      H_HPT_in = Water.specificEnthalpy_pT(P_HP*1e6, T_OTSG_initial + 273.15)/1e6;
//  K_LPB = K_LPB_initial;
//  Q_RH = 50;
    equation
//  T_OTSG = T_OTSG_initial + (if time < 2000 then 0 elseif time <= 4000 then 0.1*(time - 2000) else 200);
//  T_OTSG = T_OTSG_initial + (if time < 2000 then 0 elseif time <= 4000 then 100 else -100);
//  T_OTSG = T_OTSG_initial + (if time < 1000 then 0 else 150*sin(0.001*time));
      T_OTSG = T_OTSG_initial + 100*sin(0.001*time);
//T_OTSG = T_OTSG_initial;
      H_HPT_in = Water.specificEnthalpy_pT(P_HP*1e6, T_OTSG + 273.15)/1e6;
//======================
//VALVE AND NOZZLE CHEST EQNS:
      m_dot_v = m_dot_OTSG*K_RH;
//======================
//HIGH PRESSURE TURBINE EQNS:
      m_dot_HPT_in = m_dot_OTSG - m_dot_v;
      m_dot_HPT_bleed = K_HPB*m_dot_HPT_in;
      m_dot_HPT_out = m_dot_HPT_in - m_dot_HPT_bleed;
      H_HPT_out = Water.isentropicEnthalpy(P_LP*1e6, Water.setState_phX(P_HP*1e6, H_HPT_in*1e6))/1e6;
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
//  m_dot_vRH = 0;
//THIS IS SATURATED STEAM HEADED TO THE FWH
//  log(T_vRH_sat + 273.15) = (9.37817e-3 + 4.98951e-4*(P_OTSG/22.064) + 1.11049e-5*(P_OTSG/22.064)^2 + 3.34995e-7*(P_OTSG/22.064)^3 + 3.44102e-8*(P_OTSG/22.064)^4)^(-0.4);
//  T_vRH_sat = 327.82;
// ~~~~~~~~~~~~~
      if T_MS < 322 then
        Q_max_HP = m_dot_vRH*(Water.specificEnthalpy_pT(P_HP*1e6, T_OTSG + 273.15)/1e6 - Water.specificEnthalpy_pT(P_HP*1e6, (327.82) + 273.15)/1e6) + 1.4*((322 - T_MS)/(T_OTSG - T_MS))*(T_vRH - T_MS);
// + 1.0*((322-T_MS)/(T_OTSG-T_MS))*(327.82-T_MS);
      else
        Q_max_HP = m_dot_vRH*(Water.specificEnthalpy_pT(P_HP*1e6, T_OTSG + 273.15)/1e6 - Water.specificEnthalpy_pT(P_HP*1e6, (T_MS + 5.82) + 273.15)/1e6);
// + 0.1*(5.82);
      end if;
//  Q_max_HP = m_dot_vRH*(Water.specificEnthalpy_pT(P_HP*1e6, T_OTSG + 273.15)/1e6 - Water.specificEnthalpy_pT(P_HP*1e6, (T_MS + 10) + 273.15)/1e6);0.2
      Q_max_LP = m_dot_RH*(Water.specificEnthalpy_pT(P_LP*1e6, (T_OTSG - 10) + 273.15)/1e6 - Water.specificEnthalpy_pT(P_LP*1e6, T_MS + 273.15)/1e6);
      if Q_max_HP < Q_max_LP then
        Q_RH = Q_max_HP;
        H_RH = Water.specificEnthalpy_pT(P_LP*1e6, T_MS + 273.15)/1e6 + Q_RH/m_dot_RH;
        T_RH = Water.temperature_ph(P_LP*1e6, H_RH*1e6) - 273.15;
        H_vRH = Water.specificEnthalpy_pT(P_HP*1e6, T_OTSG + 273.15)/1e6 - Q_RH/m_dot_vRH;
        T_vRH = Water.temperature_ph(P_HP*1e6, H_vRH*1e6) - 273.15;
      else
        Q_RH = Q_max_LP;
        H_RH = Water.specificEnthalpy_pT(P_LP*1e6, T_MS + 273.15)/1e6 + Q_RH/m_dot_RH;
        T_RH = Water.temperature_ph(P_LP*1e6, H_RH*1e6) - 273.15;
        H_vRH = Water.specificEnthalpy_pT(P_HP*1e6, T_OTSG + 273.15)/1e6 - Q_RH/m_dot_vRH;
        T_vRH = Water.temperature_ph(P_HP*1e6, H_vRH*1e6) - 273.15;
      end if;
//  H_test = Water.specificEnthalpy_pT(P_HP*1e6, T_vRH + 273.15)/1e6;
// ~~~~~~~~~~~~~
//  Q_RH = 0.1183 * (T_OTSG - T_MS);
////  Q_RH = 0.26*((T_OTSG - T_RH) - (T_vRH - T_MS)) / log((T_OTSG - T_RH)/(T_vRH - T_MS));
////  Q_RH = 47.83;
////  der(Q_RH) = (Cp_RH * (m_dot_v + m_dot_RH) * (T_OTSG - T_MS))/(2*tau_RH) - Q_RH/tau_RH;
////  der(Q_RH) = (0.0021520 * (m_dot_v + m_dot_RH) * (T_OTSG - T_MS))/(2*20) - Q_RH/20;
////  der(Q_RH) = (0.001 * (m_dot_v + m_dot_RH) * (T_OTSG - T_MS))/(2*20) - Q_RH/20;
//  H_RH = Water.specificEnthalpy_pT(P_LP*1e6, T_MS + 273.15)/1e6 + Q_RH/m_dot_RH;
//  H_vRH = Water.specificEnthalpy_pT(P_HP*1e6, T_OTSG + 273.15)/1e6 - Q_RH/m_dot_vRH;
//  T_RH = Water.temperature_ph(P_LP*1e6, H_RH*1e6) - 273.15;
//  T_vRH = Water.temperature_ph(P_HP*1e6, H_vRH*1e6) - 273.15;
// ~~~~~~~~~~~~~
//  H_RH = H_MS * m_dot_MS_out/m_dot_RH + H_HPT_in * m_dot_v/m_dot_RH;
//  T_RH = Water.temperature_ph(P_LP*1e6, H_RH*1e6)-273.15;
// NEED TO MAKE A METHOD FOR Q_RH WHERE AT T_MS=T_HP_sat-5 THE LOWER ENTHALPY LIMIT ON Q_max_HP IS H_g. AT LOWER T_MS, THE EXTRA HEAT FLOW OUT OF THE HP SIDE IS: Q_extra = m_dot_RH * Cp_MS (T_HP_sat-5 - T_MS)
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
//======================
//CONDENSER EQNS:
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
//New method of determining outlet flow based on enthalpy. Hot inlet and outlet enthalpy known, along with coolant inlet. Use coolant dH_co*mdot_in as demand, which is satisfied by the difference between coolant outlet (unknown) and inlet, which can meet it up to a certain fraction of the Hot inlet. This determines the flow rate as a ratio of how well the demand is satisfied times the inlet flow rate.
//======================
//FEEDWATER HEATER EQNS:
      m_dot_fdw = (m_dot_HPT_bleed + m_dot_LPT_bleed + m_dot_liq + m_dot_vRH + m_dot_con);
      H_setpoint = Water.specificEnthalpy_pT(P_HP*1e6, T_setpoint + 273.15)/1e6;
      K_LPB_check = (m_dot_fdw/(H_LPT_bleed*m_dot_LPT_in))*(H_setpoint - ((H_HPT_bleed*m_dot_HPT_bleed/m_dot_fdw) + (H_f_MS*m_dot_liq/m_dot_fdw) + (H_vRH*m_dot_vRH/m_dot_fdw) + (H_con_out*m_dot_con/m_dot_fdw)));
      if K_LPB_check > 0 then
        K_LPB = K_LPB_check;
      elseif K_LPB_check > 1 then
        K_LPB = 1;
      else
        K_LPB = 0;
      end if;
      T_fdw = Water.temperature_ph(P_HP*1e6, H_fdw*1e6) - 273.15;
      H_fdw = (H_HPT_bleed*m_dot_HPT_bleed/m_dot_fdw) + (H_LPT_bleed*m_dot_LPT_bleed/m_dot_fdw) + (H_f_MS*m_dot_liq/m_dot_fdw) + (H_vRH*m_dot_vRH/m_dot_fdw) + (H_con_out*m_dot_con/m_dot_fdw);
      der(M_fdw) = m_dot_HPT_bleed + m_dot_LPT_bleed + m_dot_liq + m_dot_vRH + m_dot_con - m_dot_fdw;
//======================
      annotation(
        experiment(StartTime = 0, StopTime = 5000, Tolerance = 1e-06, Interval = 0.1));    end BOP
    ;
  end Rankine_Cycle;

  package Connectors
    // Creates all inlets and outlets

    connector Temp_In
      RBOP_Packages.Units.Temperature T;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, fillPattern = FillPattern.Solid, extent = {{-51, 47}, {51, -47}}), Ellipse(origin = {-8, 36}, extent = {{0, -2}, {0, 2}})}),
        Icon(graphics = {Ellipse(origin = {1, -1}, fillPattern = FillPattern.Solid, extent = {{-51, 47}, {51, -47}})}));
    end Temp_In;

    connector Temp_Out
      RBOP_Packages.Units.Temperature T;
      annotation(
        Icon(graphics = {Ellipse(origin = {1, -1}, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}),
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}));
    end Temp_Out;

    connector DecayHeat_In
      RBOP_Packages.Units.SpecificHeat Q;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineColor = {0, 225, 255}, fillColor = {0, 225, 255}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}),
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {0, 225, 255}, fillColor = {0, 225, 255}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}));
    end DecayHeat_In;

    connector DecayHeat_Out
      RBOP_Packages.Units.SpecificHeat Q;
      annotation(
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {0, 225, 255}, fillColor = {0, 225, 255}, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}),
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineColor = {0, 225, 255}, fillColor = {0, 225, 255}, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}));
    end DecayHeat_Out;

    connector NominalNeutronPopulation
      RBOP_Packages.Units.NominalNeutronPopulation n;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}),
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {20, 36, 248}, fillColor = {20, 36, 248}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}})}));
    end NominalNeutronPopulation;

    connector Reactivity
      RBOP_Packages.Units.Reactivity rho;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineColor = {78, 154, 6}, fillColor = {78, 154, 6}, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}})}),
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {78, 154, 6}, fillColor = {78, 154, 6}, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}})}));
    end Reactivity;

    connector ExternalReactivity
      RBOP_Packages.Units.Reactivity rho_Ex;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineColor = {78, 154, 6}, fillColor = {78, 154, 6}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}),
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {78, 154, 6}, fillColor = {78, 154, 6}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-51, 47}, {51, -47}}, endAngle = 360)}));
    end ExternalReactivity;

    connector FlowFraction_In
      RBOP_Packages.Units.FlowFraction FF;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, fillColor = {245, 121, 0}, fillPattern = FillPattern.Solid, extent = {{-51, 47}, {51, -47}})}),
        Icon(graphics = {Ellipse(origin = {1, -1}, fillColor = {245, 121, 0}, fillPattern = FillPattern.Solid, extent = {{-51, 47}, {51, -47}})}));
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, fillColor = {245, 121, 0}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}),
        Icon(graphics = {Ellipse(origin = {1, -1}, fillColor = {245, 121, 0}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}));
    end FlowFraction_In;

    connector FlowFraction_Out
      RBOP_Packages.Units.FlowFraction FF;
      annotation(
        Diagram(graphics = {Ellipse(origin = {1, -1}, lineColor = {245, 121, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}),
        Icon(graphics = {Ellipse(origin = {1, -1}, lineColor = {245, 121, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-51, 47}, {51, -47}})}));
    end FlowFraction_Out;
  end Connectors;

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
    type Temperature = Real(unit = "C");
    type AbsoluteTemperature = Real(unit = "K", min = 0);
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
    type Volume = Real(unit = "m^3", min = 0);
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

  replaceable package Water = Modelica.Media.Water.StandardWater;
  replaceable package MoistAir = Modelica.Media.Air.MoistAir;
end RBOP_Packages;
