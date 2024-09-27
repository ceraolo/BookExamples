package BookExamples "Package containing basic EV models"
  model FirstEV "Simulates a very basic Electric Vehicle"
    import Modelica;
    extends Modelica.Icons.Example;
    Modelica.Units.SI.Energy enP1(start = 0, fixed = true);
    Modelica.Units.SI.Energy enP2(start = 0, fixed = true);
    Modelica.Units.SI.Energy enP1Pos(start = 0, fixed = true);
    Modelica.Units.SI.Energy enP1Neg(start = 0, fixed = true);
    DragForce dragF(Cx = 0.26, S = 2.2, fc = 0.014, m = mass.m, rho(displayUnit = "kg/m3") = 1.225, v(start = 0, fixed = true)) annotation(
      Placement(transformation(origin = {80, -26}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    Modelica.Mechanics.Rotational.Components.IdealGear gear(ratio = 10, flange_b(phi(start = 0, fixed = true))) annotation(
      Placement(transformation(origin = {-22, 20}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Mechanics.Rotational.Sources.Torque torque annotation(
      Placement(visible = true, transformation(extent = {{-82, 10}, {-62, 30}}, rotation = 0)));
    PropDriver driver(CycleFileName = Modelica.Utilities.Files.loadResource("modelica://BookExamples/Resources/Sort1.txt"), extrapolation = Modelica.Blocks.Types.Extrapolation.LastTwoPoints, k = 1000, yMax = 100000.0) annotation(
      Placement(transformation(origin = {0, -4}, extent = {{-118, 14}, {-98, 34}})));
    Modelica.Mechanics.Translational.Components.Mass mass(m = 1500) annotation(
      Placement(transformation(origin = {4, 0}, extent = {{34, 10}, {54, 30}})));
    Modelica.Mechanics.Translational.Sensors.SpeedSensor velSens annotation(
      Placement(transformation(origin = {58, -8}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
    Modelica.Mechanics.Rotational.Components.IdealRollingWheel wheel(radius = 0.31) annotation(
      Placement(transformation(origin = {-12, 0}, extent = {{4, 10}, {24, 30}})));
    Modelica.Mechanics.Rotational.Components.Inertia inertia(J = 0.1) annotation(
      Placement(transformation(origin = {-4, 0}, extent = {{-54, 10}, {-34, 30}})));
    Modelica.Blocks.Nonlinear.Limiter to_mP1Pos(uMax = 1e99, uMin = 0) annotation(
      Placement(transformation(origin = {8, -6}, extent = {{18, -12}, {30, 0}})));
    Modelica.Blocks.Nonlinear.Limiter to_mP1Neg(uMax = 0, uMin = -1e99) annotation(
      Placement(transformation(origin = {8, -12}, extent = {{6, -6}, {-6, 6}})));
    Modelica.Mechanics.Translational.Sensors.PowerSensor mP1 annotation(
      Placement(transformation(origin = {26, 20}, extent = {{-8, -8}, {8, 8}})));
    Modelica.Mechanics.Translational.Sensors.PowerSensor mP2 annotation(
      Placement(transformation(origin = {80, 4}, extent = {{-8, -8}, {8, 8}}, rotation = -90)));
  equation
    der(enP1) = mP1.power;
    der(enP2) = mP2.power;
    der(enP1Pos) = to_mP1Pos.y;
    der(enP1Neg) = -to_mP1Neg.y;
    connect(torque.tau, driver.tauRef) annotation(
      Line(points = {{-84, 20}, {-97, 20}}, color = {0, 0, 127}));
    connect(driver.V, velSens.v) annotation(
      Line(points = {{-108, 9}, {-108, -36}, {58, -36}, {58, -19}}, color = {0, 0, 127}));
    connect(velSens.flange, mass.flange_b) annotation(
      Line(points = {{58, 2}, {58, 20}}, color = {0, 127, 0}));
    connect(inertia.flange_a, torque.flange) annotation(
      Line(points = {{-58, 20}, {-62, 20}}));
    connect(inertia.flange_b, gear.flange_a) annotation(
      Line(points = {{-38, 20}, {-32, 20}}));
    connect(gear.flange_b, wheel.flangeR) annotation(
      Line(points = {{-12, 20}, {-8, 20}}));
    connect(mP1.power, to_mP1Neg.u) annotation(
      Line(points = {{20, 11}, {20, -12}, {15, -12}}, color = {0, 0, 127}));
    connect(to_mP1Pos.u, mP1.power) annotation(
      Line(points = {{25, -12}, {20, -12}, {20, 11}}, color = {0, 0, 127}));
    connect(wheel.flangeT, mP1.flange_a) annotation(
      Line(points = {{12, 20}, {18, 20}}, color = {0, 127, 0}));
    connect(mass.flange_a, mP1.flange_b) annotation(
      Line(points = {{38, 20}, {34, 20}}, color = {0, 127, 0}));
    connect(mass.flange_b, mP2.flange_a) annotation(
      Line(points = {{58, 20}, {80, 20}, {80, 12}}, color = {0, 127, 0}));
    connect(mP2.flange_b, dragF.flange) annotation(
      Line(points = {{80, -4}, {80, -16}}, color = {0, 127, 0}));
    annotation(
      Documentation(info = "<html>
  <p>Very basic introductory EV model</p>
  </html>"),
      Diagram(coordinateSystem(extent = {{-120, 40}, {100, -40}}, preserveAspectRatio = false), graphics = {Rectangle(origin = {-14, 0}, lineColor = {28, 108, 200}, pattern = LinePattern.Dash, extent = {{-76, 36}, {-21, 4}}), Text(origin = {-8, 0}, textColor = {28, 108, 200}, extent = {{-82, 2}, {-26, -4}}, textString = "electric drive")}),
      Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = false, initialScale = 0.1, grid = {2, 2})),
      experiment(StartTime = 0, StopTime = 200, Tolerance = 0.0001, Interval = 0.1),
      __OpenModelica_simulationFlags(jacobian = "", s = "dassl", lv = "LOG_STATS"),
      __OpenModelica_commandLineOptions = "");
  end FirstEV;

  model SecondEV "Simulates a very basic Electric Vehicle"
    import Modelica;
    extends Modelica.Icons.Example;
    Modelica.Units.SI.Energy enBatDel(start = 0, fixed = true);
    Modelica.Units.SI.Energy enDTdel(start = 0, fixed = true);
    Modelica.Units.SI.Energy enP1(start = 0, fixed = true);
    Modelica.Units.SI.Energy enP1Pos(start = 0, fixed = true);
    Modelica.Units.SI.Energy enBattLoss(start = 0, fixed = true);
    Modelica.Units.SI.Energy enBraking(start = 0, fixed = true);
    Modelica.Mechanics.Rotational.Components.LossyGear gear(ratio = 7.94, flange_b(phi(start = 0, fixed = true)), lossTable = [0, 0.95, 0.95, 1.3, 1.3; 1088, 0.95, 0.95, 1.3, 1.3]) annotation(
      Placement(transformation(origin = {-36, 14}, extent = {{-10, -10}, {10, 10}})));
    EHPTlib.SupportModels.Miscellaneous.PropDriver driver(CycleFileName = Modelica.Utilities.Files.loadResource("modelica://BookExamples/Resources/wltc3b.txt"), extrapolation = Modelica.Blocks.Types.Extrapolation.Periodic, yMax = 100000.0, k = 2000) annotation(
      Placement(transformation(origin = {2, -6}, extent = {{-116, -10}, {-96, 10}})));
    Modelica.Mechanics.Rotational.Components.IdealRollingWheel wheel(radius = 0.316) annotation(
      Placement(transformation(origin = {14, 0}, extent = {{-4, 4}, {16, 24}})));
    EHPTlib.MapBased.OneFlange eleDrive(J = 0.0347, effTableName = "effTable", mapsFileName = Modelica.Utilities.Files.loadResource("modelica://BookExamples/Resources/NissanMaps.txt"), mapsOnFile = true, tauMax = 280, wMax = 1088, powMax = 80e3, uDcNom = 300) "Electric Drive" annotation(
      Placement(visible = true, transformation(extent = {{-74, 6}, {-54, 24}}, rotation = 0)));
    EHPTlib.SupportModels.Miscellaneous.Batt1 batt1(SOCInit = 0.7, QCellNom = 33.1*3600, ns = 2*48, C1(v(start = 0, fixed = true)), ECellMin = 2.5, ECellMax = 4.2, ICellMax = 520, np = 2) annotation(
      Placement(transformation(origin = {-2, -20}, extent = {{-112, 34}, {-92, 54}})));
    Modelica.Electrical.Analog.Basic.Ground ground annotation(
      Placement(visible = true, transformation(extent = {{-84, -20}, {-64, 0}}, rotation = 0)));
    Modelica.Mechanics.Translational.Sensors.PowerSensor mP1 annotation(
      Placement(transformation(origin = {46, 14}, extent = {{-6, -8}, {6, 8}})));
    Modelica.Mechanics.Translational.Components.Mass mass(m = 1715, v(fixed = true)) annotation(
      Placement(transformation(origin = {4, 0}, extent = {{56, 4}, {76, 24}})));
    Modelica.Mechanics.Translational.Sensors.PowerSensor mP2 annotation(
      Placement(transformation(origin = {102, 4}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Mechanics.Translational.Sensors.SpeedSensor velSens annotation(
      Placement(transformation(origin = {80, -28}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
    EHPTlib.SupportModels.Miscellaneous.DragForce dragF(Cx = 0.29, rho(displayUnit = "kg/m3") = 1.225, S = 2.276, fc = 0.013, m = mass.m, v(start = 0, fixed = true)) annotation(
      Placement(transformation(origin = {102, -24}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    Modelica.Blocks.Nonlinear.Limiter to_mP1Pos(uMax = 1e99, uMin = 0) annotation(
      Placement(transformation(origin = {28, -8}, extent = {{18, -12}, {30, 0}})));
    Modelica.Blocks.Nonlinear.Limiter to_mP1Neg(uMax = 0, uMin = -1e99) annotation(
      Placement(transformation(origin = {28, -14}, extent = {{6, -6}, {-6, 6}})));
  Modelica.Mechanics.Rotational.Components.Inertia wheels_differential(J = 3.29) annotation(
      Placement(transformation(origin = {-38, -12}, extent = {{20, 16}, {40, 36}})));
  equation
    der(enBatDel) = (batt1.p.v - batt1.n.v)*batt1.n.i;
    der(enDTdel) = eleDrive.powSensor.power;
    der(enP1) = mP1.power;
    der(enP1Pos) = to_mP1Pos.y;
    der(enBattLoss) = batt1.powerLoss;
    der(enBraking) = if mP1.power > 0 then 0 else -mP1.power;
    connect(batt1.n, eleDrive.pin_n) annotation(
      Line(points = {{-94, 18}, {-80, 18}, {-80, 10}, {-74, 10}}, color = {0, 0, 255}));
    connect(eleDrive.pin_n, ground.p) annotation(
      Line(points = {{-74, 10}, {-74, 10}, {-74, 0}, {-74, 0}}, color = {0, 0, 255}));
    connect(eleDrive.tauRef, driver.tauRef) annotation(
      Line(points = {{-74.2, 14}, {-86, 14}, {-86, -6}, {-93, -6}}, color = {0, 0, 127}));
    connect(velSens.flange, mP2.flange_a) annotation(
      Line(points = {{80, -18}, {80, 14}, {102, 14}}, color = {0, 127, 0}));
    connect(driver.V, velSens.v) annotation(
      Line(points = {{-104, -17}, {-104, -39}, {80, -39}}, color = {0, 0, 127}));
    connect(mass.flange_a, mP1.flange_b) annotation(
      Line(points = {{60, 14}, {52, 14}}, color = {0, 127, 0}));
    connect(mP2.flange_a, mass.flange_b) annotation(
      Line(points = {{102, 14}, {80, 14}}, color = {0, 127, 0}));
    connect(mP1.flange_a, wheel.flangeT) annotation(
      Line(points = {{40, 14}, {30, 14}}, color = {0, 127, 0}));
    connect(dragF.flange, mP2.flange_b) annotation(
      Line(points = {{102, -14}, {102, -6}}, color = {0, 127, 0}));
    connect(to_mP1Neg.u, to_mP1Pos.u) annotation(
      Line(points = {{35.2, -14}, {43.2, -14}}, color = {0, 0, 127}));
    connect(mP1.power, to_mP1Pos.u) annotation(
      Line(points = {{41.2, 5.2}, {40.2, 5.2}, {40.2, -13.8}, {45.2, -13.8}}, color = {0, 0, 127}));
    connect(batt1.p, eleDrive.pin_p) annotation(
      Line(points = {{-94, 30}, {-76, 30}, {-76, 18}, {-74, 18}}, color = {0, 0, 255}));
    connect(eleDrive.flange_a, gear.flange_a) annotation(
      Line(points = {{-54, 14}, {-46, 14}}));
  connect(gear.flange_b, wheels_differential.flange_a) annotation(
      Line(points = {{-26, 14}, {-18, 14}}));
  connect(wheels_differential.flange_b, wheel.flangeR) annotation(
      Line(points = {{2, 14}, {10, 14}}));
    annotation(
      Documentation(info = "<html>
  <p>Simple map-based EV model with battery.</p>
  </html>"),
      Diagram(coordinateSystem(extent = {{-120, 40}, {120, -40}}, preserveAspectRatio = false)),
      Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = false, initialScale = 0.1, grid = {2, 2})),
      experiment(StartTime = 0, StopTime = 1800, Tolerance = 0.0001, Interval = 0.128571));
  end SecondEV;
  
  block PropDriver "Simple Proportional controller driver"
    parameter String CycleFileName = "cycleName.txt" "Drive Cycle Name ex: \"sort1.txt\"";
    parameter Real k "Controller gain";
    parameter Real yMax = 1.e6 "Max output value (absolute)";
    parameter Modelica.Blocks.Types.Extrapolation extrapolation = Modelica.Blocks.Types.Extrapolation.LastTwoPoints "Extrapolation of data outside the definition range";
  protected
    parameter Boolean dummy(fixed = false);
    //Used only to render exrapolation a structural parameter.
    // This is important, otherwise the user could falsely believe that they can change it at run time,
    //which is not possible because it is passed to a CombiTable for which (at least in OpenModelica),
    //the extrapolation parameter is structural.
  public
    Modelica.Blocks.Interfaces.RealInput V annotation(
      Placement(visible = true, transformation(origin = {0, -66}, extent = {{-14, -14}, {14, 14}}, rotation = 90), iconTransformation(origin = {0, -112}, extent = {{-12, -12}, {12, 12}}, rotation = 90)));
    Modelica.Blocks.Math.UnitConversions.From_kmh from_kmh annotation(
      Placement(visible = true, transformation(extent = {{-42, -10}, {-22, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.CombiTimeTable driveCyc(columns = {2}, extrapolation = extrapolation, fileName = CycleFileName, tableName = "Cycle", tableOnFile = true) annotation(
      Placement(visible = true, transformation(extent = {{-80, -10}, {-60, 10}}, rotation = 0)));
    Modelica.Blocks.Math.Feedback feedback annotation(
      Placement(visible = true, transformation(extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Math.Gain gain(k = k) annotation(
      Placement(visible = true, transformation(extent = {{14, -10}, {34, 10}}, rotation = 0)));
    Modelica.Blocks.Nonlinear.Limiter limAcc(uMax = yMax, uMin = 0) annotation(
      Placement(visible = true, transformation(origin = {2, 40}, extent = {{52, -10}, {72, 10}}, rotation = 0)));
    Modelica.Blocks.Nonlinear.Limiter limBrak(uMax = 0, uMin = -yMax) annotation(
      Placement(visible = true, transformation(origin = {0, -40}, extent = {{52, -10}, {72, 10}}, rotation = 0)));
    Modelica.Blocks.Interfaces.RealOutput tauRef(unit = "N.m") annotation(
      Placement(visible = true, transformation(extent = {{100, -10}, {120, 10}}, rotation = 0), iconTransformation(extent = {{100, -10}, {120, 10}}, rotation = 0)));
    Modelica.Blocks.Interfaces.RealOutput accelTau(unit = "N.m") annotation(
      Placement(visible = true, transformation(origin = {110, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(extent = {{100, 52}, {120, 72}}, rotation = 0)));
    Modelica.Blocks.Interfaces.RealOutput brakeTau(unit = "N.m") annotation(
      Placement(visible = true, transformation(origin = {110, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(extent = {{100, -70}, {120, -50}}, rotation = 0)));
    Modelica.Blocks.Nonlinear.Limiter limiter1(uMax = yMax) annotation(
      Placement(visible = true, transformation(origin = {4, 0}, extent = {{52, -10}, {72, 10}}, rotation = 0)));
  initial equation
  //For the meaning of the following if see the definition of dummy
    if extrapolation == Modelica.Blocks.Types.Extrapolation.HoldLastPoint then
      dummy = true;
    else
      dummy = false;
    end if;
  equation
    connect(V, feedback.u2) annotation(
      Line(points = {{0, -66}, {0, -66}, {0, -8}, {0, -8}}, color = {0, 0, 127}));
    connect(from_kmh.u, driveCyc.y[1]) annotation(
      Line(points = {{-44, 0}, {-59, 0}}, color = {0, 0, 127}));
    connect(from_kmh.y, feedback.u1) annotation(
      Line(points = {{-21, 0}, {-8, 0}}, color = {0, 0, 127}));
    connect(feedback.y, gain.u) annotation(
      Line(points = {{9, 0}, {12, 0}}, color = {0, 0, 127}));
    connect(limBrak.y, brakeTau) annotation(
      Line(points = {{73, -40}, {104, -40}, {104, -40}, {110, -40}}, color = {0, 0, 127}));
    connect(limAcc.y, accelTau) annotation(
      Line(points = {{75, 40}, {102, 40}, {102, 40}, {110, 40}}, color = {0, 0, 127}));
    connect(limBrak.u, gain.y) annotation(
      Line(points = {{50, -40}, {40, -40}, {40, 0}, {35, 0}, {35, 0}}, color = {0, 0, 127}));
    connect(limAcc.u, gain.y) annotation(
      Line(points = {{52, 40}, {40, 40}, {40, 0}, {35, 0}, {35, 0}}, color = {0, 0, 127}));
    connect(limiter1.u, gain.y) annotation(
      Line(points = {{54, 0}, {34, 0}, {34, 0}, {35, 0}}, color = {0, 0, 127}));
    connect(limiter1.y, tauRef) annotation(
      Line(points = {{77, 0}, {102, 0}, {102, 0}, {110, 0}}, color = {0, 0, 127}));
    annotation(
      Documentation(info = "<html><head></head><body><p>Simple driver model.</p><p>It reads a reference cycle from a file then controls speed with a simple proportional feedback law.</p>
              </body></html>"),
      Icon(coordinateSystem(preserveAspectRatio = false, initialScale = 0.1), graphics = {Rectangle(fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Ellipse(fillColor = {255, 213, 170}, fillPattern = FillPattern.Solid, extent = {{-23, 22}, {-12, -4}}, endAngle = 360), Text(origin = {2, -0.1894}, lineColor = {0, 0, 255}, extent = {{-104, 142.189}, {98, 104}}, textString = "%name"), Polygon(fillColor = {215, 215, 215}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, points = {{-22, -60}, {-42, -88}, {-16, -88}, {16, -88}, {-22, -60}}), Polygon(fillColor = {135, 135, 135}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, points = {{-32, 40}, {-62, -52}, {-30, -52}, {-30, -52}, {-32, 40}}, smooth = Smooth.Bezier), Polygon(fillColor = {135, 135, 135}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, points = {{-68, -36}, {-14, -90}, {10, -50}, {0, -50}, {-68, -36}}, smooth = Smooth.Bezier), Polygon(fillColor = {175, 175, 175}, fillPattern = FillPattern.Solid, points = {{-22, 10}, {-30, 6}, {-40, -48}, {2, -46}, {2, -34}, {0, 2}, {-22, 10}}, smooth = Smooth.Bezier), Ellipse(fillColor = {255, 213, 170}, fillPattern = FillPattern.Solid, extent = {{-30, 44}, {-3, 10}}, endAngle = 360), Polygon(pattern = LinePattern.None, fillPattern = FillPattern.Solid, points = {{-38, 34}, {-16, 50}, {-2, 36}, {4, 36}, {6, 36}, {-38, 34}}, smooth = Smooth.Bezier), Polygon(fillColor = {95, 95, 95}, fillPattern = FillPattern.Solid, points = {{30, -44}, {-32, -28}, {-36, -44}, {-24, -58}, {30, -44}}, smooth = Smooth.Bezier), Polygon(fillPattern = FillPattern.Solid, points = {{42, -70}, {36, -84}, {48, -78}, {52, -72}, {50, -68}, {42, -70}}, smooth = Smooth.Bezier), Line(points = {{48, -14}, {26, 0}, {26, 0}}, thickness = 0.5), Line(points = {{20, -10}, {34, 10}, {34, 10}}, thickness = 0.5), Polygon(fillColor = {255, 213, 170}, fillPattern = FillPattern.Solid, points = {{28, 4}, {32, 8}, {28, 2}, {34, 6}, {30, 2}, {34, 4}, {30, 0}, {26, 2}, {34, 0}, {26, 0}, {26, 2}, {28, 4}, {28, 4}, {26, 2}, {26, 2}, {26, 2}, {28, 8}, {28, 6}, {28, 4}}, smooth = Smooth.Bezier), Polygon(fillColor = {175, 175, 175}, fillPattern = FillPattern.Solid, points = {{-18, 0}, {28, 6}, {26, -2}, {-16, -16}, {-20, -16}, {-24, -6}, {-18, 0}}, smooth = Smooth.Bezier), Polygon(fillColor = {215, 215, 215}, fillPattern = FillPattern.Solid, points = {{72, -6}, {48, -6}, {36, -26}, {58, -86}, {72, -86}, {72, -6}}), Polygon(fillColor = {95, 95, 95}, fillPattern = FillPattern.Solid, points = {{49, -94}, {17, -40}, {7, -44}, {-1, -50}, {49, -94}}, smooth = Smooth.Bezier), Line(points = {{-7, 31}, {-3, 29}}), Line(points = {{-9, 18}, {-5, 18}}), Line(points = {{-7, 31}, {-3, 31}}), Text(lineColor = {238, 46, 47}, extent = {{-100, 90}, {100, 58}}, textString = "%CycleFileName")}),
      Diagram(coordinateSystem(extent = {{-100, -60}, {100, 60}}, preserveAspectRatio = false, initialScale = 0.1, grid = {2, 2})));
  end PropDriver;
  
  model DragForce "Vehicle rolling and aerodinamical drag force"
    import Modelica.Constants.g_n;
    extends Modelica.Mechanics.Translational.Interfaces.PartialElementaryOneFlangeAndSupport2;
    extends Modelica.Mechanics.Translational.Interfaces.PartialFriction;
    Modelica.Units.SI.Force f "Total drag force";
    Modelica.Units.SI.Velocity v "vehicle velocity";
    Modelica.Units.SI.Acceleration a "Absolute acceleration of flange";
    Real Sign;
    parameter Modelica.Units.SI.Mass m "vehicle mass";
    parameter Modelica.Units.SI.Density rho = 1.226 "air density";
    parameter Modelica.Units.SI.Area S "vehicle cross area";
    parameter Real fc(start = 0.01) "rolling friction coefficient";
    parameter Real Cx "aerodinamic drag coefficient";
  protected
    parameter Real A = fc*m*g_n;
    parameter Real B = 1/2*rho*S*Cx;
    constant Real f_pos[:, 2] = [0, 1];
  equation
  //  s = flange.s;
    v = der(s);
    a = der(v);
  // Le seguenti definizioni seguono l'ordine e le richieste del modello "PartialFriction" di
  // Modelica.Mechanics.Translational.Interfaces"
    v_relfric = v;
    a_relfric = a;
    f0 = A "force at 0 speed 0 but with slip";
    f0_max = A "max force at 0 speed without slip";
    free = false "in principle should become true whenthe wheel loose contact with road";
  // Now the computation of f, and its attribution to the flange:
    flange.f - f = 0;
  // friction force
    if v > 0 then
      Sign = 1;
    else
      Sign = -1;
    end if;
  //The following equation equates dragForce to the force applied when locked=true, otherwise term A.
    f - B*v^2*Sign = if locked then sa*unitForce else f0*(if startForward then Modelica.Math.Vectors.interpolate(f_pos[:, 1], f_pos[:, 2], v, 1) else if startBackward then -Modelica.Math.Vectors.interpolate(f_pos[:, 1], f_pos[:, 2], -v, 1) else if pre(mode) == Forward then Modelica.Math.Vectors.interpolate(f_pos[:, 1], f_pos[:, 2], v, 1) else -Modelica.Math.Vectors.interpolate(f_pos[:, 1], f_pos[:, 2], -v, 1));
  
    annotation(
      Documentation(info = "<html>
  <p>This component models the total (rolling and aerodynamic) vehicle drag resistance: </p>
  <p>F=fc*m*g+(1/2)*rho*Cx*S*v^2 </p>
  <p>It models reliably the stuck phase. Based on Modelica-Intrerfaces.PartialFriction model </p>
  </html>"),
      Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics = {Polygon(points = {{-98, 10}, {22, 10}, {22, 41}, {92, 0}, {22, -41}, {22, -10}, {-98, -10}, {-98, 10}}, lineColor = {0, 127, 0}, fillColor = {215, 215, 215}, fillPattern = FillPattern.Solid), Line(points = {{-42, -50}, {87, -50}}, color = {0, 0, 0}), Polygon(points = {{-72, -50}, {-41, -40}, {-41, -60}, {-72, -50}}, lineColor = {0, 0, 0}, fillColor = {128, 128, 128}, fillPattern = FillPattern.Solid), Line(points = {{-90, -90}, {-70, -88}, {-50, -82}, {-30, -72}, {-10, -58}, {10, -40}, {30, -18}, {50, 8}, {70, 38}, {90, 72}, {110, 110}}, color = {0, 0, 255}, thickness = 0.5), Text(extent = {{-82, 90}, {80, 50}}, lineColor = {0, 0, 255}, textString = "%name")}),
      Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics));
  end DragForce;
  annotation(
    uses(Modelica(version = "4.0.0"), EHPTlib(version = "2.1.4")));
end BookExamples;
