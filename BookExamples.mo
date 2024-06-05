package BookExamples "Package containing basic EV models"
  model FirstEV "Simulates a very basic Electric Vehicle"
    import Modelica;
    extends Modelica.Icons.Example;
    Modelica.Units.SI.Energy enP1(start = 0, fixed = true);
    Modelica.Units.SI.Energy enP2(start = 0, fixed = true);
    Modelica.Units.SI.Energy enP1Pos(start = 0, fixed = true);
    Modelica.Units.SI.Energy enP1Neg(start = 0, fixed = true);
    EHPTlib.SupportModels.Miscellaneous.DragForce dragF(Cx = 0.26, S = 2.2, fc = 0.014, m = mass.m, rho(displayUnit = "kg/m3") = 1.225, v(start = 0, fixed = true)) annotation(
      Placement(transformation(origin = {80, -26}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    Modelica.Mechanics.Rotational.Components.IdealGear gear(ratio = 10, flange_b(phi(start = 0, fixed = true))) annotation(
      Placement(transformation(origin = {-22, 20}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Mechanics.Rotational.Sources.Torque torque annotation(
      Placement(visible = true, transformation(extent = {{-82, 10}, {-62, 30}}, rotation = 0)));
    EHPTlib.SupportModels.Miscellaneous.PropDriver driver(CycleFileName = Modelica.Utilities.Files.loadResource("modelica://BookExamples/Resources/Sort1.txt"), extrapolation = Modelica.Blocks.Types.Extrapolation.LastTwoPoints, k = 1000, yMax = 100000.0) annotation(
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
  annotation(
    uses(Modelica(version = "4.0.0"), EHPTlib(version = "2.1.4")));
end BookExamples;
