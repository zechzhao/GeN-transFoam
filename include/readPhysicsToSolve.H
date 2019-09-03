const bool solveFluidMechanics =
    runTime.controlDict().lookupOrDefault("solveFluidMechanics", false);

const bool solveEnergy =
    runTime.controlDict().lookupOrDefault("solveEnergy", false);

const bool solveNeutronics =
    runTime.controlDict().lookupOrDefault("solveNeutronics", false);

const bool solveThermalMechanics =
    runTime.controlDict().lookupOrDefault("solveThermalMechanics", false);

const bool eigenvalueNeutronics =
    runTime.controlDict().lookupOrDefault("eigenvalueNeutronics", false);

const bool tighltyCoupled =
    runTime.controlDict().lookupOrDefault("tighltyCoupled", true);

const bool implicitFuel =
    runTime.controlDict().lookupOrDefault("implicitFuel", true);

const bool powerControl =
    runTime.controlDict().lookupOrDefault("powerControl", false);

const scalar powerControlPorportional =
    runTime.controlDict().lookupOrDefault("powerControlPorportional", 1e-3);

const scalar powerControlIntegral =
    runTime.controlDict().lookupOrDefault("powerControlIntegral", 1e-3);
