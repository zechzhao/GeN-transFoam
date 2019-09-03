    fvMesh& mesh = fluidRegions[i];

    meshToMesh& fluidToTM = fluidToTMList[i];
    meshToMesh& neutroToFluid = neutroToFluidList[i];

    meshToMesh& fuelToFluid = fuelToFluidList[i];

    rhoThermo& thermo = thermoFluid[i];
    thermo.validate(args.executable(), "h", "e");

    volScalarField& rho = rhoFluid[i];
    volScalarField& rhok = rhokFluid[i];//for Boussinesq
    volVectorField& U = UFluid[i];
    volVectorField& realU = realUFluid[i];
    surfaceScalarField& phi = phiFluid[i];

    volScalarField& Re = ReFluid[i];

    compressible::turbulenceModel& turb = turbulence[i];
    Foam::porousMedium& porousMedium = porousMediumFluid[i];
    Foam::subscaleFuel& subscaleFuel = subscaleFuelFluid[i];

    const volScalarField& porosity = porousMedium.voidFraction();
    volScalarField& K = KFluid[i];
    volScalarField& dpdt = dpdtFluid[i];

	realU = U / porosity;
	realU.correctBoundaryConditions();

    volScalarField& p = thermo.p();
    const volScalarField& psi = psiFluid[i];

    volScalarField& p_rgh = p_rghFluid[i];
    const volScalarField& gh = ghFluid[i];
    const surfaceScalarField& ghf = ghfFluid[i];

    const uniformDimensionedScalarField& beta = betaFluid[i];
    const uniformDimensionedScalarField& Tref = TrefFluid[i];
    const uniformDimensionedScalarField& rhoRef = rhoRefFluid[i];

    if(compr)
    {
    	rhok = 1.0;
	}
	else
	{
		rhok = 1.0 - beta*(thermo.T() - Tref);
	}

    volScalarField& volFuelPower = volFuelPowerFluid[i];
    volScalarField& TavFuel = TavFuelFluid[i];
    volScalarField& TavClad = TavCladFluid[i];

    TavFuel = subscaleFuel.fuelAverage();
    TavClad = subscaleFuel.cladAverage();

    volScalarField& TinnerFuel = TinnerFuelFluid[i];
    volScalarField& TouterFuel = TouterFuelFluid[i];
    volScalarField& TinnerClad = TinnerCladFluid[i];
    volScalarField& TouterClad = TouterCladFluid[i];

    volScalarField& Tstructures = TstructuresFluid[i];

    radiation::radiationModel& rad = radiation[i];

    fv::IOoptionList& fvOptions = fluidFvOptions[i];

    const dimensionedScalar initialMass
    (
        "initialMass",
        dimMass,
        initialMassFluid[i]
    );

	//For incompressible
	label pRefCell = 0;
	scalar pRefValue = 0.0;
	setRefCell
	(
		p,
		p_rgh,
		mesh.solutionDict().subDict("PIMPLE"),
		pRefCell,
		pRefValue
	);
	if(!compr)
	{
		if (p_rgh.needReference())
		{
			p += dimensionedScalar
			(
				"p",
				p.dimensions(),
				pRefValue - getRefCellValue(p, pRefCell)
			);
		}
	}