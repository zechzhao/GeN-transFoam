

	// Physical fields

	volScalarField defaultFlux
	(
		IOobject
		(
			"defaultFlux",
			runTime.timeName(),
			neutroMesh,
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE
		),
		neutroMesh
	);

	Info<< "flux\n" << endl;
	PtrList<volScalarField> flux(energyGroups);
	forAll(flux,fluxI)
	{
        Info<< "    Adding to flux " << fluxI<< "\n"<< endl;
        flux.set
        (
            fluxI,
            new volScalarField
            (
                IOobject
                (
                    "flux"+Foam::name(fluxI),
                    runTime.timeName(),
                    neutroMesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                defaultFlux
            )
        );
	}

	volScalarField oneGroupFlux
	(
		IOobject
		(
			"oneGroupFlux",
			runTime.timeName(),
			neutroMesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		defaultFlux	//defaultFlux as volMesh???
	);

	oneGroupFlux *= 0.0; // not sure about this!!!
	forAll(flux,energyI)
	{
		oneGroupFlux+=flux[energyI];
	Info<< "oneGroupFlux"<< energyI<<"="<<oneGroupFlux << endl;
	}
	oneGroupFlux.write();

	Info<< "Delayed neutron precursor\n" << endl;
	PtrList<volScalarField> prec(precGroups);
	forAll(prec,precI)
	{
	if (precI>1)
	{
	Info<< " "<<prec(precI-1)<< "\n"<< endl;
	}
        Info<< "    Adding to prec"  << precI<< endl;
        prec.set
        (
            precI,
            new volScalarField
            (
                IOobject
                (
                    "prec"+Foam::name(precI),
                    runTime.timeName(),
                    neutroMesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                neutroMesh,
				dimensionedScalar("", dimensionSet(0,-3,0,0,0,0,0), 0.0),
				zeroGradientFvPatchScalarField::typeName
            )
        );
	}

	Info<< "neutroSource\n" << endl;
	volScalarField neutroSource
	(
		IOobject
		(
			"neutroSource",
			runTime.timeName(),
			neutroMesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		neutroMesh,
		dimensionedScalar("", dimensionSet(0,-3,-1,0,0,0,0), 0.0), // need to look at volScalarField's constructor!
		zeroGradientFvPatchScalarField::typeName		//look this up (D)
	);

	Info<< "delayedNeutroSource\n" << endl;
	volScalarField delayedNeutroSource
	(
		IOobject
		(
			"delayedNeutroSource",
			runTime.timeName(),
			neutroMesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		neutroMesh,
		dimensionedScalar("", dimensionSet(0,-3,-1,0,0,0,0), 0.0),
		zeroGradientFvPatchScalarField::typeName
	);

	Info<< "scatteringSourceExtra\n" << endl;
	volScalarField scatteringSourceExtra
	(
		IOobject
		(
			"scatteringSourceExtra",
			runTime.timeName(),
			neutroMesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		neutroMesh,
		dimensionedScalar("", dimensionSet(0,-3,-1,0,0,0,0), 0.0),
		zeroGradientFvPatchScalarField::typeName
	);

	Info<< "scatteringSourceIntra\n" << endl;
	volScalarField scatteringSourceIntra
	(
		IOobject
		(
			"scatteringSourceIntra",
			runTime.timeName(),
			neutroMesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		neutroMesh,
		dimensionedScalar("", dimensionSet(0,-3,-1,0,0,0,0), 0.0),
		zeroGradientFvPatchScalarField::typeName
	);

	Info<< "volFuelPowerNeutro\n" << endl;
	volScalarField volFuelPowerNeutro
	(
		IOobject
		(
			"volFuelPowerNeutro",
			runTime.timeName(),
			neutroMesh,
			IOobject::READ_IF_PRESENT,
			IOobject::AUTO_WRITE
		),
		neutroMesh,
		dimensionedScalar("", dimensionSet(1,-1,-3,0,0,0,0), 0.0),
		zeroGradientFvPatchScalarField::typeName
	);

	Info<< "TfuelNeutro\n" << endl;
	volScalarField TfuelNeutro
	(
		IOobject
		(
			"TfuelNeutro",
			runTime.timeName(),
			neutroMesh,
			IOobject::READ_IF_PRESENT,
			IOobject::AUTO_WRITE
		),
		neutroMesh,
		dimensionedScalar("", dimensionSet(0,0,0,1,0,0,0),TfuelRef),
		zeroGradientFvPatchScalarField::typeName
	);


	Info<< "TcladNeutro\n" << endl;
	volScalarField TcladNeutro
	(
		IOobject
		(
			"TcladNeutro",
			runTime.timeName(),
			neutroMesh,
			IOobject::READ_IF_PRESENT,
			IOobject::AUTO_WRITE
		),
		neutroMesh,
		dimensionedScalar("", dimensionSet(0,0,0,1,0,0,0),TcladRef),
		zeroGradientFvPatchScalarField::typeName
	);



	Info<< "radExp\n" << endl;
	volScalarField radExp
	(
		IOobject
		(
			"radExp",
			runTime.timeName(),
			neutroMesh,
			IOobject::READ_IF_PRESENT,
			IOobject::AUTO_WRITE
		),
		neutroMesh,
		dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0), 0.0),
		zeroGradientFvPatchScalarField::typeName
	);

	Info<< "axExp\n" << endl;
	volScalarField axExp
	(
		IOobject
		(
			"axExp",
			runTime.timeName(),
			neutroMesh,
			IOobject::READ_IF_PRESENT,
			IOobject::AUTO_WRITE
		),
		neutroMesh,
		dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0), 0.0),
		zeroGradientFvPatchScalarField::typeName
	);

	Info<< "rhoCoolNeutro\n" << endl;
	volScalarField rhoCoolNeutro
	(
		IOobject
		(
			"rhoCoolNeutro",
			runTime.timeName(),
			neutroMesh,
			IOobject::READ_IF_PRESENT,
			IOobject::AUTO_WRITE
		),
		neutroMesh,
		dimensionedScalar("", dimensionSet(1,-3,0,0,0,0,0), rhoCoolRef),
		zeroGradientFvPatchScalarField::typeName
	);

	//displacement field for the mesh (should come from mechanical solver)
	Info<< "DNeutro\n" << endl;
	volVectorField DNeutro
	(
		IOobject
		(
			"DNeutro",
			runTime.timeName(),
			neutroMesh,
			IOobject::READ_IF_PRESENT,
			IOobject::AUTO_WRITE
		),
		neutroMesh,
		dimensionedVector("d_zero", dimensionSet(0,1,0,0,0,0,0), vector(0,0,0)),
		zeroGradientFvPatchScalarField::typeName
	);

    // Nuclear data fields: only one set, updated before solving each equation
	Info<< "Nuclear data fields\n" << endl;
	Info<< "IV\n" << endl;
	volScalarField IV
	(
		IOobject
		(
			"IV",
			runTime.timeName(),
			neutroMesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		neutroMesh,
		dimensionedScalar("", dimensionSet(0,-1,1,0,0,0,0), 1.0),
		zeroGradientFvPatchScalarField::typeName
	);

	Info<< "D\n" << endl;
	volScalarField D	
	(
		IOobject
		(
			"D",
			runTime.timeName(),
			neutroMesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		neutroMesh,
		dimensionedScalar("", dimensionSet(0,1,0,0,0,0,0), 0.0),
		zeroGradientFvPatchScalarField::typeName
	);

	Info<< "nuSigmaEff\n" << endl;
	PtrList<volScalarField> nuSigmaEff(energyGroups);
	forAll(nuSigmaEff, i)
	{
		nuSigmaEff.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					"nuSigmaEff"+i,
					runTime.timeName(),
					neutroMesh,
					IOobject::NO_READ,
					IOobject::NO_WRITE
				),
				neutroMesh,
				dimensionedScalar("", dimensionSet(0,-1,0,0,0,0,0), 0.0),
				zeroGradientFvPatchScalarField::typeName
			)
		);
	}
	Info<< "sigmaPow\n" << endl;
	volScalarField sigmaPow
	(
		IOobject
		(
			"sigmaPow",
			runTime.timeName(),
			neutroMesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		neutroMesh,
		dimensionedScalar("", dimensionSet(1,1,-2,0,0,0,0), 0.0),
		zeroGradientFvPatchScalarField::typeName
	);

	Info<< "sigmaDisapp\n" << endl;
	volScalarField sigmaDisapp
	(
		IOobject
		(
			"sigmaDisapp",
			runTime.timeName(),
			neutroMesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		neutroMesh,
		dimensionedScalar("", dimensionSet(0,-1,0,0,0,0,0), 0.0),
		zeroGradientFvPatchScalarField::typeName
	);

	Info<< "sigmaFrom\n" << endl;
	PtrList<volScalarField > sigmaFrom(energyGroups);
	forAll(sigmaFrom,energyI)
	{
		sigmaFrom.set
		(
			energyI,
			new volScalarField
			(
				IOobject
				(
					"sigmaFrom"+Foam::name(energyI),
					runTime.timeName(),
					neutroMesh,
					IOobject::NO_READ,
					IOobject::NO_WRITE
				),
				neutroMesh,
				dimensionedScalar("", dimensionSet(0,-1,0,0,0,0,0), 0.0),
				zeroGradientFvPatchScalarField::typeName
			)
		);
	}

	Info<< "chiPrompt\n" << endl;
	volScalarField chiPrompt
	(
		IOobject
		(
			"chiPrompt",
			runTime.timeName(),
			neutroMesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		neutroMesh,
		dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0), 0.0),
		zeroGradientFvPatchScalarField::typeName
	);

	Info<< "chiDelayed\n" << endl;
	volScalarField chiDelayed
	(
		IOobject
		(
			"chiDelayed",
			runTime.timeName(),
			neutroMesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		neutroMesh,
		dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0), 0.0),
		zeroGradientFvPatchScalarField::typeName
	);

	Info<< "Beta\n" << endl;
	volScalarField Beta
	(
		IOobject
		(
			"Beta",
			runTime.timeName(),
			neutroMesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		neutroMesh,
		dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0), 0.0),
		zeroGradientFvPatchScalarField::typeName
	);

	Info<< "BetaTot\n" << endl;
	volScalarField BetaTot
	(
		IOobject
		(
			"BetaTot",
			runTime.timeName(),
			neutroMesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		neutroMesh,
		dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0), 0.0),
		zeroGradientFvPatchScalarField::typeName
	);

	Info<< "lambda\n" << endl;
	volScalarField lambda
	(
		IOobject
		(
			"lambda",
			runTime.timeName(),
			neutroMesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		neutroMesh,
		dimensionedScalar("", dimensionSet(0,0,-1,0,0,0,0), 0.0),
		zeroGradientFvPatchScalarField::typeName
	);

	Info<< "fuelFraction\n" << endl;
	volScalarField fuelFraction
	(
		IOobject
		(
			"fuelFraction",
			runTime.timeName(),
			neutroMesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		neutroMesh,
		dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0), 1.0),
		zeroGradientFvPatchScalarField::typeName
	);


	Info<< "discFactor\n" << endl;
	PtrList<volScalarField > discFactor(energyGroups);
	forAll(discFactor,energyI)
	{
		discFactor.set
		(
			energyI,
			new volScalarField
			(
				IOobject
				(
					"discFactor"+Foam::name(energyI),
					runTime.timeName(),
					neutroMesh,
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
				neutroMesh,
				dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0), 0.0),
				zeroGradientFvPatchScalarField::typeName
			)
		);
	}


	//vectorField staticPoints = neutroMesh.points();

	//stuff I need to correctly project vol Power from neutroMesh to fluidMesh

	//scalarField staticCellsVolumes = neutroMesh.V();
	scalarField displacedCellsVolumes = neutroMesh.V();
	scalarField displacedCellsVolumesOld = displacedCellsVolumes;
