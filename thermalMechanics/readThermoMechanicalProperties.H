    Info<< "Reading thermo-mechanical properties\n" << endl;

	// Read from dictionary
    IOdictionary thermoMechanicalProperties
    (
        IOobject
        (
            "thermoMechanicalProperties",
            runTime.constant(),
            mechMesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

	Switch planeStress(thermoMechanicalProperties.lookup("planeStress"));
	Switch linkedFuel(thermoMechanicalProperties.lookup("linkedFuel"));

	vector fuelOrientation(thermoMechanicalProperties.lookup("fuelOrientation"));
	//Switch thermalStress(thermoMechanicalProperties.lookup("thermalStress"));

	dimensionedScalar TrefStructures("",dimensionSet(0,0,0,1,0,0,0),thermoMechanicalProperties.lookup("TrefStructure"));

	PtrList<entry> TMEntries(thermoMechanicalProperties.lookup("zones"));

	label TMZoneNumber = TMEntries.size();

	PtrList<scalar > rhoMechList(TMZoneNumber);
	PtrList<scalar > rhoEList(TMZoneNumber);
	PtrList<scalar > nuMechList(TMZoneNumber);

	PtrList<scalar > CList(TMZoneNumber);
	PtrList<scalar > rhoKList(TMZoneNumber);
	PtrList<scalar > alphaList(TMZoneNumber);

	PtrList<scalar > TrefFuelList(TMZoneNumber);
	PtrList<scalar > alphaFuelList(TMZoneNumber);

	PtrList<scalar > TrefCRList(TMZoneNumber);
	PtrList<scalar > alphaCRList(TMZoneNumber);

	forAll(TMEntries,zoneI)
	{
		dictionary& dict = TMEntries[zoneI].dict();

		rhoMechList.set(zoneI,new scalar(dict.lookupOrDefault("rho",1.0)));
		rhoEList.set(zoneI,new scalar(dict.lookupOrDefault("E",1.0)));
		nuMechList.set(zoneI,new scalar(dict.lookupOrDefault("nu",0.0)));

		CList.set(zoneI,new scalar(dict.lookupOrDefault("C",1.0)));
		rhoKList.set(zoneI,new scalar(dict.lookupOrDefault("k",1.0)));
		alphaList.set(zoneI,new scalar(dict.lookupOrDefault("alpha",0.0)));

		TrefFuelList.set(zoneI,new scalar(dict.lookupOrDefault("TrefFuel",0.0)));
		alphaFuelList.set(zoneI,new scalar(dict.lookupOrDefault("alphaFuel",0.0)));

		TrefCRList.set(zoneI,new scalar(dict.lookupOrDefault("TrefCR",0.0)));
		alphaCRList.set(zoneI,new scalar(dict.lookupOrDefault("alphaCR",0.0)));
	}

	// Create volFields
	volScalarField rhoMech
	(
		IOobject
		(
			"rhoMech",
			runTime.timeName(),
			mechMesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mechMesh,
		dimensionedScalar("", dimensionSet(1,-3,0,0,0,0,0), 1.0),
		zeroGradientFvPatchScalarField::typeName
	);

	volScalarField rhoE
	(
		IOobject
		(
			"E",
			runTime.timeName(),
			mechMesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mechMesh,
		dimensionedScalar("", dimensionSet(1,-1,-2,0,0,0,0), 1.0),
		zeroGradientFvPatchScalarField::typeName
	);

	volScalarField nuMech
	(
		IOobject
		(
			"nuMech",
			runTime.timeName(),
			mechMesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mechMesh,
		dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0), 0.0),
		zeroGradientFvPatchScalarField::typeName
	);

	volScalarField threeKalpha
	(
		IOobject
		(
			"threeKalpha",
			runTime.timeName(),
			mechMesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mechMesh,
		dimensionedScalar("0", dimensionSet(0, 2, -2 , -1, 0), 0.0),
		zeroGradientFvPatchScalarField::typeName
	);


	volScalarField DT
	(
		IOobject
		(
			"DT",
			runTime.timeName(),
			mechMesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mechMesh,
		dimensionedScalar("0", dimensionSet(0, 2, -1 , 0, 0), 0.0),
		zeroGradientFvPatchScalarField::typeName
	);

	volScalarField C
	(
		IOobject
		(
			"C",
			runTime.timeName(),
			mechMesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mechMesh,
		dimensionedScalar("0", dimensionSet(0, 2, -2 , -1, 0), 1.0),
		zeroGradientFvPatchScalarField::typeName
	);

	volScalarField rhoK
	(
		IOobject
		(
			"rhoK",
			runTime.timeName(),
			mechMesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mechMesh,
		dimensionedScalar("0", dimensionSet(1, 1, -3 , -1, 0), 0.0),
		zeroGradientFvPatchScalarField::typeName
	);

	volScalarField alpha
	(
		IOobject
		(
			"alpha",
			runTime.timeName(),
			mechMesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mechMesh,
		dimensionedScalar("0", dimensionSet(0, 0, 0 , -1, 0), 1.0e-5),
		zeroGradientFvPatchScalarField::typeName
	);

	volScalarField TrefFuel
	(
		IOobject
		(
			"TrefFuel",
			runTime.timeName(),
			mechMesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mechMesh,
		dimensionedScalar("0", dimensionSet(0, 0, 0 , 1, 0), 0.0),
		zeroGradientFvPatchScalarField::typeName
	);

	volScalarField alphaFuel
	(
		IOobject
		(
			"alphaFuel",
			runTime.timeName(),
			mechMesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mechMesh,
		dimensionedScalar("0", dimensionSet(0, 0, 0 , -1, 0), 0.0),
		zeroGradientFvPatchScalarField::typeName
	);

	volScalarField TrefCR
	(
		IOobject
		(
			"TrefCR",
			runTime.timeName(),
			mechMesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mechMesh,
		dimensionedScalar("0", dimensionSet(0, 0, 0 , 1, 0), 0.0),
		zeroGradientFvPatchScalarField::typeName
	);

	volScalarField alphaCR
	(
		IOobject
		(
			"alphaCR",
			runTime.timeName(),
			mechMesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mechMesh,
		dimensionedScalar("0", dimensionSet(0, 0, 0 , -1, 0), 0.0),
		zeroGradientFvPatchScalarField::typeName
	);

	// Set volFields based on dictionary
	forAll(TMEntries, zoneI)
	{

		const word& name = TMEntries[zoneI].keyword();

		label zoneId = mechMesh.cellZones().findZoneID(name);

		forAll(mechMesh.cellZones()[zoneId], cellIlocal)
		{
			label cellIglobal = mechMesh.cellZones()[zoneId][cellIlocal];

			rhoMech[cellIglobal] = rhoMechList[zoneI];
			rhoE[cellIglobal] = rhoEList[zoneI];
			nuMech[cellIglobal] = nuMechList[zoneI];

			C[cellIglobal] = CList[zoneI];
			rhoK[cellIglobal] = rhoKList[zoneI];
			alpha[cellIglobal] = alphaList[zoneI];

			TrefFuel[cellIglobal] = TrefFuelList[zoneI];
			alphaFuel[cellIglobal] = alphaFuelList[zoneI];

			TrefCR[cellIglobal] = TrefCRList[zoneI];
			alphaCR[cellIglobal] = alphaCRList[zoneI];
		}

	}

	rhoMech.correctBoundaryConditions();
	rhoE.correctBoundaryConditions();
	nuMech.correctBoundaryConditions();

	C.correctBoundaryConditions();
	rhoK.correctBoundaryConditions();
	alpha.correctBoundaryConditions();

	TrefFuel.correctBoundaryConditions();
	alphaFuel.correctBoundaryConditions();

	TrefCR.correctBoundaryConditions();
	alphaCR.correctBoundaryConditions();

	Info<< "Normalising E : E/rho\n" << endl;
	volScalarField E(rhoE/rhoMech);

	Info<< "Calculating Lame's coefficients\n" << endl;

	volScalarField muMech(E/(2.0*(1.0 + nuMech)));
	volScalarField lambdaMech(nuMech*E/((1.0 + nuMech)*(1.0 - 2.0*nuMech)));
	volScalarField threeK(E/(1.0 - 2.0*nuMech));

	if (planeStress)
	{
		Info<< "Plane Stress\n" << endl;

		lambdaMech = nuMech*E/((1.0 + nuMech)*(1.0 - nuMech));
		threeK = E/(1.0 - nuMech);
	}
	else
	{
		Info<< "Plane Strain\n" << endl;
	}

    Info<< "Normalising k : k/rho\n" << endl;
    volScalarField k(rhoK/rhoMech);

    Info<< "Calculating thermal coefficients\n" << endl;

    threeKalpha = threeK*alpha;
    DT = k/C;

	Info<< "done\n" << endl;