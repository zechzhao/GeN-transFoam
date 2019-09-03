

// Physical fields
Info<< "Generating fields for fuel calcualtion"<<endl;

//input fields: mass flow rate, coolant temperature, linear generation rate
//out fields: fuel temperature(average, centraline, outer)
//			  cladding temperature (inner, average, outer)
//            Heat source to coolant
//            fuel and cladding size maybe?


	Info<< "coolantFlowFuel\n" << endl;
	volScalarField coolantFlowFuel
	(
		IOobject
		(
			"coolantFlowFuel",
			runTime.timeName(),
			fuelMesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		fuelMesh,
		dimensionedScalar("", dimensionSet(0,0,0,1,0,0,0),1.28),
		zeroGradientFvPatchScalarField::typeName
	);

	Info<< "TcoolantFuel\n" << endl;
	volScalarField TcoolantFuel
	(
		IOobject
		(
			"TcoolantFuel",
			runTime.timeName(),
			fuelMesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		fuelMesh,
		dimensionedScalar("", dimensionSet(0,0,0,1,0,0,0),0.0),
		zeroGradientFvPatchScalarField::typeName
	);

	Info<< "linearRateFuel\n" << endl;
	volScalarField linearRateFuel
	(
		IOobject
		(
			"linearRateFuel",
			runTime.timeName(),
			fuelMesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		fuelMesh,
		dimensionedScalar("", dimensionSet(0,0,0,1,0,0,0),0.0),
		zeroGradientFvPatchScalarField::typeName
	);

	// Info<< "heatSourcesFuel" << endl;
	// volScalarField heatSourcesFuel
	// (
	// 	IOobject
	// 	(
	// 		"heatSourcesFuel",
	// 		runTime.timeName(),
	// 		fuelMesh,
	// 		IOobject::NO_READ,
	// 		IOobject::AUTO_WRITE
	// 	),
	// 	fuelMesh
	// );


	Info<< "TavfuelFuel\n" << endl;
	volScalarField TavfuelFuel
	(
		IOobject
		(
			"TavfuelFuel",
			runTime.timeName(),
			fuelMesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		fuelMesh,
		dimensionedScalar("", dimensionSet(0,0,0,1,0,0,0),0.0),
		zeroGradientFvPatchScalarField::typeName
	);

	Info<< "TavcladFuel\n" << endl;
	volScalarField TavcladFuel
	(
		IOobject
		(
			"TavcladFuel",
			runTime.timeName(),
			fuelMesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		fuelMesh,
		dimensionedScalar("", dimensionSet(0,0,0,1,0,0,0),0.0),
		zeroGradientFvPatchScalarField::typeName
	);

	Info<< "TCentralFuel\n" << endl;
	volScalarField TCentralFuel
	(
		IOobject
		(
			"TCentralFuel",
			runTime.timeName(),
			fuelMesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		fuelMesh,
		dimensionedScalar("", dimensionSet(0,0,0,1,0,0,0),0.0),
		zeroGradientFvPatchScalarField::typeName
	);

	Info<< "TOuterFuel\n" << endl;
	volScalarField TOuterFuel
	(
		IOobject
		(
			"TOuterFuel",
			runTime.timeName(),
			fuelMesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		fuelMesh,
		dimensionedScalar("", dimensionSet(0,0,0,1,0,0,0),0.0),
		zeroGradientFvPatchScalarField::typeName
	);

	Info<< "TInnerClad\n" << endl;
	volScalarField TInnerClad
	(
		IOobject
		(
			"TInnerClad",
			runTime.timeName(),
			fuelMesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		fuelMesh,
		dimensionedScalar("", dimensionSet(0,0,0,1,0,0,0),0.0),
		zeroGradientFvPatchScalarField::typeName
	);

	Info<< "TOuterClad\n" << endl;
	volScalarField TOuterClad
	(
		IOobject
		(
			"TOuterClad",
			runTime.timeName(),
			fuelMesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		fuelMesh,
		dimensionedScalar("", dimensionSet(0,1,0,0,0,0,0),0.0),
		zeroGradientFvPatchScalarField::typeName
	);

	// Info<< "radiuscladFuel\n" << endl;
	// volScalarField radiuscladFuel
	// (
	// 	IOobject
	// 	(
	// 		"radiuscladFuel",
	// 		runTime.timeName(),
	// 		fuelMesh,
	// 		IOobject::NO_READ,
	// 		IOobject::AUTO_WRITE
	// 	),
	// 	fuelMesh,
		// dimensionedScalar("", dimensionSet(0,1,0,0,0,0,0),0.0),
		// zeroGradientFvPatchScalarField::typeName
	// );