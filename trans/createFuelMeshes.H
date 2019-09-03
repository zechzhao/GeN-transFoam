	const word fuelRegionName("fuelRegion");

	Foam::fvMesh fuelMesh
	(
		Foam::IOobject
		(
			fuelRegionName,
			runTime.timeName(),
			runTime,
			Foam::IOobject::MUST_READ
		)
	);

