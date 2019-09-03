
vectorField startPoints = neutroMesh.points(); //used to bring back mesh after interpolations

if(runTime.timeOutputValue()!=0)
{
	Info << "Bring back neutroMesh to underformed case to create correct interpolators" << nl << endl;

	const volPointInterpolation& neutroMeshPointInterpolationInitial = volPointInterpolation::New(neutroMesh);

	pointVectorField neutroPointsDisplacementInitial = neutroMeshPointInterpolationInitial.interpolate(DNeutro);

	vectorField displacedPointsInitial = startPoints - neutroPointsDisplacementInitial.internalField() ;

	neutroMesh.movePoints(displacedPointsInitial);
}


	scalarField staticCellsVolumes = neutroMesh.V();

	meshToMesh neutroToTM
	(
		neutroMesh,
		mechMesh,
		Foam::meshToMesh::interpolationMethod(2),
		false
	);

	PtrList<meshToMesh> fluidToTMList(fluidRegions.size());
	forAll(fluidRegions, i)
	{
		fluidToTMList.set
		(
			i,
			new meshToMesh
			(
				fluidRegions[i],
				mechMesh,
				Foam::meshToMesh::interpolationMethod(2),
				false
			)
		);
	}

	PtrList<meshToMesh> neutroToFluidList(fluidRegions.size());
	forAll(fluidRegions, i)
	{
		neutroToFluidList.set
		(
			i,
			new meshToMesh
			(
				neutroMesh,
				fluidRegions[i],
				Foam::meshToMesh::interpolationMethod(2), // imCellVolumeWeightc
				false
			)
		);
	}

	PtrList<meshToMesh> fuelToFluidList(fluidRegions.size()); // map fuel values to fluid
	forAll(fluidRegions, i)
	{
		fuelToFluidList.set
		(
			i,
			new meshToMesh
			(
				fuelMesh,
				fluidRegions[i],
				Foam::meshToMesh::interpolationMethod(2), // imCellVolumeWeightc better to use 0 (imDirect)
				false
			)
		);
	}


if(runTime.timeOutputValue()!=0)
{
	Info << "Re-derform mesh to start calculation" << nl << endl;

	neutroMesh.movePoints(startPoints);
}
