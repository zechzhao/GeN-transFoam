
{

Info << "Displace mesh" << nl << endl;

const volPointInterpolation& neutroMeshPointInterpolationOld = volPointInterpolation::New(neutroMesh);

pointVectorField neutroPointsDisplacementOld = neutroMeshPointInterpolationOld.interpolate(DNeutro);

DNeutro*=0;

volVectorField neutroDisp = (Disp-((Disp & fuelOrientation) * fuelOrientation)) + fuelDispVector ;//remove any axial variation

neutroToTM.mapTgtToSrc( neutroDisp , plusEqOp<vector>(), DNeutro.internalField());

DNeutro.correctBoundaryConditions();

const volPointInterpolation& neutroMeshPointInterpolation = volPointInterpolation::New(neutroMesh);

pointVectorField neutroPointsDisplacement = neutroMeshPointInterpolation.interpolate(DNeutro);

vectorField displacedPoints = neutroMesh.points() + neutroPointsDisplacement.internalField() - neutroPointsDisplacementOld.internalField() ;

displacedCellsVolumesOld = neutroMesh.V();

neutroMesh.movePoints(displacedPoints);

displacedCellsVolumes = neutroMesh.V();

forAll(prec,precI)
{
	prec[precI].internalField() *= (displacedCellsVolumesOld/displacedCellsVolumes);
	prec[precI].correctBoundaryConditions();
}

forAll(flux,energyI)
{
	flux[energyI].internalField() *= (displacedCellsVolumesOld/displacedCellsVolumes);
	flux[energyI].correctBoundaryConditions();
}

}

