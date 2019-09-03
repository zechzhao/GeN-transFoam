const word thermoMechanicalRegionName("thermoMechanicalRegion");

Foam::fvMesh mechMesh
(
    Foam::IOobject
    (
        thermoMechanicalRegionName,
        runTime.timeName(),
        runTime,
        Foam::IOobject::MUST_READ
    )
);
