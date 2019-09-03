    if (runTime.outputTime())
    {
        volSymmTensorField sigma
        (
            IOobject
            (
                "sigma",
                runTime.timeName(),
                mechMesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            rhoMech*sigmaD
        );


		const volScalarField& Tmech = Tptr();
		sigma = sigma - I*(rhoMech*threeKalpha*Tmech);

        volScalarField sigmaEq
        (
            IOobject
            (
                "sigmaEq",
                runTime.timeName(),
                mechMesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sqrt((3.0/2.0)*magSqr(dev(sigma)))
        );

        Info<< "Max sigmaEq = " << max(sigmaEq).value()
            << endl;

    }
