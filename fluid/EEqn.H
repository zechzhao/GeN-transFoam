{
	thermo.correct();
    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho*porosity, he) + fvm::div(phi, he)
      + fvc::ddt(rho*porosity, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            )
          : -dpdt
        )
      - fvm::laplacian(turb.alphaEff(), he)
     ==
        rad.Sh(thermo)
        + porousMedium.externalVolHeatSource()
        - pow((pow(porousMedium.externalHeatTransferCoefficient(),-1)+pow(porousMedium.heatTransferCoefficient(),-1)),-1)* porousMedium.volumetricArea()*(thermo.T()-porousMedium.externalT())
     	+ fvOptions(rho, he) //attento! no corrette per porosita'
     	+ subscaleFuel.heatSources(thermo.T(),porousMedium.heatTransferCoefficient(),volFuelPower,implicitFuel) * porousMedium.volumetricArea()
     	- porousMedium.heatTransferCoefficient()* porousMedium.volumetricArea()*(thermo.T()-Tstructures)
     	*(1.0-TavFuel/(TavFuel+dimensionedScalar("",dimensionSet(0,0,0,1,0),SMALL)))//stupid way to kill the last term where the fuel is
    );

    EEqn.relax();
	
//    Info<< "heatTransferCoe"<<porousMedium.heatTransferCoefficient()<<endl;

    fvOptions.constrain(EEqn); //attento! no corrette per porosita'

    scalar EEqnResidual(EEqn.solve(mesh.solver(he.select(finalIter))).initialResidual());
    residual = max(residual,EEqnResidual);

    fvOptions.correct(he); //attento! no corrette per porosita'

    thermo.correct();
    rad.correct();

    if(!compr)
    {
    	rhok = 1.0 - beta*(thermo.T() - Tref);
	}

    Info<< "Min/max T:" << min(thermo.T()).value() << ' '
    << max(thermo.T()).value() << endl;

     TavFuel = subscaleFuel.fuelAverage();
     TavClad = subscaleFuel.cladAverage();
     TinnerFuel = subscaleFuel.fuelIn();
     TouterFuel = subscaleFuel.fuelOut();
     TinnerClad = subscaleFuel.cladIn();
     TouterClad = subscaleFuel.cladOut();

    neutroToFluid.mapTgtToSrc( TavFuel, plusEqOp<scalar>(), TfuelNeutro.internalField());
	neutroToFluid.mapTgtToSrc( TavClad, plusEqOp<scalar>(), TcladNeutro.internalField());
	neutroToFluid.mapTgtToSrc( thermo.rho(), plusEqOp<scalar>(), rhoCoolNeutro.internalField());

	//solve for Temp of structures
    fvScalarMatrix structureEqn
    (
        fvm::ddt(porousMedium.externalRhoCp(), Tstructures)
     ==
        porousMedium.heatTransferCoefficient()* porousMedium.volumetricArea()*(thermo.T()-Tstructures)
        *(1.0-TavFuel/(TavFuel+dimensionedScalar("",dimensionSet(0,0,0,1,0),SMALL)))//stupid way to kill the last term where the fuel is
    );
	structureEqn.solve();

}


