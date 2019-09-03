    // Solve the Momentum equation
	Info<< "Solve the Momentum equation " << nl << endl;

    tmp<fvVectorMatrix> UEqn
    (
        fvm::ddt(rho, U)
      + (1/porosity)*fvm::div(phi, U)
	  + turb.divDevRhoReff(U)
      + fvOptions(rho, U) //attento! no corrette per porosita'
	  - porousMedium.semiImplicitMomentumSource(U)
    );
    UEqn().relax();

    fvOptions.constrain(UEqn());//attento! no corrette per porosita'

	//"reconstruct" stuff is used to fake a staggered grid (rhie chow)
    if (momentumPredictor)
    {
		if(solveFluidMechanics)
		{
			solve
			(
				UEqn()
			  ==
				porosity*fvc::reconstruct
				(
					(
					  - ghf*fvc::snGrad(rho*rhok)//p_rgh=p-rho*rhok*g*z
					  - fvc::snGrad(p_rgh)
					)*mesh.magSf()
				),
				mesh.solver(U.select(finalIter))

			);

			fvOptions.correct(U);//attento! no corrette per porosita'
			realU = U / porosity;
			realU.correctBoundaryConditions();
			K = 0.5*magSqr(realU);
		}
    }





