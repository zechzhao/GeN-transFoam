{
	Info<< "Solve the Pressure equation " << nl << endl;
    bool closedVolume = p_rgh.needReference();
	dimensionedScalar compressibility = fvc::domainIntegrate(psi);
	bool compressible = (compressibility.value() > SMALL);

	if(compr)
	{
	    rho = thermo.rho();
	}

    volScalarField rAU(1.0/UEqn().A());
    surfaceScalarField rhorAUf("Dp", fvc::interpolate(porosity*rho*rAU)); //qui no rhok (corr solo dove c'e' g)

    volVectorField HbyA("HbyA", U);
    HbyA = rAU*(UEqn().H());

    surfaceScalarField phig(-rhorAUf*ghf*fvc::snGrad(rho*rhok)*mesh.magSf());

	phi = phi / fvc::interpolate(porosity); //necessary because the porosity inside the divergence is already included in U
											//nn sono sicurissimomissimo sia giusto correggere phi con porosity. Se ci sono strani problemi controlla qui!

    surfaceScalarField phiHbyA
	    (
	        "phiHbyA",
	        (
	            (fvc::interpolate(rho*HbyA) & mesh.Sf())
	          + rhorAUf*fvc::ddtCorr(rho, U, phi)
	        )
	      + phig
    );

    phi = phi * fvc::interpolate(porosity);//bring phi back to normal

    //fvOptions.relativeFlux(fvc::interpolate(rho*porosity), phiHbyA);//Make the given absolute mass-flux relative: used only when MRF option is used (multi reference frame)
	fvOptions.makeRelative(fvc::interpolate(rho*porosity), phiHbyA);

    {
       	fvScalarMatrix p_rghDDtEqn
       	(
       	    (compr)*(fvc::ddt(porosity,rho) + psi*correction(fvm::ddt(porosity,p_rgh))) //correction returns the correction form of the given matrix by subtracting the matrix multiplied by the current field
       	  + fvc::div(phiHbyA)
       	);

		if(compr)
		{
        	// Thermodynamic density needs to be updated by psi*d(p) after the
        	// pressure solution - done in 2 parts. Part 1:
        	thermo.rho() -= psi*p_rgh;
		}
        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            fvScalarMatrix p_rghEqn
            (
                p_rghDDtEqn
              - fvm::laplacian(rhorAUf, p_rgh)
            );

			if(!compr)
			{
				p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));
			}

            p_rghEqn.solve
            (
                mesh.solver
                (
                    p_rgh.select
                    (
                        (
                           oCorr == nOuterCorr-1
                        && corr == nCorr-1
                        && nonOrth == nNonOrthCorr
                        )
                    )
                )
            );

            if (nonOrth == nNonOrthCorr) //velocity corrector step
            {
                phi = (phiHbyA + p_rghEqn.flux()); //porosity*rho*H/A (+rho*g*z/A*grad(gamma*rhp)) + grad(P)/A **** the term into the paranthesis is not in Jazak thesis because there is no gravity source there!
                U = HbyA
                  + rAU*porosity*fvc::reconstruct((phig + p_rghEqn.flux())/rhorAUf);
                U.correctBoundaryConditions();
                fvOptions.correct(U);
                realU = U / porosity;
                realU.correctBoundaryConditions();
                K = 0.5*magSqr(realU);
            }
        }
		if(compr)
		{
    	   	// Second part of thermodynamic density update
		  	thermo.rho() += psi*p_rgh;
		}
    }
    p = p_rgh + rho*rhok*gh;

	if(compr)
	{
		 // Update pressure time derivative if needed
		if (thermo.dpdt())
		{
			dpdt = fvc::ddt(p);
		}
    	// Solve continuity
   		#include "rhoEqn.H"
   		// Update continuity errors
    	#include "compressibleContinuityErrors.H"
    	// For closed-volume cases adjust the pressure and density levels
    	// to obey overall mass continuity
    	if (closedVolume && compressible)
		{
		    p += (initialMass - fvc::domainIntegrate(thermo.rho()))
		            /compressibility;
		    rho = thermo.rho();

    	}
	}
	else
	{
		#include "continuityErrs.H"
		if (closedVolume)
		{
			p += dimensionedScalar
			(
				"p",
				p.dimensions(),
				pRefValue - getRefCellValue(p, pRefCell)
			);
		}
	}

	p_rgh = (p - rho*rhok*gh);

}
