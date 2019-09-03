Info<< "Normalizing scalar fluxes\n" << endl;

forAll(flux, energyI)
{
    flux[energyI] *= pTarget/pTot;
	// Info<< "flux[energy I] " << flux[energyI]<<  "\n" << endl;
}

Info<< "Done\n" << endl;
