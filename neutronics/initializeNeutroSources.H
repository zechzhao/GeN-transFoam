
	delayedNeutroSource *= 0.0;
    forAll(prec,precI)
    {
		#include "setPrecConst.H"
		delayedNeutroSource += lambda*prec[precI];
	}

	neutroSource *= 0.0;
    forAll(flux,energyI)
    {
		#include "setNeutroDiffConst.H"
		neutroSource += nuSigmaEff[energyI]*flux[energyI];
	}
