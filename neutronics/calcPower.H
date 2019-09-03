	pTot*=0;
    forAll(flux,energyI)
    {
		#include "setNeutroDiffConst.H"
		pTot += fvc::domainIntegrate(flux[energyI] * sigmaPow).value();

	}
	Info << "pTot neutro: " << pTot << endl;
{
	scalar pTot2 = fvc::domainIntegrate(volFuelPowerFluid[0]).value();;
	Info << "pTot fluid: " << pTot2 << endl;
}