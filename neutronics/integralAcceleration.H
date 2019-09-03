



{


Info << "Accelerating neutronic solution" << endl ;

dimensionedScalar integralFlux = fvc::domainIntegrate(oneGroupFlux);

volScalarField oneGroupFluxOld("",oneGroupFlux);
oneGroupFluxOld *= 0.0;
forAll(flux, energyI)
{
       oneGroupFluxOld +=flux[energyI].oldTime();
}

dimensionedScalar integralFluxOld = fvc::domainIntegrate(oneGroupFluxOld);


dimensionedScalar temporalTermInt("",(dimless/dimTime),scalar(0.0));
dimensionedScalar temporalTerm("",(dimless/dimLength),scalar(0.0));

dimensionedScalar leakageTermInt("",(dimless/dimTime),scalar(0.0));
dimensionedScalar leakageTerm("",(dimless/dimLength),scalar(0.0));

dimensionedScalar disappearenceTermInt("",(dimless/dimTime),scalar(0.0));
dimensionedScalar disappearenceTerm("",(dimless/dimLength),scalar(0.0));

dimensionedScalar sourceTermInt("",(dimless/dimTime),scalar(0.0));
dimensionedScalar sourceTerm("",(dimless/dimLength),scalar(0.0));

dimensionedScalar scatteringSourceTermInt("",(dimless/dimTime),scalar(0.0));
dimensionedScalar scatteringSourceTerm("",(dimless/dimLength),scalar(0.0));

dimensionedScalar precursorTermInt = fvc::domainIntegrate(delayedNeutroSource);
dimensionedScalar precursorTerm = precursorTermInt;

for(label energyI=0; energyI<energyGroups ; energyI++)
{
	#include "setNeutroDiffConst.H"

	temporalTermInt += fvc::domainIntegrate(IV*flux[energyI]/runTime.deltaT());

	leakageTermInt += fvc::domainIntegrate(fvc::laplacian(D, flux[energyI]));

	disappearenceTermInt += fvc::domainIntegrate(sigmaDisapp*flux[energyI]);

	sourceTermInt += fvc::domainIntegrate(nuSigmaEff[energyI]*flux[energyI]);

	for(label energyJ=0; energyJ<energyGroups ; energyJ++)
	{
		if(energyJ != energyI)
		{
			scatteringSourceTermInt += fvc::domainIntegrate(sigmaFrom[energyJ] * flux[energyJ]);
		}
	}


}

temporalTerm = temporalTermInt/integralFlux;
Info<< "temporalTerm: " << temporalTerm << endl;

leakageTerm = leakageTermInt/integralFlux;
Info<< "leakageTerm: " << leakageTerm << endl;

disappearenceTerm = disappearenceTermInt/integralFlux;
Info << "disappearenceTerm: " << disappearenceTerm << endl;

sourceTerm = sourceTermInt/integralFlux;
Info<< "sourceTerm: " << sourceTerm << endl;

scatteringSourceTerm = scatteringSourceTermInt/integralFlux;
Info<< "scatteringSourceTerm: " << scatteringSourceTerm << endl;

Info<< "precursorTerm: " << precursorTerm << endl;



Info<< "Predict new integral flux " << endl;
dimensionedScalar newIntegralFlux = (precursorTerm + integralFluxOld*temporalTerm)/(temporalTerm-leakageTerm+disappearenceTerm-sourceTerm-scatteringSourceTerm);

Info<< "Old integral flux: " << integralFlux << endl;
Info<< "New predicted integral flux: " << newIntegralFlux << endl;

Info<< "Reconstruct fluxes " << endl;
forAll(flux, energyI)
{
      flux[energyI] *= newIntegralFlux/integralFlux;
}


}

