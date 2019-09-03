
	Info << "Solving neutronics" << endl ;

    neutroToFluid.mapTgtToSrc( TavFuel, plusEqOp<scalar>(), TfuelNeutro.internalField());
	neutroToFluid.mapTgtToSrc( TavClad, plusEqOp<scalar>(), TcladNeutro.internalField());
	neutroToFluid.mapTgtToSrc( thermo.rho(), plusEqOp<scalar>(), rhoCoolNeutro.internalField());
// Info<<"thermo.rho "<< thermo.rho() << endl;
// Info<<"rhoCoolNeutro" << rhoCoolNeutro << endl;
	if (eigenvalueNeutronics)
	{
		#include "calcPower.H"
		#include "normFluxes.H"
		#include "calcPower.H"
	}
    #include "initializeNeutroSources.H"

	scalar neutroResidual = 0.0;
	label neutroIter = 0;

	if(oneGroupAcceleration && couplingIter == 0 && !eigenvalueNeutronics)
	{
		#include "integralAcceleration.H"
	}

	do
	{
		Info << "Neutron iteration no:  " << neutroIter << nl << endl;
		forAll(prec,precI)
		{
			#include "setPrecConst.H"
			fvScalarMatrix precEq
			(
				fvm::ddt((1-eigenvalueNeutronics),prec[precI])
			  + fvm::Sp(lambda, prec[precI])
			  - neutroSource/keff*Beta
			);
			precEq.solve();
			delayedNeutroSource += lambda*(prec[precI]-prec[precI].oldTime());
		}

		volFuelPowerNeutro*=0.0;
		forAll(flux,energyI)
		{


			neutroSource -= nuSigmaEff[energyI]*flux[energyI]*(1-eigenvalueNeutronics);
			#include "setNeutroDiffConst.H"
			#include "calcScatteringSource.H"
			fvScalarMatrix neutroEq
			(
				fvm::ddt((1-eigenvalueNeutronics)*IV/discFactor[energyI],flux[energyI])
			  - fvm::laplacian(D/discFactor[energyI],flux[energyI])
			  - fvm::Sp((nuSigmaEff[energyI]/keff*(1.0-BetaTot)*chiPrompt*(1-eigenvalueNeutronics)-sigmaDisapp)/discFactor[energyI],flux[energyI])
			  - neutroSource/keff*(1.0-BetaTot)*chiPrompt
			  - delayedNeutroSource*chiDelayed
			  - scatteringSourceExtra
			 );
			neutroResidual = neutroEq.solve().initialResidual();
			flux[energyI] = max(flux[energyI],mag(flux[energyI]*SMALL));
			flux[energyI]/=discFactor[energyI];
			neutroSource += nuSigmaEff[energyI]*flux[energyI]*(1-eigenvalueNeutronics);
			volFuelPowerNeutro += flux[energyI] * sigmaPow;
		}
		neutroIter++;
	}while(((neutroResidual > neutronIterationResidual) && (neutroIter < maxNeutronIterations) && !eigenvalueNeutronics));

	if(oneGroupAcceleration)
	{
		oneGroupFlux *= 0.0;
		forAll(flux,energyI)
		{
			oneGroupFlux+=flux[energyI];
		}
	}

	volFuelPowerNeutro /= max(fuelFraction,SMALL) ; //determine vol power IN the fuel

	//project volPowerNeutro back to volPower on the fluid, taking into account mesh deformation

	volFuelPowerNeutro.internalField() *= displacedCellsVolumes/staticCellsVolumes;
	volFuelPowerNeutro.correctBoundaryConditions();

	volFuelPower *= 0;
	neutroToFluid.mapSrcToTgt(volFuelPowerNeutro, plusEqOp<scalar>(), volFuelPower);
	//volFuelPower.internalField() *= (1/staticCellsVolumes);
	//volFuelPower.internalField() *= (1/mechMesh.V());
	volFuelPower.correctBoundaryConditions();

	volFuelPowerNeutro.internalField() /= displacedCellsVolumes;
	volFuelPowerNeutro.correctBoundaryConditions();

	if (eigenvalueNeutronics)
	{
		#include "calcPower.H"
		Info << "pTarget = " << pTarget << nl << endl;
		Info << "keff previous = " << keff << nl << endl;
		keff *= pTot/pTarget;
		keff = min(keff,5);
		keff = max(keff,0.1);
		Info << "keff = " << keff << nl << endl;
		#include "normFluxes.H"
	}



