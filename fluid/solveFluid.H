if (finalIter)
{
    mesh.data::add("finalIteration", true);
}

if(compr)
{
	if (oCorr == 0)
	{
	    #include "rhoEqn.H"
	}
}
else
{
	rho = rhoRef;
}


#include "UEqn.H"


if(finalIter)
{
	scalar residual = 0.0;
	label couplingIter = 0;
	do
	{
		Info << "Coupling iteration no:  " << couplingIter << nl << endl;
		if(solveEnergy)
		{
			Info << "runtime value = " << runTime.value() <<"    deltaT = " << runTime.deltaTValue() << endl;
			if(runTime.value() > runTime.deltaTValue())
			{
				#include "solveFuel.H"
			}
			#include "EEqn.H"
		}
		if(solveThermalMechanics)
		{
			#include "solveThermalMechanics.H"
			#include "displaceNeutroMesh.H"
		}
		if(solveNeutronics)
		{
			#include "solveNeutronics.H"
		}
		couplingIter++;
	}while(residual > timeStepResidual && couplingIter < maxTimeStepIterations && !eigenvalueNeutronics && tighltyCoupled);

	if (!eigenvalueNeutronics && powerControl)
	{
		Info << "Power controller " << nl << endl;
		#include "calcPower.H"
		Info << "pTarget = " << pTarget << nl << endl;
		Info << "keff previous = " << keff << nl << endl;
		scalar e = (pTot-pTarget)/pTarget;
		ie += e * neutroMesh.time().deltaT().value();
		keff = keff + powerControlPorportional * e + powerControlIntegral * ie;
		Info << "keff = " << keff << nl << endl;
	}
}
else
{
	scalar residual = 0.0;
	if(solveEnergy)
	{
		#include "EEqn.H"
	}
}
if(solveFluidMechanics)
{
	for (int corr=0; corr<nCorr; corr++)
	{
		#include "pEqn.H"
	}
	turb.correct();
}


if(compr)
{
	rho = thermo.rho();
}

if (finalIter)
{
    mesh.data::remove("finalIteration");
    Re = porousMedium.ReynoldsNumber();
}
