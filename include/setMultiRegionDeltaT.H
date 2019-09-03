/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Global
    setMultiRegionDeltaT

Description
    Reset the timestep to maintain a constant maximum courant and
    diffusion Numbers. Reduction of time-step is immediate, but
    increase is damped to avoid unstable oscillations.

\*---------------------------------------------------------------------------*/

if (adjustTimeStep)
{

    if (CoNum == -GREAT)
    {
        CoNum = SMALL;
    }

    if (DiNum == -GREAT)
    {
        DiNum = SMALL;
    }

    if (DiNum == -GREAT)
    {
        DiNumFuel = SMALL;
    }

	if(!solveEnergy || implicitFuel)
	{
		DiNumFuel = SMALL;
	}

	#include "calcPower.H"

	scalar powerVariation = mag((pTot - pTotOld) / (pTotOld + SMALL)) ;


    scalar maxDeltaTFluid = maxCo/(CoNum + SMALL);
    scalar maxDeltaTSolid = maxDi/(DiNum + SMALL);
    scalar maxDeltaTFuel = maxDiFuel/(DiNumFuel + SMALL);

    scalar maxDeltaTNeutro = mag(maxPowerVariation / (powerVariation + SMALL)) ;


    scalar deltaTFluid =
        min
        (
            min(maxDeltaTFluid, 1.0 + 0.1*maxDeltaTFluid),
            1.2
        );

	if(!solveFluidMechanics)
	{
		deltaTFluid = GREAT;
	}


    runTime.setDeltaT
    (
        min
        (
            min(deltaTFluid, min(maxDeltaTSolid,min(maxDeltaTFuel,maxDeltaTNeutro)))*runTime.deltaT().value(),
            //min(deltaTFluid, min(maxDeltaTSolid,maxDeltaTFuel))*runTime.deltaT().value(),
            maxDeltaT
        )
    );

	if((maxDeltaTFuel < maxDeltaTFluid) & (maxDeltaTFuel < maxDeltaTSolid))
	{
		Info<< "Warning: time step determined by conduction in the fuel " << nl <<
		"Consider reducing the mesh size in fuel and cladding, or use implicit solver for fuel" << endl;
	}

	pTotOld = pTot;

    Info<< "deltaT = " <<  runTime.deltaT().value() << endl;
}

// ************************************************************************* //
