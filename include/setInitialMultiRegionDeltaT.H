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
    setInitialDeltaT

Description
    Set the initial timestep.

\*---------------------------------------------------------------------------*/

if (adjustTimeStep)
{
    #include "calcPower.H"
    pTotOld = pTot;

    //if ((runTime.timeIndex() == 0) && ((CoNum > SMALL) || (DiNum > SMALL) || (DiNumFuel > SMALL)))
    if (((CoNum > SMALL) || (DiNum > SMALL) || (DiNumFuel > SMALL)))
    {
        if (CoNum < SMALL)
        {
            CoNum = SMALL;
        }

        if (DiNum < SMALL)
        {
            DiNum = SMALL;
        }

        if (DiNumFuel < SMALL)
        {
            DiNumFuel = SMALL;
        }

        if(!solveEnergy || implicitFuel)
        {
			DiNumFuel = SMALL;
		}

        runTime.setDeltaT
        (
            min
            (
                min(maxCo/CoNum, min(maxDi/DiNum,maxDiFuel/DiNumFuel))*runTime.deltaT().value(),
                min(maxDeltaT,runTime.controlDict().lookupOrDefault("deltaT", maxDeltaT))
            )
        );

        Info<< "Initial deltaT = " <<  runTime.deltaT().value() << endl;
    }
}

// ************************************************************************* //
