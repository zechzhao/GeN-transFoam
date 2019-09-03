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
    readTimeControls

Description
    Read the control parameters

\*---------------------------------------------------------------------------*/

const bool adjustTimeStep =
    runTime.controlDict().lookupOrDefault("adjustTimeStep", false);

scalar maxCo =
    runTime.controlDict().lookupOrDefault<scalar>("maxCo", 1.0);

scalar maxDiFuel =
    runTime.controlDict().lookupOrDefault<scalar>("maxDiFuel", 0.5);

scalar maxDeltaT =
    runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);

const bool compr =
    runTime.controlDict().lookupOrDefault("compr", false);

scalar maxPowerVariation =
    runTime.controlDict().lookupOrDefault<scalar>("maxPowerVariation", 0.025);

scalar timeStepResidual =
    runTime.controlDict().lookupOrDefault<scalar>("timeStepResidual", 0.001);

scalar neutronIterationResidual =
    runTime.controlDict().lookupOrDefault<scalar>("neutronIterationResidual", 0.00001);

int maxTimeStepIterations =
    runTime.controlDict().lookupOrDefault<int>("maxTimeStepIterations", 5);

int maxNeutronIterations =
    runTime.controlDict().lookupOrDefault<int>("maxNeutronIterations", 10);

const bool oneGroupAcceleration =
    runTime.controlDict().lookupOrDefault("oneGroupAcceleration", false);

List<int> doNotParametrize(runTime.controlDict().lookupOrDefault<List<int> >("doNotParametrize", List<int>()));


// ************************************************************************* //
