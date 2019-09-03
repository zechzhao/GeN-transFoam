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

\*---------------------------------------------------------------------------*/

Info<< "Writing output files\n" << endl;

//powerOctave << runTime.time().value() << ", " << p_tot << ";" << endl;

powerDat << runTime.time().value() << "\t" << pTot << endl;

//fluxOctave << runTime.time().value();

fluxDat << runTime.time().value();

forAll(flux, energyI)
{
    //fluxOctave << ", " << fvc::domainIntegrate(flux[energyI]()).value();
    fluxDat << "\t" << fvc::domainIntegrate(flux[energyI]).value();
}

//fluxOctave << "; " << endl;
fluxDat << endl;

//precsOctave << runTime.time().value();

precsDat << runTime.time().value();

forAll(prec, precI)
{
    //precsOctave << ", " << fvc::domainIntegrate(prec[precI]()).value();
    precsDat << "\t" << fvc::domainIntegrate(prec[precI]).value();
}

//precsOctave << "; " << endl;
precsDat << endl;

k_eff << keff << endl ;



{

//volScalarField& T = Tptr();
dimensionedScalar Tav = fvc::domainIntegrate(thermoFluid[0].T()*(TavFuelFluid[0]/(TavFuelFluid[0]+dimensionedScalar("",dimensionSet(0,0,0,1,0),SMALL))))
						/ fvc::domainIntegrate((TavFuelFluid[0]/(TavFuelFluid[0]+dimensionedScalar("",dimensionSet(0,0,0,1,0),SMALL))));

dimensionedScalar TavAvFuel = fvc::domainIntegrate(TavFuelFluid[0])
						/ fvc::domainIntegrate((TavFuelFluid[0]/(TavFuelFluid[0]+dimensionedScalar("",dimensionSet(0,0,0,1,0),SMALL))));

dimensionedScalar TavAvClad = fvc::domainIntegrate(TavCladFluid[0])
						/ fvc::domainIntegrate((TavFuelFluid[0]/(TavFuelFluid[0]+dimensionedScalar("",dimensionSet(0,0,0,1,0),SMALL))));


scalar Tmax = max((thermoFluid[0].T().internalField())*(TavFuelFluid[0]/(TavFuelFluid[0]+dimensionedScalar("",dimensionSet(0,0,0,1,0),SMALL)))) ;

reduce(Tmax, maxOp<scalar>());

scalar Tmin = min(thermoFluid[0].T().internalField()) ;

reduce(Tmin, minOp<scalar>());

scalar TMaxFuel = max(TinnerFuelFluid[0].internalField()) ;

reduce(TMaxFuel, maxOp<scalar>());

scalar TMaxClad = max(TinnerCladFluid[0].internalField()) ;

reduce(TMaxClad, maxOp<scalar>());

//TmaxOctave << runTime.time().value() << ", " << Tmax<< ";" << endl;

TmaxDat << runTime.time().value() << "\t" << Tmax << endl;

//TavOctave << runTime.time().value() << ", " << Tav.value() << ";" << endl;

TavDat << runTime.time().value() << "\t" << Tav.value() << endl;

TavFuelDat << runTime.time().value() << "\t" << TavAvFuel.value() << endl;

TavCladDat << runTime.time().value() << "\t" << TavAvClad.value() << endl;

TMaxFuelDat << runTime.time().value() << "\t" <<  TMaxFuel << endl;

TMaxCladDat << runTime.time().value() << "\t" <<  TMaxClad << endl;

//TminOctave << runTime.time().value() << ", " << Tmin<< ";" << endl;

TminDat << runTime.time().value() << "\t" << Tmin << endl;

//totalDisplacedCellsVolumeOctave << runTime.time().value() << ", " << fvc::domainIntegrate((TavFuel/(TavFuel+dimensionedScalar("",dimensionSet(0,0,0,1,0),SMALL)))) << ";" << endl;

totalDisplacedCellsVolumeDat << runTime.time().value() << "\t" << (fvc::domainIntegrate((sigmaPow/(sigmaPow+dimensionedScalar("",dimensionSet(1,1,-2,0,0,0,0),SMALL))))).value() << endl;

}