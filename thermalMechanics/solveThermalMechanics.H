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

Description
    Transient segregated finite-volume solver of linear-elastic,
    small-strain deformation of a solid body, with optional thermal
    diffusion and thermal stresses.

    Simple linear elasticity structural analysis code.
    Solves for the displacement vector field D, also generating the
    stress tensor field sigma.

    Derived from solidDisplacementFoam, modified for multi-zone
    applications through use of a dictionary specifying the properties in each
    mesh region

Reference publication (for the modification of solidDisplacementFoam)
	Carlo Fiorina, Ivor Clifford, Manuele Aufiero, Konstantin Mikityuk,
	"GeN-Foam: a novel OpenFOAM® based multi-physics solver for 2D/3D transient
	analysis of nuclear reactors", Nuclear Engineering and Design, 2015

Author (for the modification of solidDisplacementFoam)
    Carlo Fiorina <carlo.fiorina@hotmail.it; carlo.fiorina@psi.ch; carlo.fiorina@polimi.it>

\*---------------------------------------------------------------------------*/





    Info<< "\nCalculating displacement field\n" << endl;

        int iCorrMech = 0;
        scalar initialResidual = 0;

        do
        {
			Info << "ThermoMechanical iteration no:  " << iCorrMech << nl << endl;
            {
                volScalarField& Tmech = Tptr();
                //copy from porous solver instead of solving

                /*
                solve
                (
                    fvm::ddt(Tmech) == fvm::laplacian(DT, Tmech)
                );
                */

				Tmech *= 0;
				TavFuelMech *= 0;


				fluidToTM.mapSrcToTgt(Tstructures, plusEqOp<scalar>(), Tmech);
				if(linkedFuel)
				{
					fluidToTM.mapSrcToTgt(TavClad, plusEqOp<scalar>(), TavFuelMech);
				}
				else
				{
					fluidToTM.mapSrcToTgt(TavFuel, plusEqOp<scalar>(), TavFuelMech);
				}


				Tmech.correctBoundaryConditions();
				TavFuelMech.correctBoundaryConditions();
            }


            {
                fvVectorMatrix DEqn
                (
                    fvm::d2dt2(Disp)
                 ==
                    fvm::laplacian(2*muMech + lambdaMech, Disp, "laplacian(DD,D)")
                  + divSigmaExp

                );

                //if (thermalStress)
                {
                    const volScalarField& Tmech = Tptr();
                    DEqn += fvc::grad(threeKalpha*(Tmech-TrefStructures));
                }

                initialResidual = DEqn.solve().initialResidual();

                if (!compactNormalStress)
                {
                    divSigmaExp = fvc::div(DEqn.flux());
                }
            }

            {
                volTensorField gradD(fvc::grad(Disp));
                sigmaD = muMech*twoSymm(gradD) + (lambdaMech*I)*tr(gradD);

                if (compactNormalStress)
                {
                    divSigmaExp = fvc::div
                    (
                        sigmaD - (2*muMech + lambdaMech)*gradD,
                        "div(sigmaD)"
                    );
                }
                else
                {
                    divSigmaExp += fvc::div(sigmaD);
                }
            }

        } while (initialResidual > convergenceTolerance && ++iCorrMech < nCorrMech);

		residual = max(residual,initialResidual);

        #include "calculateStress.H"

        // calculating axial fuel expansion
		{

			volScalarField& Tmech = Tptr();
			Tmech *= 0;
			fluidToTM.mapSrcToTgt(Tstructures, plusEqOp<scalar>(), Tmech);
			Tmech.correctBoundaryConditions();


			surfaceScalarField faceArea(fuelOrientation & mechMesh.Sf());

			fvScalarMatrix fuelDispEqn
			(
				fvm::div(faceArea,fuelDisp,"div(fuelDisp)")
			 ==
				alphaFuel*(TavFuelMech-TrefFuel)
			);

			vector CROrientation = -fuelOrientation;
			surfaceScalarField faceAreaCR(CROrientation & mechMesh.Sf());

			fvScalarMatrix CRDispEqn
			(
				fvm::div(faceAreaCR,CRDisp,"div(CRDisp)")
			 ==
				-alphaCR*(Tmech-TrefCR)
			);

			fuelDispEqn.solve();
			CRDispEqn.solve();

			fuelDispVector = (fuelDisp+CRDisp)*fuelOrientation;
		//Info << "faceArea: " << faceArea << endl;
		//Info << "mechMesh.Sf(): " << mechMesh.Sf() << endl;
		//Info << "alphaFuel*(TavFuelMech-TrefFuel): " << alphaFuel*(TavFuelMech-TrefFuel) << endl;
		}
// ************************************************************************* //
