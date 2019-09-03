/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Application
    GeN-transFoam

Description
    Multi-physics solver for nuclear reactor analysis. It couples together
    a multi-scale fine/coarse mesh (porous medium) sub-solver for thermal-hydraulics,
    a multi-group diffusion sub-solver for neutronics, a displacement-based
    sub-solver for thermal-mechanics and a finite-difference model for the
    temperature field in the fuel. It is targeted towards the analysis of
    pin-based reactors (e.g., liquid metal fast reactors or light water reactors)
    or homogeneous reactors (e.g., fast-spectrum molten salt reactors).
    Derived from chtMultiRegionFoam

    Transuranus extension of GeN-Foam solver.

Reference publications
	Carlo Fiorina, Ivor Clifford, Manuele Aufiero, Konstantin Mikityuk,
	"GeN-Foam: a novel OpenFOAM® based multi-physics solver for 2D/3D transient
	analysis of nuclear reactors", Nuclear Engineering and Design, submitted

	Carlo Fiorina, Konstantin Mikityuk, " Application of the new GeN-Foam multi-physics
	solver to the European Sodium Fast Reactor and verification against available codes",
	Proceedings of ICAPP 2015, May 03-06, 2015 - Nice (France), Paper 15226
	
	Xinyu Zhao, Eugene Shwageraus, "Extending the capabilities of fast react transient
	analysis tool to account for fuel behaviour", Proceedings of ICAPP 2018, April 08-11, 
	2018 - Charlotte, Nc, US.

	Xinyu Zhao, Eugene Shwageraus, "Couplign between GeN-Foam and multiphysics solver and 
	Transuranus fuel performce code", Proceedings of PHYSOR 2018, April 22-26, 2018
	-Cancun Mexico. 


Coupling made by:
    Xinyu Zhao <xz323@cam.ac.uk> 


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "rhoThermo.H"
#include "turbulenceModel.H"
#include "fixedGradientFvPatchFields.H"
#include "regionProperties.H"
#include "compressibleCourantNo.H"
#include "solidRegionDiffNo.H"
#include "solidThermo.H"
#include "radiationModel.H"
#include "fvIOoptionList.H"
#include "coordinateSystem.H"
#include "porousMedium.H"
#include "subscaleFuel.H"
#include "FieldFields.H"
#include "FieldField.H"
#include "volPointInterpolation.H"
#include "meshToMesh.H"

#include "transtruct.H"

extern "C" {
void passtotrans(transtruct*, int *icase);
void transinitialiser(transtruct*);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "readPhysicsToSolve.H"


    regionProperties rp(runTime);

    #include "createFluidMeshes.H"
    #include "createSolidMeshes.H"
    #include "createNeutroMesh.H"
    #include "createThermoMechanicalMesh.H"
    #include "createFuelMeshes.H"

    #include "readNuclearDataExt.H"
    #include "readThermoMechanicalProperties.H"
    #include "readFuelData.H"

    #include "createFluidFields.H"
    #include "createSolidFields.H"
    #include "createNeutroFields.H"
    #include "createFuelFields.H"

    #include "initContinuityErrs.H"
    #include "readTimeControls.H"
    #include "readSolidTimeControls.H"
    #include "readSolidDisplacementFoamControls.H"

    #include "createThermoMechanicalFields.H"

    #include "compressibleMultiRegionCourantNo.H"
    #include "solidRegionDiffusionNo.H"
    #include "multiRegionFuelDiffNo.H"
    #include "setInitialMultiRegionDeltaT.H"

    #include "createMeshInterpolators.H"

    #include "openOutputFiles.H"

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "readSolidTimeControls.H"
        #include "readPIMPLEControls.H"
        #include "readSolidDisplacementFoamControls.H"

        #include "compressibleMultiRegionCourantNo.H"
        #include "solidRegionDiffusionNo.H"
        #include "multiRegionFuelDiffNo.H"
        if((runTime.timeIndex()-runTime.startTimeIndex())>0)
        {
        	#include "setMultiRegionDeltaT.H"
		}
        runTime++;

        Info << "Time = " << runTime.timeName() << nl << endl;


        if (nOuterCorr != 1)
        {
            forAll(fluidRegions, i)
            {
                #include "setRegionFluidFields.H"
                #include "storeOldFluidFields.H"
            }
        }


        // --- PIMPLE loop
        for (int oCorr=0; oCorr<nOuterCorr; oCorr++)
        {
			Info << "PIMPLE iteration no:  " << oCorr << nl << endl;
            bool finalIter = oCorr == nOuterCorr-1;

            forAll(fluidRegions, i)
            {
                Info<< "\nSolving for fluid region "
                    << fluidRegions[i].name() << endl;
                #include "setRegionFluidFields.H"
                #include "readFluidMultiRegionPIMPLEControls.H"
                #include "solveFluid.H"
            }

            forAll(solidRegions, i)
            {
                Info<< "\nSolving for solid region "
                    << solidRegions[i].name() << endl;
                #include "setRegionSolidFields.H"
                #include "readSolidMultiRegionPIMPLEControls.H"
                #include "solveSolid.H"
            }

        }

		#include "writeOutputs.H"

		runTime.write();


        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
