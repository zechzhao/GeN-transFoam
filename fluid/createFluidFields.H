    // Initialise fluid field pointer lists

    PtrList<rhoThermo> thermoFluid(fluidRegions.size());
    PtrList<volScalarField> rhoFluid(fluidRegions.size());
    PtrList<volScalarField> psiFluid(fluidRegions.size());
    PtrList<volScalarField> rhokFluid(fluidRegions.size());
    PtrList<uniformDimensionedScalarField> betaFluid(fluidRegions.size());
    PtrList<uniformDimensionedScalarField> TrefFluid(fluidRegions.size());
    PtrList<uniformDimensionedScalarField> rhoRefFluid(fluidRegions.size());
    PtrList<volVectorField> realUFluid(fluidRegions.size());
    PtrList<volVectorField> UFluid(fluidRegions.size());
    PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
    PtrList<uniformDimensionedVectorField> gFluid(fluidRegions.size());
    PtrList<compressible::turbulenceModel> turbulence(fluidRegions.size());
    PtrList<Foam::porousMedium> porousMediumFluid(fluidRegions.size());
    PtrList<Foam::subscaleFuel> subscaleFuelFluid(fluidRegions.size());
    PtrList<volScalarField> p_rghFluid(fluidRegions.size());
    PtrList<volScalarField> ghFluid(fluidRegions.size());
    PtrList<surfaceScalarField> ghfFluid(fluidRegions.size());
    PtrList<radiation::radiationModel> radiation(fluidRegions.size());
    PtrList<volScalarField> KFluid(fluidRegions.size());
    PtrList<volScalarField> dpdtFluid(fluidRegions.size());
    PtrList<volScalarField> ReFluid(fluidRegions.size());
    List<scalar> initialMassFluid(fluidRegions.size());
    PtrList<fv::IOoptionList> fluidFvOptions(fluidRegions.size());
	PtrList<volScalarField> volFuelPowerFluid(fluidRegions.size());
    PtrList<volScalarField> TavFuelFluid(fluidRegions.size());
    PtrList<volScalarField> TavCladFluid(fluidRegions.size());
    PtrList<volScalarField> TinnerFuelFluid(fluidRegions.size());
    PtrList<volScalarField> TouterFuelFluid(fluidRegions.size());
    PtrList<volScalarField> TinnerCladFluid(fluidRegions.size());
    PtrList<volScalarField> TouterCladFluid(fluidRegions.size());

    PtrList<volScalarField> TstructuresFluid(fluidRegions.size());
    // Populate fluid field pointer lists
	forAll(fluidRegions, i)
	{

        Info<< "*** Reading fluid mesh thermophysical properties for region "
            << fluidRegions[i].name() << nl << endl;

        Info<< "    Adding to thermoFluid\n" << endl;

        thermoFluid.set
        (
            i,
            rhoThermo::New(fluidRegions[i]).ptr()
        );

        Info<< "    Adding to betaFluid\n" << endl;
        betaFluid.set
        (
            i,
            new uniformDimensionedScalarField
            (
                IOobject
                (
                    "beta",
                    runTime.constant(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );

        Info<< "    Adding to TrefFluid\n" << endl;
        TrefFluid.set
        (
            i,
            new uniformDimensionedScalarField
            (
                IOobject
                (
                    "Tref",
                    runTime.constant(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );

        Info<< "    Adding to rhoRefFluid\n" << endl;
        rhoRefFluid.set
        (
            i,
            new uniformDimensionedScalarField
            (
                IOobject
                (
                    "rhoRef",
                    runTime.constant(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );


		Info<< "    Adding to rhoFluid\n" << endl;
		rhoFluid.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					"rho",
					runTime.timeName(),
					fluidRegions[i],
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
				thermoFluid[i].rho()
			)
		);


        Info<< "    Adding to psiFluid\n" << endl;
        psiFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "psi",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                thermoFluid[i].psi()
            )
        );

        Info<< "    Adding to rhokFluid\n" << endl; //for Boussinesq
        rhokFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "rhok",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
        	fluidRegions[i]
            )
        );

        Info<< "    Adding to UFluid\n" << endl;
        UFluid.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    "U",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i]
            )
        );

		Info<< "    Adding to realUFluid\n" << endl;
        realUFluid.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    "realU",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                UFluid[i]
            )
        );

        Info<< "    Adding to phiFluid\n" << endl;
        phiFluid.set
        (
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    "phi",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                linearInterpolate(rhoFluid[i]*realUFluid[i])
                    & fluidRegions[i].Sf()
            )
        );

        Info<< "    Adding to gFluid\n" << endl;
        gFluid.set
        (
            i,
            new uniformDimensionedVectorField
            (
                IOobject
                (
                    "g",
                    runTime.constant(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );

        Info<< "    Adding porousMedium \n" << endl;
        porousMediumFluid.set
        (
            i,
            Foam::porousMedium::New
            (
                rhoFluid[i],
                realUFluid[i],
                thermoFluid[i]
            ).ptr()
        );

        Info<< "    Adding subscaleFuel \n" << endl;
        subscaleFuelFluid.set
        (
            i,
            Foam::subscaleFuel::New
            (
                fluidRegions[i]
            ).ptr()
        );

        Info<< "    Adding turbulence \n" << endl;
        turbulence.set
        (
            i,
            compressible::turbulenceModel::New
            (
                rhoFluid[i],
                UFluid[i],
                phiFluid[i],
                thermoFluid[i]
            ).ptr()
        );

        Info<< "    Adding to ghFluid\n" << endl;
        ghFluid.set
        (
            i,
            new volScalarField("gh", gFluid[i] & fluidRegions[i].C())
        );

        Info<< "    Adding to ghfFluid\n" << endl;
        ghfFluid.set
        (
            i,
            new surfaceScalarField("ghf", gFluid[i] & fluidRegions[i].Cf())
        );

        p_rghFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "p_rgh",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i]
            )
        );

        // Force p_rgh to be consistent with p
        p_rghFluid[i] = thermoFluid[i].p() - rhoFluid[i]*ghFluid[i];

        radiation.set
        (
            i,
            radiation::radiationModel::New(thermoFluid[i].T())
        );

        initialMassFluid[i] = fvc::domainIntegrate(rhoFluid[i]).value();

        Info<< "    Adding to KFluid\n" << endl;
        KFluid.set
        (
            i,
            new volScalarField
            (
                "K",
                0.5*magSqr(UFluid[i])
            )
        );

        Info<< "    Adding to dpdtFluid\n" << endl;
        dpdtFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "dpdt",
                    runTime.timeName(),
                    fluidRegions[i]
                ),
                fluidRegions[i],
                dimensionedScalar
                (
                    "dpdt",
                    thermoFluid[i].p().dimensions()/dimTime,
                    0
                )
            )
        );

        Info<< "    Adding fvOptions\n" << endl;
        fluidFvOptions.set
        (
            i,
            new fv::IOoptionList(fluidRegions[i])
        );

        Info<< "    Adding to ReFluid\n" << endl;
        ReFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "Re",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i],
				dimensionedScalar("", dimless, 0.0),
				zeroGradientFvPatchScalarField::typeName
            )
        );

        Info<< "    Adding to volFuelPowerFluid\n" << endl;
        volFuelPowerFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "volFuelPower",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i],
				dimensionedScalar("", dimensionSet(1,-1,-3,0,0,0,0), 0.0),
				zeroGradientFvPatchScalarField::typeName
            )
        );

        Info<< "    Adding to TavFuelFluid\n" << endl;
        TavFuelFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "TavFuel",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i],
				dimensionedScalar("", dimensionSet(0,0,0,1,0,0,0), 0.0),
				zeroGradientFvPatchScalarField::typeName
            )
        );

         Info<< "    Adding to TavCladFluid\n" << endl;
         TavCladFluid.set
         (
             i,
             new volScalarField
             (
                 IOobject
                 (
                     "TavClad",
                     runTime.timeName(),
                     fluidRegions[i],
                     IOobject::READ_IF_PRESENT,
                     IOobject::AUTO_WRITE
                 ),
                 fluidRegions[i],
 				dimensionedScalar("", dimensionSet(0,0,0,1,0,0,0), 0.0),
 				zeroGradientFvPatchScalarField::typeName
             )
        );

        Info<< "    Adding to TinnerFuelFluid\n" << endl;
        TinnerFuelFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "TinnerFuel",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i],
				dimensionedScalar("", dimensionSet(0,0,0,1,0,0,0), 0.0),
				zeroGradientFvPatchScalarField::typeName
            )
        );

        Info<< "    Adding to TouterFuelFluid\n" << endl;
        TouterFuelFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "TouterFuel",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i],
				dimensionedScalar("", dimensionSet(0,0,0,1,0,0,0), 0.0),
				zeroGradientFvPatchScalarField::typeName
            )
        );

        Info<< "    Adding to TinnerCladFluid\n" << endl;
        TinnerCladFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "TinnerClad",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i],
				dimensionedScalar("", dimensionSet(0,0,0,1,0,0,0), 0.0),
				zeroGradientFvPatchScalarField::typeName
            )
        );

        Info<< "    Adding to TouterCladFluid\n" << endl;
        TouterCladFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "TouterClad",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i],
				dimensionedScalar("", dimensionSet(0,0,0,1,0,0,0), 0.0),
				zeroGradientFvPatchScalarField::typeName
            )
        );

        Info<< "    Adding to TstructuresFluid\n" << endl;
        TstructuresFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "Tstructures",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                thermoFluid[i].T()
            )
        );
    }


