
	volScalarField DNeutroRad = DNeutro & radialOrientation;
	volScalarField DNeutroAx = DNeutro & axialOrientation;

	DNeutroRad.correctBoundaryConditions();
	DNeutroAx.correctBoundaryConditions();

	radExp = fvc::grad(DNeutroRad)& radialOrientation;
	axExp = fvc::grad(DNeutroAx)& axialOrientation;

	volScalarField logT = log(max(0.01,(TfuelNeutro/dimensionedScalar("",dimensionSet(0,0,0,1,0,0,0),TfuelRef))));

	bool parametrize = true;
	forAll(doNotParametrize, groupI)
	{
		if(energyI == doNotParametrize[groupI])
		{
			parametrize = false;
			Info << "Not Parametrizing xs for group: " << doNotParametrize[groupI] << endl;
		}
	}


	forAll(CRentries,CRzoneI)
	{

		if(runTime.value()>CRstart[CRzoneI] && runTime.value()<CRfinish[CRzoneI])
		{
			CRposition.set(CRzoneI,new scalar(CRinitialPosition[CRzoneI]-CRspeed[CRzoneI]*(runTime.value()-CRstart[CRzoneI])));
			Info << "Control rod position " << CRposition[CRzoneI] << endl;
		}


	}



	forAll(entries, zoneI)
	{

		label zone = zoneI;

		const word& name = entries[zoneI].keyword();

		label zoneId = neutroMesh.cellZones().findZoneID(name);
		if(zoneId == -1)
		{
		  Info << "Mesh zone: " << name << " does not exists" << endl;
		}


		//check if you are in a follower or CR
		bool isFollower = false;
		bool isCR = false;
		bool CRisInserting = true;
		label CRzone = 0;

		forAll(CRentries,CRzoneI)
		{

			const word& CRname = CRentries[CRzoneI].keyword();
			if(name == CRname)
			{
				isCR = true ;
				CRzone = CRzoneI;
			}
			if(name == CRFollowerName[CRzoneI])
			{
				isFollower = true ;
				CRzone = CRzoneI;
			}

		}


		if(CRNumber>0)
		{
			if(CRposition[CRzone]>CRinitialPosition[CRzone])
			{
				CRisInserting = false;
			}
		}

		forAll(neutroMesh.cellZones()[zoneId], cellIlocal)
		{

			label cellIglobal = neutroMesh.cellZones()[zoneId][cellIlocal];
			zone = zoneI;

			if(CRNumber>0)
			{
				if(isFollower && CRisInserting && neutroMesh.C()[cellIglobal].z()>CRposition[CRzone])
				{
					forAll(entries, zoneII)
					{
						if(entries[zoneII].keyword() == CRentries[CRzone].keyword())
						{
							zone = zoneII;
						}
					}
				}

				if(isCR && !CRisInserting && neutroMesh.C()[cellIglobal].z()<CRposition[CRzone])
				{
					forAll(entries, zoneII)
					{
						if(entries[zoneII].keyword() == CRFollowerName[CRzone])
						{
							zone = zoneII;
						}
					}
				}
			}



			//constant
			IV[cellIglobal] = IVList[zone][energyI];
			chiPrompt[cellIglobal] = chiPromptList[zone][energyI];
			chiDelayed[cellIglobal] = chiDelayedList[zone][energyI];
			BetaTot[cellIglobal] = BetaTotList[zone];

			//fields dependent
			D[cellIglobal] = DList[zone][energyI] + parametrize * (
								  fuelTempDList[zone][energyI]*logT[cellIglobal]
								+ axialExpDList[zone][energyI]*(axExp[cellIglobal])
								+ radialExpDList[zone][energyI]*(radExp[cellIglobal])
								+ rhoCoolDList[zone][energyI]*(rhoCoolNeutro[cellIglobal] - rhoCoolRef)
								+ cladExpDList[zone][energyI]*(TcladNeutro[cellIglobal] - TcladRef) );

			nuSigmaEff[energyI][cellIglobal] = nuSigmaEffList[zone][energyI] + parametrize * (
								  fuelTempNuSigmaEffList[zone][energyI]*logT[cellIglobal]
								+ axialExpNuSigmaEffList[zone][energyI]*(axExp[cellIglobal])
								+ radialExpNuSigmaEffList[zone][energyI]*(radExp[cellIglobal])
								+ rhoCoolNuSigmaEffList[zone][energyI]*(rhoCoolNeutro[cellIglobal] - rhoCoolRef)
								+ cladExpNuSigmaEffList[zone][energyI]*(TcladNeutro[cellIglobal] - TcladRef) );

			sigmaPow[cellIglobal] = sigmaPowList[zone][energyI] + parametrize * (
								  fuelTempSigmaPowList[zone][energyI]*logT[cellIglobal]
								+ axialExpSigmaPowList[zone][energyI]*(axExp[cellIglobal])
								+ radialExpSigmaPowList[zone][energyI]*(radExp[cellIglobal])
								+ rhoCoolSigmaPowList[zone][energyI]*(rhoCoolNeutro[cellIglobal] - rhoCoolRef)
								+ cladExpSigmaPowList[zone][energyI]*(TcladNeutro[cellIglobal] - TcladRef) );

			sigmaDisapp[cellIglobal] = sigmaDisappList[zone][energyI] + parametrize * (
								  fuelTempSigmaDisappList[zone][energyI]*logT[cellIglobal]
								+ axialExpSigmaDisappList[zone][energyI]*(axExp[cellIglobal])
								+ radialExpSigmaDisappList[zone][energyI]*(radExp[cellIglobal])
								+ rhoCoolSigmaDisappList[zone][energyI]*(rhoCoolNeutro[cellIglobal] - rhoCoolRef)
								+ cladExpSigmaDisappList[zone][energyI]*(TcladNeutro[cellIglobal] - TcladRef) );

			forAll(sigmaFrom,energyJ)
			{
					sigmaFrom[energyJ][cellIglobal] = sigmaFromToList[zone][energyJ][energyI] + parametrize * (
									  fuelTempSigmaFromToList[zone][energyJ][energyI]*logT[cellIglobal]
									+ axialExpSigmaFromToList[zone][energyJ][energyI]*(axExp[cellIglobal])
									+ radialExpSigmaFromToList[zone][energyJ][energyI]*(radExp[cellIglobal])
									+ rhoCoolSigmaFromToList[zone][energyJ][energyI]*(rhoCoolNeutro[cellIglobal] - rhoCoolRef)
									+ cladExpSigmaFromToList[zone][energyJ][energyI]*(TcladNeutro[cellIglobal] - TcladRef) );

			}

			fuelFraction[cellIglobal] = fuelFractionList[zone];

			discFactor[energyI][cellIglobal] = discFactorList[zone][energyI];

		}

		IV.correctBoundaryConditions();
		D.correctBoundaryConditions();
		nuSigmaEff[energyI].correctBoundaryConditions();
		sigmaPow.correctBoundaryConditions();
		sigmaDisapp.correctBoundaryConditions();
		forAll(sigmaFrom,energyJ)
		{
			sigmaFrom[energyJ].correctBoundaryConditions();
		}
		chiPrompt.correctBoundaryConditions();
		chiDelayed.correctBoundaryConditions();
		BetaTot.correctBoundaryConditions();
		fuelFraction.correctBoundaryConditions();
		discFactor[energyI].correctBoundaryConditions();
	}




