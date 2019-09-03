	Info<< "\nReading nuclear data\n" << endl;

	IOdictionary nuclearData
	(
		IOobject
		(
			"nuclearData",
			runTime.constant(),
			neutroMesh,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		)
	);

	IOdictionary nuclearDataRadialExp
	(
		IOobject
		(
			"nuclearDataRadialExp",
			runTime.constant(),
			neutroMesh,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		)
	);

	IOdictionary nuclearDataAxialExp
		(
			IOobject
			(
				"nuclearDataAxialExp",
				runTime.constant(),
				neutroMesh,
				IOobject::MUST_READ,
				IOobject::NO_WRITE
			)
	);

	IOdictionary nuclearDataFuelTemp
		(
			IOobject
			(
				"nuclearDataFuelTemp",
				runTime.constant(),
				neutroMesh,
				IOobject::MUST_READ,
				IOobject::NO_WRITE
			)
	);

	IOdictionary nuclearDataRhoCool
		(
			IOobject
			(
				"nuclearDataRhoCool",
				runTime.constant(),
				neutroMesh,
				IOobject::MUST_READ,
				IOobject::NO_WRITE
			)
	);

	IOdictionary nuclearDataCladExp
		(
			IOobject
			(
				"nuclearDataCladExp",
				runTime.constant(),
				neutroMesh,
				IOobject::MUST_READ,
				IOobject::NO_WRITE
			)
	);

	label energyGroups(nuclearData.lookupOrDefault("energyGroups",1));
	label precGroups(nuclearData.lookupOrDefault("precGroups",1));
	scalar pTarget(nuclearData.lookupOrDefault("pTarget",1.0));
	scalar pTot(1.0);
	scalar pTotOld(1.0);
	scalar ie(0.0);//intgral of power error used in the PI controller
	scalar keff(nuclearData.lookupOrDefault("keff",1.0));


	PtrList<entry> entries(nuclearData.lookup("zones"));
	label zoneNumber = entries.size();
	PtrList<entry> entriesRadialExp(nuclearDataRadialExp.lookup("zones"));
	vector axialOrientation(nuclearDataRadialExp.lookupOrDefault("axialOrientation",Foam::vector(0,1,0)));
	vector radialOrientation(nuclearDataRadialExp.lookupOrDefault("radialOrientation",Foam::vector(1,0,0)));
	PtrList<entry> entriesAxialExp(nuclearDataAxialExp.lookup("zones"));
	PtrList<entry> entriesFuelTemp(nuclearDataFuelTemp.lookup("zones"));
	PtrList<entry> entriesRhoCool(nuclearDataRhoCool.lookup("zones"));
	PtrList<entry> entriesCladExp(nuclearDataCladExp.lookup("zones"));

	// fuel fraction used to determine vol power IN the fuel starting from vol power from sigmaPow
	PtrList<scalar > fuelFractionList(zoneNumber);

	//discontinuity factors
	PtrList<scalarField > discFactorList(zoneNumber);

	//Nuclear data (assumed to be) not affected by feedbacks
	PtrList<scalarField > IVList(zoneNumber);
	PtrList<scalarField > chiPromptList(zoneNumber);
	PtrList<scalarField > chiDelayedList(zoneNumber);
	PtrList<scalarField > BetaList(zoneNumber);
	PtrList<scalar > BetaTotList(zoneNumber);
	PtrList<scalarField > lambdaList(zoneNumber);

	//Nuclear data affected by feedbacks
		//Nominal
		PtrList<scalarField > DList(zoneNumber);
		PtrList<scalarField > nuSigmaEffList(zoneNumber);
		PtrList<scalarField > sigmaPowList(zoneNumber);
		PtrList<scalarField > sigmaDisappList(zoneNumber);
		PtrList<scalarSquareMatrix > sigmaFromToList(zoneNumber);	//scattering matrix
		//RadExp
		scalar RadExp(nuclearDataRadialExp.lookupOrDefault("expansionFromNominal",1.0));
		PtrList<scalarField > DRadialExpList(zoneNumber);
		PtrList<scalarField > nuSigmaEffRadialExpList(zoneNumber);
		PtrList<scalarField > sigmaPowRadialExpList(zoneNumber);
		PtrList<scalarField > sigmaDisappRadialExpList(zoneNumber);
		PtrList<scalarSquareMatrix > sigmaFromToRadialExpList(zoneNumber);	//scattering matrix
		//AxExp
		scalar AxExp(nuclearDataAxialExp.lookupOrDefault("expansionFromNominal",1.0));
		PtrList<scalarField > DAxialExpList(zoneNumber);
		PtrList<scalarField > nuSigmaEffAxialExpList(zoneNumber);
		PtrList<scalarField > sigmaPowAxialExpList(zoneNumber);
		PtrList<scalarField > sigmaDisappAxialExpList(zoneNumber);
		PtrList<scalarSquareMatrix > sigmaFromToAxialExpList(zoneNumber);	//scattering matrix
		//fuelTemp
		scalar TfuelRef(nuclearDataFuelTemp.lookupOrDefault("TfuelRef",900.0));
		scalar TfuelPerturbed(nuclearDataFuelTemp.lookupOrDefault("TfuelPerturbed",1200.0));
		PtrList<scalarField > DFuelTempList(zoneNumber);
		PtrList<scalarField > nuSigmaEffFuelTempList(zoneNumber);
		PtrList<scalarField > sigmaPowFuelTempList(zoneNumber);
		PtrList<scalarField > sigmaDisappFuelTempList(zoneNumber);
		PtrList<scalarSquareMatrix > sigmaFromToFuelTempList(zoneNumber);	//scattering matrix
		//RhoCool
		scalar rhoCoolRef(nuclearDataRhoCool.lookupOrDefault("rhoCoolRef",860.0));
		scalar rhoCoolPerturbed(nuclearDataRhoCool.lookupOrDefault("rhoCoolPerturbed",1.0));
		PtrList<scalarField > DRhoCoolList(zoneNumber);
		PtrList<scalarField > nuSigmaEffRhoCoolList(zoneNumber);
		PtrList<scalarField > sigmaPowRhoCoolList(zoneNumber);
		PtrList<scalarField > sigmaDisappRhoCoolList(zoneNumber);
		PtrList<scalarSquareMatrix > sigmaFromToRhoCoolList(zoneNumber);	//scattering matrix
		//cladExp
		scalar TcladRef(nuclearDataCladExp.lookupOrDefault("TcladRef",900.0));
		scalar TcladPerturbed(nuclearDataCladExp.lookupOrDefault("TcladPerturbed",1200.0));
		PtrList<scalarField > DCladExpList(zoneNumber);
		PtrList<scalarField > nuSigmaEffCladExpList(zoneNumber);
		PtrList<scalarField > sigmaPowCladExpList(zoneNumber);
		PtrList<scalarField > sigmaDisappCladExpList(zoneNumber);
		PtrList<scalarSquareMatrix > sigmaFromToCladExpList(zoneNumber);	//scattering matrix

	//Feedbacks Coefficients
		// fuel temperature coefficient
		PtrList<scalarField > fuelTempDList(zoneNumber);
		PtrList<scalarField > fuelTempNuSigmaEffList(zoneNumber);
		PtrList<scalarField > fuelTempSigmaPowList(zoneNumber);
		PtrList<scalarField > fuelTempSigmaDisappList(zoneNumber);
		PtrList<scalarSquareMatrix > fuelTempSigmaFromToList(zoneNumber);

		// axial expansion coefficient
		PtrList<scalarField > axialExpDList(zoneNumber);
		PtrList<scalarField > axialExpNuSigmaEffList(zoneNumber);
		PtrList<scalarField > axialExpSigmaPowList(zoneNumber);
		PtrList<scalarField > axialExpSigmaDisappList(zoneNumber);
		PtrList<scalarSquareMatrix > axialExpSigmaFromToList(zoneNumber);

		// radial expansion coefficient
		PtrList<scalarField > radialExpDList(zoneNumber);
		PtrList<scalarField > radialExpNuSigmaEffList(zoneNumber);
		PtrList<scalarField > radialExpSigmaPowList(zoneNumber);
		PtrList<scalarField > radialExpSigmaDisappList(zoneNumber);
		PtrList<scalarSquareMatrix > radialExpSigmaFromToList(zoneNumber);

		// coolant expansion coefficient
		PtrList<scalarField > rhoCoolDList(zoneNumber);
		PtrList<scalarField > rhoCoolNuSigmaEffList(zoneNumber);
		PtrList<scalarField > rhoCoolSigmaPowList(zoneNumber);
		PtrList<scalarField > rhoCoolSigmaDisappList(zoneNumber);
		PtrList<scalarSquareMatrix > rhoCoolSigmaFromToList(zoneNumber);

		// cladding expansion coefficient
		PtrList<scalarField > cladExpDList(zoneNumber);
		PtrList<scalarField > cladExpNuSigmaEffList(zoneNumber);
		PtrList<scalarField > cladExpSigmaPowList(zoneNumber);
		PtrList<scalarField > cladExpSigmaDisappList(zoneNumber);
		PtrList<scalarSquareMatrix > cladExpSigmaFromToList(zoneNumber);

	//Fill the data
	forAll(entries,zoneI)
	{
		dictionary& dict = entries[zoneI].dict();


		const word& name = entries[zoneI].keyword();

		Info << "name: " << name << endl;

		fuelFractionList.set(zoneI,new scalar(dict.lookupOrDefault("fuelFraction",0.0)));


		//Nuclear data (assumed to be) not affected by feedbacks
		IVList.set(zoneI,new scalarField("IV",dict,energyGroups));
		chiPromptList.set(zoneI,new scalarField("chiPrompt",dict,energyGroups));
		chiDelayedList.set(zoneI,new scalarField("chiDelayed",dict,energyGroups));
		discFactorList.set(zoneI,new scalarField("discFactor",dict,energyGroups));
		BetaList.set(zoneI,new scalarField("Beta",dict,precGroups));
		BetaTotList.set(zoneI,new scalar(sum(BetaList[zoneI])));
		lambdaList.set(zoneI,new scalarField("lambda",dict,precGroups));
		Info << "Reading nominal " << endl;
		//Nominal
		DList.set(zoneI,new scalarField("D",dict,energyGroups));
		nuSigmaEffList.set(zoneI,new scalarField("nuSigmaEff",dict,energyGroups));
		sigmaPowList.set(zoneI,new scalarField("sigmaPow",dict,energyGroups));
		sigmaDisappList.set(zoneI,new scalarField("sigmaDisapp",dict,energyGroups));
		sigmaFromToList.set(zoneI,new scalarSquareMatrix(dict.lookup("scatteringMatrix")));
		Info << "Finishing reading nominal " << endl;
		//RadExp
		bool read = false;
		int perturbedEntry = 0;
		forAll(entriesRadialExp,entryI)
		{
		  if(entriesRadialExp[entryI].keyword() == name)
		  {
		    read = true;
		    perturbedEntry = entryI;
		  }
		}
		if(read)
		{
		  dictionary& dictRadialExp = entriesRadialExp[perturbedEntry].dict();
		  DRadialExpList.set(zoneI,new scalarField("D",dictRadialExp,energyGroups));
		  nuSigmaEffRadialExpList.set(zoneI,new scalarField("nuSigmaEff",dictRadialExp,energyGroups));
		  sigmaPowRadialExpList.set(zoneI,new scalarField("sigmaPow",dictRadialExp,energyGroups));
		  sigmaDisappRadialExpList.set(zoneI,new scalarField("sigmaDisapp",dictRadialExp,energyGroups));
		  sigmaFromToRadialExpList.set(zoneI,new scalarSquareMatrix(dictRadialExp.lookup("scatteringMatrix")));
		}
		else
		{
		  DRadialExpList.set(zoneI,new scalarField(DList[zoneI]));
		  nuSigmaEffRadialExpList.set(zoneI,new scalarField(nuSigmaEffList[zoneI]));
		  sigmaPowRadialExpList.set(zoneI,new scalarField(sigmaPowList[zoneI]));
		  sigmaDisappRadialExpList.set(zoneI,new scalarField(sigmaDisappList[zoneI]));
		  sigmaFromToRadialExpList.set(zoneI,new scalarSquareMatrix(sigmaFromToList[zoneI]));
		}
		//AxExp
		read = false;
		perturbedEntry = 0;
		forAll(entriesAxialExp,entryI)
		{
		  if(entriesAxialExp[entryI].keyword() == name)
		  {
		    read = true;
		    perturbedEntry = entryI;
		  }
		}
		if(read)
		{
		  dictionary& dictAxialExp = entriesAxialExp[perturbedEntry].dict();
		  DAxialExpList.set(zoneI,new scalarField("D",dictAxialExp,energyGroups));
		  nuSigmaEffAxialExpList.set(zoneI,new scalarField("nuSigmaEff",dictAxialExp,energyGroups));
		  sigmaPowAxialExpList.set(zoneI,new scalarField("sigmaPow",dictAxialExp,energyGroups));
		  sigmaDisappAxialExpList.set(zoneI,new scalarField("sigmaDisapp",dictAxialExp,energyGroups));
		  sigmaFromToAxialExpList.set(zoneI,new scalarSquareMatrix(dictAxialExp.lookup("scatteringMatrix")));
		}
		else
		{
		  DAxialExpList.set(zoneI,new scalarField(DList[zoneI]));
		  nuSigmaEffAxialExpList.set(zoneI,new scalarField(nuSigmaEffList[zoneI]));
		  sigmaPowAxialExpList.set(zoneI,new scalarField(sigmaPowList[zoneI]));
		  sigmaDisappAxialExpList.set(zoneI,new scalarField(sigmaDisappList[zoneI]));
		  sigmaFromToAxialExpList.set(zoneI,new scalarSquareMatrix(sigmaFromToList[zoneI]));
		}
		//fuelTemp
		read = false;
		perturbedEntry = 0;
		forAll(entriesFuelTemp,entryI)
		{
		  if(entriesFuelTemp[entryI].keyword() == name)
		  {
		    read = true;
		    perturbedEntry = entryI;
		  }
		}
		if(read)
		{
		  dictionary& dictFuelTemp = entriesFuelTemp[perturbedEntry].dict();
		  DFuelTempList.set(zoneI,new scalarField("D",dictFuelTemp,energyGroups));
		  nuSigmaEffFuelTempList.set(zoneI,new scalarField("nuSigmaEff",dictFuelTemp,energyGroups));
		  sigmaPowFuelTempList.set(zoneI,new scalarField("sigmaPow",dictFuelTemp,energyGroups));
		  sigmaDisappFuelTempList.set(zoneI,new scalarField("sigmaDisapp",dictFuelTemp,energyGroups));
		  sigmaFromToFuelTempList.set(zoneI,new scalarSquareMatrix(dictFuelTemp.lookup("scatteringMatrix")));
		}
		else
		{
		  DFuelTempList.set(zoneI,new scalarField(DList[zoneI]));
		  nuSigmaEffFuelTempList.set(zoneI,new scalarField(nuSigmaEffList[zoneI]));
		  sigmaPowFuelTempList.set(zoneI,new scalarField(sigmaPowList[zoneI]));
		  sigmaDisappFuelTempList.set(zoneI,new scalarField(sigmaDisappList[zoneI]));
		  sigmaFromToFuelTempList.set(zoneI,new scalarSquareMatrix(sigmaFromToList[zoneI]));
		}
		//RhoCool
		read = false;
		perturbedEntry = 0;
		forAll(entriesRhoCool,entryI)
		{
		  if(entriesRhoCool[entryI].keyword() == name)
		  {
		    read = true;
		    perturbedEntry = entryI;
		  }
		}
		if(read)
		{
		  dictionary& dictRhoCool = entriesRhoCool[perturbedEntry].dict();
		  DRhoCoolList.set(zoneI,new scalarField("D",dictRhoCool,energyGroups));
		  nuSigmaEffRhoCoolList.set(zoneI,new scalarField("nuSigmaEff",dictRhoCool,energyGroups));
		  sigmaPowRhoCoolList.set(zoneI,new scalarField("sigmaPow",dictRhoCool,energyGroups));
		  sigmaDisappRhoCoolList.set(zoneI,new scalarField("sigmaDisapp",dictRhoCool,energyGroups));
		  sigmaFromToRhoCoolList.set(zoneI,new scalarSquareMatrix(dictRhoCool.lookup("scatteringMatrix")));
		}
		else
		{
		  DRhoCoolList.set(zoneI,new scalarField(DList[zoneI]));
		  nuSigmaEffRhoCoolList.set(zoneI,new scalarField(nuSigmaEffList[zoneI]));
		  sigmaPowRhoCoolList.set(zoneI,new scalarField(sigmaPowList[zoneI]));
		  sigmaDisappRhoCoolList.set(zoneI,new scalarField(sigmaDisappList[zoneI]));
		  sigmaFromToRhoCoolList.set(zoneI,new scalarSquareMatrix(sigmaFromToList[zoneI]));
		}
		//cladExp
		read = false;
		perturbedEntry = 0;
		forAll(entriesCladExp,entryI)
		{
		  if(entriesCladExp[entryI].keyword() == name)
		  {
		    read = true;
		    perturbedEntry = entryI;
		  }
		}
		if(read)
		{
		  dictionary& dictCladExp = entriesCladExp[perturbedEntry].dict();
		  DCladExpList.set(zoneI,new scalarField("D",dictCladExp,energyGroups));
		  nuSigmaEffCladExpList.set(zoneI,new scalarField("nuSigmaEff",dictCladExp,energyGroups));
		  sigmaPowCladExpList.set(zoneI,new scalarField("sigmaPow",dictCladExp,energyGroups));
		  sigmaDisappCladExpList.set(zoneI,new scalarField("sigmaDisapp",dictCladExp,energyGroups));
		  sigmaFromToCladExpList.set(zoneI,new scalarSquareMatrix(dictCladExp.lookup("scatteringMatrix")));
		}
		else
		{
		  DCladExpList.set(zoneI,new scalarField(DList[zoneI]));
		  nuSigmaEffCladExpList.set(zoneI,new scalarField(nuSigmaEffList[zoneI]));
		  sigmaPowCladExpList.set(zoneI,new scalarField(sigmaPowList[zoneI]));
		  sigmaDisappCladExpList.set(zoneI,new scalarField(sigmaDisappList[zoneI]));
		  sigmaFromToCladExpList.set(zoneI,new scalarSquareMatrix(sigmaFromToList[zoneI]));
		}

		//Feedbacks Coefficients
		dimensionedScalar TrefRatio("",dimensionSet(0,0,0,0,0,0,0),(TfuelRef/TfuelPerturbed));
		fuelTempDList.set(zoneI,new scalarField((DList[zoneI]-DFuelTempList[zoneI])/log(TrefRatio).value()));
		fuelTempNuSigmaEffList.set(zoneI,new scalarField((nuSigmaEffList[zoneI]-nuSigmaEffFuelTempList[zoneI])/log(TrefRatio).value()));
		fuelTempSigmaPowList.set(zoneI,new scalarField((sigmaPowList[zoneI]-sigmaPowFuelTempList[zoneI])/log(TrefRatio).value()));
		fuelTempSigmaDisappList.set(zoneI,new scalarField((sigmaDisappList[zoneI]-sigmaDisappFuelTempList[zoneI])/log(TrefRatio).value()));
		fuelTempSigmaFromToList.set(zoneI,new scalarSquareMatrix(sigmaFromToList[zoneI]));
		for(label energyJ=0; energyJ<energyGroups ; energyJ++)
		{
			for(label energyI=0; energyI<energyGroups ; energyI++)
			{
				fuelTempSigmaFromToList[zoneI][energyJ][energyI] = (sigmaFromToList[zoneI][energyJ][energyI]-sigmaFromToFuelTempList[zoneI][energyJ][energyI])/log(TrefRatio).value();
			}
		}
		axialExpDList.set(zoneI,new scalarField(-(DList[zoneI]-DAxialExpList[zoneI])/AxExp));
		axialExpNuSigmaEffList.set(zoneI,new scalarField(-(nuSigmaEffList[zoneI]-nuSigmaEffAxialExpList[zoneI])/AxExp));
		axialExpSigmaPowList.set(zoneI,new scalarField(-(sigmaPowList[zoneI]-sigmaPowAxialExpList[zoneI])/AxExp));
		axialExpSigmaDisappList.set(zoneI,new scalarField(-(sigmaDisappList[zoneI]-sigmaDisappAxialExpList[zoneI])/AxExp));
		axialExpSigmaFromToList.set(zoneI,new scalarSquareMatrix(sigmaFromToList[zoneI]));
		for(label energyJ=0; energyJ<energyGroups ; energyJ++)
		{
			for(label energyI=0; energyI<energyGroups ; energyI++)
			{
				axialExpSigmaFromToList[zoneI][energyJ][energyI] = -(sigmaFromToList[zoneI][energyJ][energyI]-sigmaFromToAxialExpList[zoneI][energyJ][energyI])/AxExp;
			}
		}

		radialExpDList.set(zoneI,new scalarField(-(DList[zoneI]-DRadialExpList[zoneI])/RadExp));
		radialExpNuSigmaEffList.set(zoneI,new scalarField(-(nuSigmaEffList[zoneI]-nuSigmaEffRadialExpList[zoneI])/RadExp));
		radialExpSigmaPowList.set(zoneI,new scalarField(-(sigmaPowList[zoneI]-sigmaPowRadialExpList[zoneI])/RadExp));
		radialExpSigmaDisappList.set(zoneI,new scalarField(-(sigmaDisappList[zoneI]-sigmaDisappRadialExpList[zoneI])/RadExp));
		radialExpSigmaFromToList.set(zoneI,new scalarSquareMatrix(sigmaFromToList[zoneI]));
		for(label energyJ=0; energyJ<energyGroups ; energyJ++)
		{
			for(label energyI=0; energyI<energyGroups ; energyI++)
			{
				radialExpSigmaFromToList[zoneI][energyJ][energyI] = -(sigmaFromToList[zoneI][energyJ][energyI]-sigmaFromToRadialExpList[zoneI][energyJ][energyI])/RadExp;
			}
		}

		scalar deltaRho(rhoCoolRef - rhoCoolPerturbed);
		rhoCoolDList.set(zoneI,new scalarField((DList[zoneI]-DRhoCoolList[zoneI])/deltaRho));
		rhoCoolNuSigmaEffList.set(zoneI,new scalarField((nuSigmaEffList[zoneI]-nuSigmaEffRhoCoolList[zoneI])/deltaRho));
		rhoCoolSigmaPowList.set(zoneI,new scalarField((sigmaPowList[zoneI]-sigmaPowRhoCoolList[zoneI])/deltaRho));
		rhoCoolSigmaDisappList.set(zoneI,new scalarField((sigmaDisappList[zoneI]-sigmaDisappRhoCoolList[zoneI])/deltaRho));
		rhoCoolSigmaFromToList.set(zoneI,new scalarSquareMatrix(sigmaFromToList[zoneI]));
		for(label energyJ=0; energyJ<energyGroups ; energyJ++)
		{
			for(label energyI=0; energyI<energyGroups ; energyI++)
			{
				rhoCoolSigmaFromToList[zoneI][energyJ][energyI] = (sigmaFromToList[zoneI][energyJ][energyI]-sigmaFromToRhoCoolList[zoneI][energyJ][energyI])/deltaRho;
			}
		}

		scalar deltaTClad(TcladRef - TcladPerturbed);
		cladExpDList.set(zoneI,new scalarField((DList[zoneI]-DCladExpList[zoneI])/deltaTClad));
		cladExpNuSigmaEffList.set(zoneI,new scalarField((nuSigmaEffList[zoneI]-nuSigmaEffCladExpList[zoneI])/deltaTClad));
		cladExpSigmaPowList.set(zoneI,new scalarField((sigmaPowList[zoneI]-sigmaPowCladExpList[zoneI])/deltaTClad));
		cladExpSigmaDisappList.set(zoneI,new scalarField((sigmaDisappList[zoneI]-sigmaDisappCladExpList[zoneI])/deltaTClad));
		cladExpSigmaFromToList.set(zoneI,new scalarSquareMatrix(sigmaFromToList[zoneI]));
		for(label energyJ=0; energyJ<energyGroups ; energyJ++)
		{
			for(label energyI=0; energyI<energyGroups ; energyI++)
			{
				cladExpSigmaFromToList[zoneI][energyJ][energyI] = (sigmaFromToList[zoneI][energyJ][energyI]-sigmaFromToCladExpList[zoneI][energyJ][energyI])/deltaTClad;
			}
		}
	}


	Info << "Reading CRmove " << endl;




	IOdictionary CRmove
		(
			IOobject
			(
				"CRmove",
				runTime.constant(),
				neutroMesh,
				IOobject::MUST_READ,
				IOobject::NO_WRITE
			)
	);


	PtrList<entry> CRentries(CRmove.lookup("zones"));

	label CRNumber = CRentries.size();
	PtrList<scalar > CRstart(CRNumber);
	PtrList<scalar > CRfinish(CRNumber);
	PtrList<scalar > CRspeed(CRNumber);//positive for extraction
	PtrList<word > CRFollowerName(CRNumber);
	PtrList<scalar > CRinitialPosition(CRNumber);
	PtrList<scalar > CRposition(CRNumber);

	forAll(CRentries,CRzoneI)
	{

		dictionary& CRdict = CRentries[CRzoneI].dict();
		const word& CRname = CRentries[CRzoneI].keyword();
		CRstart.set(CRzoneI,new scalar(CRdict.lookupOrDefault("startTime",0.0)));
		CRfinish.set(CRzoneI,new scalar(CRdict.lookupOrDefault("endTime",0.0)));
		CRspeed.set(CRzoneI,new scalar(CRdict.lookupOrDefault("speed",0.0)));
		CRFollowerName.set(CRzoneI,new word(CRdict.lookup("followerName")));


		label CRzoneId = neutroMesh.cellZones().findZoneID(CRname);

		if(CRzoneId == -1)
		{
		  Info << "Control rod: " << CRname << " does not exists" << endl;
		}

		scalar cellHeight = 100.0 ;
		scalar minCRHeight = cellHeight ;

		volScalarField centersVSF(neutroMesh.C() & vector(0,0,1));

		Field<scalar> centersSF(centersVSF.internalField());

		forAll(neutroMesh.cellZones()[CRzoneId], cellIlocal)
		{

			label cellIglobal = neutroMesh.cellZones()[CRzoneId][cellIlocal];

			minCRHeight = min(cellHeight,centersSF[cellIglobal]) ;

			cellHeight = minCRHeight;

		}
		reduce(minCRHeight, minOp<scalar>());
		CRinitialPosition.set(CRzoneI, new scalar(minCRHeight));
		CRposition.set(CRzoneI,new scalar(minCRHeight));
		Info << "Control rod position " << CRname << ": " << CRposition[CRzoneI] << endl;

	}

	Info << "Finished reading nuclear data " << endl;