//one assembly in a zone 
   	Info<< "\nReading assembly data\n" << endl;

	IOdictionary assemblyData
	(
		IOobject
		(
			"assemblyData",
			runTime.constant(),
			fuelMesh,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		)
	);

    PtrList<entry> assemblies(assemblyData.lookup("zones"));
    label assemblySize = assemblies.size();
	
	PtrList<transtruct> fuelDataList;

	double transTime = 9840.1;

	for(label i = 0; i<assemblySize; ++i)
	{
		transtruct * newtranstruct = new transtruct;
		fuelDataList.append(newtranstruct);
	}

	int icase;

	forAll(fuelDataList,zoneI)
	{
	transinitialiser(&(fuelDataList[zoneI]));


// 2: linear rod power (kW/m)
// 4: coolant flow rate (g/h)
// 9: coolant temperature (C)
// 0 iwert; 1 iaxvar; 2 zeit; 3 wert; 4 dwert; 5-65  axvar

    // Info<< fuelDataList(zoneI)->start_loca <<endl;
    fuelDataList[zoneI].messageToTrans[0][0] = 1;fuelDataList[zoneI].messageToTrans[0][1] = 1;
    fuelDataList[zoneI].messageToTrans[1][0] = 2; fuelDataList[zoneI].messageToTrans[1][1] = 1; fuelDataList[zoneI].messageToTrans[1][4] = 0;
    fuelDataList[zoneI].messageToTrans[1][5] = 0.718;fuelDataList[zoneI].messageToTrans[1][6] = 0.80;fuelDataList[zoneI].messageToTrans[1][7] = 0.85;fuelDataList[zoneI].messageToTrans[1][8] = 0.90;fuelDataList[zoneI].messageToTrans[1][9] = 0.936;
    fuelDataList[zoneI].messageToTrans[1][10] = 0.96;fuelDataList[zoneI].messageToTrans[1][11] = 0.977;fuelDataList[zoneI].messageToTrans[1][12] = 0.989;fuelDataList[zoneI].messageToTrans[1][13] = 0.996;fuelDataList[zoneI].messageToTrans[1][14] = 1;
    fuelDataList[zoneI].messageToTrans[1][15] = 1;fuelDataList[zoneI].messageToTrans[1][16] = 0.999;fuelDataList[zoneI].messageToTrans[1][17] = 0.995;fuelDataList[zoneI].messageToTrans[1][18] = 0.987;fuelDataList[zoneI].messageToTrans[1][19] = 0.974;
    fuelDataList[zoneI].messageToTrans[1][20] = 0.954;fuelDataList[zoneI].messageToTrans[1][21] = 0.925;fuelDataList[zoneI].messageToTrans[1][22] = 0.8865;fuelDataList[zoneI].messageToTrans[1][23] = 0.8332;fuelDataList[zoneI].messageToTrans[1][24] = 0.76;
    fuelDataList[zoneI].messageToTrans[2][0] = 4; fuelDataList[zoneI].messageToTrans[2][3] = 0.5728e6;
    fuelDataList[zoneI].messageToTrans[3][0] = 9; fuelDataList[zoneI].messageToTrans[3][3] = 0.395e3;
    fuelDataList[zoneI].messageToTrans[4][0] = 30;
    fuelDataList[zoneI].messageToTrans[5][0] = 1; fuelDataList[zoneI].messageToTrans[5][2] = 2.7;
    fuelDataList[zoneI].messageToTrans[6][0] = 2; fuelDataList[zoneI].messageToTrans[6][2] = 2.7;fuelDataList[zoneI].messageToTrans[6][3] = 5;
    fuelDataList[zoneI].messageToTrans[7][0] = 4; fuelDataList[zoneI].messageToTrans[7][2] = 2.7;fuelDataList[zoneI].messageToTrans[7][3] = 0.58e6;
    fuelDataList[zoneI].messageToTrans[8][0] = 9; fuelDataList[zoneI].messageToTrans[8][2] = 2.7;fuelDataList[zoneI].messageToTrans[8][3] = 20;
    fuelDataList[zoneI].messageToTrans[9][2] = 2.71;
    fuelDataList[zoneI].messageToTrans[10][2] = 2.71;
    fuelDataList[zoneI].macroi = 0;
    fuelDataList[zoneI].zoneNo = zoneI;

    icase = 0;

    passtotrans (&(fuelDataList[zoneI]),&icase);

    cout << "\n initialising for zone number " << zoneI << '\t'  << assemblies[zoneI].keyword() << "\n transTime = 2.7h " << endl;
	

    for (int j = 0; j<13; ++j)
            for (int i = 0; i<65; ++i){
                fuelDataList[zoneI].messageToTrans[j][i] = 0;
        }

    fuelDataList[zoneI].messageToTrans[0][0] = 1;fuelDataList[zoneI].messageToTrans[0][2] = 30;
    fuelDataList[zoneI].messageToTrans[1][0] = 2;fuelDataList[zoneI].messageToTrans[1][1] = 1;fuelDataList[zoneI].messageToTrans[1][2] = 30;fuelDataList[zoneI].messageToTrans[1][3] = 0.2932780000000e2;
    fuelDataList[zoneI].messageToTrans[1][5] = 0.718;fuelDataList[zoneI].messageToTrans[1][6] = 0.80;fuelDataList[zoneI].messageToTrans[1][7] = 0.85;fuelDataList[zoneI].messageToTrans[1][8] = 0.90;fuelDataList[zoneI].messageToTrans[1][9] = 0.936;
    fuelDataList[zoneI].messageToTrans[1][10] = 0.96;fuelDataList[zoneI].messageToTrans[1][11] = 0.977;fuelDataList[zoneI].messageToTrans[1][12] = 0.989;fuelDataList[zoneI].messageToTrans[1][13] = 0.996;fuelDataList[zoneI].messageToTrans[1][14] = 1;
    fuelDataList[zoneI].messageToTrans[1][15] = 1;fuelDataList[zoneI].messageToTrans[1][16] = 0.999;fuelDataList[zoneI].messageToTrans[1][17] = 0.995;fuelDataList[zoneI].messageToTrans[1][18] = 0.987;fuelDataList[zoneI].messageToTrans[1][19] = 0.974;
    fuelDataList[zoneI].messageToTrans[1][20] = 0.954;fuelDataList[zoneI].messageToTrans[1][21] = 0.925;fuelDataList[zoneI].messageToTrans[1][22] = 0.8865;fuelDataList[zoneI].messageToTrans[1][23] = 0.8332;fuelDataList[zoneI].messageToTrans[1][24] = 0.76;

    fuelDataList[zoneI].messageToTrans[2][0] = 3;fuelDataList[zoneI].messageToTrans[2][1] = 0;fuelDataList[zoneI].messageToTrans[2][2] = 30;fuelDataList[zoneI].messageToTrans[2][3] = 0.531e15;
    fuelDataList[zoneI].messageToTrans[3][0] = 9;fuelDataList[zoneI].messageToTrans[3][1] = 0;fuelDataList[zoneI].messageToTrans[3][2] = 30;fuelDataList[zoneI].messageToTrans[3][3] = 395;
    fuelDataList[zoneI].messageToTrans[4][2] = 30 + 0.1;
    fuelDataList[zoneI].messageToTrans[5][2] = 30 + 0.1;
    fuelDataList[zoneI].macroi = 0;
    fuelDataList[zoneI].zoneNo = zoneI; 

    icase = 1;

    passtotrans (&(fuelDataList[zoneI]),&icase);

    cout << "\n transTime = 30h " << endl;

    }

    
