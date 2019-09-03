
Info << "Solving fuel" << endl ;

    fuelToFluid.mapTgtToSrc( thermo.T(), plusEqOp<scalar>(), TcoolantFuel.internalField());
	fuelToFluid.mapTgtToSrc( volScalarField(U.component(2)), plusEqOp<scalar>(), coolantFlowFuel.internalField()); //U.component(2) or mag(U)
	fuelToFluid.mapTgtToSrc( volFuelPower, plusEqOp<scalar>(), linearRateFuel.internalField()); //need to divide the pin height
	
	transTime  += 1;
	
//core nodes only for now (only generate the core mesh for fuelMesh)
	

	forAll(assemblies,zoneI)
	{
		const word& name = assemblies[zoneI].keyword();
        Info << "\nCalculating " << name << endl;
		label zoneId = fuelMesh.cellZones().findZoneID(name);

		if(zoneId == -1)
		{
		  Info << "Mesh zone: " << name << " does not exists" << endl;
		}

//TcoolantFuel,coolantFlowFuel,linearRateFuel -> Trans  //mendmax = 60

		scalar height = SMALL;

		List<scalar> avgFlow;
		List<scalar> avgRate;
		List<scalar> avgTemp;

		List<scalar> heightNum;

		List<scalar> cellTemp;
		List<scalar> cellRate;
		List<scalar> cellFlow;

		label meshNum = 20;

		forAll(fuelMesh.cellZones()[zoneId], cellIlocal)
		{
		//	Info<< "size" << fuelMesh.cellZones()[zoneId].size()  << endl;

			label cellIglobal = fuelMesh.cellZones()[zoneId][cellIlocal];

			// Info<< fuelMesh.C()[cellIglobal].component(2)  << endl;

			label i = int(cellIlocal/meshNum);
			label r = cellIlocal%meshNum;

			if (i<1)
			{
				heightNum.append(fuelMesh.C()[cellIglobal].component(2));
				avgTemp.append(TcoolantFuel[cellIglobal]);
				avgFlow.append(coolantFlowFuel[cellIglobal]);
				avgRate.append(linearRateFuel[cellIglobal]);
			} else {
				avgTemp[r] = (avgTemp[r] * i + TcoolantFuel[cellIglobal]) / (i+1);
				avgFlow[r] = (avgFlow[r] * i + coolantFlowFuel[cellIglobal]) / (i+1);
				avgRate[r] = (avgRate[r] * i + linearRateFuel[cellIglobal]) / (i+1);
			} 
		}

		double maxRate = max(avgRate);
		double maxflow = max(avgFlow);
		double maxtemp = max(avgTemp);
		if( sum(avgRate)/avgRate.size() < 1) {maxRate = 1;}
		if( maxflow < 1) {maxflow = 1; avgFlow = 1.283295;}
		

// 2: linear rod power (kW/m)
// 4: coolant flow rate (g/h)
// 9: coolant temperature (C)

//		label transTime = fuelDataList[zoneI].messageToTrans[12][2]+0.06;
//
		for (int j = 0; j<13; ++j)
			for (int i = 0; i<65; ++i){
				fuelDataList[zoneI].messageToTrans[j][i] = 0;
		}		
					
		forAll(avgTemp, tempId)
		{
			fuelDataList[zoneI].messageToTrans[1][tempId+5] = avgRate[tempId]/maxRate;
			// fuelDataList[zoneI].messageToTrans[2][tempId+5] = avgFlow[tempId]/maxflow;
			fuelDataList[zoneI].messageToTrans[3][tempId+5] = (avgTemp[tempId]-273)/(maxtemp-273);
		}
	

		fuelDataList[zoneI].messageToTrans[0][0] = 1;fuelDataList[zoneI].messageToTrans[0][2] = transTime;
		fuelDataList[zoneI].messageToTrans[1][0] = 2;fuelDataList[zoneI].messageToTrans[1][1] = 1;fuelDataList[zoneI].messageToTrans[1][2] = transTime;fuelDataList[zoneI].messageToTrans[1][3] = maxRate*6.8128E-8;
		fuelDataList[zoneI].messageToTrans[2][0] = 4;fuelDataList[zoneI].messageToTrans[2][1] = 0;fuelDataList[zoneI].messageToTrans[2][2] = transTime;fuelDataList[zoneI].messageToTrans[2][3] = maxflow*446170.02;
		fuelDataList[zoneI].messageToTrans[3][0] = 9;fuelDataList[zoneI].messageToTrans[3][1] = 1;fuelDataList[zoneI].messageToTrans[3][2] = transTime;fuelDataList[zoneI].messageToTrans[3][3] = maxtemp-273;
		fuelDataList[zoneI].messageToTrans[4][2] = transTime + 0.01;
		fuelDataList[zoneI].messageToTrans[5][2] = transTime + 0.01;
		fuelDataList[zoneI].macroi = 0;
		fuelDataList[zoneI].zoneNo = zoneI; 

		Info << "message passing to trans:" << endl;
		for (int j = 0; j<5; ++j)
		{
			for(int i = 0; i<65; ++i) 
				Info<< fuelDataList[zoneI].messageToTrans[j][i] << '\t';
				// Info<< ((fuelDataList[zoneI].messageToTrans[j][i] == 0 ) ?  : fuelDataList[zoneI].messageToTrans[j][i] ) << '\t';
			Info << endl;
		}

		Info << "t 	= " << fuelDataList[zoneI].T << "\tt1 	= "  << fuelDataList[zoneI].T1<< "\tt2 	= "  << fuelDataList[zoneI].T2 << endl;  
		icase = 1;
		passtotrans(&(fuelDataList[zoneI]), &icase);

		List<scalar> transFuelCentralTemp;
		List<scalar> transFuelOuterTemp;
		List<scalar> transCladinnerTemp;
		List<scalar> transCladOuterTemp;
		List<scalar> transFuelAvgTemp;
		List<scalar> transCladAvgTemp;

		forAll(avgTemp, tempId)
		{			
			transFuelCentralTemp.append(fuelDataList[zoneI].TEBI[tempId]+273);
			transFuelOuterTemp.append(fuelDataList[zoneI].TEBA[tempId]+273);
			transCladinnerTemp.append(fuelDataList[zoneI].TEHI[tempId]+273);
			transCladOuterTemp.append(fuelDataList[zoneI].TEHA[tempId]+273);
			transFuelAvgTemp.append(fuelDataList[zoneI].AvTempFu[0][tempId]+273);
			transCladAvgTemp.append(fuelDataList[zoneI].AvTempCl[0][tempId]+273);
		}
 		
 		avgTemp.clear();
 		avgFlow.clear();
 		avgRate.clear();

		height = SMALL;
		label h = SMALL;

		forAll(fuelMesh.cellZones()[zoneId], cellIlocal)
		{
			label cellIglobal = fuelMesh.cellZones()[zoneId][cellIlocal];

			h = cellIlocal%meshNum;

			TCentralFuel[cellIglobal] = transFuelCentralTemp[h];
			TOuterFuel[cellIglobal] = transFuelOuterTemp[h];
			TInnerClad[cellIglobal] = transCladinnerTemp[h];
			TOuterClad[cellIglobal] = transCladOuterTemp[h];
			TavfuelFuel[cellIglobal] = transFuelAvgTemp[h];
			TavcladFuel[cellIglobal] = transCladAvgTemp[h];

			// Info<< "graptag h" << h << '\t' << "cellIglocal" << cellIglobal << '\t' << "cellIlocal" << cellIlocal << endl;
		}
		transFuelCentralTemp.clear();
		transFuelOuterTemp.clear();
		transCladinnerTemp.clear();
		transCladOuterTemp.clear();
		transFuelAvgTemp.clear();
		transCladAvgTemp.clear();
	}
 
// mapping back

TinnerFuel *= 0;

// Info<< "TinnerFuel" << TinnerFuel << endl;

fuelToFluid.mapSrcToTgt(TCentralFuel, plusEqOp<scalar>(), TinnerFuel);
// Info<< "TinnerFuel" << TinnerFuel << endl;

TinnerFuel.correctBoundaryConditions();
// Info<< "TinnerFuel" << TinnerFuel << endl;

TouterFuel *= 0;
fuelToFluid.mapSrcToTgt(TOuterFuel, plusEqOp<scalar>(), TouterFuel);
TouterFuel.correctBoundaryConditions();

TinnerClad *= 0;
fuelToFluid.mapSrcToTgt(TInnerClad, plusEqOp<scalar>(), TinnerClad);
TinnerClad.correctBoundaryConditions();

TouterClad *= 0;
fuelToFluid.mapSrcToTgt(TOuterClad, plusEqOp<scalar>(), TouterClad);
TouterClad.correctBoundaryConditions();

TavFuel *= 0;
fuelToFluid.mapSrcToTgt(TavfuelFuel, plusEqOp<scalar>(), TavFuel);
TavFuel.correctBoundaryConditions();

TavClad *= 0;
fuelToFluid.mapSrcToTgt(TavcladFuel, plusEqOp<scalar>(), TavClad);
TavClad.correctBoundaryConditions();
