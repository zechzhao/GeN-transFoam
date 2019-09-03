

	forAll(entries, zoneI)
	{
		const word& name = entries[zoneI].keyword();

		label zoneId = neutroMesh.cellZones().findZoneID(name);// verify what is returned

		const unallocLabelList& addr = neutroMesh.cellZones()[zoneId]; //same as |^

		IndirectList<scalar>(Beta.internalField(), addr)
			= scalar(BetaList[zoneI][precI]);

		IndirectList<scalar>(lambda.internalField(), addr)
			= scalar(lambdaList[zoneI][precI]);

	}

    Beta.correctBoundaryConditions();// why correct boundaryCondition? what has been corrected?
    lambda.correctBoundaryConditions();
