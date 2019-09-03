
scatteringSourceExtra *= 0.0;
scatteringSourceIntra *= 0.0;
forAll(flux,energyJ)
{
	if(energyJ != energyI)
	{
		scatteringSourceExtra += sigmaFrom[energyJ] * flux[energyJ];
	}
	else
	{
		scatteringSourceIntra += sigmaFrom[energyJ] * flux[energyJ];
	}
}
