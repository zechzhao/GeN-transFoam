/*
OFstream powerOctave( runTime.path()/"power.m" );
powerOctave << "% time - power " << endl;
powerOctave << "power = [" << endl;
*/
OFstream powerDat( runTime.path()/"power.dat" );
powerDat << "# time - power" << endl;
/*
OFstream fluxOctave( runTime.path()/"flux.m" );
fluxOctave << "% time - flux " << endl;
fluxOctave << "flux = [" << endl;
*/
OFstream fluxDat( runTime.path()/"flux.dat" );
fluxDat << "# time - flux" << endl;
/*
OFstream precsOctave( runTime.path()/"precs.m" );
precsOctave << "% time - precs " << endl;
precsOctave << "precs = [" << endl;
*/
OFstream precsDat( runTime.path()/"precs.dat" );
precsDat << "# time - precs" << endl;
/*
OFstream TmaxOctave( runTime.path()/"Tmax.m" );
TmaxOctave << "% time - Tmax " << endl;
TmaxOctave << "Tmax = [" << endl;
*/
OFstream TmaxDat( runTime.path()/"Tmax.dat" );
TmaxDat << "# time - Tmax" << endl;
/*
OFstream TavOctave( runTime.path()/"Tav.m" );
TavOctave << "% time - Tav " << endl;
TavOctave << "Tav = [" << endl;
*/
OFstream TavDat( runTime.path()/"Tav.dat" );
TavDat << "# time - Tav" << endl;

OFstream TavFuelDat( runTime.path()/"TavFuel.dat" );
TavFuelDat << "# time - TavFuel" << endl;

OFstream TavCladDat( runTime.path()/"TavClad.dat" );
TavCladDat << "# time - TavClad" << endl;

OFstream TMaxFuelDat( runTime.path()/"TMaxFuel.dat" );
TMaxFuelDat << "# time - TMaxFuel" << endl;

OFstream TMaxCladDat( runTime.path()/"TMaxClad.dat" );
TMaxCladDat << "# time - TMaxClad" << endl;

/*
OFstream TminOctave( runTime.path()/"Tmin.m" );
TminOctave << "% time - Tmin " << endl;
TminOctave << "Tmin = [" << endl;
*/
OFstream TminDat( runTime.path()/"Tmin.dat" );
TminDat << "# time - Tmin" << endl;
/*
OFstream totalDisplacedCellsVolumeOctave( runTime.path()/"totalDisplacedCellsVolume.m" );
totalDisplacedCellsVolumeOctave << "% time - totalDisplacedCellsVolume " << endl;
totalDisplacedCellsVolumeOctave << "totalDisplacedCellsVolume = [" << endl;
*/
OFstream totalDisplacedCellsVolumeDat( runTime.path()/"totalDisplacedCellsVolume.dat" );
totalDisplacedCellsVolumeDat << "# time - totalDisplacedCellsVolume" << endl;

OFstream k_eff( runTime.path()/"k_eff.txt" );
