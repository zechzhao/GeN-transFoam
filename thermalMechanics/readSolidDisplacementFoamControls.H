const dictionary& stressControl = mechMesh.solutionDict().subDict("stressAnalysis");

const int nCorrMech = stressControl.lookupOrDefault<int>("nCorrectors", 1);

scalar convergenceTolerance(readScalar(stressControl.lookup("D")));
Switch compactNormalStress(stressControl.lookup("compactNormalStress"));
