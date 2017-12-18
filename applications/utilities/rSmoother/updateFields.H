Info << "averaging fields for time = " << runTime.value() <<", time index = " << index << endl;

U_smooth = recurrenceBase.recM().exportAveragedVolVectorField("UMean", threshold, index)();
U_smooth.write();
