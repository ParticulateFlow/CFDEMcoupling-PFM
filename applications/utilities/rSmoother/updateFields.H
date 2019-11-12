Info << "averaging fields for time = " << runTime.value() << ", time index = " << index << endl;

recurrenceBase.recM().exportAveragedVolVectorField(U_smooth, "UMean", threshold, index);
U_smooth.write();
