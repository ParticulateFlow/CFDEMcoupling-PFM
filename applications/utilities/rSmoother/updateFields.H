Info << "averaging fields for time = " << runTime.value() <<", time index = " << index << endl;

U_smooth = recurrenceBase.recM().exportAveragedVolVectorField("U", threshold, index);
U_smooth.write();
