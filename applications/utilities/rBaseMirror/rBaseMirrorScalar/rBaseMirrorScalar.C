/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    rBaseMirror


Description
    Read time series and extend it by mirrored fields if geometry possesses
    the same symmetry

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    timeSelector::addOptions();

    #include "setRootCase.H"
    #include "createTime.H"

   // read in start and end time from controlDict

    scalar startTime=runTime.startTime().value();
    scalar endTime=runTime.endTime().value();
    scalar origTimeRange = endTime - startTime;

    Info << "start time = " << runTime.startTime() << endl;
    Info << "end time = " << runTime.endTime() << endl;

    // check which time directories are present
    // instantList timeDirs = timeSelector::select0(runTime, args);
    // runTime.setTime(timeDirs[0], 0);

    #include "createMesh.H"

    #include "createFields.H"

    Foam::Time recTime(fileName(dataBaseName), "", "../system", "../constant", false);
    instantList timeDirs(recTime.times());
    recTime.setTime(timeDirs[0],0);

    #include "readFields.H"

    Info << fieldName << endl;

    volScalarField transformedField = origField;

    scalar t;

    label shiftedTimeI = 0;

    // check number of time directories
    label shift = 0;
    forAll(timeDirs, timeI)
    {

       if (recTime.timeName() == "constant") continue;
       recTime.setTime(timeDirs[timeI], timeI);
       t = recTime.value();
       if(t < startTime) continue;
       if(t > endTime) continue;
       shift++;
    }

    scalar dt = origTimeRange / (shift - 1.0);
    recTime.setEndTime(startTime + 2 * origTimeRange + dt);

    label cellI_transformed = -1;
    forAll(timeDirs, timeI)
    {

       recTime.setTime(timeDirs[timeI], timeI);
       t = recTime.value();
       if(t < startTime) continue;
       if(t > endTime) continue;
       Info << "time = " << t << ", time index = " << timeI << endl;

       #include "readFields.H"

       forAll(transformedField, cellI)
       {
           vector position = mesh.C()[cellI];
           vector transformedPosition = 2 * ((refPoint - position) & refDirection) * refDirection / (refDirection & refDirection) + position;
           cellI_transformed = mesh.findCell(transformedPosition);
           if(cellI_transformed < 0)
           {
               Info << "Couldn't find transformed cell. Stopping." << endl;
               return 0;
           }

           scalar value = origField[cellI_transformed];
           scalar transformedValue = value;

           transformedField[cellI] = transformedValue;
       }

       shiftedTimeI = timeI + shift;
       t = recTime.value() + origTimeRange + dt;
       runTime.setTime(t, shiftedTimeI);
       Info << "creating transformed fields for time = " << t << ", time index = " << shiftedTimeI << endl;
       transformedField.write();
    }

    Info << "\nEnd" << endl;

    return 0;
}

// ************************************************************************* //
