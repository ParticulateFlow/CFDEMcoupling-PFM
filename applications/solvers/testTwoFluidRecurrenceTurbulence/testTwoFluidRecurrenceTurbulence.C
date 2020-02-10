/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    testTwoFluidRecurrenceTurbulence

Description
    A modified variant of the two-fluid, recurrence model A solver
    with the extension of recurrence-based, multi-phase turbulence modelling.
    This application is used to test whether turbulent fields can be provided
    by the recurrence-based turbulence models.
    
    Run this test application in a recurrence case, with turbulence enabled and 
    the necessary turbulent field quantities present in the data base.
    
    Note the initialisation in checkTurbulenceModels.H
    Updating the turbulence model is done by calling phaseX.correctTurbulence()
    in the file readFields.H

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "twoPhaseSystem.H"
#include "phaseCompressibleTurbulenceModel.H"
#include "pimpleControl.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

#include "recBase.H"
#include "recModel.H"

#include "recurrenceTurbulenceModel.H"

/* // uncomment for OpenFOAM-5.0
namespace Foam
{
    tmp<volScalarField> byDt(const volScalarField& vf)
    {
        if (fv::localEulerDdt::enabled(vf.mesh()))
        {
            return fv::localEulerDdt::localRDeltaT(vf.mesh())*vf;
        }
        else
        {
            return vf/vf.mesh().time().deltaT();
        }
    }

    tmp<surfaceScalarField> byDt(const surfaceScalarField& sf)
    {
        if (fv::localEulerDdt::enabled(sf.mesh()))
        {
            return fv::localEulerDdt::localRDeltaTf(sf.mesh())*sf;
        }
        else
        {
            return sf/sf.mesh().time().deltaT();
        }
    }
}
*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createRDeltaT.H" // remove for OpenFOAM-5.0
    #include "createFields.H"
    #include "createFieldRefs.H"
    
    #include "createTransportFields.H"

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }
    
    Switch faceMomentum
    (
        pimple.dict().lookupOrDefault<Switch>("faceMomentum", false)
    );

    
    recBase recurrenceBase(mesh);
    
    #include "checkTurbulenceModels.H"

    #include "pUf/createDDtU.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting recurrence-based time loop\n" << endl;
    
    label recTimeIndex(0);
    scalar recTimeStep_=recurrenceBase.recM().recTimeStep();

    while (runTime.run())
    {
        #include "readTimeControls.H"

        #include "CourantNos.H"
        #include "setDeltaT.H"

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        
        #include "CEqn.H"
        
        
        if ( runTime.timeOutputValue() - (recTimeIndex+1)*recTimeStep_ + 1.0e-5 > 0.0 )
        {
	        Info << "Updating fields at run time " << runTime.timeOutputValue()
	            << " corresponding to recurrence time " << (recTimeIndex+1)*recTimeStep_ << ".\n" << endl;
	        
            recurrenceBase.updateRecFields();
	        #include "readFields.H"
            recTimeIndex++;
        }
        

        runTime.write();
        
        #include "writeCField.H"

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n\n" << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
