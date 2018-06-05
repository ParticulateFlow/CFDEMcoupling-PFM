//===============================
// recurrence turbulence
//===============================


    // check both phases for turbulence models
    forAllIter(PtrListDictionary<phaseModel>, fluid.phases(), iter)
    {
        phaseModel& phase = iter();
        
        Info << "Checking phase " << phase.name() << "'s turbulence model: "
            << phase.turbulence().type() << endl;
        
        /*
            Check for laminar turbulence. This works with OpenFOAM-4.0 and OpenFOAM-5.0,
            as the laminar, multi-phase turbulence model is named "laminar" in OF-4.0
            and "Stokes" in OF-5.0
        */
        if (phase.turbulence().type() == "laminar" || phase.turbulence().type() == "Stokes")
        {
            // do nothing
        }
        else if (isA<Foam::recurrenceTurbulenceModel>(phase.turbulence()))
        {
            /*
                create a reference of the type recurrenceTurbulenceModel
                register the recurrence model with the recurrenceTurbulenceModel
            */
            
            // get const-reference to the turbulence model
            const phaseCompressibleTurbulenceModel& turbConstRef = phase.turbulence();
            
            // cast away const-ness, the underlying turbulence model is not a const object, so this is bad but fine
            phaseCompressibleTurbulenceModel& turbRef = const_cast<phaseCompressibleTurbulenceModel&>(turbConstRef);
            
            // cast away the wrapper class, to get a reference to the turbulence models' base class
            PhaseCompressibleTurbulenceModel<phaseModel>& baseTurbRef
            (
                static_cast<PhaseCompressibleTurbulenceModel<phaseModel>&>(turbRef)
            );
            
            // casting down the family tree
            Foam::recurrenceTurbulenceModel& recTurbRef
            (
                dynamic_cast<Foam::recurrenceTurbulenceModel&>(baseTurbRef)
            );
            
            
            // set recurrenceBase pointer
            recTurbRef.setRecurrenceBasePtr(&recurrenceBase);
            
            
            // check model settings
            turbRef.validate();
        }
        else
        {
            /*
                In a recurrence run, we do not compute any turbulence as we do not solve the fluid flow
                At this point, the phase is not laminar (i.e. not using turbulence) or 
                    using recurrenceTurbulence (i.e. taking turbulent quantities from the data base).
                Hence, abort!
            */
            FatalError
                << "Wrong turbulence model type "
                << phase.turbulence().type() << " for phase " << phase.name() << nl << nl
                << "Valid turbulence model types are types derived from recurrenceTurbulenceModel or laminar" << endl
                << exit(FatalError);
        }
    }
