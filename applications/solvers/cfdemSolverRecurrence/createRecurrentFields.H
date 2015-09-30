
    
    // to be moved into cfdemCloudRec
	
	
	Foam::Info<< "Create recurrent-time\n" << Foam::endl;
    
    // the documentation is in the code ;-)
    /* Time
        (
            const fileName& rootPath,
            const fileName& caseName,
            const word& systemName = "system",
            const word& constantName = "constant",
            const bool enableFunctionObjects = true
        );
    */
    Foam::Time recTime("dataBase", "", "../system", "../constant", false);
    
    if (verbose)
    {
    	// be informative on properties of the "recTime" Time-object
    	Info << "recTime.rootPath() " << recTime.rootPath() << endl;
    	Info << "recTime.caseName() " << recTime.caseName() << endl;
    	Info << "recTime.path() " << recTime.path() << endl;
    	Info << "recTime.timePath() " << recTime.timePath() << endl;
	Info << "recTime.timeName() " << recTime.timeName() << endl;
    }
    
    // all time directories including constant
    instantList timeDirs = recTime.times();
    
    if (verbose)
    {
    	Info << "timeDirs " << timeDirs << endl;
    }
    
    // iterator for timeDirs
    instantList::iterator it;
    
    // create a data structure for the time indices
    //	constant will not be contained
    //	runTimeIndex -> continuousIndex
    HashTable<label,word> timeIndexList(label(timeDirs.size())-1);
    label contTimeIndex(0);
    
    // create a data structure for the time values
    //	constant will not be contained
    //	continuousIndex -> time.value()
    HashTable<label,scalar> timeValueList(label(timeDirs.size())-1);
    
    // fill the data structure for the time indices
    for (it=timeDirs.begin(); it != timeDirs.end(); ++it)
    {
    	// set run-time
    	recTime.setTime(*it, it->value());
    	
    	
    	// skip constant
    	if (recTime.timeName() == "constant")
    	{
        	continue;
        }
        
        
        // insert the time name into the hash-table with a continuous second index
        timeIndexList.insert(recTime.timeName(), contTimeIndex);
        
        
        if (verbose)
    	{
    		Info << "current time " << recTime.timeName() << endl;
    		Info << "insert " << recTime.timeName() << " , " << contTimeIndex << endl;
    	}
        
        
        // insert the time value
        timeValueList.insert(contTimeIndex, recTime.timeOutputValue());
        
        // increment continuousIndex
        contTimeIndex++;
        
        if (verbose)
    	{
        	Info << "contTimeIndex " << contTimeIndex << endl;
        }
    }
    
    if (verbose)
    {
    	Info << endl;
    	Info << "Found " << label(timeDirs.size()) << " time folders" << endl;
    	Info << "Used " << (label(timeDirs.size())-1) << " time folders" << endl;
	    Info << "Found " << label(timeIndexList.size()) << " time steps" << endl;
    }
    
    // check time step of provided data
    scalar dtCur(0.0);
    scalar dtOld(0.0);
    
    if (verbose)
    {
    	Info << "timeValueList : " << timeValueList << endl;
    }
    
    forAll(timeValueList, i)
    {
    	// compute time step
    	if (timeDirs[i].value() == timeDirs.last().value())
    	{
    		if (verbose)
    		{
    			Info << ".. leaving loop at " << timeDirs[i] << endl;
    		}
    		// leave loop
    		break;
    	}
    	
    	if (verbose)
    	{
    		Info << "timeDirs.fcIndex(i)].value(),  timeDirs[i].value() : " 
    			<< timeDirs[timeDirs.fcIndex(i)].value() << "   " << timeDirs[i].value()
    			<< endl;
    	}
    	
    	// the documentation is in the code ;-)
    	//	fcIndex() - return forward circular index, i.e. the next index
        dtCur = timeDirs[timeDirs.fcIndex(i)].value() - timeDirs[i].value();
        
        if (dtOld == 0.0)
        {
        	dtOld = dtCur;
        }
        
        // do not use dtOld != dtCur as test, you will run into numerical trouble
        // e.g. dtOld = 0.1    dtCur = 0.1    delta = 1.41553e-15
        if (abs(dtOld - dtCur) > SMALL)
        {
        	Info << "dtCur, dtOld = " << dtCur << "   " << dtOld << endl;
        	FatalError << "    in setting up data" << nl
				<< "    non-constant time-step of provided simulation data" 
				<< abort(FatalError);
        }
    }
    
    
    // set deltaT
	recTime.setDeltaT(dtCur, false);
	
	if (verbose)
    {
		Info << "Setting deltaT to " << dtCur << endl;
		Info << "Actual recTime.deltaT = " << recTime.deltaTValue() << endl;
		Info << "Actual runTime.deltaT = " << runTime.deltaTValue() << endl;
    }
    
    // allocate space for flow fields
    //	this is not the final stage
    //	see TEqn.H or other corresponding files for the fields that actually needed
    PtrList<volScalarField> alpha1pl(label(timeIndexList.size()));
    PtrList<volScalarField> alpha2pl(label(timeIndexList.size()));
    PtrList<volVectorField> U1pl(label(timeIndexList.size()));
    PtrList<volVectorField> U2pl(label(timeIndexList.size()));
    
    Info << "Reading fields" << endl;
    
    for (it=timeDirs.begin(); it != timeDirs.end(); ++it)
    {
        // set time
        recTime.setTime(*it, it->value());
        
        // skip constant
        if (recTime.timeName() == "constant")
        {
        	continue;
        }
        
        if (verbose)
    	{
        	Info << "Reading at t = " << recTime.timeName() << endl;
        }
        
        alpha1pl.set
        (
            timeIndexList(recTime.timeName()),
            new volScalarField
            (
                IOobject
                (
                    "voidfraction",
                    recTime.timePath(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            )
        );
        
        IOobject headerAlpha2
        (
            "particlefraction",
            recTime.timePath(),
            mesh,
            IOobject::NO_READ
        );
        
        if (headerAlpha2.headerOk())
        {
	        alpha2pl.set
    	    (
    	        timeIndexList(recTime.timeName()),
    	        new volScalarField
    	        (
    	            IOobject
    	            (
    	                "particlefraction",
                    	recTime.timePath(),
    	                mesh,
    	                IOobject::MUST_READ,
    	                IOobject::NO_WRITE
    	            ),
    	            mesh
    	        )
    	    );
    	}
    	else
    	{
    		alpha2pl.set
    	    (
    	        timeIndexList(recTime.timeName()),
    	        new volScalarField
    	        (
    	            IOobject
    	            (
    	                "particlefraction",
                    	recTime.timePath(),
    	                mesh,
    	                IOobject::NO_READ,
    	                IOobject::NO_WRITE
    	            ),
    	            scalar(1.0) - alpha1pl[timeIndexList(recTime.timeName())]
    	        )
    	    );
    	}
        
        U1pl.set
        (
        	timeIndexList(recTime.timeName()),
            new volVectorField
            (
                IOobject
                (
                    "U",
                    recTime.timePath(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            )
        );
        U2pl.set
        (
        	timeIndexList(recTime.timeName()),
            new volVectorField
            (
                IOobject
                (
                    "Us",
                    recTime.timePath(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            )
        );
        
    }
    
    
    // go to the earliest time step
    it=timeDirs.begin();
    recTime.setTime(*it, it->value());
    
    // skip constant
    if (recTime.timeName() == "constant")
    {
    	it++;
    	recTime.setTime(*it, it->value());
    }
    
    
    
    Info<< "\nComputing recurrence matrix\n" << endl;
    
    // create output file
    std::ostringstream str_pid;
    str_pid << pid();
    const std::string my_pid(str_pid.str());
    OFstream matrixFile("recurrenceMatrix."+my_pid);
    
    // create recurrence matrix
    SymmetricSquareMatrix<scalar> recurrenceMatrix(timeIndexList.size(), timeIndexList.size(), scalar(-666));
    scalar maxElemVal(0.0);
    
    // compute recurrence matrix elements
    forAll(timeIndexList, ti)
    {
    	forAll(timeIndexList, tj)
    	{
    		// main diagonal
    		if (ti == tj)
    		{
    			recurrenceMatrix[ti][tj] = 1;
    			continue;
    		}
    		
    		// skip one half of the matrix
    		if (ti > tj)
    		{
    			continue;
    		}
    		
    		// compute elements
    		recurrenceMatrix[ti][tj]
    			= sumSqr(alpha1pl[ti].internalField() - alpha1pl[tj].internalField())
    			+ sum(magSqr(U1pl[ti].internalField() - U1pl[tj].internalField()));
    		
    		recurrenceMatrix[tj][ti] = recurrenceMatrix[ti][tj];
    		
    		if (maxElemVal < recurrenceMatrix[ti][tj])
    		{
    			maxElemVal = recurrenceMatrix[ti][tj];
    		}
    	}
    }
    
    // normalize matrix elements
    forAll(timeIndexList, ti)
    {
    	forAll(timeIndexList, tj)
    	{
    		recurrenceMatrix[ti][tj] /= maxElemVal;
    	}
    }
    
    
    // write recurrence matrix
    matrixFile << recurrenceMatrix;
    
    
    
    // create time-sequence output file
    OFstream timeSequenceFile("timeSequence."+my_pid);
    
    
    
    // create working fields, these will change according to the recurrent rules
    volScalarField alpha1 = alpha1pl[0];
    volScalarField alpha2 = alpha2pl[0];
    
    volVectorField U1 = U1pl[0];
    volVectorField U2 = U2pl[0];
    
    
    // write stuff ?!
    //runTime.write();
    // runTime.write does not work somehow. As of now, I have no idea why.
    alpha2.write();
    U2.write();
    
    
