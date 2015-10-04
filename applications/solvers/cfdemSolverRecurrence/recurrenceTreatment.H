// if still in current sequence
        //	advance time = increment virtualTimeIndex
        //	do nothing else
        // else
        //	generate new sequence
        //		random sequence length
        //		reset virtualTimeIndex
        //		jump to beginning of sequence
        
        // if (in sequence)
        if (virtualTimeIndex < sequenceStart + sequenceLength)
        {
        	// do nothing        	
        }
        else
        {
        	
        	timeSequenceFile << "new sequence : ";
        	
        	int startLoop = 0;
        	int endLoop = 0;
        	
        	// search the other half of the recurrence matrix for 
        	// the new starting point of the next sequence
        	if (virtualTimeIndex < recurrenceMatrix.n()/2)
        	{
        		startLoop = recurrenceMatrix.n()/2;
        		endLoop = recurrenceMatrix.n();
        	}
        	else
        	{
        		startLoop = 0;
        		endLoop = recurrenceMatrix.n()/2-1;
        	}
        	
        	timeSequenceFile << "vti: " << virtualTimeIndex << " st " << startLoop << " end " << endLoop;
        	
        	for (label j = startLoop; j < endLoop; j++)
        	{
        		if (recurrenceMatrix[j][virtualTimeIndex] < nextBestMinimum)
        		{
        			secondBestMinimum = nextBestMinimum;
        			nextBestMinimum = recurrenceMatrix[j][virtualTimeIndex];
        			sequenceStart2 = sequenceStart;
        			sequenceStart = j;
        			continue;
        		}
        		if (recurrenceMatrix[j][virtualTimeIndex] < secondBestMinimum)
        		{
        			secondBestMinimum = recurrenceMatrix[j][virtualTimeIndex];
        			sequenceStart2 = j;
        		}
        	}
        	
        	
        	
        	
        	// use next best or second best minimum
        	if (ranGen.scalar01() <= 0.5)
        	{
        		if (verbose)
    			{
        			Info << " ... using second best minimum" << endl;
        		}
        		
        		sequenceStart = sequenceStartOld;
        		
        		timeSequenceFile << " second best minimum : ";
        	}
        	
        	if (sequenceStartOld == sequenceStart)
        	{
        		if (verbose)
    			{
    				Info << " ... using second best minimum to prevent starting again at the same time" << endl;
    			}
        		sequenceStart = sequenceStart2;
        	}
        	
        	sequenceStartOld = sequenceStart;
        	
        	// write to timeSequenceFile
        	timeSequenceFile << sequenceStart << endl;
        	
        	
        	
        	
        	// set new starting point
        	nextBestMinimum = GREAT;
        	secondBestMinimum = GREAT;
        	
        	// generate new sequence length
        	sequenceLength = ranGen.integer(lowerSeqLim, upperSeqLim);
        	
        	// be wordy?
        	if (verbose)
    		{
    			Info << "Found new starting point at i = " << sequenceStart << endl;
    			Info << "new sequence length = " << sequenceLength << endl;
    		}
        	
        	// trim new sequence length
        	if (sequenceStart + sequenceLength > timeIndexList.size()-1)
        	{
        		//sequenceLength = sequenceLength - (sequenceStart + sequenceLength - timeIndexList.size());
        		sequenceLength = timeIndexList.size()-1 - sequenceStart; // simplified version of the above line
        	}
        	
        	// reset virtual time index
        	virtualTimeIndex = sequenceStart;
        	
        } // end else if(in sequence)
