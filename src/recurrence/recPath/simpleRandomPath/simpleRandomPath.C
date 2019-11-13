/*---------------------------------------------------------------------------*\
    CFDEMcoupling academic - Open Source CFD-DEM coupling

    Contributing authors:
    Thomas Lichtenegger
    Copyright (C) 2015- Johannes Kepler University, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling academic.

    CFDEMcoupling academic is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CFDEMcoupling academic is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling academic.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "error.H"
#include "simpleRandomPath.H"
#include "Random.H"
#include "recModel.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(simpleRandomPath, 0);

addToRunTimeSelectionTable
(
    recPath,
    simpleRandomPath,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
simpleRandomPath::simpleRandomPath
(
    const dictionary& dict,
    recBase& base
)
:
    recPath(dict, base),
    propsDict_(dict.subDict(typeName + "Props"))
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

simpleRandomPath::~simpleRandomPath()
{}

// * * * * * * * * * * * * * protected Member Functions  * * * * * * * * * * * * //

void simpleRandomPath::computeRecPath()
{
    Info << "\nComputing recurrence path\n" << endl;

    const label seed = 0;
    Random ranGen(seed);

    label virtualTimeIndex = 0;
    label recSteps = 0;
    label seqStart = 0;
    label lowerSeqLim( base_.recM().lowerSeqLim() );
    label upperSeqLim( base_.recM().upperSeqLim() );
#if OPENFOAM_VERSION_MAJOR < 6
    label seqLength = ranGen.integer(lowerSeqLim, upperSeqLim);
#else
    label seqLength = ranGen.sampleAB(lowerSeqLim, upperSeqLim);
#endif

    virtualTimeIndex = seqEnd(seqStart,seqLength);
    labelPair seqStartEnd(seqStart,virtualTimeIndex);
    virtualTimeIndexList_.append(seqStartEnd);
    recSteps += seqLength;

    SymmetricSquareMatrix<scalar>& recurrenceMatrix( base_.recM().recurrenceMatrix() );
    label numRecFields( base_.recM().numRecFields() );

    if(base_.recM().totRecSteps() == 1)
    {
        Info << "\nPrimitive recurrence path with one element.\n" << endl;
        return;
    }

    while(recSteps <= base_.recM().totRecSteps() )
    {
        label startLoop = 0;
        label endLoop = 0;

        // search the other half of the recurrence matrix for
        // the new starting point of the next sequence
        if (virtualTimeIndex < numRecFields/2)
        {
            startLoop = numRecFields/2;
            endLoop = numRecFields-1 - upperSeqLim/2;  // avoid running into end of recurrence database too often
            endLoop--; // start of next sequence one snapshot AFTER minimum position
        }
        else
        {
            startLoop = 0;
            endLoop = numRecFields/2-1;
        }

        scalar nextMinimum(GREAT);
        for (label j = startLoop; j < endLoop; j++)
        {
            if (recurrenceMatrix[j][virtualTimeIndex] < nextMinimum)
            {
                nextMinimum = recurrenceMatrix[j][virtualTimeIndex];
                seqStart = j+1;
                continue;
            }
        }

#if OPENFOAM_VERSION_MAJOR < 6
        seqLength = ranGen.integer(lowerSeqLim, upperSeqLim);
#else
        seqLength = ranGen.sampleAB(lowerSeqLim, upperSeqLim);
#endif
        virtualTimeIndex = seqEnd(seqStart,seqLength);
        labelPair seqStartEnd(seqStart,virtualTimeIndex);
        virtualTimeIndexList_.append(seqStartEnd);
        recSteps += seqLength;
    }

    Info << "\nComputing recurrence path done\n" << endl;

    if (verbose_)
    {
        Info << " virtualTimeIndexList_ : " << virtualTimeIndexList_ << endl;
    }

    computeJumpVector();
}

label simpleRandomPath::seqEnd(label seqStart, label & seqLength)
{
    label numRecFields( base_.recM().numRecFields() );
    if(seqStart+seqLength > numRecFields-1)
    {
        seqLength=numRecFields-1-seqStart;
    }

    return seqStart+seqLength;
}

void simpleRandomPath::computeJumpVector()
{
    Info << "\nComputing recurrence jump vector\n" << endl;

    OFstream jumpvec("rec_jump.dat");
    label numRecFields( base_.recM().numRecFields() );
    label startLoop = 0;
    label endLoop = 0;
    SymmetricSquareMatrix<scalar>& recurrenceMatrix( base_.recM().recurrenceMatrix() );

    for (label i = 0; i < numRecFields; i++)
    {
        if (i < numRecFields/2)
        {
            startLoop = numRecFields/2;
            endLoop = numRecFields-1;
        }
        else
        {
            startLoop = 0;
            endLoop = numRecFields/2-1;
        }

        scalar nextMinimum(GREAT);
        label jumpdest = 0;

        for (label j = startLoop; j < endLoop; j++)
        {
            if (recurrenceMatrix[j][i] < nextMinimum)
            {
                nextMinimum = recurrenceMatrix[j][i];
                jumpdest = j+1;

                continue;
            }
        }

        jumpvec << jumpdest << endl;
    }

    Info << "\nComputing recurrence jump vector done\n" << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
