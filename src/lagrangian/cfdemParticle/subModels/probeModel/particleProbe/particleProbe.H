/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).

Class
    particleProbe

SourceFiles
    particleProbe.C

\*---------------------------------------------------------------------------*/

#ifndef particleProbe_H
#define particleProbe_H

#include "probeModel.H"
#include "fvCFD.H"
#include "polyMesh.H"
#include "cfdemCloud.H"
#include "OFstream.H"
#include <vector>
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class particleProbe Declaration
\*---------------------------------------------------------------------------*/

class particleProbe
:
    public probeModel
{

private:

    dictionary propsDict_;

    word name_;

    cfdemCloud& particleCloud_;

    const bool verbose_;

    const bool verboseToFile_;

    const int writePrecision_;

    word dirName_;

    int rank_;

    mutable OFstream* sPtr;

    const int printEvery_;

    const bool sampleAll_;

    const bool probeDebug_;

    const bool includePosition_;

    const labelList particleIDsToSample_;

    wordList itemsToSample_;

    List<OFstream*> sPtrList_;

    int itemCounter_;

    int currItemId_;

    int printCounter_;

    bool printNow_;

    std::vector<std::string> probeName_;

    int probeIndex_;

    void setCounter();

public:

    //- Runtime type information
    TypeName("particleProbe");


    // Constructors

        //- Construct from components
        particleProbe
        (
            const dictionary& dict,
            cfdemCloud& sm,
            const word& typeName,
            const char* logFileName
        );

    // Destructor

        ~particleProbe();

    // Member Functions
    void initialize(const word& modelName, const word& logFileName);
    void setOutputFile(const word& logFileName);
    void writeHeader() const;
    void writeProbe(int index, Field<scalar> sValues, Field<vector> vValues);
    bool checkIDForPrint(int) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
