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
\*---------------------------------------------------------------------------*/

#include "error.H"

#include "oneWayVTK.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(oneWayVTK, 0);

addToRunTimeSelectionTable
(
    dataExchangeModel,
    oneWayVTK,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
oneWayVTK::oneWayVTK
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    dataExchangeModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    filename_(propsDict_.lookup("couplingFilename")),
    relativePath_(propsDict_.lookup("relativePath"))
{
    readDEMtsfromDict(propsDict_);

    // set max nr of particles from dict
    maxNumberOfParticles_ = readScalar(propsDict_.lookup("maxNumberOfParticles"));
    setNumberOfParticles(maxNumberOfParticles_);

    Info << "relativePath_" << relativePath_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

oneWayVTK::~oneWayVTK()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void oneWayVTK::getData
(
    word name,
    word type,
    double ** const& field,
    label step
) const
{
    if (type == "scalar-atom")
    {
        // get path to particle VTK files
        char index[100];
        sprintf(index, filename_.c_str(), step);
        fileName paricleFilePath(particleCloud_.mesh().time().path()/relativePath_/index);
        Info << "opening file: " << paricleFilePath << endl;

        // open file
        ifstream input(paricleFilePath.c_str());

        if (!input.is_open()) cerr << "File not found!, " << paricleFilePath << endl;

        if (name == "radius")
        {
            // read data
            string just_read;
            while (just_read.compare(name) != 0)  input >> just_read; //read until we read "name"
            input >> just_read;                                    // skip text for dataType
            input >> just_read;                                    // skip text for "1"
            input >> just_read;                                    // skip text for "LookUp"
            input >> just_read;                                    // skip text for "default"
            for (int index = 0; index<particleCloud_.numberOfParticles(); ++index)
            {
                input >> field[index][0];
            }
        }
        else
        {
            // read data
            string just_read;
            while (just_read.compare(name) != 0)  input >> just_read; //read until we read "name"
            input >> just_read;                                    // skip text for dataType
            for (int index = 0; index<particleCloud_.numberOfParticles(); ++index)
            {
                input >> field[index][0];
            }
        }

        // close inputStream
        input.close();
    }
    else if (type == "vector-atom")
    {
        // get path to particle VTK files
        char index[100];
        sprintf(index, filename_.c_str(), step);
        fileName paricleFilePath(particleCloud_.mesh().time().path()/relativePath_/index);
        Info << "opening file: " << paricleFilePath << endl;

        // open file
        ifstream input(paricleFilePath.c_str());

        if (!input.is_open()) cerr << "File not found!, " << paricleFilePath << endl;

        // read position data from VTK file
        //NP: special case as position data has no "name" in the vtk file
        if (name == "x")
        {
            int numberOfParticles;  // remove this?
            string just_read;
            while (just_read.compare("POINTS") != 0)  input >> just_read; // read until we read "POINTS"
            input >> numberOfParticles;                             // this is now the number of points in the file
            input >> just_read;                                     // skip text for dataType

            // give nr of particles to cloud
            setNumberOfParticles(numberOfParticles);

            // re-allocate arrays of cloud
            particleCloud_.reAllocArrays();

            for (int index = 0; index<numberOfParticles; ++index)
            {
                input >> field[index][0] >> field[index][1] >> field[index][2];
            }
        }
        else
        {
            string just_read;
            while (just_read.compare(name) != 0)  input >> just_read;  // read until we read "name"
            input >> just_read;                                        // skip text for dataType
            // need to distinguish between file formats from dump custom/vtk and from LPP
            if(just_read.compare("3") == 0)
            {
                input >> just_read;
                input >> just_read;
            }
            for (int index = 0; index<particleCloud_.numberOfParticles(); ++index)
            {
                input >> field[index][0] >> field[index][1] >> field[index][2];
            }
        }

        // close inputStream
        input.close();
    }
    else
    {
        Info << "unknown type in getData!!!" << endl;
    }
}

void oneWayVTK::giveData
(
    word name,
    word type,
    double ** const& field,
    const char* datatype
) const
{
    // do nothing
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
