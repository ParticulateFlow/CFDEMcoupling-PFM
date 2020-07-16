/*---------------------------------------------------------------------------*\
License

    This is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This code is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with this code.  If not, see <http://www.gnu.org/licenses/>.

    Copyright (C) 2015- Thomas Lichtenegger, JKU Linz, Austria

    Description
    Model to initiate multilayer radii from dump files

\*---------------------------------------------------------------------------*/

#ifndef initMultiLayers_H
#define initMultiLayers_H

#include <sstream>
#include <string>

#include "fvCFD.H"
#include "cfdemCloudEnergy.H"
#include "chemistryModel.H"
#include "vectorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class species Declaration
\*---------------------------------------------------------------------------*/

// for future use:
// + communicate every N steps

class initMultiLayers
:
    public chemistryModel
{
private:

    dictionary propsDict_;

    bool interpolation_;

    const fvMesh& mesh_;

    bool verbose_;

    label maxNumLayers_;

    labelList numLayers_;

    labelList partTypes_;

    mutable double** relRadii_;

    std::string filepath_;

    void allocateMyArrays() const;

    bool initialized_;

    void init();

    void readDump(std::string, labelList &, vectorList &, vectorList &);

    label findNearestCellWithValue(label, volScalarField &) const;

    label getListIndex(label) const;

public:


    //- Runtime type information
    TypeName("initMultiLayers");

    // Constructors

        //- Construct from components
        initMultiLayers
        (
            const dictionary& dict,
            cfdemCloudEnergy& sm
        );

    // Destructor

        virtual ~initMultiLayers();

    // Member Functions
        void execute();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
