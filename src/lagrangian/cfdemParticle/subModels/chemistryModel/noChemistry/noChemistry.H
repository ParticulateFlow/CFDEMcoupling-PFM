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
                        M. Efe Kinaci, JKU Linz, Austria

    Description
    Chemistry turned off

\*---------------------------------------------------------------------------*/

#ifndef noChemistry_H
#define noChemistry_H

#include "fvCFD.H"
#include "cfdemCloudEnergy.H"
#include "chemistryModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class noChemistry Declaration
\*---------------------------------------------------------------------------*/


class noChemistry
:
    public chemistryModel
{
public:

    //- Runtime type information
    TypeName("off");

    // Constructors

        //- Construct from components
        noChemistry
        (
            const dictionary& dict,
            cfdemCloudEnergy& sm
        );


    // Destructor

        ~noChemistry();


    // Member Functions

        void execute();

        tmp <volScalarField> Smi(const label i) const;

        tmp <volScalarField> Sm() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
