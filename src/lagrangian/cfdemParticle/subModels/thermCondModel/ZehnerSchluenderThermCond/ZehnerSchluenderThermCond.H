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
    thermal conductivity of PARTICLE phase in presence of a fluid according to
    Zehner and Schluender (1970)

Class
    ZehnerSchluenderThermCond

SourceFiles
    ZehnerSchluenderThermCond.C

\*---------------------------------------------------------------------------*/


#ifndef ZehnerSchluenderThermCond_H
#define ZehnerSchluenderThermCond_H

#include "thermCondModel.H"
#include "scalarList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ZehnerSchluenderThermCond Declaration
\*---------------------------------------------------------------------------*/

class ZehnerSchluenderThermCond
:
    public thermCondModel
{

private:

    dictionary propsDict_;

    word partKsFieldName_;

    volScalarField& partKsField_;

    word voidfractionFieldName_;

    const volScalarField& voidfraction_;

    scalarList typeKs_;

    mutable double **partKs_;

    void allocateMyArrays() const;

    void calcPartKsField() const;

public:

    //- Runtime type information
    TypeName("ZehnerSchluenderThermCond");


    // Constructors

        //- Construct from components
        ZehnerSchluenderThermCond
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~ZehnerSchluenderThermCond();


    // Member Functions

        void calcThermCond();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
