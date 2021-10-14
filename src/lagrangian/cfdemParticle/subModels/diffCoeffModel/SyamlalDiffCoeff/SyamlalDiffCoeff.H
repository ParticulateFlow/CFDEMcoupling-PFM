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
    diffusion coefficient of species in fluid phase in presence of particles
    according to Syamlal, M. and Gidaspow, M. AIChE Journal 31.1 (1985)

Class
    SyamlalDiffCoeff

SourceFiles
    SyamlalDiffCoeff.C

\*---------------------------------------------------------------------------*/


#ifndef SyamlalDiffCoeff_H
#define SyamlalDiffCoeff_H

#include "diffCoeffModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class SyamlalDiffCoeff Declaration
\*---------------------------------------------------------------------------*/

class SyamlalDiffCoeff
:
    public diffCoeffModel
{

private:

    dictionary propsDict_;

    word voidfractionFieldName_;

    const volScalarField& voidfraction_;

public:

    //- Runtime type information
    TypeName("SyamlalDiffCoeff");


    // Constructors

        //- Construct from components
        SyamlalDiffCoeff
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~SyamlalDiffCoeff();


    // Member Functions

        tmp<volScalarField> diffCoeff() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
