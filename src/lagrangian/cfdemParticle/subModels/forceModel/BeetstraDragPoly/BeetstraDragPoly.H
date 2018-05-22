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
    drag law for monodisperse systems according to
    Beetstra et al. AIChE J 53.2 (2007)

SourceFiles
    BeetstraDragPoly.C
\*---------------------------------------------------------------------------*/

#ifndef BeetstraDragPoly_H
#define BeetstraDragPoly_H

#include "BeetstraDrag.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class BeetstraDragPoly Declaration
\*---------------------------------------------------------------------------*/

class BeetstraDragPoly
:
    public BeetstraDrag
{
protected:

    const bool fines_;
    
    scalar dFine_;

    autoPtr<volScalarField> alphaP_;
    
    autoPtr<volScalarField> alphaSt_;

    autoPtr<volScalarField> dSauter_;

    void adaptVoidfraction(double&, label) const;

    scalar effDiameter(double, label, label) const;
    
    scalar meanSauterDiameter(double, label) const;

public:

    //- Runtime type information
    TypeName("BeetstraDragPoly");


    // Constructors

        //- Construct from components
        BeetstraDragPoly
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~BeetstraDragPoly();


    // Member Functions
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
