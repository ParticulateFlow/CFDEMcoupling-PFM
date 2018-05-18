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
    BeetstraDrag.C
\*---------------------------------------------------------------------------*/

#ifndef BeetstraDrag_H
#define BeetstraDrag_H

#include "forceModel.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class BeetstraDrag Declaration
\*---------------------------------------------------------------------------*/

class BeetstraDrag
:
    public forceModel
{
protected:
    dictionary propsDict_;

    bool multiTypes_;

    word velFieldName_;

    const volVectorField& U_;

    word voidfractionFieldName_;

    const volScalarField& voidfraction_;

    const scalar minVoidfraction_;

    word UsFieldName_;

    const volVectorField& UsField_;

    mutable scalar scaleDia_;

    scalarList typeCG_;

    mutable scalar scaleDrag_;

    mutable scalar rhoP_;

    mutable scalar rho_;

    mutable scalar Lc2_;

    mutable scalar dPrim_;

    mutable scalar nuf_;

    mutable scalar g_;

    mutable scalar k_;

    bool useGC_;

    bool usePC_;

    virtual void adaptVoidfraction(double&, label) const {}

    virtual scalar effDiameter(double d, label cellI, label index) const {return d;}
    
    virtual scalar meanSauterDiameter(double d, label cellI) const {return d;}

    double F(double, double) const;

    double terminalVelocity(double, double, double, double, double, double) const; 

    double a(double) const;

    double h(double) const;



public:

    //- Runtime type information
    TypeName("BeetstraDrag");


    // Constructors

        //- Construct from components
        BeetstraDrag
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~BeetstraDrag();


    // Member Functions
        void setForce() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
