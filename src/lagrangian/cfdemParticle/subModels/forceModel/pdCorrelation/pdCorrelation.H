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


Description
    calculates the correlation between particle diameter and momentum

SourceFiles
    pdCorrelation.C

Contributing Author
    2018    Paul Kieckhefen, TUHH
\*---------------------------------------------------------------------------*/

#ifndef pdCorrelation_H
#define pdCorrelation_H

#include "forceModel.H"
#include "averagingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class pdCorrelation Declaration
\*---------------------------------------------------------------------------*/

class pdCorrelation
:
    public forceModel
{
private:

    dictionary propsDict_;

    mutable double **d_;
    mutable double **p_;
    mutable double **d2_;
    mutable double **pd_;
    mutable double **cg3_;

    mutable volScalarField dField_;
    mutable volVectorField pField_;
    mutable volScalarField d2Field_;
    mutable volVectorField pdField_;

    mutable volScalarField cg3Field_;

    const scalarList typeCG_;
    const scalarList particleDensities_;
    const Switch     constantCG_;
    const Switch     CG_;
    const Switch     runOnWriteOnly_;

    void allocateMyArrays() const;

public:

    //- Runtime type information
    TypeName("pdCorrelation");

    // Constructors

        //- Construct from components
        pdCorrelation
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~pdCorrelation();


    // Member Functions
        void setForce() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
