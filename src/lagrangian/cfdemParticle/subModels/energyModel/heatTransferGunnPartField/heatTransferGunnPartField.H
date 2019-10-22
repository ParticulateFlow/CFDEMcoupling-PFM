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
    Correlation for Nusselt number according to
    Gunn, D. J. International Journal of Heat and Mass Transfer 21.4 (1978)

\*---------------------------------------------------------------------------*/

#ifndef heatTransferGunnPartField_H
#define heatTransferGunnPartField_H

#include "fvCFD.H"
#include "cfdemCloudEnergy.H"
#include "heatTransferGunn.H"
#include "fvOptions.H"
#include "scalarList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
class thermCondModel;
/*---------------------------------------------------------------------------*\
                           Class heatTransferGunnPartField Declaration
\*---------------------------------------------------------------------------*/

class heatTransferGunnPartField
:
    public heatTransferGunn
{
private:

    volScalarField partCpField_;

    const volScalarField& partRhoField_;

    scalarList typeCp_;

    mutable double **partCp_;

    dimensionedScalar pTMax_;

    dimensionedScalar pTMin_;

    autoPtr<thermCondModel> thermCondModel_;

    fv::options& fvOptions;

    void allocateMyArrays() const;

    void giveData();

public:

    //- Runtime type information
    TypeName("heatTransferGunnPartField");

    // Constructors

        //- Construct from components
        heatTransferGunnPartField
        (
            const dictionary& dict,
            cfdemCloudEnergy& sm
        );


    // Destructor

        virtual ~heatTransferGunnPartField();


    // Member Functions
        void addEnergyContribution(volScalarField&) const;

        void calcEnergyContribution();

        void postFlow();

        void solve();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
