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
    Correlation for wall-to-bed heat transfer according to
    Yagi, S. A.I.Ch.E. Journal 5.1 (1959)

\*---------------------------------------------------------------------------*/

#ifndef wallHeatTransferYagi_H
#define wallHeatTransferYagi_H

#include "fvCFD.H"
#include "cfdemCloudEnergy.H"
#include "energyModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class wallHeatTransferYagi Declaration
\*---------------------------------------------------------------------------*/

class wallHeatTransferYagi
:
    public energyModel
{
protected:

    dictionary propsDict_;

    bool verbose_;

    bool implicit_;

    word QWallFluidName_;

    volScalarField QWallFluid_;

    word QWallFluidCoeffName_;

    volScalarField QWallFluidCoeff_;

    word wallTempName_;

    volScalarField wallTemp_;

    volScalarField dpField_;

    volVectorField GField_;

    volScalarField magGField_;

    volScalarField ReField_;

    volScalarField PrField_;

    word tempFieldName_;

    const volScalarField& tempField_;

    word voidfractionFieldName_;

    const volScalarField& voidfraction_;

    scalar voidfractionMax_;

    scalar maxSource_;                           // max (limited) value of src field

    word velFieldName_;

    const volVectorField& U_;

    word densityFieldName_;

    const volScalarField& rho_;

    const word dpArrayRegName_;

    virtual void heatFlux(label, scalar, scalar, scalar, scalar);

public:

    //- Runtime type information
    TypeName("wallHeatTransferYagi");

    // Constructors

        //- Construct from components
        wallHeatTransferYagi
        (
            const dictionary& dict,
            cfdemCloudEnergy& sm
        );


    // Destructor

        virtual ~wallHeatTransferYagi();


    // Member Functions

    void addEnergyContribution(volScalarField&) const;

    void addEnergyCoefficient(volScalarField&) const;

    void calcEnergyContribution();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
