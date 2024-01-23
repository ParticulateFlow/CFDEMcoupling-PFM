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

#ifndef heatTransferGunn_H
#define heatTransferGunn_H

#include "fvCFD.H"
#include "cfdemCloudEnergy.H"
#include "energyModel.H"
#include "autoPtr.H"
#include "cellSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class heatTransferGunn Declaration
\*---------------------------------------------------------------------------*/

class heatTransferGunn
:
    public energyModel
{
protected:

    dictionary propsDict_;

    bool multiTypes_;

    bool expNusselt_;

    bool interpolation_;

    bool verbose_;

    bool implicit_;

    bool calcTotalHeatFlux_;

    bool initPartTemp_;

    word scaleNuCellsName_;

    autoPtr<cellSet> scaleNuCells_;

    bool allScaleNuCells_;

    bool scaleNu_;

    scalar Tmin_;

    scalar Tmax_;

    scalar totalHeatFlux_;

    scalar NusseltScalingFactor_;

    const word QPartFluidName_;

    volScalarField QPartFluid_;

    const word QPartFluidCoeffName_;

    volScalarField QPartFluidCoeff_;

    mutable volScalarField partTempField_;

    volScalarField partRelTempField_;

    volScalarField ReField_;

    volScalarField NuField_;

    dimensionedScalar partRefTemp_;

    bool calcPartTempField_;

    dimensionedScalar partTempAve_;

    const word tempFieldName_;

    const volScalarField& tempField_;           // ref to temperature field

    const word voidfractionFieldName_;

    const volScalarField& voidfraction_;        // ref to voidfraction field

    scalar maxSource_;                          // max (limited) value of src field

    const word velFieldName_;

    const volVectorField& U_;

    const word densityFieldName_;

    const volScalarField& rho_;

    const word partTempName_;

    const word partHeatFluxName_;

    const word partHeatFluxCoeffRegName_;

    const word partReRegName_;

    const word partNuRegName_;

    mutable scalar scaleDia_;

    scalarList typeCG_;

    const label maxTypeCG_;

    void partTempField();

    scalar Nusselt(scalar, scalar, scalar) const;

    bool scaleNuCell(label) const;

    virtual void giveData();

    virtual void heatFlux(label, scalar, scalar, scalar, scalar cg3 = 1.0);

    virtual void initPartTemp();

public:

    //- Runtime type information
    TypeName("heatTransferGunn");

    // Constructors

        //- Construct from components
        heatTransferGunn
        (
            const dictionary& dict,
            cfdemCloudEnergy& sm
        );


    // Destructor

        virtual ~heatTransferGunn();


    // Member Functions

        void addEnergyContribution(volScalarField&) const;

        void addEnergyCoefficient(volScalarField&) const;

        void calcEnergyContribution();

        void postFlow();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

