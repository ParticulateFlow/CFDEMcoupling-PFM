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
    heat transfer from granular conduction according to the
    Effective Medium Theory discussed e.g. by
    Carson et al. Int. J. Heat Mass Transfer 48 (2005) 2150–2158
    optional radiative contribution according to
    G. Breitbach, H. Barthels. Nucl. Technol. 49 (1980) 392–399
    as discussed in
    Qian et al. Int. J. Heat Mass Transf. 127 (2018) 573–584

\*---------------------------------------------------------------------------*/

#ifndef heatTransferInterGrain_H
#define heatTransferInterGrain_H

#include "fvCFD.H"
#include "cfdemCloudEnergy.H"
#include "energyModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class heatTransferInterGrain Declaration
\*---------------------------------------------------------------------------*/

class heatTransferInterGrain
:
    public energyModel
{
protected:

    dictionary propsDict_;

    bool multiTypes_;

    bool verbose_;

    bool implicit_;

    bool calcTotalHeatFlux_;

    bool radiativeHeatTransfer_;

    scalar totalHeatFlux_;

    word partTempName_;

    const word QPartPartName_;

    volScalarField QPartPart_;

    volScalarField partEffThermCondField_;

    volScalarField partThermCondField_;

    volScalarField partThermCapField_;

    volScalarField partThermRadField_;

    const volScalarField& partTempField_;

    const word prescribedVoidfractionFieldName_;

    const volScalarField& prescribedVoidfraction_;

    const word voidfractionFieldName_;

    const volScalarField& voidfraction_;

    scalarList typeCG_;

    const label maxTypeCG_;

    const word partHeatFluxName_;

    scalarList typePartThermCond_;

    const word partThermCondRegName_;

    scalarList typePartThermCap_;

    const word partThermCapRegName_;

    const word partThermRadRegName_;

    const scalar TMin = 10.0;

    const scalar voidfracMax = 0.95;

    const scalar voidfracMin = 0.05;

    const scalar StefanBoltzmannConst_;

    scalarList typePartEmissivity_;

    scalar kMax_;

    void calcPartEffThermCond();

    void calcPartThermCond();

    void calcPartThermRad();

    scalar FE(scalar, scalar, scalar);

    virtual void giveData();

    virtual void heatFlux(label, scalar, scalar, scalar);



public:

    //- Runtime type information
    TypeName("heatTransferInterGrain");

    // Constructors

        //- Construct from components
        heatTransferInterGrain
        (
            const dictionary& dict,
            cfdemCloudEnergy& sm
        );


    // Destructor

        virtual ~heatTransferInterGrain();


    // Member Functions

        void addEnergyContribution(volScalarField&) const {}

        void addEnergyCoefficient(volScalarField&) const {}

        void calcEnergyContribution();

        void postFlow();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

