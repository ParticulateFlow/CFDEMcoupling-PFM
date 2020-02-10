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
    Carson et al. Int. J. Heat Mass Transfer 48 (2005)

\*---------------------------------------------------------------------------*/

#ifndef heatTransferGranConduction_H
#define heatTransferGranConduction_H

#include "fvCFD.H"
#include "cfdemCloudEnergy.H"
#include "energyModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class heatTransferGranConduction Declaration
\*---------------------------------------------------------------------------*/

class heatTransferGranConduction
:
    public energyModel
{
protected:

    dictionary propsDict_;

    bool multiTypes_;

    bool verbose_;

    bool calcTotalHeatFlux_;

    scalar totalHeatFlux_;

    word QPartPartName_;

    volScalarField QPartPart_;

    volScalarField partEffThermCondField_;

    volScalarField partThermCondField_;

    const volScalarField& partTempField_;

    word prescribedVoidfractionFieldName_;

    const volScalarField& prescribedVoidfraction_;

    word voidfractionFieldName_;

    const volScalarField& voidfraction_;

    word partHeatFluxName_;

    mutable double **partHeatFlux_;

    scalarList typePartThermCond_;

    mutable double **partThermCond_;

    void allocateMyArrays() const;

    void calcPartEffThermCond();

    void calcPartThermCond();

    virtual void giveData();

    virtual void heatFlux(label, scalar, scalar, scalar);



public:

    //- Runtime type information
    TypeName("heatTransferGranConduction");

    // Constructors

        //- Construct from components
        heatTransferGranConduction
        (
            const dictionary& dict,
            cfdemCloudEnergy& sm
        );


    // Destructor

        virtual ~heatTransferGranConduction();


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

