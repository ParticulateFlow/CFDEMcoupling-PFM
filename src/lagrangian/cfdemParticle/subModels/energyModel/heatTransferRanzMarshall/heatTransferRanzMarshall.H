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
    Ranz and Marshall, Chem. Eng. Prog. 48.3 (1952)

\*---------------------------------------------------------------------------*/

#ifndef heatTransferRanzMarshall_H
#define heatTransferRanzMarshall_H

#include "fvCFD.H"
#include "cfdemCloudEnergy.H"
#include "energyModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class heatTransferRanzMarshall Declaration
\*---------------------------------------------------------------------------*/

class heatTransferRanzMarshall
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

    scalar Tmin_;

    scalar Tmax_;

    scalar totalHeatFlux_;

    scalar NusseltScalingFactor_;

    word QPartFluidName_;

    volScalarField QPartFluid_;

    word QPartFluidCoeffName_;

    volScalarField QPartFluidCoeff_;

    mutable volScalarField partTempField_;

    volScalarField partRelTempField_;

    volScalarField ReField_;

    volScalarField NuField_;

    dimensionedScalar partRefTemp_;

    bool calcPartTempField_;

    dimensionedScalar partTempAve_;

    word tempFieldName_;

    const volScalarField& tempField_;            // ref to temperature field

    word voidfractionFieldName_;

    const volScalarField& voidfraction_;    // ref to voidfraction field

    scalar maxSource_;                           // max (limited) value of src field

    word velFieldName_;

    const volVectorField& U_;

    word densityFieldName_;

    const volScalarField& rho_;

    word partTempName_;

    mutable double **partTemp_;

    word partHeatFluxName_;

    mutable double **partHeatFlux_;

    mutable double **partHeatFluxCoeff_;

    mutable double **partRe_;

    mutable double **partNu_;

    mutable scalar scaleDia_;

    scalarList typeCG_;

    const label maxTypeCG_;

    void allocateMyArrays() const;

    void partTempField();

    scalar Nusselt(scalar, scalar, scalar) const;

    virtual void giveData();

    virtual void heatFlux(label, scalar, scalar, scalar, scalar cg3 = 1.0);

    virtual void initPartTemp();

public:

    //- Runtime type information
    TypeName("heatTransferRanzMarshall");

    // Constructors

        //- Construct from components
        heatTransferRanzMarshall
        (
            const dictionary& dict,
            cfdemCloudEnergy& sm
        );


    // Destructor

        virtual ~heatTransferRanzMarshall();


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

