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
    Correlation for Sherwood number according to
    Gunn, D. J. International Journal of Heat and Mass Transfer 21.4 (1978)

\*---------------------------------------------------------------------------*/

#ifndef massTransferGunn_H
#define massTransferGunn_H

#include "fvCFD.H"
#include "cfdemCloudEnergy.H"
#include "massTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class massTransferGunn Declaration
\*---------------------------------------------------------------------------*/

class massTransferGunn
:
    public massTransferModel
{
protected:

    dictionary propsDict_;

    bool multiTypes_;

    bool expSherwood_;

    bool interpolation_;

    bool verbose_;

    bool implicit_;

    word SPartFluidName_;

    volScalarField SPartFluid_;

    word SPartFluidCoeffName_;

    volScalarField SPartFluidCoeff_;

    //mutable volScalarField partTempField_;

    //volScalarField partRelTempField_;

    volScalarField ReField_;

    volScalarField ShField_;

    //dimensionedScalar partRefTemp_;

    //bool calcPartTempField_;

    //bool calcPartTempAve_;

    //scalar partTempAve_;

    word concFieldName_;

    const volScalarField& concField_;           // ref to concentration field

	word satConcFieldName_;

    const volScalarField& satConcField_;        // ref to saturation concentration field

    word voidfractionFieldName_;

    const volScalarField& voidfraction_;        // ref to voidfraction field

    scalar maxSource_;                          // max (limited) value of src field

    word velFieldName_;

    const volVectorField& U_;

    word densityFieldName_;

    const volScalarField& rho_;

    //word partTempName_;

    //mutable double **partTemp_;

	bool coupleDEM_;

    word partMassFluxName_;

    mutable double **partMassFlux_;

    mutable double **partMassFluxCoeff_;

    mutable double **partRe_;

    mutable double **partSh_;

    mutable scalar scaleDia_;

    scalarList typeCG_;

    void allocateMyArrays() const;

    //void partTempField();

    scalar Sherwood(scalar, scalar, scalar) const;

    virtual void giveData();

    virtual void massFlux(label, scalar, scalar, scalar, scalar, scalar cg3 = 1.0);

public:

    //- Runtime type information
    TypeName("massTransferGunn");

    // Constructors

        //- Construct from components
        massTransferGunn
        (
            const dictionary& dict,
            cfdemCloudEnergy& sm
        );


    // Destructor

        virtual ~massTransferGunn();


    // Member Functions

        void addMassContribution(volScalarField&) const;

        void addMassCoefficient(volScalarField&) const;

        void calcMassContribution();

        void postFlow();

        //scalar aveTpart() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

