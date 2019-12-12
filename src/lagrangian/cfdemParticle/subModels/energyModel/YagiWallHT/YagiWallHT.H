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

#ifndef YagiWallHT_H
#define YagiWallHT_H

#include "fvCFD.H"
#include "cfdemCloudEnergy.H"
#include "energyModel.H"
#include "fixedGradientFvPatchFields.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class heatTransferGunn Declaration
\*---------------------------------------------------------------------------*/

class YagiWallHT
:
    public energyModel
{
protected:

    dictionary propsDict_;

    bool interpolation_;

    bool verbose_;

    word QWallFluidName_;

    volScalarField QWallFluid_;

    volScalarField WallTempField_;

    volScalarField ReField_;

    volScalarField PrField_;

    word tempFieldName_;

    volScalarField tempField_;

    word voidfractionFieldName_;

    const volScalarField& voidfraction_;

    scalar maxSource_;

    word velFieldName_;

    const volVectorField& U_;

    word UsFieldName_;

    const volVectorField& Us_;

    word densityFieldName_;

    const volScalarField& rho_;

    word WallTempName_;

    const volScalarField& WallTemp_;

    mutable double **partRe_;

    bool multiphase_;

    word kfFieldName_;

    const volScalarField& kfField_;

    word CpFieldName_;

    const volScalarField& CpField_;

    void allocateMyArrays() const;


public:

    //- Runtime type information
    TypeName("YagiWallHT");

    // Constructors

        //- Construct from components
        YagiWallHT
        (
            const dictionary& dict,
            cfdemCloudEnergy& sm
        );


    // Destructor

        virtual ~YagiWallHT();


    // Member Functions

        void addEnergyContribution(volScalarField&) const;

	void addEnergyCoefficient(volScalarField&) const {}

	void calcEnergyContribution();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
