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
                        M. Efe Kinaci, JKU Linz, Austria

    Description
    Model to communicate the diffusion coefficients that are necessary for particle
    reduction.

\*---------------------------------------------------------------------------*/

#ifndef diffusionCoefficients_H
#define diffusionCoefficients_H

#include "fvCFD.H"
#include "cfdemCloudEnergy.H"
#include "chemistryModel.H"
#include "species.H"

#include "HashPtrTable.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class diffusionCoefficient Declaration
\*---------------------------------------------------------------------------*/

// for future use:
// + communicate every N steps

class diffusionCoefficient
:
    public chemistryModel
{
private:

    dictionary propsDict_;

    bool verbose_;

    bool interpolation_;

    const fvMesh& mesh_;

    // Species dictionary - To be read where the species are from
    dictionary specDict_;

    wordList speciesNames_;

    word tempFieldName_;

    const volScalarField& tempField_;              // ref to gas temperature field

    // gas pressure at particle location
    word pressureFieldName_;

    const volScalarField& P_;

    word partPressureName_;

    mutable double **partPressure_;

    UPtrList<volScalarField> X_;

    wordList diffusantGasNames_;

    mutable List<double**> diffusionCoefficients_;

    UPtrList<volScalarField> Xdiffusant_;

    HashTable<scalar, word> coeffs;

    HashTable<scalar, word> molWeight;

    void createCoeffs();

    void molWeightTable();

    // calculate numerator part molar mass addition
    double calcMolNum(int, int);

    // calculate denominator part diffusion volume equation
    double calcDiffVol(int, int);

    void allocateMyArrays() const;

    bool initialized_;

    void init();

 public:

    //- Runtime type information
    TypeName("diffusionCoefficients");

    // Constructors

        //- Construct from components
        diffusionCoefficient
        (
            const dictionary& dict,
            cfdemCloudEnergy& sm
        );

    // Destructor

        virtual ~diffusionCoefficient();

    // Member Functions
        void execute();

        void reAllocMyArrays() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
