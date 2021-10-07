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
    cloud class managing DEM data for CFD-DEM coupling with energy models,
    e.g. for compressible, reacting flows

Class
    Foam::cfdemCloudEnergy derived from cfdemCloud

SourceFiles
    cfdemCloudEnergy.C

\*---------------------------------------------------------------------------*/

#ifndef cfdemCloudEnergy_H
#define cfdemCloudEnergy_H

#include "cfdemCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// forward declarations
class energyModel;
class massTransferModel;
class thermCondModel;
class diffCoeffModel;  
class chemistryModel;
/*---------------------------------------------------------------------------*\
                           Class cfdemCloudEnergy Declaration
\*---------------------------------------------------------------------------*/

class cfdemCloudEnergy
:
    public cfdemCloud
{
protected:

    const wordList energyModels_;

    const wordList massTransferModels_;

    bool implicitEnergyModel_;

    bool implicitMassTransferModel_;

    const wordList chemistryModels_;

    PtrList<energyModel> energyModel_;

    PtrList<massTransferModel> massTransferModel_;

    autoPtr<thermCondModel> thermCondModel_;

    autoPtr<diffCoeffModel> diffCoeffModel_;

    PtrList<chemistryModel> chemistryModel_;

    void calcEnergyContributions();

    void calcMassContributions();

    void speciesExecute();

public:

    friend class energyModel;

    friend class massTransferModel;

    friend class chemistryModel;

    // Constructors

        //- Construct from components
        cfdemCloudEnergy
        (
            const fvMesh& mesh
        );

    // Destructor

        ~cfdemCloudEnergy();


    // Member Functions


       const energyModel& energyM(int);

       const massTransferModel& massTransferM(int);

       const thermCondModel& thermCondM();

       const diffCoeffModel& diffCoeffM();

       const chemistryModel& chemistryM(int);

       label nrEnergyModels() const;

       label nrMassTransferModels() const;

       bool implicitEnergyModel() const;

       bool implicitMassTransferModel() const;

       int nrChemistryModels();

       inline const wordList& energyModels() const;

       inline const wordList& massTransferModels() const;

       inline const wordList& chemistryModels() const;

       void energyContributions(volScalarField&);

       void energyCoefficients(volScalarField&);

       void massContributions(volScalarField&);

       void massCoefficients(volScalarField&);

       bool evolve(volScalarField&,volVectorField&,volVectorField&);

       void postFlow();

       void solve();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cfdemCloudEnergyI.H"

#endif

// ************************************************************************* //
