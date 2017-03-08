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
class thermCondModel;
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
    
    bool implicitEnergyModel_;
    
    const wordList chemistryModels_;

    autoPtr<energyModel>* energyModel_;
    
    autoPtr<thermCondModel> thermCondModel_;
    
    autoPtr<chemistryModel> chemistryModel_;
    
    void calcEnergyContributions();

    void speciesExecute();
    
public:

    friend class energyModel;

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
       
       const thermCondModel& thermCondM();
       
       const chemistryModel& chemistryM();

       int nrEnergyModels();

       inline const wordList& energyModels();
       
       bool& implicitEnergyModel();
       
       void energyContributions(volScalarField&);
       
       void energyCoefficients(volScalarField&);
       
       bool evolve(volScalarField&,volVectorField&,volVectorField&);
       
       void postFlow();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cfdemCloudEnergyI.H"

#endif

// ************************************************************************* //
