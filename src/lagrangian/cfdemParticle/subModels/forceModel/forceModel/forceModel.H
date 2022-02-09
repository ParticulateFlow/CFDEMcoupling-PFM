/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).

Class
    forceModel

SourceFiles
    forceModel.C

\*---------------------------------------------------------------------------*/

#ifndef forceModel_H
#define forceModel_H

#include "fvCFD.H"
#include "cfdemCloud.H"
#include "probeModel.H"
#include "forceSubModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class forceModel Declaration
\*---------------------------------------------------------------------------*/

class forceModel
{

protected:

    // Protected data
        const dictionary& dict_;

        cfdemCloud& particleCloud_;

        mutable volVectorField impParticleForces_;  // sum of implicit particle forces [N]

        mutable volVectorField expParticleForces_;  // sum of explicit particle forces [N]

        bool coupleForce_;

        const word modelType_;

        bool probeIt_;

        bool requiresEx_;           //requires a orientation vector

        wordList forceSubModels_;

        mutable PtrList<forceSubModel> forceSubModel_;

public:

    //- Runtime type information
    TypeName("forceModel");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            forceModel,
            dictionary,
            (
                const dictionary& dict,
                cfdemCloud& sm
            ),
            (dict,sm)
        );


    // Constructors

        //- Construct from components
        forceModel
        (
            const dictionary& dict,
            cfdemCloud& sm
        );


    // Destructor

        virtual ~forceModel();


    // Selector

        static autoPtr<forceModel> New
        (
            const dictionary& dict,
            cfdemCloud& sm,
            word forceType
        );


    // Member Functions
        virtual void setForce() const = 0;

        //tmp<volScalarField> provideScalarField();

        virtual void manipulateScalarField(volScalarField&) const;

    // Access
        const word& modelType() const { return modelType_; }

        inline volVectorField& impParticleForces() const { return impParticleForces_; }

        inline volVectorField& expParticleForces() const { return expParticleForces_; }






        inline bool coupleForce() const { return coupleForce_; }

        virtual inline bool& requiresEx() { return requiresEx_; }

        void repartitionImExForces() const; //Repartition Implixit/Explicit forces

        void treatVoidCells() const;

        inline const wordList& forceSubModels() { return forceSubModels_; }

        inline forceSubModel& forceSubM(int i) const { return forceSubModel_[i]; }

        inline int nrForceSubModels() const { return forceSubModels_.size(); }

        void setForceSubModels(dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
