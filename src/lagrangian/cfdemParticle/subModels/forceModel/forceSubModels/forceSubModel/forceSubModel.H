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
    forceSubModel

SourceFiles
    forceSubModel.C

\*---------------------------------------------------------------------------*/

#ifndef forceSubModel_H
#define forceSubModel_H

#include "fvCFD.H"
#include "cfdemCloud.H"
#include "probeModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

enum {
    SW_TREAT_FORCE_EXPLICIT = 0,
    SW_TREAT_FORCE_DEM,
    SW_IMPL_FORCE_DEM,
    SW_VERBOSE,
    SW_INTERPOLATION,
    SW_FILTERED_DRAG_MODEL,
    SW_PARCEL_SIZE_DEPENDENT_FILTERED_DRAG,
    SW_IMPL_FORCE_DEM_ACCUMULATED,
    SW_SCALAR_VISCOSITY,
    SW_MAX
};

/*---------------------------------------------------------------------------*\
                           Class forceSubModel Declaration
\*---------------------------------------------------------------------------*/

class forceSubModel
{

protected:

    // Protected data
        const dictionary& dict_;

        cfdemCloud& particleCloud_;

        forceModel& forceModel_;

        label nrDefaultSwitches_;                   // nr of switches defined in mother class

        wordList switchesNameList_;                 // names of switches available

        List<Switch> switchesList_;         // switches which are requested from dict

        List<Switch> switches_;             // switch status

        mutable volScalarField nu_;

        mutable volVectorField divTau_;

        mutable volVectorField IBDragPerV_;

        word densityFieldName_;

        const volScalarField& rho_;

public:

    //- Runtime type information
    TypeName("forceSubModel");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            forceSubModel,
            dictionary,
            (
                const dictionary& dict,
                cfdemCloud& sm,
                forceModel& fm
            ),
            (dict,sm,fm)
        );


    // Constructors

        //- Construct from components
        forceSubModel
        (
            const dictionary& dict,
            cfdemCloud& sm,
            forceModel& fm
        );


    // Destructor

        virtual ~forceSubModel();


    // Selector

        static autoPtr<forceSubModel> New
        (
            const dictionary& dict,
            cfdemCloud& sm,
            forceModel& fm,
            word forceType
        );


    // Member Functions
        virtual void partToArray(label, vector&, const vector&, const vector& Ufluid=vector::zero, scalar Cd=scalar(0)) const;

        virtual void explicitCorr(vector&, vector&, scalar&, vector&, const vector&, vector&, const vector&, bool,label index=100) const;

    // Access

        inline bool verbose() const      { return switches_[SW_VERBOSE]; }

        inline bool interpolation() const { return switches_[SW_INTERPOLATION]; }

        inline bool useFilteredDragModel() const { return switches_[SW_FILTERED_DRAG_MODEL]; }

        inline bool useParcelSizeDependentFilteredDrag() const { return switches_[SW_PARCEL_SIZE_DEPENDENT_FILTERED_DRAG]; }

        virtual word myType() const = 0;

        inline forceModel& myForceM() const { return forceModel_; }

        inline const List<Switch>& switches() const { return switches_; }

        inline const wordList& switchesNameList() const { return switchesNameList_; }

        void setSwitchesList(label i, Switch v) { switchesList_[i] = v; }

        void setSwitches(label i, Switch v) { switches_[i] = v; }

        virtual void readSwitches();

        const volScalarField& nuField() const;

        const volScalarField& muField() const;

        const volScalarField& rhoField() const;

        const volVectorField& divTauField(const volVectorField&) const;

        const volVectorField& IBDragPerV(const volVectorField&,const volScalarField&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
