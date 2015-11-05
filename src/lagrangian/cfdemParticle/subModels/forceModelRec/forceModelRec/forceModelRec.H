/*---------------------------------------------------------------------------*\
    CFDEMcoupling academic - Open Source CFD-DEM coupling
    
    Contributing authors:
    Thomas Lichtenegger, Gerhard Holzinger
    Copyright (C) 2015- Johannes Kepler University, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling academic.

    CFDEMcoupling academic is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CFDEMcoupling academic is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling academic.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#ifndef forceModelRec_H
#define forceModelRec_H

#include "cfdemCloudRec.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class forceModelRec Declaration
\*---------------------------------------------------------------------------*/

class forceModelRec
{

protected:

    // Protected data
        const dictionary& dict_;

        cfdemCloudRec& particleCloud_;

        void coupleRecForce();
public:

    //- Runtime type information
    TypeName("forceModelRec");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            forceModelRec,
            dictionary,
            (
                const dictionary& dict,
                cfdemCloudRec& sm
            ),
            (dict,sm)
        );


    // Constructors

        //- Construct from components
        forceModelRec
        (
            const dictionary& dict,
            cfdemCloudRec& sm
        );


    // Destructor

        virtual ~forceModelRec();


    // Selector

        static autoPtr<forceModelRec> New
        (
            const dictionary& dict,
            cfdemCloudRec& sm,
            word forceType
        );


    // Member Functions
        virtual void setForce() const = 0;


    // Access
        inline double ** DEMForces() const { return particleCloud_.DEMForces_; }

        inline double ** fluidVel() const { return particleCloud_.fluidVel_; }
        
        inline double ** particleVel() const { return particleCloud_.velocities_; }
        
        void partToArray(label&, const vector&) const;
	
	void partToArrayU(label&, const vector&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
