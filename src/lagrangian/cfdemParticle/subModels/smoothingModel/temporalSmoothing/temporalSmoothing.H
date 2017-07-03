/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
                                Copyright (C) 2013-     Graz University of  
                                                        Technology, IPPT
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
    temporalSmoothing

SourceFiles
    temporalSmoothing.C

\*---------------------------------------------------------------------------*/

#ifndef temporalSmoothing_H
#define temporalSmoothing_H

#include "smoothingModel.H"
//#include "wallFvPatch.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class temporalSmoothing Declaration
\*---------------------------------------------------------------------------*/

class temporalSmoothing
:
    public smoothingModel
{

private:

    dictionary propsDict_;
    scalar lowerLimit_;
    scalar upperLimit_;
    mutable volScalarField DT_;
    bool verbose_;
    word refFieldName_;
    const volScalarField& refField_;
    scalar alpha_;

public:

    //- Runtime type information
    TypeName("temporalSmoothing");


    // Constructors

        //- Construct from components
        temporalSmoothing
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~temporalSmoothing();


    // Member Functions
        bool doSmoothing() const;

        //void dSmoothing(volScalarField&) const;
        
        void smoothen(volScalarField&) const;

        void smoothen(volVectorField&) const;

        void smoothenReferenceField(volVectorField&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
