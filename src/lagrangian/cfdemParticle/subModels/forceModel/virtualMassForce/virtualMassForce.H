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

Reference
    Crowe, Clayton T.; Sommerfeld, Martin; Tsuji, Yutaka (1998).
    "Multiphase flows with droplets and particles"

Class
    virtualMassForce

SourceFiles
    virtualMassForce.C

\*---------------------------------------------------------------------------*/

#ifndef virtualMassForce_H
#define virtualMassForce_H

#include "forceModel.H"
#include "dataExchangeModel.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class virtualMassForce Declaration
\*---------------------------------------------------------------------------*/

class virtualMassForce
:
    public forceModel
{
private:
    dictionary propsDict_;

    word velFieldName_;

    const volVectorField& U_;

    word voidfractionFieldName_;

    const volScalarField& voidfraction_;

    word UsFieldName_;

    const volVectorField& Us_;

    word phiFieldName_;

    const surfaceScalarField& phi_;

    const word UrelOldRegName_;

    const bool splitUrelCalculation_;      //indicator to split calculation of Urel between CFDEM and LIGGGHTS
                                           //requires the integration fix to take dv/dt into account!

    mutable bool useUs_;                   //indicator to use averaged Us instead of v_p

    mutable bool useFelderhof_;            //indicator to use Cadd values from Felderhof

    mutable double Cadd_;                  //custom Cadd value

    mutable volVectorField DDtUrel_;

    autoPtr<smoothingModel> smoothingModel_;

public:

    //- Runtime type information
    TypeName("virtualMassForce");


    // Constructors

        //- Construct from components
        virtualMassForce
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~virtualMassForce();


    // Member Functions
        void setForce() const;

        inline const smoothingModel& smoothingM() const
        {
            return smoothingModel_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
