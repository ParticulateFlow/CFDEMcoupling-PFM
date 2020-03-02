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
    ParmarBassetForce

SourceFiles
    ParmarBassetForce.C

\*---------------------------------------------------------------------------*/

#ifndef ParmarBassetForce_H
#define ParmarBassetForce_H

#include "forceModel.H"
#include "dataExchangeModel.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ParmarBassetForce Declaration
\*---------------------------------------------------------------------------*/

class ParmarBassetForce
:
    public forceModel
{
private:
    dictionary propsDict_;

    word velFieldName_;

    const volVectorField& U_;

    word UsFieldName_;

    const volVectorField& Us_;

    scalar nInt_;                          //no of timesteps to solve full integral

    scalar discOrder_;                     //ODE discretisation order

    scalar nHist_;                         //no of timesteps to save ddtUrel history for

    mutable List<double**> ddtUrelHist_;

    mutable List<double**> rHist_;

    mutable List<List<double**>> FHist_;

    mutable double** gH0_;

    mutable double** tRef_;

    mutable double** mRef_;

    mutable double** lRef_;

    mutable volVectorField Urel_;

    mutable volVectorField ddtUrel_;

    autoPtr<smoothingModel> smoothingModel_;


public:

    //- Runtime type information
    TypeName("ParmarBassetForce");


    // Constructors

        //- Construct from components
        ParmarBassetForce
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~ParmarBassetForce();


    // Member Functions
        void setForce() const;


        void reAllocArrays() const;

        inline const smoothingModel& smoothingM() const
        {
            return smoothingModel_;
        }

        scalar calculateK0(scalar r, scalar dt) const;

        scalar trapWeight(int i, int n) const;

        void update_ddtUrelHist(vector ddtUrel, int index) const;

        void update_rHist(scalar r, int index) const;

        void update_FHist(vector F1, vector F2, int index) const;

        void calculateCoeffs(int k, scalar t0, scalar r, double c[2][4][3], double chi[2][4][2], double C[4]) const;

        void solveFlongODE(int k, double C[4], scalar dt, int index) const;

        void rescaleHist(scalar tScale, scalar mScale, scalar lScale, scalar rScale, int index) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
