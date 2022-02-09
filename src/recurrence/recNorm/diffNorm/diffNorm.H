/*---------------------------------------------------------------------------*\
    CFDEMcoupling academic - Open Source CFD-DEM coupling

    Contributing authors:
    Thomas Lichtenegger
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

#ifndef diffNorm_H
#define diffNorm_H

#include "recNorm.H"
#include "cellSet.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class diffNorm Declaration
\*---------------------------------------------------------------------------*/

class diffNorm
:
    public recNorm
{
protected:

    // Protected data

    dictionary propsDict_;

    const scalar normConstant_;

    word fieldType_;

    word fieldName_;

    word integrationDomainName_;

    autoPtr<cellSet> integrationDomain_;

    bool existIntegrationDomain_;

    void computeRecMatrix();

    void restrictToIntegrationDomain(volScalarField &vsf);

    void restrictToIntegrationDomain(volVectorField &vvf);

    void restrictToIntegrationDomain(surfaceScalarField &ssf);

    virtual scalar normVSF(label, label) {return 0.0;}

    virtual scalar normVVF(label, label) {return 0.0;}

    virtual scalar normSSF(label, label) {return 0.0;}



public:

    //- Runtime type information
    TypeName("diffNorm");

    // Constructors

        //- Construct from components
        diffNorm
        (
            const dictionary& dict,
            recBase& base,
            word type = "diffNorm"
        );


    // Destructor

        virtual ~diffNorm();


    // Member Functions


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
