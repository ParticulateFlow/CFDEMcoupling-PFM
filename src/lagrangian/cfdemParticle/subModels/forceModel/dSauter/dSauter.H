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
    calculates the Sauter mean diameter \sum d_i^3 / \sum d_i^2

SourceFiles
    dSauter.C
\*---------------------------------------------------------------------------*/

#ifndef dSauter_H
#define dSauter_H

#include "forceModel.H"
#include "averagingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dSauter Declaration
\*---------------------------------------------------------------------------*/

class dSauter
:
    public forceModel
{
private:

    dictionary propsDict_;

    bool multiTypes_;

    mutable double **d2_;
    
    mutable double **d3_;

    scalarList typeCG_;
    
    mutable volScalarField d2Field_;
    
    mutable volScalarField d3Field_;
    
    mutable volScalarField dSauter_;

    void allocateMyArrays() const;

public:

    //- Runtime type information
    TypeName("dSauter");

    // Constructors

        //- Construct from components
        dSauter
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~dSauter();


    // Member Functions
        void setForce() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
