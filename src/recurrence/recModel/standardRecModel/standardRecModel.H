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

#ifndef standardRecModel_H
#define standardRecModel_H

#include "recModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class standardRecModel Declaration
\*---------------------------------------------------------------------------*/

class standardRecModel
:
    public recModel
{
protected:
  
    void readFieldSeries();

public:

    //- Runtime type information
    TypeName("standardRecModel");


    // Constructors

        //- Construct from components
        standardRecModel
        (
            const dictionary& dict,
            recBase& base
        );

    // Destructor

        ~standardRecModel();


    void updateRecFields();
    
    void exportVolScalarField(word, volScalarField&) const;
    void exportVolVectorField(word, volVectorField&) const;
    void exportSurfaceScalarField(word, surfaceScalarField&) const;
    
    const volScalarField& exportVolScalarField(word, label) const;
    const volVectorField& exportVolVectorField(word, label) const;
    const surfaceScalarField& exportSurfaceScalarField(word, label) const;
    
    inline label getVolScalarFieldIndex(word, label) const;
    inline label getVolVectorFieldIndex(word, label) const;
    inline label getSurfaceScalarFieldIndex(word, label) const;
	
private:
  
    dictionary propsDict_;
    
    List<PtrList<volScalarField> > volScalarFieldList_;
    List<PtrList<volVectorField> > volVectorFieldList_;
    List<PtrList<surfaceScalarField> > surfaceScalarFieldList_;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "standardRecModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
