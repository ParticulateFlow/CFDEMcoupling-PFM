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
  
    void computeRecMatrix();
    void computeRecPath();
    void readFieldSeries();
    scalar norm(label,label);

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
//    void setRecFields();
    
//     const volScalarField* voidfraction() const;
//     const volVectorField* U() const;
//     const volVectorField* Us() const;
//     const surfaceScalarField* phi() const;
//     const surfaceScalarField* phiS() const;
    
    void exportVolScalarField(word, volScalarField&) const;
    void exportVolVectorField(word, volVectorField&) const;
    void exportSurfaceScalarField(word, surfaceScalarField&) const;
    
	
private:
  
    dictionary propsDict_;
  
//     word voidfractionFieldName_;
//     word UFieldName_;
//     word UsFieldName_;
//     word phiFieldName_;
//     word phiSFieldName_;
//     bool readPhi_,readPhiS_;
  
//     PtrList<volScalarField> voidfractionRecpl;
//     PtrList<volVectorField> URecpl;
//     PtrList<volVectorField> UsRecpl;
//     PtrList<surfaceScalarField> phiRecpl;
//     PtrList<surfaceScalarField> phiSRecpl;
    
    List<PtrList<volScalarField> volScalarFieldList_;
    List<PtrList<volVectorField> volVectorFieldList_;
    List<PtrList<surfaceScalarField> surfaceScalarFieldList_;
    
//     const volScalarField *voidfractionRec_;
//     const volVectorField *URec_;
//     const volVectorField *UsRec_;
//     const surfaceScalarField *phiRec_;
//     const surfaceScalarField *phiSRec_;
    
    word normType_;
    scalar refVol_;
    scalar refVel_;
    
    label seqEnd(label, label&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
