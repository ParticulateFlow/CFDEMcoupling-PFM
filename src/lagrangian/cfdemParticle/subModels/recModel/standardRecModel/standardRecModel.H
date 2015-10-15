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

public:

    //- Runtime type information
    TypeName("standardRecModel");


    // Constructors

        //- Construct from components
        standardRecModel
        (
            const dictionary& dict,
            cfdemCloudRec& sm
        );

    // Destructor

        ~standardRecModel();


    void initRecFields();
    void updateRecFields();
    void writeRecFields() const;
    
    autoPtr<const volScalarField> voidfraction() const;
    autoPtr<const volVectorField> U() const;
    autoPtr<const volVectorField> Us() const;
    
    tmp<volScalarField> tvoidfraction() const;
    tmp<volVectorField> tU() const;
    tmp<volVectorField> tUs() const;
	
private:
  
    dictionary propsDict_;
  
    word UFieldName_;
    word voidfractionFieldName_;
    word UsFieldName_;
  
    PtrList<volScalarField> voidfractionRecpl;
    PtrList<volVectorField> URecpl;
    PtrList<volVectorField> UsRecpl;
    
    autoPtr<const volScalarField> voidfractionRec_;
    autoPtr<const volVectorField> URec_;
    autoPtr<const volVectorField> UsRec_;
    
    label seqEnd(label, label&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
