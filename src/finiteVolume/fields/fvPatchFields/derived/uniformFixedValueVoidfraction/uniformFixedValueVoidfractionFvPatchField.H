/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright (C) 2011 OpenFOAM Foundation
                                Copyright (C) 2012-     DCS Computing GmbH,Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::uniformFixedValueVoidfractionFvPatchField

Description
    Foam::uniformFixedValueVoidfractionFvPatchField is an evolution of
    Foam::uniformFixedValueFvPatchField

SourceFiles
    uniformFixedValueVoidfractionFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef uniformFixedValueVoidfractionFvPatchField_H
#define uniformFixedValueVoidfractionFvPatchField_H

#include "Random.H"
#include "fixedValueFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class uniformFixedValueVoidfractionFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class uniformFixedValueVoidfractionFvPatchField
:
    public fixedValueFvPatchField<Type>
{
    // Private data

        autoPtr<Function1<Type> > uniformValue_;

        word voidfractionFieldName_;

public:

    //- Runtime type information
    TypeName("uniformFixedValueVoidfraction");


    // Constructors

        //- Construct from patch and internal field
        uniformFixedValueVoidfractionFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        uniformFixedValueVoidfractionFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given uniformFixedValueVoidfractionFvPatchField
        //  onto a new patch
        uniformFixedValueVoidfractionFvPatchField
        (
            const uniformFixedValueVoidfractionFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        uniformFixedValueVoidfractionFvPatchField
        (
            const uniformFixedValueVoidfractionFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new uniformFixedValueVoidfractionFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        uniformFixedValueVoidfractionFvPatchField
        (
            const uniformFixedValueVoidfractionFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new uniformFixedValueVoidfractionFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "uniformFixedValueVoidfractionFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
