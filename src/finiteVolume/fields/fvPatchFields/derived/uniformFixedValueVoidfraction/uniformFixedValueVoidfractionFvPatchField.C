/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "uniformFixedValueVoidfractionFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
uniformFixedValueVoidfractionFvPatchField<Type>::uniformFixedValueVoidfractionFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    uniformValue_(),
    voidfractionFieldName_("voidfraction")
{}


template<class Type>
uniformFixedValueVoidfractionFvPatchField<Type>::uniformFixedValueVoidfractionFvPatchField
(
    const uniformFixedValueVoidfractionFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper&
)
:
    fixedValueFvPatchField<Type>(p, iF),
    uniformValue_(ptf.uniformValue_, false),
    voidfractionFieldName_("voidfraction")
{
    const scalar t = this->db().time().timeOutputValue();
    fvPatchField<Type>::operator==(uniformValue_->value(t));
}


template<class Type>
uniformFixedValueVoidfractionFvPatchField<Type>::uniformFixedValueVoidfractionFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF),
    uniformValue_(Function1<Type>::New("uniformValue", dict)),
    voidfractionFieldName_("voidfraction")
{
    const scalar t = this->db().time().timeOutputValue();
    fvPatchField<Type>::operator==(uniformValue_->value(t));
}


template<class Type>
uniformFixedValueVoidfractionFvPatchField<Type>::uniformFixedValueVoidfractionFvPatchField
(
    const uniformFixedValueVoidfractionFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    uniformValue_(ptf.uniformValue_, false),
    voidfractionFieldName_("voidfraction")
{
    const scalar t = this->db().time().timeOutputValue();
    fvPatchField<Type>::operator==(uniformValue_->value(t));
}


template<class Type>
uniformFixedValueVoidfractionFvPatchField<Type>::uniformFixedValueVoidfractionFvPatchField
(
    const uniformFixedValueVoidfractionFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    uniformValue_(ptf.uniformValue_, false),
    voidfractionFieldName_("voidfraction")
{
    const scalar t = this->db().time().timeOutputValue();
    fvPatchField<Type>::operator==(uniformValue_->value(t));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void uniformFixedValueVoidfractionFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    this->setSize(m.size());
    const scalar t = this->db().time().timeOutputValue();
    fvPatchField<Type>::operator==(uniformValue_->value(t));
}


template<class Type>
void uniformFixedValueVoidfractionFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const fvPatchField<scalar>& voidfraction = this->patch().template lookupPatchField<volScalarField, scalar>(voidfractionFieldName_);
    const scalar t = this->db().time().timeOutputValue();
    fvPatchField<Type>::operator==(1./voidfraction*uniformValue_->value(t));

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void uniformFixedValueVoidfractionFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    this->writeEntry("value", os);
    uniformValue_->writeData(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
