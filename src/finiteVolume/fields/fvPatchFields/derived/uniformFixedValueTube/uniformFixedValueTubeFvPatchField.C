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

#include "uniformFixedValueTubeFvPatchField.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
uniformFixedValueTubeFvPatchField<Type>::uniformFixedValueTubeFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    uniformValue_(),
	pName_("p"), //JOKER
	phiName_("phi"), //JOKER
    velocityFieldName_("U"),
    densityFieldName_("rho"),
    tubeLength_(-1.),
    tubeDiameter_(-1.),
    zeta_(-1)
{
    Info << "error, wrong constructor1!" << endl;
}


template<class Type>
uniformFixedValueTubeFvPatchField<Type>::uniformFixedValueTubeFvPatchField
(
    const uniformFixedValueTubeFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper&
)
:
    fixedValueFvPatchField<Type>(p, iF),
    uniformValue_(ptf.uniformValue_().clone().ptr()),
	pName_("p"), //JOKER
	phiName_("phi"), //JOKER
    velocityFieldName_("U"),
    densityFieldName_("rho"),
    tubeLength_(ptf.tubeLength_),
    tubeDiameter_(ptf.tubeDiameter_),
    zeta_(ptf.zeta_)
{
    const scalar t = this->db().time().timeOutputValue();
    fvPatchField<Type>::operator==(uniformValue_->value(t));
    Info << "error, wrong constructor2!" << endl;
}


template<class Type>
uniformFixedValueTubeFvPatchField<Type>::uniformFixedValueTubeFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF),
    uniformValue_(Function1<Type>::New("uniformValue", dict)),
	pName_("p"), //JOKER
	phiName_("phi"), //JOKER
    velocityFieldName_("U"),
    densityFieldName_("rho"),
    tubeLength_(readScalar(dict.lookup("tubeLength"))),
    tubeDiameter_(readScalar(dict.lookup("tubeDiameter"))),
    zeta_(readScalar(dict.lookup("zeta")))
{
    const scalar t = this->db().time().timeOutputValue();
    fvPatchField<Type>::operator==(uniformValue_->value(t));
}


template<class Type>
uniformFixedValueTubeFvPatchField<Type>::uniformFixedValueTubeFvPatchField
(
    const uniformFixedValueTubeFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    uniformValue_(ptf.uniformValue_().clone().ptr()),
	pName_("p"), //JOKER
	phiName_("phi"), //JOKER
    velocityFieldName_("U"),
    densityFieldName_("rho"),
    tubeLength_(ptf.tubeLength_),
    tubeDiameter_(ptf.tubeDiameter_),
    zeta_(ptf.zeta_)
{
    const scalar t = this->db().time().timeOutputValue();
    fvPatchField<Type>::operator==(uniformValue_->value(t));
}


template<class Type>
uniformFixedValueTubeFvPatchField<Type>::uniformFixedValueTubeFvPatchField
(
    const uniformFixedValueTubeFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    uniformValue_(ptf.uniformValue_().clone().ptr()),
	pName_("p"), //JOKER
	phiName_("phi"), //JOKER
    velocityFieldName_("U"),
    densityFieldName_("rho"),
    tubeLength_(ptf.tubeLength_),
    tubeDiameter_(ptf.tubeDiameter_),
    zeta_(ptf.zeta_)
{
    const scalar t = this->db().time().timeOutputValue();
    fvPatchField<Type>::operator==(uniformValue_->value(t));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void uniformFixedValueTubeFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    this->setSize(m.size());
    const scalar t = this->db().time().timeOutputValue();
    fvPatchField<Type>::operator==(uniformValue_->value(t));
}


template<class Type>
void uniformFixedValueTubeFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

//    const fvPatchField<vector>& velocity = this->patch().template lookupPatchField<volVectorField, vector>(velocityFieldName_);
    const fvPatchField<scalar>& density = this->patch().template lookupPatchField<volScalarField, scalar>(densityFieldName_);
	const fvPatchField<scalar>& pressure = this->patch().template lookupPatchField<volScalarField, scalar>(pName_);
    const scalar t = this->db().time().timeOutputValue();

	const fvsPatchField<scalar>& phip = this->patch().template lookupPatchField<surfaceScalarField, scalar>(phiName_);

	const scalar uRelFact = 1.e-3;

	//scalar zeta = 0.1; //must be read from dict

	//calc cell velocity from flux phip
	//scalar vel = phip/this->patch().magSf();
    // some relaxation might be useful?
    
    // calc pressure drop
//    fvPatchField<Type>::operator==(0.5*zeta*density*mag(velocity)*mag(velocity)*uniformValue_->value(t));
	//dP = zeta * l/d * rho u^2 / 2
	fvPatchField<Type>::operator==((pressure+(0.5*zeta_*tubeLength_/tubeDiameter_*density*phip/this->patch().magSf()*phip/this->patch().magSf()-pressure)*uRelFact)*uniformValue_->value(t));
	//fvPatchField<Type>::operator==((100.*uniformValue_->value(t)));
    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void uniformFixedValueTubeFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    this->writeEntry("value", os);
    uniformValue_->writeData(os);
    os.writeKeyword("tubeLength") << tubeLength_ << token::END_STATEMENT << nl;
    os.writeKeyword("tubeDiameter") << tubeDiameter_ << token::END_STATEMENT << nl;
    os.writeKeyword("zeta") << zeta_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
