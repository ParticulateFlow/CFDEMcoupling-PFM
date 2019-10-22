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

    Copyright (C) 2018- Mathias Vångö, JKU Linz, Austria

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "temporalSmoothing.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(temporalSmoothing, 0);

addToRunTimeSelectionTable
(
    smoothingModel,
    temporalSmoothing,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
temporalSmoothing::temporalSmoothing
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    smoothingModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    lowerLimit_(readScalar(propsDict_.lookup("lowerLimit"))),
    upperLimit_(readScalar(propsDict_.lookup("upperLimit"))),
    verbose_(false),
    refFieldName_(propsDict_.lookup("refField")),
    gamma_(readScalar(propsDict_.lookup("gamma")))
{

    if(propsDict_.found("verbose"))
        verbose_ = true;

    checkFields(sSmoothField_);
    checkFields(vSmoothField_);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

temporalSmoothing::~temporalSmoothing()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
bool temporalSmoothing::doSmoothing() const
{
    return true;
}


void Foam::temporalSmoothing::smoothen(volScalarField& fieldSrc) const
{
    // Create scalar smooth field from virgin scalar smooth field template
    volScalarField sSmoothField = sSmoothField_;

    sSmoothField.dimensions().reset(fieldSrc.dimensions());
    sSmoothField.ref()=fieldSrc.internalField();
    sSmoothField.correctBoundaryConditions();
    sSmoothField.oldTime().dimensions().reset(fieldSrc.dimensions());
    sSmoothField.oldTime()=fieldSrc;
    sSmoothField.oldTime().correctBoundaryConditions();

    volScalarField refField = particleCloud_.mesh().lookupObject<volScalarField>(refFieldName_);

    // do smoothing
    dimensionedScalar deltaT = sSmoothField.mesh().time().deltaT();
    solve
    (
        fvm::ddt(sSmoothField)
        -
        gamma_/deltaT * (refField - fvm::Sp(1.0,sSmoothField))
    );

    // bound sSmoothField_
    forAll(sSmoothField,cellI)
    {
        sSmoothField[cellI]=max(lowerLimit_,min(upperLimit_,sSmoothField[cellI]));
    }

    // get data from working sSmoothField - will copy only values at new time
    fieldSrc=sSmoothField;
    fieldSrc.correctBoundaryConditions();

    if(verbose_)
    {
        Info << "min/max(fieldoldTime) (unsmoothed): " << min(sSmoothField.oldTime()) << tab << max(sSmoothField.oldTime()) << endl;
        Info << "min/max(fieldSrc): " << min(fieldSrc) << tab << max(fieldSrc) << endl;
        Info << "min/max(fieldSrc.oldTime): " << min(fieldSrc.oldTime()) << tab << max(fieldSrc.oldTime()) << endl;
    }

}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void Foam::temporalSmoothing::smoothen(volVectorField& fieldSrc) const
{
    // Create scalar smooth field from virgin scalar smooth field template
    volVectorField vSmoothField = vSmoothField_;

    vSmoothField.dimensions().reset(fieldSrc.dimensions());
    vSmoothField.ref()=fieldSrc.internalField();
    vSmoothField.correctBoundaryConditions();
    vSmoothField.oldTime().dimensions().reset(fieldSrc.dimensions());
    vSmoothField.oldTime()=fieldSrc;
    vSmoothField.oldTime().correctBoundaryConditions();

    volVectorField refField = particleCloud_.mesh().lookupObject<volVectorField>(refFieldName_);

    dimensionedScalar deltaT = vSmoothField.mesh().time().deltaT();
    solve
    (
        fvm::ddt(vSmoothField)
        -
        gamma_/deltaT * (refField - fvm::Sp(1.0,vSmoothField))
    );

    // get data from working vSmoothField
    fieldSrc=vSmoothField;
    fieldSrc.correctBoundaryConditions();

    if(verbose_)
    {
        Info << "min/max(fieldoldTime) (unsmoothed): " << min(vSmoothField.oldTime()) << tab << max(vSmoothField.oldTime()) << endl;
        Info << "min/max(fieldSrc): " << min(fieldSrc) << tab << max(fieldSrc) << endl;
        Info << "min/max(fieldSrc.oldTime): " << min(fieldSrc.oldTime()) << tab << max(fieldSrc.oldTime()) << endl;
    }

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void Foam::temporalSmoothing::smoothenReferenceField(volVectorField& fieldSrc) const
{
    FatalError << "Smoothen reference field is not implemented for this smoothing model!" << abort(FatalError);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
