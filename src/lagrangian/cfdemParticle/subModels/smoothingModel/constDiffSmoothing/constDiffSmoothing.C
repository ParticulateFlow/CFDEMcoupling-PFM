/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
                                Copyright (C) 2013-     Graz University of
                                                        Technology, IPPT
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
\*---------------------------------------------------------------------------*/

#include "error.H"

#include "constDiffSmoothing.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(constDiffSmoothing, 0);

addToRunTimeSelectionTable
(
    smoothingModel,
    constDiffSmoothing,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
constDiffSmoothing::constDiffSmoothing
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    smoothingModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    lowerLimit_(readScalar(propsDict_.lookup("lowerLimit"))),
    upperLimit_(readScalar(propsDict_.lookup("upperLimit"))),
    smoothingLength_(propsDict_.lookupOrDefault<scalar>("smoothingLength", -1.0)),
    smoothingLengthReference_(propsDict_.lookupOrDefault<scalar>("smoothingLengthReference",smoothingLength_)),
    smoothingLengthFieldName_(propsDict_.lookupOrDefault<word>("smoothingLengthFieldName","smoothingLengthField")),
    smoothingLengthField_
    (   IOobject
        (
            smoothingLengthFieldName_,
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("smoothingLength", dimensionSet(0,1,0,0,0,0,0), smoothingLength_),
        "zeroGradient"
    ),
    smoothingLengthReferenceFieldName_(propsDict_.lookupOrDefault<word>("smoothingLengthReferenceFieldName","smoothingLengthReferenceField")),
    smoothingLengthReferenceField_
    (   IOobject
        (
            smoothingLengthReferenceFieldName_,
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("smoothingLengthReference", dimensionSet(0,1,0,0,0,0,0), smoothingLengthReference_),
        "zeroGradient"
    ),
    DT_
    (   IOobject
        (
            "DT",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("DT", dimensionSet(0,2,-1,0,0,0,0), 0.0),
        "zeroGradient"
    ),
    verbose_(propsDict_.found("verbose"))
{
    // either use scalar or field parameters for smoothing
    if (smoothingLength_ > 0.0 || smoothingLengthReference_ > 0.0)
    {
        if (smoothingLengthField_.headerOk() || smoothingLengthReferenceField_.headerOk())
        {
            FatalError <<"constDiffSmoothing: Either use scalar or field parameter for smoothing.\n" << abort(FatalError);
        }
    }

    if (smoothingLength_ < 0.0 && !smoothingLengthField_.headerOk())
    {
        FatalError <<"constDiffSmoothing: Provide scalar or field parameter for smoothing.\n" << abort(FatalError);
    }
  
    // if no scalar length for smoothing wrt reference field is provided and no
    // such field, use smoothingLengthField
    if (smoothingLengthReference_ < 0.0 && !smoothingLengthReferenceField_.headerOk())
    {
        smoothingLengthReferenceField_ = smoothingLengthField_;
    }

    checkFields(sSmoothField_);
    checkFields(vSmoothField_);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

constDiffSmoothing::~constDiffSmoothing()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
bool constDiffSmoothing::doSmoothing() const
{
    return true;
}


void constDiffSmoothing::smoothen(volScalarField& fieldSrc) const
{
    // Create scalar smooth field from virgin scalar smooth field template
    volScalarField sSmoothField = sSmoothField_;

    sSmoothField.dimensions().reset(fieldSrc.dimensions());
    sSmoothField.ref()=fieldSrc.internalField();
    sSmoothField.correctBoundaryConditions();
    sSmoothField.oldTime().dimensions().reset(fieldSrc.dimensions());
    sSmoothField.oldTime()=fieldSrc;
    sSmoothField.oldTime().correctBoundaryConditions();

    dimensionedScalar deltaT = sSmoothField.mesh().time().deltaT();
    DT_ = smoothingLengthField_ * smoothingLengthField_ / deltaT;

    // do smoothing
    solve
    (
        fvm::ddt(sSmoothField)
       -fvm::laplacian(DT_, sSmoothField)
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
void constDiffSmoothing::smoothen(volVectorField& fieldSrc) const
{
    // Create scalar smooth field from virgin scalar smooth field template
    volVectorField vSmoothField = vSmoothField_;

    vSmoothField.dimensions().reset(fieldSrc.dimensions());
    vSmoothField.ref()=fieldSrc.internalField();
    vSmoothField.correctBoundaryConditions();
    vSmoothField.oldTime().dimensions().reset(fieldSrc.dimensions());
    vSmoothField.oldTime()=fieldSrc;
    vSmoothField.oldTime().correctBoundaryConditions();

    dimensionedScalar deltaT = vSmoothField.mesh().time().deltaT();
    DT_ = smoothingLengthField_ * smoothingLengthField_ / deltaT;

    // do smoothing
    solve
    (
        fvm::ddt(vSmoothField)
       -fvm::laplacian(DT_, vSmoothField)
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
void constDiffSmoothing::smoothenReferenceField(volVectorField& fieldSrc) const
{
    // Create scalar smooth field from virgin scalar smooth field template
    volVectorField vSmoothField = vSmoothField_;

    vSmoothField.dimensions().reset(fieldSrc.dimensions());
    vSmoothField.ref()=fieldSrc.internalField();
    vSmoothField.correctBoundaryConditions();
    vSmoothField.oldTime().dimensions().reset(fieldSrc.dimensions());
    vSmoothField.oldTime()=fieldSrc;
    vSmoothField.oldTime().correctBoundaryConditions();

    double sourceStrength = 1e5; //large number to keep reference values constant

    dimensionedScalar deltaT = vSmoothField.mesh().time().deltaT();
    DT_ = smoothingLengthReferenceField_ * smoothingLengthReferenceField_ / deltaT;
        
    tmp<volScalarField> NLarge
    (
        new volScalarField
        (
            IOobject
            (
                "NLarge",
                particleCloud_.mesh().time().timeName(),
                particleCloud_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            particleCloud_.mesh(),
            0.0
        )
    );


    //loop over particles and map max particle diameter to Euler Grid
    forAll(vSmoothField,cellI)
    {
        if ( mag(vSmoothField.oldTime().internalField()[cellI]) > 0.0f)  // have a vector in the OLD vSmoothField, so keep it!
            NLarge.ref()[cellI] = sourceStrength;
    }

    // do the smoothing
    solve
    (
        fvm::ddt(vSmoothField)
       -fvm::laplacian( DT_, vSmoothField)
       ==
        NLarge() / deltaT * vSmoothField.oldTime()  //add source to keep cell values constant
       -fvm::Sp( NLarge() / deltaT, vSmoothField)   //add sink to keep cell values constant
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

} // End namespace Foam

// ************************************************************************* //
