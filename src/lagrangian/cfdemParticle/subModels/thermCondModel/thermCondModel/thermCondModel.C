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

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "thermCondModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(thermCondModel, 0);

defineRunTimeSelectionTable(thermCondModel, dictionary);

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
thermCondModel::thermCondModel
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    dict_(dict),
    particleCloud_(sm),
    mesh_(particleCloud_.mesh()),
    transportProperties_
    (
        IOobject
        (
            "transportProperties",
            sm.mesh().time().constant(),
            sm.mesh(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    kf0_(transportProperties_.lookup("kf")),
    thermCondField_(const_cast<volScalarField&>(sm.mesh().lookupObject<volScalarField> ("thCond"))),
    wallQFactor_
    (   IOobject
        (
            "wallQFactor",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0,0,0), 1.0)
    ),
    hasWallQFactor_(false),
    wallBoundaryLayerThickness_
    (   IOobject
        (
            "wallBoundaryLayerThickness",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0,1,0,0,0,0,0), 1.0)
    ),
    hasWallBoundaryLayerThickness_(false),
    wallHeatLoss_
    (   IOobject
        (
            "wallHeatLoss",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(1,0,-3,-1,0,0,0), 0.0)
    ),
    hasWallHeatLoss_(false)
{
    if (wallQFactor_.headerOk())
    {
        Info << "Found field for scaling wall heat flux.\n" << endl;
        hasWallQFactor_ = true;
        wallQFactor_.writeOpt() = IOobject::AUTO_WRITE;
    }
    if (wallBoundaryLayerThickness_.headerOk())
    {
        Info << "Found field near-wall boundary layer thickness.\n" << endl;
        hasWallBoundaryLayerThickness_ = true;
        wallBoundaryLayerThickness_.writeOpt() = IOobject::AUTO_WRITE;
    }
    if (wallHeatLoss_.headerOk())
    {
        Info << "Found field for wall heat loss.\n" << endl;
        hasWallHeatLoss_ = true;
        wallHeatLoss_.writeOpt() = IOobject::AUTO_WRITE;
    }

    if (hasWallQFactor_ + hasWallBoundaryLayerThickness_ + hasWallHeatLoss_ > 1)
    {
        FatalError << "thermCondModel: cannot use more than one option of wallQFactor, wallBoundaryLayerThickness and wallHeatLoss." << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

thermCondModel::~thermCondModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void thermCondModel::calcBoundaryCorrections()
{
    if(!hasWallQFactor_ && !hasWallHeatLoss_ && !hasWallBoundaryLayerThickness_) return;

    const fvPatchList& patches = mesh_.boundary();
    // if a wallQFactor field is present, use it to scale heat transport through a patch
    if (hasWallQFactor_)
    {
        wallQFactor_.correctBoundaryConditions();
        forAll(patches, patchi)
        {
            const fvPatch& curPatch = patches[patchi];
            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];
                thermCondField_.boundaryFieldRef()[patchi][facei] = thermCondField_[faceCelli]*wallQFactor_.boundaryField()[patchi][facei];
            }

        }
    }
    else if (hasWallHeatLoss_)
    {
        wallHeatLoss_.correctBoundaryConditions();
        forAll(patches, patchi)
        {
            // since not explicitly looped over all faces, need to use == to force assignment
            thermCondField_.boundaryFieldRef()[patchi] == wallHeatLoss_.boundaryField()[patchi] / mesh_.deltaCoeffs().boundaryField()[patchi];
        }
    }
    else if (hasWallBoundaryLayerThickness_)
    {
        wallBoundaryLayerThickness_.correctBoundaryConditions();
        forAll(patches, patchi)
        {
            const fvPatch& curPatch = patches[patchi];
            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];
                thermCondField_.boundaryFieldRef()[patchi][facei] = thermCondField_[faceCelli] /
                    (wallBoundaryLayerThickness_.boundaryField()[patchi][facei] * mesh_.deltaCoeffs().boundaryField()[patchi][facei]);
            }
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam

// ************************************************************************* //
