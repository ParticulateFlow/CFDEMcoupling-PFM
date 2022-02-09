/*---------------------------------------------------------------------------*\
    CFDEMcoupling academic - Open Source CFD-DEM coupling

    Contributing authors:
    Thomas Lichtenegger
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

#include "error.H"

#include "turbulentDispersion.H"
#include "addToRunTimeSelectionTable.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(turbulentDispersion, 0);

addToRunTimeSelectionTable
(
    forceModel,
    turbulentDispersion,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
turbulentDispersion::turbulentDispersion
(
    const dictionary& dict,
    cfdemCloud& sm,
    word type
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(type + "Props")),
    interpolate_(propsDict_.lookupOrDefault<bool>("interpolation", false)),
    mesh_(sm.mesh()),
    ignoreCellsName_(propsDict_.lookupOrDefault<word>("ignoreCellsName","none")),
    ignoreCells_(),
    existIgnoreCells_(true),
    usePreCalcNut_(propsDict_.lookupOrDefault<bool>("usePreCalcNut",false)),
    nutName_(propsDict_.lookupOrDefault<word>("nutName","nut")),
    usePreCalcK_(propsDict_.lookupOrDefault<bool>("usePreCalcK",false)),
    kName_(propsDict_.lookupOrDefault<word>("kName","k")),
    usePreCalcDispField_(propsDict_.lookupOrDefault<bool>("usePreCalcDispField",false)),
    dispFieldName_(propsDict_.lookupOrDefault<word>("dispFieldName","dispVarField")),
    nut_
    (   IOobject
        (
            "nuTurb",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0,2,-1,0,0,0,0), 0.0)
    ),
    wallIndicatorField_
    (   IOobject
        (
            "wallIndicator",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0,0,0), 0.0)
    ),
    delta_
    (   IOobject
        (
            "delta",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("delta", dimLength, 1.0)
    ),
    turbulentSchmidtNumber_(readScalar(propsDict_.lookup("turbulentSchmidtNumber"))),
    voidfractionFieldName_(propsDict_.lookupOrDefault<word>("voidfractionFieldName","voidfraction")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    critVoidfraction_(propsDict_.lookupOrDefault<scalar>("critVoidfraction", 0.9)),
    Ck_(propsDict_.lookupOrDefault<scalar>("turbKineticEnergyCoeff", 0.094)),
    ranGen_(clock::getTime()+pid())
{
    if (ignoreCellsName_ != "none")
    {
       ignoreCells_.set(new cellSet(particleCloud_.mesh(),ignoreCellsName_));
       Info << type << ": ignoring fluctuations in cellSet " << ignoreCells_().name() <<
        " with " << ignoreCells_().size() << " cells." << endl;
    }
    else existIgnoreCells_ = false;

    if ((usePreCalcNut_ && usePreCalcK_) || (usePreCalcNut_ && usePreCalcDispField_) || (usePreCalcDispField_ && usePreCalcK_))
    {
        FatalError<< "Cannot use more than one precalculated nut, k and displacement fluctuations at the same time. Choose one." << abort(FatalError);
    }

    if (usePreCalcK_)
    {
        delta_.primitiveFieldRef()=Foam::pow(mesh_.V(),1.0/3.0);
    }

    // define a field to indicate if a cell is next to boundary
     label cellI = -1;
     forAll (mesh_.boundary(),patchI)
     {
         word patchName = mesh_.boundary()[patchI].name();
         if (patchName.rfind("procB",0) == 0) continue;

         forAll(mesh_.boundary()[patchI], faceI)
         {
             cellI = mesh_.boundary()[patchI].faceCells()[faceI];
             wallIndicatorField_[cellI] = 1.0;
         }
     }

    scalar dtCFD = voidfraction_.mesh().time().deltaTValue();
    scalar dtDEM = particleCloud_.dataExchangeM().DEMts();
    // if CFD step is larger than DEM step, a corresponding number of DEM steps is taken to reach the CFD step;
    // if DEM step is larger than CFD step, no update occurs so that full DEM step needs to be taken
    dt_ = max(dtCFD, dtDEM);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

turbulentDispersion::~turbulentDispersion()
{}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //

bool turbulentDispersion::ignoreCell(label cell) const
{
    if (!existIgnoreCells_) return false;
    else return ignoreCells_()[cell];
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void turbulentDispersion::setForce() const
{
    if (!usePreCalcNut_ && !usePreCalcK_ && !usePreCalcDispField_)
    {
        nut_ = particleCloud_.turbulence().nut()();
    }
    else if (usePreCalcNut_)
    {
        volScalarField& nutRef (const_cast<volScalarField&>(mesh_.lookupObject<volScalarField> (nutName_)));
        nut_ = nutRef;
    }
    else if (usePreCalcK_)
    {
        volScalarField& kRef (const_cast<volScalarField&>(mesh_.lookupObject<volScalarField> (kName_)));
        nut_ = Ck_ * delta_ * sqrt(kRef);
    }
    else
    {
        volVectorField& dispFieldRef (const_cast<volVectorField&>(mesh_.lookupObject<volVectorField> (dispFieldName_)));
        dispField_ = &dispFieldRef;
    }

    label cellI = -1;
    label patchID = -1;
    label faceIGlobal = -1;
    scalar flucProjection = 0.0;
    scalar D = 0.0;
    scalar randScalar = 0.0;
    vector faceINormal = vector::zero;
    vector flucU = vector::zero;
    vector position = vector::zero;
    word patchName("");

    interpolationCellPoint<scalar> nutInterpolator_(nut_);

    for(int index = 0; index < particleCloud_.numberOfParticles(); ++index)
    {
        cellI = particleCloud_.cellIDs()[index][0];
        if (cellI > -1 && !ignoreCell(cellI))
        {
            if (usePreCalcDispField_)
            {
                for (label comp=0; comp<3; comp++)
                {
#if OPENFOAM_VERSION_MAJOR < 6
                    randScalar = ranGen_.GaussNormal();
#else
                    randScalar = ranGen_.scalarNormal();
#endif
                    flucU.component(comp) = randScalar * (*dispField_)[cellI].component(comp);
                }
            }
            else
            {
                if (interpolate_)
                {
                    position = particleCloud_.position(index);
                    D = nutInterpolator_.interpolate(position,cellI) / turbulentSchmidtNumber_;
                }
                else
                {
                    D = nut_[cellI] / turbulentSchmidtNumber_;
                }

                // include concentration dependence on the diffusivity at this point if necessary

                flucU=unitFlucDir()*Foam::sqrt(6.0*D/dt_);
            }

            // prevent particles being pushed through walls by regulating velocity fluctuations
            // check if cell is adjacent to wall and remove corresponding components
            if (wallIndicatorField_[cellI] > 0.5)
            {
                const cell& faces = mesh_.cells()[cellI];
                forAll (faces, faceI)
                {
                    faceIGlobal = faces[faceI];
                    patchID = mesh_.boundaryMesh().whichPatch(faceIGlobal);
                    if (patchID < 0) continue;
                    patchName = mesh_.boundary()[patchID].name();

                    if (patchName.rfind("procB",0) == 0) continue;

                    faceINormal = mesh_.Sf()[faceIGlobal];
                    faceINormal /= mag(faceINormal);
                    flucProjection = faceINormal&flucU;
                    if (flucProjection > 0.0) flucU -= flucProjection*faceINormal;
                }
            }

            for(int j=0;j<3;j++)
            {
                particleCloud_.particleFlucVels()[index][j] += flucU[j];
            }
        }
    }
}

vector turbulentDispersion::unitFlucDir() const
{
    // unit random vector
    // algorithm according to:
    // Marsaglia. "Choosing a point from the surface of a sphere." The Annals of Mathematical Statistics 43.2 (1972): 645-646.
    scalar v1(0.0);
    scalar v2(0.0);
    scalar s(10.0);
    scalar s2(0.0);
    vector rvec(0,0,0);
    while(s>1.0)
    {
        v1=2*(ranGen_.scalar01()-0.5);
        v2=2*(ranGen_.scalar01()-0.5);
        s=v1*v1+v2*v2;
    }
    s2=Foam::sqrt(1-s);
    rvec[0]=2*v1*s2;
    rvec[1]=2*v2*s2;
    rvec[2]=1-2*s;
    return rvec;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
