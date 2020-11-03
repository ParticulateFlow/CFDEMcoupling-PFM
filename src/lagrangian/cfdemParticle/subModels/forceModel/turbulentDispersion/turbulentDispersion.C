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
    minTurbKinetcEnergy_(propsDict_.lookupOrDefault<scalar>("minTurbKinetcEnergy", 0.0)),
    voidfractionFieldName_(propsDict_.lookupOrDefault<word>("voidfractionFieldName","voidfraction")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    critVoidfraction_(propsDict_.lookupOrDefault<scalar>("critVoidfraction", 0.9)),
    ranGen_(clock::getTime()+pid())
{
    if(ignoreCellsName_ != "none")
    {
       ignoreCells_.set(new cellSet(particleCloud_.mesh(),ignoreCellsName_));
       Info << type << ": ignoring fluctuations in cellSet " << ignoreCells_().name() <<
        " with " << ignoreCells_().size() << " cells." << endl;
    }
    else existIgnoreCells_ = false;

    // define a field to indicate if a cell is next to boundary
     label cellI = -1;
    // set wall indicator field
Info << "Setting wall indicator field." << endl;
     forAll(mesh_.boundary(),patchI)
     {
         word patchName = mesh_.boundary()[patchI].name();
         if (patchName.rfind("procBoundary",0) == 0) continue;
Info << "patch = " << mesh_.boundary()[patchI].name() << endl;
         forAll(mesh_.boundary()[patchI], faceI)
         {
             cellI = mesh_.boundary()[patchI].faceCells()[faceI];
             wallIndicatorField_[cellI] = 1.0;

// testing
                    label patchID = -1;
                    label faceind = -1;
                    const cell& faces = mesh_.cells()[cellI];
Pout << "cellI = " << cellI << " has faces = " << faces << endl;
                    forAll (faces, faceJ)        // loop over all faces in cellI
                    {
                        faceind = faces[faceJ];
                        patchID = mesh_.boundaryMesh().whichPatch(faceind);
                        if (patchID < 0) continue;
                        vector faceINormal = mesh_.Sf()[faceind] / mesh_.magSf()[faceind] ;
      //                  Pout << " faceind = " << faceind << "at " << mesh_.C()[cellI] << ", faceINormal = " << faceINormal << endl ;
                    }

// testing done

         }
     }

wallIndicatorField_.write();

    // make sure this is the last force model in list so that fluid velocity does not get overwritten
    label numLastForceModel = sm.nrForceModels();
    word lastForceModel = sm.forceModels()[numLastForceModel-1];
    if (lastForceModel != "turbulentDispersion")
    {
        FatalError <<"Force model 'turbulentDispersion' needs to be last in list!\n" << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

turbulentDispersion::~turbulentDispersion()
{
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //

bool turbulentDispersion::ignoreCell(label cell) const
{
    if (!existIgnoreCells_) return false;
    else return ignoreCells_()[cell];
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void turbulentDispersion::setForce() const
{
    const volScalarField turbKinetcEnergy(particleCloud_.turbulence().k());

    vector position(0,0,0);
    scalar k = 0.0;

    vector flucU(0,0,0);
    label cellI = 0;

    interpolationCellPoint<scalar> turbKinetcEnergyInterpolator_(turbKinetcEnergy);

    for(int index = 0;index <  particleCloud_.numberOfParticles(); ++index)
    {
        cellI = particleCloud_.cellIDs()[index][0];
        flucU = vector(0,0,0);
        k=0.0;

        if (cellI > -1 && !ignoreCell(cellI))
        {
            // particles in dilute regions follow fluid without fluctuations
            if (voidfraction_[cellI] < critVoidfraction_)
            {
                if (interpolate_)
                {
                    position = particleCloud_.position(index);
                    k = turbKinetcEnergyInterpolator_.interpolate(position,cellI);
                }
                else
                {
                    k = turbKinetcEnergy[cellI];
                }

                if (k < minTurbKinetcEnergy_) k = minTurbKinetcEnergy_;

                flucU=unitFlucDir()*Foam::sqrt(2.0*k);

                // if particles are pushed through walls, the velocity fluctuations may be regulated here
                // check if cell is adjacent to wall and invert corresponding component
                label patchID = -1;
                if (wallIndicatorField_[cellI] > 0.5)
                {
                    const cell& faces = mesh_.cells()[cellI];
                    forAll (faces, faceI)        // loop over all faces in cellI
                    {
                        patchID = mesh_.boundaryMesh().whichPatch(faceI);
                        if (patchID < 0) continue;

                        vector faceINormal = mesh_.Sf()[faceI] / mesh_.magSf()[faceI] ;
                        Info << " faceI = " << faceI << "at " << mesh_.C()[cellI] << ", faceINormal = " << faceINormal << endl ;
                    }
                }

                for(int j=0;j<3;j++)
                {
                    particleCloud_.fluidVels()[index][j] += flucU[j];
                }
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
