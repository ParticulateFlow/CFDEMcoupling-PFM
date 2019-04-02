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

#include "isotropicFluctuations.H"
#include "addToRunTimeSelectionTable.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(isotropicFluctuations, 0);

addToRunTimeSelectionTable
(
    forceModel,
    isotropicFluctuations,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
isotropicFluctuations::isotropicFluctuations
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    interpolate_(propsDict_.lookupOrDefault<bool>("interpolation", false)),
    measureDiff_(propsDict_.lookupOrDefault<bool>("measureDiff", false)),
    recErrorFile_("recurrenceError"),
    ignoreCellsName_(propsDict_.lookupOrDefault<word>("ignoreCellsName","none")),
    ignoreCells_(),
    existIgnoreCells_(true),
    voidfractionFieldName_(propsDict_.lookupOrDefault<word>("voidfractionFieldName","voidfraction")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    voidfractionRecFieldName_(propsDict_.lookupOrDefault<word>("voidfractionRecFieldName","voidfractionRec")),
    voidfractionRec_(sm.mesh().lookupObject<volScalarField> (voidfractionRecFieldName_)),
    critVoidfraction_(propsDict_.lookupOrDefault<scalar>("critVoidfraction", 1.0)),
    D0_(propsDict_.lookupOrDefault<scalar>("D0",0.0)),
    D0FieldName_(propsDict_.lookupOrDefault<word>("D0FieldName","D0Field")),
    D0Field_
    (   IOobject
        (
            D0FieldName_,
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("D0", dimensionSet(0,0,0,0,0,0,0), D0_)
    ),
    dtDEM_(particleCloud_.dataExchangeM().DEMts()),
    ranGen_(osRandomInteger())
{
    if(ignoreCellsName_ != "none")
    {
       ignoreCells_.set(new cellSet(particleCloud_.mesh(),ignoreCellsName_));
       Info << "isotropicFluctuations: ignoring fluctuations in cellSet " << ignoreCells_().name() <<
        " with " << ignoreCells_().size() << " cells." << endl;
    }
    else existIgnoreCells_ = false;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

isotropicFluctuations::~isotropicFluctuations()
{
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //

bool isotropicFluctuations::ignoreCell(label cell) const
{
    if (!existIgnoreCells_) return false;
    else return ignoreCells_()[cell];
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void isotropicFluctuations::setForce() const
{
    vector position(0,0,0);
    scalar voidfraction(0.0);
    scalar voidfractionRec(0.0);
    scalar deltaVoidfrac(0.0);
    scalar D = 0.0;

    vector flucU(0,0,0);
    label cellI=0;
    scalar relVolfractionExcess(0.0);
   
    interpolationCellPoint<scalar> voidfractionInterpolator_(voidfraction_);
    interpolationCellPoint<scalar> voidfractionRecInterpolator_(voidfractionRec_);
    
    for(int index = 0;index <  particleCloud_.numberOfParticles(); ++index)
    {
            cellI = particleCloud_.cellIDs()[index][0];
            flucU =vector(0,0,0);
            voidfraction=0.0;
            voidfractionRec=0.0;
            deltaVoidfrac=0.0;

            if (cellI > -1 && !ignoreCell(cellI))
            {
                // particles in empty regions follow trajectories subject to gravity
                if(voidfractionRec_[cellI] < critVoidfraction_)
                {
                    if( interpolate_ )
                    {
                        position = particleCloud_.position(index);
                        voidfraction = voidfractionInterpolator_.interpolate(position,cellI);
                        voidfractionRec = voidfractionRecInterpolator_.interpolate(position,cellI);
                    }
                    else
                    {
                        voidfraction = voidfraction_[cellI];
                        voidfractionRec = voidfractionRec_[cellI];
                    }
                    // write particle based data to global array
                    
                    deltaVoidfrac=voidfractionRec-voidfraction;
                    relVolfractionExcess=deltaVoidfrac/(1-voidfraction+SMALL);
                    if(deltaVoidfrac>0)
                    {
                        D = D0Field_[cellI];
                        flucU=unitRndVec()*fluctuationMag(relVolfractionExcess,D);
                    }              

                    // write particle based data to global array
                    for(int j=0;j<3;j++)
                    {
                        particleCloud_.particleFlucVels()[index][j] += flucU[j];
                    }
                }
            }
    }
    
    if (measureDiff_)
    {
        dimensionedScalar diff( fvc::domainIntegrate( sqr( voidfraction_ - voidfractionRec_ ) ) );
        scalar t = particleCloud_.mesh().time().timeOutputValue(); 
        recErrorFile_ << t << "\t" << diff.value() << endl;
    }
}

scalar isotropicFluctuations::fluctuationMag(const scalar relVolfractionExcess, const scalar D) const
{
    // magnitude of steps dr = sqrt(6*D0*relVolfractionExcess * dt)
    // fluctuation velocity v = dr / dt
    scalar fluctuation;
    if(relVolfractionExcess<0.0 || D < 1e-10)
    {
        return 0.0;
    }
    else
    {
        fluctuation=Foam::sqrt(6*D*relVolfractionExcess/dtDEM_);
        return fluctuation;
    }
}

vector isotropicFluctuations::unitRndVec() const
{
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
