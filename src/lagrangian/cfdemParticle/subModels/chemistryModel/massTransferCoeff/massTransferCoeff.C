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
                        M.Efe Kinaci, JKU Linz, Austria

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "massTransferCoeff.H"
#include "addToRunTimeSelectionTable.H"

#include "dataExchangeModel.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(massTransferCoeff, 0);

addToRunTimeSelectionTable
(
        chemistryModel,
        massTransferCoeff,
        dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
massTransferCoeff::massTransferCoeff
(
    const dictionary& dict,
    cfdemCloudEnergy& sm
)
:
    chemistryModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    interpolation_(propsDict_.lookupOrDefault<bool>("interpolation",false)),
    mesh_(sm.mesh()),
    velFieldName_(propsDict_.lookup("velFieldName")),
    U_(mesh_.lookupObject<volVectorField>(velFieldName_)),
    voidfractionFieldName_(propsDict_.lookup("voidfractionFieldName")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    partNuName_(propsDict_.lookup("partNuName")),
    partNu_(NULL),
    partReynolds_(propsDict_.lookup("partReynolds")),
    Rep_(NULL)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

massTransferCoeff::~massTransferCoeff()
{
    int nP_ = particleCloud_.numberOfParticles();

    particleCloud_.dataExchangeM().destroy(partNu_,nP_);
    particleCloud_.dataExchangeM().destroy(Rep_,nP_);
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //
void massTransferCoeff::allocateMyArrays() const
{
    double initVal=0.0;
    if (particleCloud_.dataExchangeM().maxNumberOfParticles() > 0)
    {
        // get memory for 2d arrays
        particleCloud_.dataExchangeM().allocateArray(partNu_,initVal,1,"nparticles");
        particleCloud_.dataExchangeM().allocateArray(Rep_,initVal,1,"nparticles");
    }
}

void massTransferCoeff::reAllocMyArrays() const
{
    if (particleCloud_.numberOfParticlesChanged())
    {
        double initVal=0.0;
        particleCloud_.dataExchangeM().allocateArray(partNu_,initVal,1);
        particleCloud_.dataExchangeM().allocateArray(Rep_,initVal,1);
    }
}

// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //

void massTransferCoeff::execute()
{
    // realloc the arrays
    reAllocMyArrays();

    const volScalarField& nufField_   =   particleCloud_.turbulence().nu();

    label  cellI=0;

    // use to calculate Rep
    scalar voidfraction(1);
    vector Ufluid(0,0,0);
    vector Us(0,0,0);
    vector Ur(0,0,0);
    scalar ds(0);
    scalar magUr(0);

    // give Rep & kin. visc. to DEM
    scalar Rep(0);
    scalar nuf(0);

    // defining interpolators for U, voidfraction
    interpolationCellPoint <vector> UfluidInterpolator_(U_);
    interpolationCellPoint <scalar> voidfractionInterpolator_(voidfraction_);

    for (int index=0; index<particleCloud_.numberOfParticles(); index++)
    {
        cellI=particleCloud_.cellIDs()[index][0];
        if (cellI >=0)
        {
            if(interpolation_)
            {
                 vector position    =   particleCloud_.position(index);
                 Ufluid             =   UfluidInterpolator_.interpolate(position,cellI);
                 voidfraction       =   voidfractionInterpolator_.interpolate(position,cellI);
            }
            else
            {
                Ufluid          =   U_[cellI];
                voidfraction    =   voidfraction_[cellI];
            }

            nuf =   nufField_[cellI];
            Us  =   particleCloud_.velocity(index);
            Ur  =   Ufluid  -   Us;
            ds  =   2*particleCloud_.radius(index);
            magUr   =   mag(Ur);

            // calculate particle Reynolds number
            Rep =   ds*voidfraction*magUr/(nuf+SMALL);

            // fill arrays
            partNu_[index][0]   =   nuf;
            Rep_[index][0]      =   Rep;
        }

        //if (particleCloud_.verbose() && index >=0 && index < 2)
        //{
            Info << "Nufield = " << nuf << endl;
            Info << "Rep = " << Rep << endl;
        //}
    }

    Info << "give data done" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
