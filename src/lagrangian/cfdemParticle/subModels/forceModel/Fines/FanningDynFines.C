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

#include "FanningDynFines.H"
#include "addToRunTimeSelectionTable.H"
#include "averagingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(FanningDynFines, 0);

addToRunTimeSelectionTable
(
    forceModel,
    FanningDynFines,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
FanningDynFines::FanningDynFines
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    velFieldName_(propsDict_.lookup("velFieldName")),
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
    voidfractionFieldName_(propsDict_.lookup("voidfractionFieldName")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    UsFieldName_(propsDict_.lookup("granVelFieldName")),
    UsField_(sm.mesh().lookupObject<volVectorField> (UsFieldName_)),
    UDyn_(sm.mesh().lookupObject<volVectorField> ("uDyn")),
    FanningCoeff_(sm.mesh().lookupObject<volScalarField> ("FanningCoeff")),
    alphaP_(sm.mesh().lookupObject<volScalarField> ("alphaP")),
    dSauter_(sm.mesh().lookupObject<volScalarField> ("dSauter")),
    scaleDia_(1.),
    scaleDrag_(1.)
{
    // init force sub model
    setForceSubModels(propsDict_);
    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(SW_TREAT_FORCE_EXPLICIT,true); // activate treatExplicit switch
    forceSubM(0).setSwitchesList(SW_VERBOSE,true); // activate search for verbose switch
    forceSubM(0).readSwitches();
    forceSubM(0).setSwitches(SW_TREAT_FORCE_EXPLICIT,true);

    particleCloud_.checkCG(true);
    if (propsDict_.found("scale"))
        scaleDia_ = scalar(readScalar(propsDict_.lookup("scale")));
    if (propsDict_.found("scaleDrag"))
        scaleDrag_ = scalar(readScalar(propsDict_.lookup("scaleDrag")));

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

FanningDynFines::~FanningDynFines()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void FanningDynFines::setForce() const
{
    if(forceSubM(0).verbose())
        Info << "Entering force loop of FanningDynFines.\n" << endl;

    if (scaleDia_ > 1)
    {
        Info << "FanningDynFines using scale = " << scaleDia_ << endl;
    }
    else if (particleCloud_.cg() > 1)
    {
        scaleDia_=particleCloud_.cg();
        Info << "FanningDynFines using scale from liggghts cg = " << scaleDia_ << endl;
    }

    vector UDyn(0,0,0);
    vector drag(0,0,0);
    label cellI = 0;

    vector Us(0,0,0);
    vector Ur(0,0,0);
    scalar ds(0);
    scalar ds_scaled(0);
    scalar scaleDia3 = scaleDia_*scaleDia_*scaleDia_;

    scalar dragCoefficient(0);


    #include "setupProbeModel.H"

    for(int index = 0; index < particleCloud_.numberOfParticles(); ++index)
    {
            cellI = particleCloud_.cellIDs()[index][0];
            drag = vector::zero;
            UDyn = vector::zero;
            dragCoefficient = 0;

            if (cellI > -1) // particle found
            {
                UDyn = UDyn_[cellI];
                Us = UsField_[cellI];
                Ur = UDyn-Us;
                ds = 2 * particleCloud_.radius(index);
                ds_scaled = ds/scaleDia_;

                dragCoefficient = FanningCoeff_[cellI];

                // calc particle's drag
                dragCoefficient *= M_PI/6 * ds_scaled * ds_scaled / alphaP_[cellI] * dSauter_[cellI] * scaleDia3 * scaleDrag_;
                if (modelType_ == "B")
                    dragCoefficient /= voidfraction_[cellI];

                drag = dragCoefficient * Ur;
            }

            // write particle based data to global array
            forceSubM(0).partToArray(index,drag,vector::zero);
    }

    if (forceSubM(0).verbose())
        Info << "Leaving force loop of FanningDynFines.\n" << endl;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
