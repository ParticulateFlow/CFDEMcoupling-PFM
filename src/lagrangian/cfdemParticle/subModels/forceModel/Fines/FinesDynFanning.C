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

#include "FinesDynFanning.H"
#include "addToRunTimeSelectionTable.H"
#include "averagingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(FinesDynFanning, 0);

addToRunTimeSelectionTable
(
    forceModel,
    FinesDynFanning,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
FinesDynFanning::FinesDynFanning
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
    alphaDyn_(sm.mesh().lookupObject<volScalarField> ("alphaDyn")),
    alphaP_(sm.mesh().lookupObject<volScalarField> ("alphaP")),
    dHydMix_(sm.mesh().lookupObject<volScalarField> ("dHydMix")),
    dSauter_(sm.mesh().lookupObject<volScalarField> ("dSauter")),
    Froude_(sm.mesh().lookupObject<volScalarField> ("Froude")),
    rhoDyn_(readScalar(propsDict_.lookup ("rhoDyn"))),
    prefactor_(readScalar(propsDict_.lookup ("prefactor"))),
    exponent_(readScalar(propsDict_.lookup ("exponent"))),
    scaleDia_(1.),
    scaleDrag_(1.)
{
    // init force sub model
    setForceSubModels(propsDict_);
    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(0,true); // activate treatExplicit switch
    forceSubM(0).setSwitchesList(2,true); // activate implDEM switch
    forceSubM(0).setSwitchesList(3,true); // activate search for verbose switch
    forceSubM(0).readSwitches();
    forceSubM(0).setSwitches(0,true);

    particleCloud_.checkCG(true);
    if (propsDict_.found("scale"))
        scaleDia_=scalar(readScalar(propsDict_.lookup("scale")));
    if (propsDict_.found("scaleDrag"))
        scaleDrag_=scalar(readScalar(propsDict_.lookup("scaleDrag")));

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

FinesDynFanning::~FinesDynFanning()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void FinesDynFanning::setForce() const
{
    vector UDyn(0,0,0);
    vector drag(0,0,0);
    label cellI=0;

    vector Us(0,0,0);
    vector Ur(0,0,0);
    scalar ds(0);
    scalar ds_scaled(0);
    scalar scaleDia3 = scaleDia_*scaleDia_*scaleDia_;
   
    scalar magUr(0);
    scalar Fk(0);

    vector dragExplicit(0,0,0);
    scalar dragCoefficient(0);
    

    #include "setupProbeModel.H"

    for(int index = 0;index <  particleCloud_.numberOfParticles(); ++index)
    {
            cellI = particleCloud_.cellIDs()[index][0];
            drag = vector(0,0,0);
            dragExplicit = vector(0,0,0);
            UDyn = vector(0,0,0);
            dragCoefficient = 0;

            if (cellI > -1) // particle found
            {
                UDyn = UDyn_[cellI];
                Us = UsField_[cellI];
                Ur = UDyn-Us;
                magUr = mag(Ur);
                ds = 2*particleCloud_.radius(index);
		ds_scaled = ds/scaleDia_;
		
		Fk = prefactor_*Foam::pow(Froude_[cellI], exponent_);
		
                dragCoefficient = alphaDyn_[cellI] * rhoDyn_ * magUr * Fk / (2 * dHydMix_[cellI] )

                // calc particle's drag
                dragCoefficient *= M_PI/6 * ds_scaled * ds_scaled / alphaP_[cellI] * dSauter_[cellI] *scaleDia3*scaleDrag_;
                if (modelType_=="B")
                    dragCoefficient /= voidfraction_[cellI];

                drag = dragCoefficient * Ur;
            }

            // write particle based data to global array
            forceSubM(0).partToArray(index,drag,dragExplicit);
    }
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
