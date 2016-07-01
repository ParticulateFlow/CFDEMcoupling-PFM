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

#include "Fines.H"
#include "addToRunTimeSelectionTable.H"
#include "averagingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Fines, 0);

addToRunTimeSelectionTable
(
    forceModel,
    Fines,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Fines::Fines
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    finesFields_(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    velFieldName_(propsDict_.lookup("velFieldName")),
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
    voidfractionFieldName_(propsDict_.lookup("voidfractionFieldName")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    UsFieldName_(propsDict_.lookup("granVelFieldName")),
    UsField_(sm.mesh().lookupObject<volVectorField> (UsFieldName_)),
    scaleDia_(1.),
    scaleDrag_(1.)
{
    //Append the field names to be probed
    particleCloud_.probeM().initialize(typeName, "Fines.logDat");
    particleCloud_.probeM().vectorFields_.append("dragForce"); //first entry must  be the force
    particleCloud_.probeM().vectorFields_.append("Urel");
    particleCloud_.probeM().scalarFields_.append("Rep");
    particleCloud_.probeM().scalarFields_.append("betaP");
    particleCloud_.probeM().scalarFields_.append("voidfraction");
    particleCloud_.probeM().writeHeader();

    // init force sub model
    setForceSubModels(propsDict_);
    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(0,true); // activate treatExplicit switch
    forceSubM(0).setSwitchesList(2,true); // activate implDEM switch
    forceSubM(0).setSwitchesList(3,true); // activate search for verbose switch
    forceSubM(0).setSwitchesList(4,true); // activate search for interpolate switch
    forceSubM(0).setSwitchesList(8,true); // activate scalarViscosity switch
    forceSubM(0).readSwitches();

    particleCloud_.checkCG(true);
    if (propsDict_.found("scale"))
        scaleDia_=scalar(readScalar(propsDict_.lookup("scale")));
    if (propsDict_.found("scaleDrag"))
        scaleDrag_=scalar(readScalar(propsDict_.lookup("scaleDrag")));

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Fines::~Fines()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Fines::setForce() const
{
    finesFields_.update();

    

    vector position(0,0,0);
    scalar voidfraction(1);
    vector UDyn(0,0,0);
    vector drag(0,0,0);
    label cellI=0;

    vector Us(0,0,0);
    vector Ur(0,0,0);
    scalar ds(0);
    scalar ds_scaled(0);
    scalar scaleDia3 = scaleDia_*scaleDia_*scaleDia_;
   
    scalar magUr(0);
    scalar Froude(0);


    vector dragExplicit(0,0,0);
    scalar dragCoefficient(0);
    

    #include "setupProbeModel.H"

    for(int index = 0;index <  particleCloud_.numberOfParticles(); ++index)
    {
            cellI = particleCloud_.cellIDs()[index][0];
            drag = vector(0,0,0);
            dragExplicit = vector(0,0,0);
            UDyn = vector(0,0,0);
            voidfraction=0;
            dragCoefficient = 0;

            if (cellI > -1) // particle found
            {
                voidfraction = voidfraction_[cellI];
                UDyn = finesFields_.uDyn(cellI);
                Us = UsField_.[cellI];
                Ur = UDyn-Us;
                magUr = mag(Ur);
                ds = 2*particleCloud_.radius(index);
		ds_scaled = ds/scaleDia_;
               
		// scalinge!!
                Froude=finesFields_.Froude(cellI);
                ds_scaled*voidfraction*magUr/nuf+SMALL;
                dragCoefficient = 

                // calc particle's drag
                dragCoefficient *= 3*M_PI*ds_scaled*nuf*rho*voidfraction*scaleDia3*scaleDrag_;
                if (modelType_=="B")
                    dragCoefficient /= voidfraction;

                drag = dragCoefficient * Ur;

                // explicitCorr
                forceSubM(0).explicitCorr(drag,dragExplicit,dragCoefficient,Ufluid,U_[cellI],Us,UsField_[cellI],forceSubM(0).verbose());
            }

            // write particle based data to global array
            forceSubM(0).partToArray(index,drag,dragExplicit,Ufluid,dragCoefficient);
    }
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
