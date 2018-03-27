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

#include "ErgunStatFines.H"
#include "addToRunTimeSelectionTable.H"
#include "averagingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(ErgunStatFines, 0);

addToRunTimeSelectionTable
(
    forceModel,
    ErgunStatFines,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
ErgunStatFines::ErgunStatFines
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
    dSauter_(sm.mesh().lookupObject<volScalarField> ("dSauter")),
    dSauterMix_(sm.mesh().lookupObject<volScalarField> ("dSauterMix")),
    alphaP_(sm.mesh().lookupObject<volScalarField> ("alphaP")),
    alphaSt_(sm.mesh().lookupObject<volScalarField> ("alphaSt")),
    phi_(readScalar(propsDict_.lookup("phi"))),
    UsFieldName_(propsDict_.lookup("granVelFieldName")),
    UsField_(sm.mesh().lookupObject<volVectorField> (UsFieldName_)),
    scaleDia_(1.),
    scaleDist_(1.),
    scaleDrag_(1.),
    switchingVoidfraction_(0.8)
{
    //Append the field names to be probed
    particleCloud_.probeM().initialize(typeName, typeName+".logDat");
    particleCloud_.probeM().vectorFields_.append("dragForce"); //first entry must  be the force
    particleCloud_.probeM().vectorFields_.append("Urel");
    particleCloud_.probeM().scalarFields_.append("Rep");
    particleCloud_.probeM().scalarFields_.append("betaP");
    particleCloud_.probeM().scalarFields_.append("voidfraction");
    particleCloud_.probeM().writeHeader();

    // init force sub model
    setForceSubModels(propsDict_);
    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(SW_TREAT_FORCE_EXPLICIT,true); // activate treatExplicit switch
    forceSubM(0).setSwitchesList(SW_IMPL_FORCE_DEM,true); // activate implDEM switch
    forceSubM(0).setSwitchesList(SW_VERBOSE,true); // activate search for verbose switch
    forceSubM(0).setSwitchesList(SW_INTERPOLATION,true); // activate search for interpolate switch
    forceSubM(0).setSwitchesList(SW_SCALAR_VISCOSITY,true); // activate scalarViscosity switch
    forceSubM(0).readSwitches();

    particleCloud_.checkCG(true);
    if (propsDict_.found("scale"))
        scaleDia_ = scalar(readScalar(propsDict_.lookup("scale")));
    if (propsDict_.found("scaleDrag"))
        scaleDrag_ = scalar(readScalar(propsDict_.lookup("scaleDrag")));

    if (propsDict_.found("switchingVoidfraction"))
        switchingVoidfraction_ = readScalar(propsDict_.lookup("switchingVoidfraction"));

    dictionary SauterDict(dict.subDict("dSauterProps"));
    if (SauterDict.found("scaleDist"))
        scaleDist_ = scalar(readScalar(SauterDict.lookup("scaleDist")));

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

ErgunStatFines::~ErgunStatFines()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
scalar ErgunStatFines::dSauter(label cellI) const
{
    // Sauter mean diameter without influence of medium-scale fines
    scalar dS = dSauter_[cellI] / scaleDist_;
    return dS;
}

void ErgunStatFines::setForce() const
{
    if (scaleDia_ > 1)
    {
        Info << "ErgunStatFines using scale = " << scaleDia_ << endl;
    }
    else if (particleCloud_.cg() > 1)
    {
        scaleDia_=particleCloud_.cg();
        Info << "ErgunStatFines using scale from liggghts cg = " << scaleDia_ << endl;
    }

    const volScalarField& nufField = forceSubM(0).nuField();
    const volScalarField& rhoField = forceSubM(0).rhoField();

    vector position(0,0,0);
    scalar voidfraction(1);
    vector Ufluid(0,0,0);
    vector drag(0,0,0);
    label cellI = 0;

    vector Us(0,0,0);
    vector Ur(0,0,0);
    scalar ds(0);
    scalar dSauterMix(0.0);
    scalar nuf(0);
    scalar rho(0);
    scalar magUr(0);
    scalar Rep(0);
    scalar alphaPartEff(0);

    scalar CdMagUrLag(0);       //Cd of the very particle
    scalar betaP(0);            //momentum exchange of the very particle

    vector dragExplicit(0,0,0);
    scalar dragCoefficient(0);

    scalar scaleDia3 = scaleDia_*scaleDia_*scaleDia_;

    interpolationCellPoint<scalar> voidfractionInterpolator_(voidfraction_);
    interpolationCellPoint<vector> UInterpolator_(U_);

    #include "setupProbeModel.H"

    if(forceSubM(0).verbose())
                Info << "Entering force loop of ErgunStatFines.\n" << endl;

    for(int index = 0; index < particleCloud_.numberOfParticles(); ++index)
    {
            cellI = particleCloud_.cellIDs()[index][0];
            drag = vector::zero;
            dragExplicit = vector::zero;
            betaP = 0;
            Ufluid = vector::zero;
            voidfraction = 0;
            dragCoefficient = 0;

            if (cellI > -1) // particle found
            {

                if( forceSubM(0).interpolation() )
                {
                    position     = particleCloud_.position(index);
                    voidfraction = voidfractionInterpolator_.interpolate(position,cellI);
                    Ufluid       = UInterpolator_.interpolate(position,cellI);
                }
                else
                {
                    voidfraction = voidfraction_[cellI];
                    Ufluid = U_[cellI];
                }

                // ensure voidfraction to be meaningful
                // problems could arise from interpolation or empty cells

                if(voidfraction > 0.999)
                    voidfraction = 0.999;
                else if(voidfraction < 0.05)
                    voidfraction = 0.05;

                Us = particleCloud_.velocity(index);
                Ur = Ufluid-Us;
                magUr = mag(Ur);
                dSauterMix = dSauterMix_[cellI];
                ds = 2*particleCloud_.radius(index);
                rho = rhoField[cellI];
                nuf = nufField[cellI];

                Rep = 0.0;
                alphaPartEff = 1.0 - voidfraction + alphaSt_[cellI] + SMALL;

                // calc particle's drag coefficient (i.e., Force per unit slip velocity and per m³ PARTICLE)
                if(voidfraction > switchingVoidfraction_) //dilute, no static hold-up present
                {
                    Rep=dSauterMix*voidfraction*magUr/nuf;
                    CdMagUrLag = (24.0*nuf/(dSauterMix*voidfraction)) //1/magUr missing here, but compensated in expression for betaP!
                                 *(scalar(1.0)+0.15*Foam::pow(Rep, 0.687));

                    betaP = 0.75* alphaPartEff * (
                                            rho*voidfraction*CdMagUrLag
                                          /
                                            (dSauterMix*Foam::pow(voidfraction,2.65))
                                          );
                }
                else  //dense
                {
                    betaP = (150 * alphaPartEff * alphaPartEff *nuf*rho)
                             /  ((1-alphaPartEff) * dSauterMix*phi_*dSauterMix*phi_)
                            +
                              (1.75 * magUr * rho * alphaPartEff)
                             /((dSauterMix*phi_));
                }

                // calc particle's drag
                betaP /= (1-alphaPartEff);
                dragCoefficient = M_PI/6 * ds/scaleDia_ * ds/scaleDia_ * dSauter(cellI) * voidfraction / (1 - voidfraction) * betaP * scaleDrag_;
                dragCoefficient *= scaleDia3;
                if (modelType_ == "B")
                    dragCoefficient /= voidfraction;

                drag = dragCoefficient * Ur;

                // explicitCorr
                forceSubM(0).explicitCorr(drag,dragExplicit,dragCoefficient,Ufluid,U_[cellI],Us,UsField_[cellI],forceSubM(0).verbose());

                if(forceSubM(0).verbose() && index >= 0 && index < 2)
                {
                    Pout << "cellI = " << cellI << endl;
                    Pout << "index = " << index << endl;
                    Pout << "Us = " << Us << endl;
                    Pout << "Ur = " << Ur << endl;
                    Pout << "ds = " << ds << endl;
                    Pout << "ds/scale = " << ds/scaleDia_ << endl;
                    Pout << "phi = " << phi_ << endl;
                    Pout << "rho = " << rho << endl;
                    Pout << "nuf = " << nuf << endl;
                    Pout << "voidfraction = " << voidfraction << endl;
                    Pout << "Rep = " << Rep << endl;
                    Pout << "betaP = " << betaP << endl;
                    Pout << "drag = " << drag << endl;
                }

                //Set value fields and write the probe
                if(probeIt_)
                {
                    #include "setupProbeModelfields.H"
                    vValues.append(drag);   //first entry must the be the force
                    vValues.append(Ur);
                    sValues.append(Rep);
                    sValues.append(betaP);
                    sValues.append(voidfraction);
                    particleCloud_.probeM().writeProbe(index, sValues, vValues);
                }
            }

            // write particle based data to global array
            forceSubM(0).partToArray(index,drag,dragExplicit,Ufluid,dragCoefficient);

    }// end loop particles

    if(forceSubM(0).verbose())
        Pout << "Leaving force loop of ErgunStatFines.\n" << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
