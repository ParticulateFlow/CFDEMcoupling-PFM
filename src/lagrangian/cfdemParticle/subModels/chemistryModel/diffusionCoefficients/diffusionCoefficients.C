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
#include "diffusionCoefficients.H"
#include "addToRunTimeSelectionTable.H"

#include "dataExchangeModel.H"
#include "IFstream.H"

#define SMALL 1e-7
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(diffusionCoefficient, 0);

addToRunTimeSelectionTable
(
        chemistryModel,
        diffusionCoefficient,
        dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
diffusionCoefficient::diffusionCoefficient
(
    const dictionary& dict,
    cfdemCloudEnergy& sm
)
:
    chemistryModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    verbose_(propsDict_.lookupOrDefault<bool>("verbose",false)),
    interpolation_(propsDict_.lookupOrDefault<bool>("interpolation",false)),
    mesh_(sm.mesh()),
    // define a file name in the coupling properties that contains the species
    specDict_
    (
        IFstream
        (
            fileName(propsDict_.lookup("ChemistryFile")).expand()
        )()
    ),
    // create a list from the Species table in the specified diffusionCoefficient dictionary
    speciesNames_(specDict_.lookup("species")),
    tempFieldName_(propsDict_.lookup("tempFieldName")),
    tempField_(sm.mesh().lookupObject<volScalarField> (tempFieldName_)),
    pressureFieldName_(propsDict_.lookup("pressureFieldName")),
    P_(sm.mesh().lookupObject<volScalarField>(pressureFieldName_)),
    densityFieldName_(propsDict_.lookup("densityFieldName")),
    rho_(sm.mesh().lookupObject<volScalarField> (densityFieldName_)),
    totalMoleFieldName_(propsDict_.lookup("totalMoleFieldName")),
    // needed to calculate the mixture diffusion coefficient
    // dcoeff is dependent on molar fraction not mass fraction
    N_(sm.mesh().lookupObject<volScalarField>(totalMoleFieldName_)),
    Y_(speciesNames_.size()),
    // X_(speciesNames_.size()),
    diffusantGasNames_(propsDict_.lookup("diffusantGasNames")),
    diffusionCoefficients_(diffusantGasNames_.size(),NULL)
 /*   diffusionCoefficientNames_(propsDict_.lookup("diffusionCoefficientNames")),
                               //volumeScalarFields created in the ts folders
    diffusionCoefficients_(diffusionCoefficientNames_.size(),NULL),         //the value of diffusionCoefficient concentration for every diffusionCoefficient 
    */
{
    Info << " Reading diffusionCoefficient list: " << diffusantGasNames_ << endl;
    for (int i = 0; i < diffusantGasNames_.size(); i++)
        Info << " Diffusant names: " << diffusantGasNames_[i] << endl;

    for (int i=0; i<speciesNames_.size(); i++)
    {
        volScalarField& Y = const_cast<volScalarField&>
                (sm.mesh().lookupObject<volScalarField>(speciesNames_[i]));
        Y_.set(i, &Y);
        particleCloud_.checkCG(false);
    }

    // if mole fractions field should be generated ??
    /*forAll(Y, i)
    {
        X_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "X_" + Y_[i].name(),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("X", dimless, 0)
            )
        );
    } */

    allocateMyArrays();
    createCoeffs();
    molWeightTable();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

diffusionCoefficient::~diffusionCoefficient()
{
    coeffs.clearStorage();
    molWeight.clearStorage();

    int nP_ = particleCloud_.numberOfParticles();
    for (int i=0; i<diffusantGasNames_.size(); i++)
    {
        particleCloud_.dataExchangeM().destroy(diffusionCoefficients_[i],nP_);
    }
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //
 void diffusionCoefficient::allocateMyArrays() const
{
    double initVal=0.0;
    if (particleCloud_.dataExchangeM().maxNumberOfParticles() > 0)
    {
        for (int i=0; i<diffusantGasNames_.size(); i++)
        {
            particleCloud_.dataExchangeM().allocateArray(diffusionCoefficients_[i],initVal,1,"nparticles");

        }
    }
}

void diffusionCoefficient::reAllocMyArrays() const
{
    if (particleCloud_.numberOfParticlesChanged())
    {
        double initVal=0.0;

        for (int i=0; i<diffusantGasNames_.size(); i++)
        {
            particleCloud_.dataExchangeM().allocateArray(diffusionCoefficients_[i],initVal,1);
        }
    }
}

// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //

void diffusionCoefficient::execute()
{
    // realloc the arrays
    reAllocMyArrays();

    label  cellI=0;
    scalar Tfluid(0);
    scalar rhofluid(0);
    List<scalar> Yfluid_;
    Yfluid_.setSize(speciesNames_.size());
    scalar Pfluid(0);
    scalar Nfluid(0);
    scalar Texp(0);
    List<scalar> Xfluid_;
    Xfluid_.setSize(speciesNames_.size());
    List<scalar> dBinary_;
    dBinary_.setSize(diffusantGasNames_.size());
    dCoeff_.setSize(diffusantGasNames_.size());

    double **molNum_ = new double*[diffusantGasNames_.size()];
    double **volDiff_ = new double*[diffusantGasNames_.size()];

    // defining interpolators for T, rho, voidfraction, N
    interpolationCellPoint <scalar> TInterpolator_(tempField_);
    interpolationCellPoint <scalar> rhoInterpolator_(rho_);
    interpolationCellPoint <scalar> PInterpolator_(P_);
    interpolationCellPoint <scalar> NInterpolator_(N_);

    for (int index=0; index<particleCloud_.numberOfParticles(); index++)
    {
        cellI=particleCloud_.cellIDs()[index][0];
        if (cellI >=0)
        {
            if(interpolation_)
            {
                vector position     =   particleCloud_.position(index);
                Tfluid              =   TInterpolator_.interpolate(position,cellI);
                rhofluid            =   rhoInterpolator_.interpolate(position,cellI);
                Pfluid              =   PInterpolator_.interpolate(position,cellI);
                Nfluid              =   NInterpolator_.interpolate(position,cellI);
            }
            else
            {
                Tfluid          =   tempField_[cellI];
                Texp            =   pow(Tfluid,1.75);

                rhofluid        =   rho_[cellI];
                Pfluid          =   P_[cellI];
                Nfluid          =   N_[cellI];

                for (int i = 0; i<speciesNames_.size();i++)
                {
                    Yfluid_[i] = Y_[i][cellI];
                }
            }

            for (int i=0; i<diffusantGasNames_.size();i++)
            {
                molNum_[i]  =   new double [speciesNames_.size()];
                volDiff_[i]  =   new double [speciesNames_.size()];
                for (int j=0; j < speciesNames_.size();j++)
                {
                    if (diffusantGasNames_[i] != speciesNames_[j])
                    {
                        Info << "molar weights diffuser gases: " << molWeight(speciesNames_[j]) << nl << endl;
                        Info << "N fluid" << Nfluid << nl << endl;
                        Info << "rho fluid" << rhofluid << nl << endl;
                        Info << "Y fluid" << Yfluid_[j] << nl << endl;
                        Info << "molar weights diffusant gases: " << molWeight(diffusantGasNames_[i]) << nl << endl;

                        // Nfluid is 0 for first ts, division by zero is no good...
                        if (Nfluid == 0)
                        {
                            Xfluid_[j]  =   0.0;
                            Xfluid_[i]  =   0.0;
                        } else
                        {
                            // convert mass to molar fractions
                            Xfluid_[j]  =   Yfluid_[j]*rhofluid/(Nfluid*molWeight(speciesNames_[j]));
                            Xfluid_[i]  =   Yfluid_[i]*rhofluid/(Nfluid*molWeight(diffusantGasNames_[i]));
                            Info << "molar fraction diffuser gases:" << Xfluid_[j] << nl << endl;
                            Info << "molar fraction diffusant gases:" << Xfluid_[i] << nl << endl;
                        }

                        if (Xfluid_[j] > 0.0)
                        {
                            calcMolNum(i,j,molNum_);
                            calcDiffVol(i,j,volDiff_);

                            // convert mass to mole fraction
                            // if ( i != j) but checks speciesPairs anyways so not needed.
                            if(coeffs.found(diffusantGasNames_[i]) && coeffs.found(speciesNames_[j]))
                            {
                                dBinary_[i] = 0.001*Texp*molNum_[i][j]/(Pfluid*volDiff_[i][j]);
                                Info << "dBinary: "  << dBinary_[i] << nl << endl;
                                // According to literature i.e Valipour 2006, Elnashaie et al. 1993, Taylor and Krishna (1993), Natsui et al.
                                // dCoeff = 1/(1-X[j])*sum(X[i]/D_[i,j])^-1
                                // X is molar fraction / Dij binary diff coeff.
                                dCoeff_[i]  +=  (Xfluid_[j]/dBinary_[i]);
                                dCoeff_[i]  =   (1-Xfluid_[i])*(1/dCoeff_[i]);
                                // According to Maier (who referred to wilke but in Wilke's paper its written as the previous eq.)
                                // and Nietros
                                // dCoeff_[i]   +=   Xfluid_[i]/dCoeff_[i];
                                // dCoeff_[i]   =   1/dCoeff_[i];
                                Info << "dCoeff: " << dCoeff_[i] << nl << endl;
                            }else
                            {
                                FatalError
                                        << "check tables for species diffusion volume"
                                        << endl
                                        << abort(FatalError);
                            }
                        }
                    }
                }
                diffusionCoefficients_[i][index][0]= dCoeff_[i];
            }
        }

        if(verbose_ && index >=0 && index < 2)
        {
            for(int i =0; i<diffusantGasNames_.size();i++)
            {
                Info << "effective diffusionCoefficient of species " << diffusantGasNames_[i] << " = " << diffusionCoefficients_[i][index][0] << endl;
            }
        }

    }

    for (int i=0; i<diffusantGasNames_.size();i++)
    {
        word pushName = diffusantGasNames_[i] + "_diffCoeff";
        particleCloud_.dataExchangeM().giveData(pushName,"scalar-atom",diffusionCoefficients_[i]);
    };

    Info << "give data done" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void diffusionCoefficient::createCoeffs()
{
    // diffusion volume of species i^(0.33333)
    coeffs.insert("CO", 2.6636);
    coeffs.insert("CO2", 2.9962);
    coeffs.insert("O2", 2.5509);
    coeffs.insert("N2", 2.6158);
    coeffs.insert("H2", 1.919);
    coeffs.insert("H2O", 2.3331);
    coeffs.insert("Air", 2.71893);
    coeffs.insert("Ne", 1.774750);
    coeffs.insert("N2O", 3.29886);
    coeffs.insert("NH3", 2.4607);
    coeffs.insert("H", 1.255707);
    coeffs.insert("O", 1.76303);
    coeffs.insert("C", 2.54582);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void diffusionCoefficient::calcDiffVol(int i, int j, double **volDiff_)
{
    volDiff_[i][j]  =   coeffs(diffusantGasNames_[i])+coeffs(speciesNames_[j]);
    volDiff_[i][j]  *=   volDiff_[i][j];
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void diffusionCoefficient::molWeightTable()
{
    // table for molecular weights
    molWeight.insert("CO", 28.0101);
    molWeight.insert("CO2", 44.01);
    molWeight.insert("O2", 32.00);
    molWeight.insert("N2", 28.01);
    molWeight.insert("H2", 2.02);
    molWeight.insert("H2O", 2.02);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// either calculate Molecular Weight addition (eq. D_ij) or consturct hashtable with diffusant and fifuser species
void diffusionCoefficient::calcMolNum(int i, int j, double **molNum_)
{
    double& W1 = molWeight(diffusantGasNames_[i]);
    double& W2 = molWeight(speciesNames_[j]);

    molNum_[i][j]   =   (W1 * W2) / (W1 + W2);
    molNum_[i][j]   =   sqrt(molNum_[i][j]);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
