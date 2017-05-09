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
    totalMoleFieldName_(propsDict_.lookup("totalMoleFieldName")),
    // needed to calculate the mixture diffusion coefficient
    // dcoeff is dependent on molar fraction not mass fraction
    N_(sm.mesh().lookupObject<volScalarField>(totalMoleFieldName_)),
    Y_(speciesNames_.size()),
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

    allocateMyArrays();
    createCoeffs();
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
    List<scalar> Yfluid_;
    Yfluid_.setSize(speciesNames_.size());
    scalar Pfluid(0);
    scalar Nfluid(0);
    
    scalar dCoeff(0.0);
    
    // word speciesPair("none");

    // defining interpolators for T, rho, voidfraction, N
    interpolationCellPoint <scalar> TInterpolator_(tempField_);
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
                Pfluid              =   PInterpolator_.interpolate(position,cellI);
                Nfluid              =   NInterpolator_.interpolate(position,cellI);
            }
            else
            {
                Tfluid          =   tempField_[cellI];
                Pfluid          =   P_[cellI];
                Nfluid          =   N_[cellI];

                for (int i = 0; i<speciesNames_.size();i++)
                {
                    Yfluid_[i] = Y_[i][cellI];
                }
            }

            /*for (int i=0; i<diffusantGasNames_.size();i++)
            {
                // do the calculation
                dCoeff = 0.0;
                for (int j=0; j < speciesNames_.size();j++)
                {
                    // speciesPair = diffusantGasNames_[i] + "_" + speciesNames_[j];
                    // According to literature i.e Valipour 2006, Elnashaie et al. 1993, Taylor and Krishna (1993), Natsui et al.
                    // dCoeff = (1-X[j])*sum(X[i]/D_[i,j])
                    // X is molar fraction / Dij binary diff coeff.

                     if(coeffs.found(speciesPair))
                    {
                        dCoeff += Y[j] / coeffs.find(speciesPair)();
                    }
                }
                // diffusionCoefficients_[i][index][0]= *1.0/dCoeff;
            } */
        }

        if(particleCloud_.verbose() && index >=0 && index < 2)
        {
            for(int i =0; i<diffusantGasNames_.size();i++)
            {
                Info << "effective diffusionCoefficient of species " << diffusantGasNames_[i] << " = " << diffusionCoefficients_[i][index][0] << endl;
            }
        }
    }

    /*for (int i=0; i<diffusionCoefficientNames_.size();i++)
    {
        word pushName = diffusantGasNames_[i] + "_diffCoeff";
        particleCloud_.dataExchangeM().giveData(pushName,"scalar-atom",diffusionCoefficients_[i]);
    };*/

    Info << "give data done" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void diffusionCoefficient::createCoeffs()
{
    // add all relevant combinations
    coeffs.insert("CO", 18.9);
    coeffs.insert("CO2", 26.9);
    coeffs.insert("O2", 16.6);
    coeffs.insert("N2", 17.9);
    coeffs.insert("H2", 7.07);
    coeffs.insert("H2O", 12.7);

    // coeffs for pairs (Va^(1/3)+Vb^(1/3))
    coeffs.insert("CO_CO2", 5.66);
    coeffs.insert("H2_H2O", 4.25239);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void diffusionCoefficient::molWeightTable()
{
    // table for molecular weights
    molWeight.insert("CO", 28.01);
    molWeight.insert("CO2", 44.01);
    molWeight.insert("O2", 32.00);
    molWeight.insert("N2", 28.01);
    molWeight.insert("H2", 2.02);
    molWeight.insert("H2O", 2.02);

    //Molecular Weight eq. solution in D_ij for species pairs (reactant-product)
    molWeight.insert("CO_CO2",0.2417);
    molWeight.insert("H2_H2O",0.74198);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// either calculate Molecular Weight addition (eq. D_ij) or consturct hashtable with diffusant and fifuser species
/*void diffusionCoefficient::calcMolNum(int i, int j, double *molNum_)
{
    molNum_ = (1/molWeight.find(diffusantGasNames_[i])+1/molWeight.find(speciesNames_[j]));
    molNum_ = pow(molNum_,0.5);
} */

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// add dummy volScalarFields, used in YEqn as Smi
tmp <volScalarField> diffusionCoefficient::Smi(const label i) const
{
    PtrList<volScalarField> dummy_;
    dummy_.set
            (
            i,
            new volScalarField
            (
                IOobject
                (
                    "empty1",
                    particleCloud_.mesh().time().timeName(),
                    particleCloud_.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                particleCloud_.mesh(),
                dimensionedScalar("zero",dimless,0.0)
            )
        );
    return tmp<volScalarField> (dummy_[i]);
}
// add dummy volScalarFields, used in YEqn as Sm
tmp <volScalarField> diffusionCoefficient::Sm() const
{
    tmp<volScalarField> dummy
    (
            new volScalarField
            (
                IOobject
                (
                    "empty2",
                    particleCloud_.mesh().time().timeName(),
                    particleCloud_.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                    particleCloud_.mesh(),
                    dimensionedScalar("zero",dimless,0.0)
            )
        );
    return dummy;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
