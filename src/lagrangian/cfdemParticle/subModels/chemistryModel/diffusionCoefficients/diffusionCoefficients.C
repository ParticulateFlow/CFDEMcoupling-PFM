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

#define VSMALL 1e-15
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
    specDict_
    (
        IFstream
        (
            fileName(propsDict_.lookup("ChemistryFile")).expand()
        )()
    ),
    // create a list from the Species table in the specified diffusionCoefficient dictionary
    speciesNames_(specDict_.lookup("species")),
    tempFieldName_(propsDict_.lookupOrDefault<word>("tempFieldName","T")),
    tempField_(sm.mesh().lookupObject<volScalarField> (tempFieldName_)),
    pressureFieldName_(propsDict_.lookupOrDefault<word>("pressureFieldName","p")),
    P_(sm.mesh().lookupObject<volScalarField>(pressureFieldName_)),
    partPressureName_(propsDict_.lookupOrDefault<word>("partPressureName","partP")),
    partPressure_(NULL),
    X_(speciesNames_.size()),
    diffusantGasNames_(propsDict_.lookup("diffusantGasNames")),
    diffusionCoefficients_(diffusantGasNames_.size(),NULL),
    Xdiffusant_(diffusantGasNames_.size()),
    initialized_(false)
{
    particleCloud_.checkCG(false);
    allocateMyArrays();
    createCoeffs();
    molWeightTable();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

diffusionCoefficient::~diffusionCoefficient()
{
    particleCloud_.dataExchangeM().destroy(partPressure_,1);
    for (int i=0; i<diffusantGasNames_.size(); i++) particleCloud_.dataExchangeM().destroy(diffusionCoefficients_[i],1);

    coeffs.clearStorage();
    molWeight.clearStorage();
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //

void diffusionCoefficient::allocateMyArrays() const
{
    double initVal=0.0;
    if (particleCloud_.dataExchangeM().maxNumberOfParticles() > 0)
    {
        particleCloud_.dataExchangeM().allocateArray(partPressure_,initVal,1,"nparticles");
        for (int i=0; i<diffusantGasNames_.size(); i++)
        {
            particleCloud_.dataExchangeM().allocateArray(diffusionCoefficients_[i],initVal,1,"nparticles");
        }
    }
}

void diffusionCoefficient::reAllocMyArrays() const
{
    double initVal=0.0;
    particleCloud_.dataExchangeM().allocateArray(partPressure_,initVal,1,"nparticles");
    for (int i=0; i<diffusantGasNames_.size(); i++)
    {
        particleCloud_.dataExchangeM().allocateArray(diffusionCoefficients_[i],initVal,1);
    }
}

void diffusionCoefficient::init()
{
    for (int i=0; i<speciesNames_.size(); i++)
    {
        // Defining the Species volume scalar fields
        volScalarField& X = const_cast<volScalarField&>
                (mesh_.lookupObject<volScalarField>("X_"+speciesNames_[i]));
        X_.set(i, &X);
    }

    for (int j = 0; j < diffusantGasNames_.size(); j++)
    {
        volScalarField& Xdiffusant = const_cast<volScalarField&>
                (mesh_.lookupObject<volScalarField>("X_"+diffusantGasNames_[j]));
        Xdiffusant_.set(j, &Xdiffusant);

        if (verbose_)
        {
            Info << " Reading diffusing gas species list: " << diffusantGasNames_ << endl;
            Info << " Looking up diffusin gas species fields: " << "X_"+diffusantGasNames_[j] << endl;
            Info << " The molar fraction fields (Xdiffusant_i): " << Xdiffusant_[j].name() << nl << endl;
        }
    }

    initialized_ = true;
}

// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //

void diffusionCoefficient::execute()
{
    if (!initialized_)
    {
        init();
    }

    // realloc the arrays
    reAllocMyArrays();

    label  cellI=0;
    scalar Tfluid(0);
    scalar Pfluid(0);
    scalar Texp(0);
    scalar dBinary_(0);
    scalar Xnegative(0);

    List<scalar> TotalFraction_(diffusantGasNames_.size(),Zero);

    // defining interpolators for T and Pressure
    interpolationCellPoint <scalar> TInterpolator_(tempField_);
    interpolationCellPoint <scalar> PInterpolator_(P_);

    for (int index=0; index<particleCloud_.numberOfParticles(); ++index)
    {
        cellI = particleCloud_.cellIDs()[index][0];
        if (cellI >=0)
        {
            if(interpolation_)
            {
                vector position = particleCloud_.position(index);
                Tfluid          = TInterpolator_.interpolate(position,cellI);
                Pfluid          = PInterpolator_.interpolate(position,cellI);
            }
            else
            {
                Tfluid = tempField_[cellI];
                Pfluid = P_[cellI];

                for (int i = 0; i<speciesNames_.size(); i++)
                {
                    // total amount of negative molar fractions in the domain
                    // check it and then delete it
                    if (X_[i][cellI] < 0.0)
                    {
                        Xnegative += X_[i][cellI] * mesh_.time().deltaTValue();
                        Info << "total negative molar fractions = " << Xnegative << endl;
                    }
                }
            }

            partPressure_[index][0] = Pfluid;
            // change fluid pressure to 1 bar instead of Pa
            Pfluid = Pfluid / 100000.0;
            Texp   = Tfluid * sqrt(sqrt(Tfluid*Tfluid*Tfluid));

            for (int j=0; j<diffusantGasNames_.size(); j++)
            {
                TotalFraction_[j] = 0.0;
                dBinary_ = 0.0;

                for (int i=0; i < speciesNames_.size(); i++)
                {
                    // get molecular diffusion coefficients if diffusant gas and reactant gas are not equal
                    if (diffusantGasNames_[j] != speciesNames_[i])
                    {
                        if (verbose_)
                        {
                            Info << "molar weights diffuser gases: "  << molWeight(speciesNames_[i])      << nl << endl;
                            Info << "molar weights diffusant gases: " << molWeight(diffusantGasNames_[j]) << nl << endl;
                            Info << "Pressure: " << Pfluid << nl << endl;
                        }

                        if (coeffs.found(diffusantGasNames_[j]) && coeffs.found(speciesNames_[i]))
                        {
                            //  Fuller-Schettler-Giddings Equation
                            //  Unit of dBinary is [m^2/s]
                            //  INFO:: Normally unit of dBinary is cm^2/s, but the 1st term in RHS is 10^-3 instead
                            //  So here it is already converted
                            dBinary_ = 1e-7 * Texp * calcMolNum(j,i) / (Pfluid * calcDiffVol(j,i));

                            if (verbose_)
                            {
                                Info << "Xfluid for "     << speciesNames_[i]      << " : " << X_[i][cellI]          << nl << endl;
                                Info << "Xdiffusant for " << diffusantGasNames_[j] << " : " << Xdiffusant_[j][cellI] << nl << endl;
                            }

                            TotalFraction_[j] += X_[i][cellI] / dBinary_;

                            if (verbose_)
                                Info << "Total Fraction = " << TotalFraction_[j] << nl << endl;

                            // pass on dCoeff values to array
                            if (TotalFraction_[j] < VSMALL)
                                diffusionCoefficients_[j][index][0] = VSMALL;
                            else
                                diffusionCoefficients_[j][index][0] = (1.0 - Xdiffusant_[j][cellI]) / TotalFraction_[j];
                        }
                        else
                        {
                            FatalError
                                << "check tables for species diffusion volume"
                                << endl
                                << abort(FatalError);
                        }
                    }
                }

                if(verbose_)
                    Info << "diffusionCoefficient of species " << diffusantGasNames_[j] << " = " << diffusionCoefficients_[j][index][0] << endl;
            }
        }
    }

    particleCloud_.dataExchangeM().giveData(partPressureName_,"scalar-atom",partPressure_);

    for (int j=0; j<diffusantGasNames_.size(); j++)
    {
        word pushName = diffusantGasNames_[j] + "_diffCoeff";
        particleCloud_.dataExchangeM().giveData(pushName,"scalar-atom",diffusionCoefficients_[j]);
    }

    Info << "give data done" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void diffusionCoefficient::createCoeffs()
{
    // diffusion volume of species i^(0.33333)
    // unit = [cm^3/g-mole]
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

double diffusionCoefficient::calcDiffVol(int j, int i)
{
    double sqrtvolDiff = coeffs(diffusantGasNames_[j])+coeffs(speciesNames_[i]);
    return sqrtvolDiff * sqrtvolDiff;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void diffusionCoefficient::molWeightTable()
{
    // table for molecular weights
    // molecular weights here are in [g/mole]
    molWeight.insert("CO", 28.0101);
    molWeight.insert("CO2", 44.01);
    molWeight.insert("O2", 32.00);
    molWeight.insert("N2", 28.01);
    molWeight.insert("H2", 2.02);
    molWeight.insert("H2O", 18.01);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// either calculate Molecular Weight addition (eq. D_ij) or consturct hashtable with diffusant and diffuser species
double diffusionCoefficient::calcMolNum(int j, int i)
{
    double molNum_ = 0.0;

    double W1 = molWeight(diffusantGasNames_[j]);
    double W2 = molWeight(speciesNames_[i]);

    molNum_ = (W1 + W2) / (W1 * W2);
    molNum_ = sqrt(molNum_);

    return molNum_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
