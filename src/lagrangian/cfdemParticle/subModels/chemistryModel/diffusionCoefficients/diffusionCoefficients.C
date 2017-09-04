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
    tempFieldName_(propsDict_.lookupOrDefault<word>("tempFieldName","T")),
    tempField_(sm.mesh().lookupObject<volScalarField> (tempFieldName_)),
    pressureFieldName_(propsDict_.lookupOrDefault<word>("pressureFieldName","p")),
    P_(sm.mesh().lookupObject<volScalarField>(pressureFieldName_)),
    densityFieldName_(propsDict_.lookupOrDefault<word>("densityFieldName","rho")),
    rho_(sm.mesh().lookupObject<volScalarField> (densityFieldName_)),
    molarConcFieldName_(propsDict_.lookupOrDefault<word>("totalMoleFieldName","molarConc")),
    molarConc_(sm.mesh().lookupObject<volScalarField>(molarConcFieldName_)),
    X_(speciesNames_.size()),
    diffusantGasNames_(propsDict_.lookup("diffusantGasNames")),
    diffusionCoefficients_(diffusantGasNames_.size(),NULL),
    initialized_(false)
{
    Info << " Reading diffusionCoefficient list: " << diffusantGasNames_ << endl;
    for (int i = 0; i < diffusantGasNames_.size(); i++)
    {
        Info << "Diffusant names: " << diffusantGasNames_[i] << endl;
    }

    particleCloud_.checkCG(false);
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

void diffusionCoefficient::init()
{
    // look-up of molar fraction fields can't happen in constructor because functionObject
    // has not been created at that time

    Info << " Read species list from: " << specDict_.name() << endl;
    Info << " Reading species list: " << speciesNames_ << endl;

    for (int i=0; i<speciesNames_.size(); i++)
    {
        // Defining the Species volume scalar fields
        Info << " Looking up species fields \n " << "X_"+speciesNames_[i] << endl;
        volScalarField& X = const_cast<volScalarField&>
                (mesh_.lookupObject<volScalarField>("X_"+speciesNames_[i]));
        X_.set(i, &X);

         Info << "The molar fraction fields (X_i): \n" << X_[i].name() << endl;
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
    scalar rhofluid(0);
    List<scalar> Xfluid_;
    Xfluid_.setSize(speciesNames_.size());
    scalar Pfluid(0);
    scalar molarConcfluid(0);
    scalar Texp(0);
    List<scalar> dBinarytot_;
    dBinarytot_.setSize(diffusantGasNames_.size());
    List<scalar> molecularDiffusion_;
    molecularDiffusion_.setSize(diffusantGasNames_.size());
    List<scalar> Xstag_tot_;
    Xstag_tot_.setSize(diffusantGasNames_.size());

    double  **dBinary_   = new double*[diffusantGasNames_.size()];
    double  **molNum_    = new double*[diffusantGasNames_.size()];
    double  **volDiff_   = new double*[diffusantGasNames_.size()];

    // defining interpolators for T, rho, voidfraction, N
    interpolationCellPoint <scalar> TInterpolator_(tempField_);
    interpolationCellPoint <scalar> rhoInterpolator_(rho_);
    interpolationCellPoint <scalar> PInterpolator_(P_);
    interpolationCellPoint <scalar> molarConcInterpolator_(molarConc_);

    for (int index=0; index<particleCloud_.numberOfParticles(); ++index)
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
                molarConcfluid      =   molarConcInterpolator_.interpolate(position,cellI);
            }
            else
            {
                Tfluid          =   tempField_[cellI];


                rhofluid        =   rho_[cellI];
                Pfluid          =   P_[cellI];
                molarConcfluid  =   molarConc_[cellI];

                for (int i = 0; i<speciesNames_.size();i++)
                {
                    Xfluid_[i] = X_[i][cellI];
                    if (Xfluid_[i] <= 0.) Xfluid_[i] = 0;
                    Info << "X fluid for species " << speciesNames_[i] << " : " << Xfluid_[i] << nl << endl;
                }
            }

            Texp            =   pow(Tfluid,1.75);

            for (int i=0; i<diffusantGasNames_.size();i++)
            {
                dBinary_[i]     =   new double [speciesNames_.size()];
                molNum_[i]      =   new double [speciesNames_.size()];
                volDiff_[i]     =   new double [speciesNames_.size()];

                dBinarytot_[i]  =   0.0;
                molecularDiffusion_[i]      =   0.0;
                Xstag_tot_[i]   =   0.0;
                for (int j=0; j < speciesNames_.size();j++)
                {
                    // get molecular diffusion coefficients if diffusant gas and reactant gas are not equal
                    if (diffusantGasNames_[i] != speciesNames_[j])
                    {
                        Info << "molar weights diffuser gases: " << molWeight(speciesNames_[j]) << nl << endl;
                        Info << "molarConc fluid: " << molarConcfluid << nl << endl;
                        Info << "rho fluid:  " << rhofluid << nl << endl;
                        Info << "molar weights diffusant gases: " << molWeight(diffusantGasNames_[i]) << nl << endl;
                        Info << "Pressure: " << Pfluid << nl << endl;
                        Info << "Temperature: " << Tfluid << nl << endl;

                        calcMolNum(i,j,molNum_);
                        calcDiffVol(i,j,volDiff_);

                        if(coeffs.found(diffusantGasNames_[i]) && coeffs.found(speciesNames_[j]))
                        {
                            //  Fuller-Schettler-Giddings Equation
                            //  Unit of dBinary is [m^2/s]
                            //  INFO:: Normally unit of dBinary is cm^2/s, but the 1st term in RHS is 10^-3 instead
                            //  So here it is already converted
                            dBinary_[i][j] = pow(10,-7)*Texp*molNum_[i][j]/(Pfluid*volDiff_[i][j]);
                            Info << "Molecular diffusion for species " << diffusantGasNames_[i] << " in "
                                 << speciesNames_[j] << " is : " << dBinary_[i][j] << nl << endl;

                            // sum of all binary diffusion coefficients except i = j
                            dBinarytot_[i]  +=  dBinary_[i][j];
                            Info << "Sum of binary mol. diffusion coefficients for diffusing gas " << diffusantGasNames_[i]
                                 << " : " << dBinarytot_[i] << nl << endl;

                            Info << "Molar fraction of species " << speciesNames_[j] << " : " << Xfluid_[j] << nl << endl;
                            Info << "What is value of xstag? " << Xstag_tot_[i] << nl << endl;

                            // sum of all molar fractions except the diffusant gas
                            Xstag_tot_[i]   +=  Xfluid_[j];
                            Info << "Sum of stagnant gases molar fractions according to diffusing gas " << diffusantGasNames_[i]
                                 << " : " << Xstag_tot_[i] << nl << endl;

                            // dCoeff -- diffusion component of diffusant gas
                            if (Xstag_tot_[i] == 0.0) //(Xstag_tot_[i] == 0. || (1-Xfluid_[i]) == 0.)
                            {
                                molecularDiffusion_[i] = 0.0;
                            } else
                            {
                                molecularDiffusion_[i]  =   (1-Xfluid_[i])*dBinarytot_[i]/Xstag_tot_[i];
                            }

                            Info << "Multicomp. mix diffusion for species " << diffusantGasNames_[i]
                                 << " is: " << molecularDiffusion_[i] << nl << endl;

                            // pass on dCoeff values to array
                            diffusionCoefficients_[i][index][0]= molecularDiffusion_[i];

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

void diffusionCoefficient::calcDiffVol(int i, int j, double **volDiff_)
{
    volDiff_[i][j]  =   coeffs(diffusantGasNames_[i])+coeffs(speciesNames_[j]);
    volDiff_[i][j]  *=   volDiff_[i][j];
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
    molWeight.insert("H2O", 2.02);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// either calculate Molecular Weight addition (eq. D_ij) or consturct hashtable with diffusant and diffuser species
void diffusionCoefficient::calcMolNum(int i, int j, double **molNum_)
{
    double& W1 = molWeight(diffusantGasNames_[i]);
    double& W2 = molWeight(speciesNames_[j]);

    molNum_[i][j]   =   (W1 + W2) / (W1 * W2);
    molNum_[i][j]   =   sqrt(molNum_[i][j]);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
