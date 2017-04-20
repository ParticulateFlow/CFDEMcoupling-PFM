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
#include "diffusionCoefficient.H"
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
    // define a file name in the coupling properties that contains the diffusionCoefficient
    specDict_
    (
        IFstream
        (
            fileName(propsDict_.lookup("ChemistryFile")).expand()
        )()
    ),
    // create a list from the Species table in the specified diffusionCoefficient dictionary
    speciesNames_(specDict_.lookup("species")),
    diffusionCoefficientNames_(propsDict_.lookup("diffusionCoefficientNames")),
    Y_(speciesNames_.size()),                           //volumeScalarFields created in the ts folders
    diffusionCoefficients_(diffusionCoefficientNames_.size(),NULL),         //the value of diffusionCoefficient concentration for every diffusionCoefficient 
    tempFieldName_(propsDict_.lookup("tempFieldName")),
    tempField_(sm.mesh().lookupObject<volScalarField> (tempFieldName_)),
    pressureFieldName_(propsDict_.lookup("pressureFieldName")),
    P_(sm.mesh().lookupObject<volScalarField>(pressureFieldName_))
{
    Info << " Read species list from: " << specDict_.name() << endl;
    Info << " Reading diffusionCoefficient list: " << diffusionCoefficientNames_ << endl;

    for (int i=0; i<speciesNames_.size(); i++)
    {
        // Defining the Species volume scalar fields
        Info << " Looking up diffusionCoefficient fields \n " << speciesNames_[i] << endl;
        volScalarField& Y = const_cast<volScalarField&>
                (sm.mesh().lookupObject<volScalarField>(speciesNames_[i]));
        Y_.set(i, &Y);

         Info << "The concentration fields (Y_i): \n" << Y_[i].name() << endl;
        // define the modified diffusionCoefficient names

        particleCloud_.checkCG(false);
    }

    allocateMyArrays();
    createCoeffs();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

diffusionCoefficient::~diffusionCoefficient()
{
    int nP_ = particleCloud_.numberOfParticles();
    for (int i=0; i<diffusionCoefficientNames_.size(); i++)
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
        for (int i=0; i<diffusionCoefficientNames_.size(); i++)
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

        for (int i=0; i<diffusionCoefficientNames_.size(); i++)
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
    
    scalar dCoeff(0.0);
    
    word speciesPair("none");


    // defining interpolators for T, rho, voidfraction, N
    interpolationCellPoint <scalar> TInterpolator_(tempField_);
    interpolationCellPoint <scalar> PInterpolator_(P_);


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
            }
            else
            {
                Tfluid          =   tempField_[cellI];
                Pfluid          =   P_[cellI];

                for (int i = 0; i<speciesNames_.size();i++)
                {
                    Yfluid_[i] = Y_[i][cellI];
                }
            }


            for (int i=0; i<diffusionCoefficientNames_.size();i++)
            {
	      // do the calculation
	        dCoeff = 0.0;
		for (int j=0; j < speciesNames_.size();j++)
		{
		    speciesPair = diffusionCoefficientNames_[i] + speciesNames_[j];
		    if(coeffs.find(speciesPair))
		    {
		        dCoeff += Y[j] * 
		    }
		}
                diffusionCoefficients_[i][index][0]=
            }
        }

        if(particleCloud_.verbose() && index >=0 && index < 2)
        {
            for(int i =0; i<diffusionCoefficientNames_.size();i++)
            {
                Info << "effective diffusionCoefficient of species " << diffusionCoefficientNames_[i] << " = " << diffusionCoefficients_[i][index][0] << endl;
            }
        }
    }

    for (int i=0; i<diffusionCoefficientNames_.size();i++)
    {
        word pushName = diffusionCoefficientNames_[i] + "_diffCoeff";
        particleCloud_.dataExchangeM().giveData(pushName,"scalar-atom",diffusionCoefficients_[i]);
    };

    Info << "give data done" << endl;

       
}

void diffusionCoefficient::createCoeffs()
{
  // add all relevant combinations
//   coeffs.insert("CO_CO2", );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
