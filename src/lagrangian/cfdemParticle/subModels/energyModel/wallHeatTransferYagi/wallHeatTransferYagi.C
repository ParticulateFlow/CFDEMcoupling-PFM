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
#include "wallHeatTransferYagi.H"
#include "addToRunTimeSelectionTable.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

  // * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(wallHeatTransferYagi, 0);

addToRunTimeSelectionTable(energyModel, wallHeatTransferYagi, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


wallHeatTransferYagi::wallHeatTransferYagi
(
    const dictionary& dict,
    cfdemCloudEnergy& sm
)
:
    energyModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    verbose_(propsDict_.lookupOrDefault<bool>("verbose",false)),
    implicit_(propsDict_.lookupOrDefault<bool>("implicit",true)),
    QWallFluidName_(propsDict_.lookupOrDefault<word>("QWallFluidName","QWallFluid")),
    QWallFluid_
    (   IOobject
        (
            QWallFluidName_,
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(1,-1,-3,0,0,0,0), 0.0)
    ),
    QWallFluidCoeffName_(propsDict_.lookupOrDefault<word>("QWallFluidCoeffName","QWallFluidCoeff")),
    QWallFluidCoeff_
    (   IOobject
        (
            QWallFluidCoeffName_,
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(1,-1,-3,-1,0,0,0), 0.0)
    ),
    wallTempName_(propsDict_.lookup("wallTempName")),
    wallTemp_
    (   IOobject
        (
        wallTempName_,
        sm.mesh().time().timeName(),
        sm.mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0,0,0,1,0,0,0), 0.0)
    ),
	dpField_
    (   IOobject
        (
        "dpField",
        sm.mesh().time().timeName(),
        sm.mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0,1,0,0,0,0,0), 0.0)
    ),
	GField_
    (   IOobject
        (
        "GField",
        sm.mesh().time().timeName(),
        sm.mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedVector("zero", dimensionSet(1,-2,-1,0,0,0,0), vector(0.0, 0.0, 0.0))
    ),
	magGField_
    (   IOobject
        (
        "magGField",
        sm.mesh().time().timeName(),
        sm.mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(1,-2,-1,0,0,0,0), 0.0)
    ),
    ReField_
    (   IOobject
        (
        "ReField",
        sm.mesh().time().timeName(),
        sm.mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0,0,0), 0.0)
    ),
    PrField_
    (   IOobject
        (
        "PrField",
        sm.mesh().time().timeName(),
        sm.mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0,0,0), 0.0)
    ),
    tempFieldName_(propsDict_.lookupOrDefault<word>("tempFieldName","T")),
    tempField_(sm.mesh().lookupObject<volScalarField> (tempFieldName_)),
    voidfractionFieldName_(propsDict_.lookupOrDefault<word>("voidfractionFieldName","voidfraction")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    voidfractionMax_(readScalar(propsDict_.lookup("voidfractionMax"))),
    maxSource_(1e30),
    velFieldName_(propsDict_.lookupOrDefault<word>("velFieldName","U")),
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
    densityFieldName_(propsDict_.lookupOrDefault<word>("densityFieldName","rho")),
    rho_(sm.mesh().lookupObject<volScalarField> (densityFieldName_)),
    dpArray_(NULL)
{
    allocateMyArrays();

    if (propsDict_.found("maxSource"))
    {
        maxSource_=readScalar(propsDict_.lookup ("maxSource"));
        Info << "limiting wall source field to: " << maxSource_ << endl;
    }

    if (!implicit_)
    {
        QWallFluidCoeff_.writeOpt() = IOobject::NO_WRITE;
    }

    if (verbose_)
    {
		dpField_.writeOpt() = IOobject::AUTO_WRITE;
        GField_.writeOpt() = IOobject::AUTO_WRITE;
		magGField_.writeOpt() = IOobject::AUTO_WRITE;
        ReField_.writeOpt() = IOobject::AUTO_WRITE;
        PrField_.writeOpt() = IOobject::AUTO_WRITE;

		dpField_.write();
        GField_.write();
		magGField_.write();
        ReField_.write();
        PrField_.write();
    }
    
    // currently it is detected if field was auto generated or defined
    // improvement would be changing the type here automatically
    forAll(wallTemp_.boundaryField(),patchI)
        if(wallTemp_.boundaryField()[patchI].type()=="calculated")
            FatalError <<"Scalar field:"<< wallTemp_.name() << " must be defined.\n" << abort(FatalError);

    wallTemp_.writeOpt() = IOobject::AUTO_WRITE;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

wallHeatTransferYagi::~wallHeatTransferYagi()
{
    particleCloud_.dataExchangeM().destroy(dpArray_,1);
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //
void wallHeatTransferYagi::allocateMyArrays() const
{
    // get memory for 2d arrays
    double initVal=0.0;

    if(verbose_)
    {
        particleCloud_.dataExchangeM().allocateArray(dpArray_,initVal,1);
    }
}

  // * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //

void wallHeatTransferYagi::calcEnergyContribution()
{
    allocateMyArrays();

    // reset Scalar field
    QWallFluid_.primitiveFieldRef() = 0.0;
    if (implicit_)
        QWallFluidCoeff_.primitiveFieldRef() = 0.0;

    #ifdef compre
    const volScalarField mufField = particleCloud_.turbulence().mu();
    #else
    const volScalarField mufField = particleCloud_.turbulence().nu()*rho_;
    #endif

	const volScalarField& CpField_  = CpField();
	const volScalarField& kf0Field_ = kf0Field();

    // calculate mean dp
    for(int index = 0;index < particleCloud_.numberOfParticles(); ++index)
    {
        label cellI = particleCloud_.cellIDs()[index][0];
        if(cellI >= 0)
        {
            scalar ds = 2.*particleCloud_.radius(index);
            dpArray_[index][0] = ds;
        }
    }

    // calculate mean dp field
    particleCloud_.averagingM().resetWeightFields();
    particleCloud_.averagingM().setScalarAverage
    (
        dpField_,
        dpArray_,
        particleCloud_.particleWeights(),
        particleCloud_.averagingM().UsWeightField(),
        NULL
    );

    // calculate G field (superficial mass velocity)
	GField_ = U_*voidfraction_*rho_;
	magGField_ = mag(GField_);

	// calculate Re field
	ReField_ = dpField_ * magGField_ / mufField;

	// calculate Pr field
	PrField_ = CpField_ * mufField / kf0Field_;         

    const fvPatchList& patches = U_.mesh().boundary();

    // calculate flux
    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (wallTemp_.boundaryField().types()[patchi] == "fixedValue")
        {
            if(tempField_.boundaryField().types()[patchi] == "zeroGradient")
            {
                forAll(curPatch, facei)
                {							
                    label faceCelli = curPatch.faceCells()[facei];

                    // calculate H
                    scalar JH;
                    if (voidfraction_[faceCelli]<=voidfractionMax_)
                        JH = 0.20 * pow(ReField_[faceCelli]+SMALL,-0.20); // Yagi eq. 6a
                    else
                        JH = 0;

					scalar h = JH * CpField_[faceCelli] * magGField_[faceCelli] / (pow(PrField_[faceCelli],0.666) + SMALL);

                    // get delta T (wall-fluid)
                    scalar Twall  = wallTemp_.boundaryField()[patchi][facei];
                    scalar Tfluid = tempField_[faceCelli];

                    // get area
                    scalar area = curPatch.magSf()[facei];

                    // calculate heat flux
                    heatFlux(faceCelli, h, area, Twall, Tfluid);

                    if(verbose_ && facei >=0 && facei <2)
                    {
						scalar deltaT = Twall - Tfluid;
                        Info << "####################" << endl;
                        Info << "cellID: " << faceCelli << endl;
						Info << "kf: " << kf0Field_[faceCelli] << " J/msK" << endl;
						Info << "Cp: " << CpField_[faceCelli] << " J/kgK" << endl;
						Info << "ro: " << rho_[faceCelli] << " kg/m3" << endl;
						Info << "mu: " << mufField[faceCelli] << " Pa s" << endl;
						Info << "dp: " << dpField_[faceCelli] << " m" << endl;
						Info << "ep: " << voidfraction_[faceCelli] << endl;
						Info << "U : " << U_[faceCelli] << " m/s" << endl;
						Info << "G : " << GField_[faceCelli] << " kg/m2s" << endl;
                        Info << "mG: " << magGField_[faceCelli] << " kg/m2s" << endl;
                        Info << "Re: " << ReField_[faceCelli] << endl;
						Info << "Pr: " << PrField_[faceCelli] << endl;
						Info << "JH: " << JH << endl;
                        Info << "h : " << h << " J/m2sK" << endl;
                        Info << "Tw: " << Twall << " K" << endl;
                        Info << "Tf: " << Tfluid << " K" << endl;
                        Info << "dT: " << deltaT << " K" << endl;
                        Info << "q : " << h*deltaT << " J/m2s" << endl;
                        Info << "A : " << area << " m2" << endl;
                        Info << "Q : " << h*deltaT*area << " J/s" << endl;
                    }
                }		    
            }
            else
            {
                FatalError << "wallHeatTransferYagi requires zeroGradient BC for temperature field" << endl;
            }		    
        }	    
    }

    QWallFluid_.primitiveFieldRef() /= QWallFluid_.mesh().V();
    if(implicit_)
        QWallFluidCoeff_.primitiveFieldRef() /= QWallFluidCoeff_.mesh().V();

    // limit source term in explicit treatment
    if(!implicit_)
    {
        forAll(QWallFluid_,cellI)
        {
            scalar EuFieldInCell = QWallFluid_[cellI];

            if(mag(EuFieldInCell) > maxSource_ )
            {
                Pout << "limiting source term"  << endl  ;
                QWallFluid_[cellI] = sign(EuFieldInCell) * maxSource_;
            }
        }
    }

    QWallFluid_.correctBoundaryConditions();
  }

void wallHeatTransferYagi::addEnergyContribution(volScalarField& Qsource) const
{
    Qsource += QWallFluid_;
}

void wallHeatTransferYagi::addEnergyCoefficient(volScalarField& Qcoeff) const
{
    if(implicit_)
    {
        Qcoeff += QWallFluidCoeff_;
    }
}

void wallHeatTransferYagi::heatFlux(label faceCelli, scalar h, scalar area, scalar Twall, scalar Tfluid)
{
    if(!implicit_)
    {
        QWallFluid_[faceCelli] += h * area * (Twall - Tfluid);
    }
    else
    {
        QWallFluid_[faceCelli]      += h * area * Twall;
        QWallFluidCoeff_[faceCelli] -= h * area;
    }
}

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
