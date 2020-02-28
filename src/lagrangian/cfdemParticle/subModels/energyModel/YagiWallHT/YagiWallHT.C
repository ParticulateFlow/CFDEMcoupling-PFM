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
#include "YagiWallHT.H"
#include "addToRunTimeSelectionTable.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

  // * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

  defineTypeNameAndDebug(YagiWallHT, 0);

  addToRunTimeSelectionTable(energyModel, YagiWallHT, dictionary);

  // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


  YagiWallHT::YagiWallHT
  (
    const dictionary& dict,
    cfdemCloudEnergy& sm
  )
  :
  energyModel(dict,sm),
  propsDict_(dict.subDict(typeName + "Props")),
  interpolation_(propsDict_.lookupOrDefault<bool>("interpolation",false)),
  verbose_(propsDict_.lookupOrDefault<bool>("verbose",false)),
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
  UsFieldName_(propsDict_.lookup("granVelFieldName")),
  Us_(sm.mesh().lookupObject<volVectorField> (UsFieldName_)),
  densityFieldName_(propsDict_.lookupOrDefault<word>("densityFieldName","rho")),
  rho_(sm.mesh().lookupObject<volScalarField> (densityFieldName_)),
  partRe_(NULL),
  multiphase_(propsDict_.lookupOrDefault<bool>("multiphase",false)),
  kfFieldName_(propsDict_.lookupOrDefault<word>("kfFieldName",voidfractionFieldName_)), // use voidfractionField as dummy to prevent lookup error when not using multiphase
  kfField_(sm.mesh().lookupObject<volScalarField> (kfFieldName_)),
  CpFieldName_(propsDict_.lookupOrDefault<word>("CpFieldName",voidfractionFieldName_)), // use voidfractionField as dummy to prevent lookup error when not using multiphase
  CpField_(sm.mesh().lookupObject<volScalarField> (CpFieldName_))
  {
    allocateMyArrays();

    if (propsDict_.found("maxSource"))
    {
        maxSource_=readScalar(propsDict_.lookup ("maxSource"));
        Info << "limiting wall source field to: " << maxSource_ << endl;
    }

    if (verbose_)
    {
      ReField_.writeOpt() = IOobject::AUTO_WRITE;
      PrField_.writeOpt() = IOobject::AUTO_WRITE;
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

  YagiWallHT::~YagiWallHT()
  {
    particleCloud_.dataExchangeM().destroy(partRe_,1);
  }

  // * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //
  void YagiWallHT::allocateMyArrays() const
  {
    // get memory for 2d arrays
    double initVal=0.0;

    if(verbose_)
    {
      particleCloud_.dataExchangeM().allocateArray(partRe_,initVal,1);
    }
  }

  // * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //

  void YagiWallHT::calcEnergyContribution()
  {
    allocateMyArrays();

    // reset Scalar field
    QWallFluid_.primitiveFieldRef() = 0.0;

    #ifdef compre
    const volScalarField mufField = particleCloud_.turbulence().mu();
    #else
    const volScalarField mufField = particleCloud_.turbulence().nu()*rho_;
    #endif

    interpolationCellPoint<scalar> voidfractionInterpolator_(voidfraction_);
    interpolationCellPoint<vector> UInterpolator_(U_);
    interpolationCellPoint<scalar> TInterpolator_(tempField_);

    // calculate Rep
    for(int index = 0;index < particleCloud_.numberOfParticles(); ++index)
    {
      label cellI = particleCloud_.cellIDs()[index][0];
      if(cellI >= 0)
      {
        scalar voidfraction;
        vector Ufluid;

        if(interpolation_)
        {
          vector position = particleCloud_.position(index);
          voidfraction = voidfractionInterpolator_.interpolate(position,cellI);
          Ufluid = UInterpolator_.interpolate(position,cellI);
        }
        else
        {
          voidfraction = voidfraction_[cellI];
          Ufluid = U_[cellI];
        }

        if (voidfraction < 0.01)
          voidfraction = 0.01;

        vector Us = particleCloud_.velocity(index);
        scalar magUr = mag(Ufluid - Us);
        scalar ds = 2.*particleCloud_.radius(index);
        scalar muf = mufField[cellI];

        scalar Rep = ds * magUr * voidfraction * rho_[cellI]/ muf;
        partRe_[index][0] = Rep;
      }
    }

    // calculate Rep field
    particleCloud_.averagingM().resetWeightFields();
    particleCloud_.averagingM().setScalarAverage
    (
      ReField_,
      partRe_,
      particleCloud_.particleWeights(),
      particleCloud_.averagingM().UsWeightField(),
      NULL
    );

    // calculate Pr field
    PrField_ = CpField_ * mufField / kfField_;

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

            // calculate Urel
            scalar magG = mag(U_[faceCelli]-Us_[faceCelli])*voidfraction_[faceCelli]*rho_[faceCelli];

            // calculate H
            scalar H;
            if (voidfraction_[faceCelli]<=voidfractionMax_)
              H = 0.2087 * (pow(ReField_[faceCelli]+SMALL,-0.20)) * CpField_[faceCelli] * magG / (pow(PrField_[faceCelli],2/3) + SMALL);
            else
              H = 0;

            // get delta T (wall-fluid)
            scalar Twall  = wallTemp_.boundaryField()[patchi][facei];
            scalar Tfluid = tempField_[faceCelli];
            scalar deltaT = Twall - Tfluid;

            // get area
            scalar area = curPatch.magSf()[facei];

            // calculate heat flux
            heatFlux(faceCelli, H, area, Twall, Tfluid);
            heatFluxCoeff(faceCelli, H, area);

            if(verbose_ && facei >=0 && facei <2)
            {
              Info << "####################" << endl;
              Info << "cellID: " << faceCelli << endl;
              Info << "G : " << magG << endl;
              Info << "Re: " << ReField_[faceCelli] << endl;
              Info << "Pr: " << PrField_[faceCelli] << endl;
              Info << "Cp: " << CpField_[faceCelli] << endl;
              Info << "kf: " << kfField_[faceCelli] << endl;
              Info << "H : " << H << endl;
              Info << "Twall: " << Twall << endl;
              Info << "Tfluid: " << Tfluid << endl;
              Info << "dT: " << deltaT << endl;
              Info << "q: " << H*deltaT << endl;
              Info << "area: " << area << endl;
              Info << "Q:" << H*deltaT*area << endl;
            }
          }		    
        }
        else
        {
            FatalError << "YagiWallHT requires zeroGradient BC for temperature field" << endl;
        }		    
      }	    
    }

    QWallFluid_.primitiveFieldRef() /= QWallFluid_.mesh().V();

    // limit source term
    forAll(QWallFluid_,cellI)
    {
        scalar EuFieldInCell = QWallFluid_[cellI];

        if(mag(EuFieldInCell) > maxSource_ )
        {
             Pout << "limiting source term"  << endl  ;
             QWallFluid_[cellI] = sign(EuFieldInCell) * maxSource_;
        }
    }

    QWallFluid_.correctBoundaryConditions();

  }

  void YagiWallHT::addEnergyContribution(volScalarField& Qsource) const
  {
    Qsource += QWallFluid_;
  }

  void YagiWallHT::heatFlux(label faceCelli, scalar H, scalar area, scalar Twall, scalar Tfluid)
  {
    QWallFluid_[faceCelli] += H * area * (Twall - Tfluid);
  }

  void YagiWallHT::heatFluxCoeff(label faceCelli, scalar H, scalar area)
  {
    //no heat transfer coefficient in explicit model
  }

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
