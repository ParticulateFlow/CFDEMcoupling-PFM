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
  WallTempField_
  (   IOobject
    (
      "WallTemp",
      sm.mesh().time().timeName(),
      sm.mesh(),
      IOobject::READ_IF_PRESENT,
      IOobject::NO_WRITE
    ),
    sm.mesh(),
    dimensionedScalar("zero", dimensionSet(0,0,0,1,0,0,0), 0.0)      // should i fix the wall temp here if its constant or fix it in BC?
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
  maxSource_(1e30),
  velFieldName_(propsDict_.lookupOrDefault<word>("velFieldName","U")),
  U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
  UsFieldName_(propsDict_.lookup("granVelFieldName")),
	Us_(sm.mesh().lookupObject<volVectorField> (UsFieldName_)),
  densityFieldName_(propsDict_.lookupOrDefault<word>("densityFieldName","rho")),
  rho_(sm.mesh().lookupObject<volScalarField> (densityFieldName_)),
  WallTempName_(propsDict_.lookup("WallTempName")),
  WallTemp_(sm.mesh().lookupObject<volScalarField> (WallTempName_)),
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
      Info << "limiting eulerian source field to: " << maxSource_ << endl;
    }

    if (verbose_)
    {
      ReField_.writeOpt() = IOobject::AUTO_WRITE;
      PrField_.writeOpt() = IOobject::AUTO_WRITE;
      ReField_.write();
      PrField_.write();
    }
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

    QWallFluid_.primitiveFieldRef() = 0.0;



    #ifdef compre
    const volScalarField mufField = particleCloud_.turbulence().mu();
    #else
    const volScalarField mufField = particleCloud_.turbulence().nu()*rho_;
    #endif



    // calc La based heat flux
    scalar voidfraction(1);
    vector Ufluid(0,0,0);
    label cellI=0;
    vector Us(0,0,0);
    scalar ds(0);
    scalar muf(0);
    scalar magUr(0);
    scalar Rep(0);
    scalar H(0);

    interpolationCellPoint<scalar> voidfractionInterpolator_(voidfraction_);
    interpolationCellPoint<vector> UInterpolator_(U_);
    interpolationCellPoint<scalar> TInterpolator_(tempField_);

    for(int index = 0;index < particleCloud_.numberOfParticles(); ++index)
    {
      cellI = particleCloud_.cellIDs()[index][0];
      if(cellI >= 0)
      {
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

        Us = particleCloud_.velocity(index);
        magUr = mag(Ufluid - Us);
        ds = 2.*particleCloud_.radius(index);
        muf = mufField[cellI];

        Rep = ds * magUr * voidfraction * rho_[cellI]/ muf;
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
    //    mass_flow_rate = magUr * rho_[cellI];

    const fvPatchList& patches = U_.mesh().boundary();

        forAll(patches, patchi)
        {
            const fvPatch& curPatch = patches[patchi];

            if (WallTemp_.boundaryField().types()[patchi] == "fixedValue")
            {
                forAll(curPatch, facei)
                {
                    label faceCelli = curPatch.faceCells()[facei];
                    // calculate Urel
                    scalar magG = mag(U_[faceCelli]-Us_[faceCelli])*voidfraction_[faceCelli]*rho_[faceCelli];
                    // calculate H
                    //H = 0.2087 * (pow(ReField_[faceCelli] , -0.20)) * CpField_[faceCelli] * magG / (pow(PrField_[faceCelli] , (2/3)));
                    H = 1;
                    // get delta T
                    //scalar deltaT = WallTemp_.boundaryField()[patchi][facei] - tempField_.boundaryField()[patchi][facei];
                    scalar deltaT = 100;
                    // get facei.area
                    scalar area = curPatch.magSf()[facei];
                    scalar volume = U_.mesh().V()[faceCelli];

                    //Info << "faceID: " << facei << " area: " << area << " volume: " << volume << endl;

                    QWallFluid_[faceCelli] = H*deltaT*area/volume;
                }
            }
        }

//    h = 0.2087 * (pow(partRe_ , -0.20)) * CpField_ * magUr * rho_ / (pow(PrField_ , (2/3)));

    // limit source term
    forAll(QWallFluid_,cellI)
    {
      scalar EuFieldInCell = QWallFluid_[cellI];

      if(mag(EuFieldInCell) > maxSource_ )
      {
        Pout << "limiting source term\n"  << endl  ;
        QWallFluid_[cellI] = sign(EuFieldInCell) * maxSource_;
      }
    }

    QWallFluid_.correctBoundaryConditions();


  }

  void YagiWallHT::addEnergyContribution(volScalarField& Qsource) const
  {
    Qsource += QWallFluid_;
  }


  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
