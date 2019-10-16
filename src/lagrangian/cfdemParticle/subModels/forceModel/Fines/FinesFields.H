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

Description
    class managing fines fields

SourceFiles
    finesFields.C
\*---------------------------------------------------------------------------*/

#ifndef FinesFields_H
#define FinesFields_H

#include "fvCFD.H"
#include "cfdemCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class FinesFields Declaration
\*---------------------------------------------------------------------------*/

class FinesFields
{
private:
    cfdemCloud& particleCloud_;

    dictionary propsDict_;

    bool smoothing_;

    bool verbose_;

    bool clogKin_;

    bool clogStick_;

    bool movingBed_;

    word fluxFieldName_;

    const surfaceScalarField& phi_;

    word velFieldName_;

    const volVectorField& U_;

    word voidfractionFieldName_;

    volScalarField& voidfraction_;

    word UsFieldName_;

    const volVectorField& UsField_;

    word pFieldName_;

    const volScalarField& p_;

    word rhoGFieldName_;

    const volScalarField& rhoG_;

    const volScalarField& dSauter_;

    volScalarField alphaG_;

    volScalarField alphaDyn_;

    volScalarField alphaP_;

    volScalarField alphaSt_;

    volScalarField deltaAlpha_;

    volScalarField dHydMix_;

    volScalarField DragCoeff_;

    volScalarField dSauterMix_;

    volScalarField FanningCoeff_;

    volScalarField Froude_;

    volScalarField Sds_;

    surfaceScalarField massFluxDyn_;

    volVectorField uDyn_;

    autoPtr<volVectorField> uReconstructed_; // velocity field can show oscillations at porosity steps

    dimensionedScalar dFine_;

    dimensionedScalar diffCoeff_;

    dimensionedScalar nuAve_;

    dimensionedScalar rhoFine_;

    const dimensionedVector g_;

    scalar alphaDynMax_;

    scalar alphaMax_;

    scalar alphaMinClog_;

    scalar critVoidfraction_;

    scalar deltaT_;

    scalar depositionLength_;

    scalar exponent_;

    scalar nCrit_;

    scalar poresizeWidth_;

    scalar prefactor_;

    scalar ratioHydraulicPore_;

    scalar tauRelease_;

    scalar uBindHigh_;

    scalar uBindLow_;

    scalar uMin_;

    void calcSource();

    void integrateFields();

    void updateAlphaG();

    void updateAlphaP();

    void updateDHydMix();

    void updateDragCoeff();

    void updateDSauter();

    void updateFanningCoeff();

    void updateFroude();

    void updateUDyn();

    void updateUReconstructed();

public:

    //- Runtime type information
//    TypeName("FinesFields");


    // Constructors

        //- Construct from components
        FinesFields
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        virtual ~FinesFields();


    // Member Functions


    void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

