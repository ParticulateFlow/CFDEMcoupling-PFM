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
    Correlation for Nusselt number according to
    Gunn, D. J. International Journal of Heat and Mass Transfer 21.4 (1978)

\*---------------------------------------------------------------------------*/

#ifndef heatTransferGunn_H
#define heatTransferGunn_H

#include "fvCFD.H"
#include "cfdemCloudEnergy.H"
#include "energyModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class heatTransferGunn Declaration
\*---------------------------------------------------------------------------*/

class heatTransferGunn
:
    public energyModel
{
protected:

    dictionary propsDict_;
    
    bool interpolation_;
    
    bool verbose_;
    
    word QPartFluidName_;
        
    volScalarField QPartFluid_;
    
    volScalarField partTempField_;
    
    volScalarField partRelTempField_;
    
    volScalarField ReField_;
    
    volScalarField NuField_;
    
    dimensionedScalar partRefTemp_;
    
    bool calcPartTempField_;

    word tempFieldName_;

    const volScalarField& tempField_;            // ref to temperature field

    word voidfractionFieldName_;

    const volScalarField& voidfraction_;    // ref to voidfraction field

    scalar maxSource_;                           // max (limited) value of src field

    word velFieldName_;

    const volVectorField& U_;
    
    word densityFieldName_;
    
    const volScalarField& rho_;

    word partTempName_;

    mutable double **partTemp_;                          // Lagrangian array

    word partHeatFluxName_;

    mutable double **partHeatFlux_;              // Lagrangian array
    
    mutable double **partRe_;
    
    mutable double **partNu_;

    bool multiphase_;

    word kfFieldName_;

    const volScalarField& kfField_;

    word CpFieldName_;

    const volScalarField& CpField_;

    void allocateMyArrays() const;
    
    scalar Nusselt(scalar, scalar, scalar) const;
    
    virtual void giveData(int);
    
    virtual void heatFlux(label, scalar, scalar, scalar);
    
    virtual void heatFluxCoeff(label, scalar, scalar);
       
public:

    //- Runtime type information
    TypeName("heatTransferGunn");

    // Constructors

        //- Construct from components
        heatTransferGunn
        (
            const dictionary& dict,
            cfdemCloudEnergy& sm
        );


    // Destructor

        virtual ~heatTransferGunn();


    // Member Functions

        void addEnergyContribution(volScalarField&) const;
	
	void addEnergyCoefficient(volScalarField&) const {}
	
	void calcEnergyContribution();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

