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

#ifndef heatTransferGunnImplicit_H
#define heatTransferGunnImplicit_H

#include "fvCFD.H"
#include "cfdemCloudEnergy.H"
#include "heatTransferGunn.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class heatTransferGunnImplicit Declaration
\*---------------------------------------------------------------------------*/

class heatTransferGunnImplicit
:
    public heatTransferGunn
{
private:
    
    word QPartFluidCoeffName_;
        
    volScalarField QPartFluidCoeff_;

    mutable double **partHeatFluxCoeff_;              // Lagrangian array

    void allocateMyArrays() const;
    
    void giveData(int);
    
    void heatFlux(label, scalar, scalar, scalar);
    
    void heatFluxCoeff(label, scalar, scalar);
       
public:

    //- Runtime type information
    TypeName("heatTransferGunnImplicit");

    // Constructors

        //- Construct from components
        heatTransferGunnImplicit
        (
            const dictionary& dict,
            cfdemCloudEnergy& sm
        );


    // Destructor

        virtual ~heatTransferGunnImplicit();


    // Member Functions
	
	void addEnergyCoefficient(volScalarField&) const;
	
	void calcEnergyContribution();
	
	void postFlow();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

