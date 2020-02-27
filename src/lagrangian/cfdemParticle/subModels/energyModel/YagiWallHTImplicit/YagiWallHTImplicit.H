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
    Correlation for wall-to-bed heat transfer according to
    Yagi, S. A.I.Ch.E. Journal 5.1 (1959)

\*---------------------------------------------------------------------------*/

#ifndef YagiWallHTImplicit_H
#define YagiWallHTImplicit_H

#include "fvCFD.H"
#include "cfdemCloudEnergy.H"
#include "YagiWallHT.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class YagiWallHTImplicit Declaration
\*---------------------------------------------------------------------------*/

class YagiWallHTImplicit
:
    public YagiWallHT
{
private:
    
    word QWallFluidCoeffName_;
        
    volScalarField QWallFluidCoeff_;

    mutable double **partHeatFluxCoeff_;              // Lagrangian array

    void allocateMyArrays() const;

    virtual void heatFlux(label, scalar, scalar, scalar, scalar);

    virtual void heatFluxCoeff(label, scalar, scalar);
       
public:

    //- Runtime type information
    TypeName("YagiWallHTImplicit");

    // Constructors

        //- Construct from components
        YagiWallHTImplicit
        (
            const dictionary& dict,
            cfdemCloudEnergy& sm
        );


    // Destructor

        virtual ~YagiWallHTImplicit();


    // Member Functions
	
	void addEnergyCoefficient(volScalarField&) const;
	
	void calcEnergyContribution();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
