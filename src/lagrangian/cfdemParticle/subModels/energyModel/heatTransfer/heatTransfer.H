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

#ifndef heatTransfer_H
#define heatTransfer_H

#include "fvCFD.H"
#include "cfdemCloudEnergy.H"
#include "energyModel.H"
#include "LaEuScalarTemp.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class heatTransfer Declaration
\*---------------------------------------------------------------------------*/

class heatTransfer
:
    public energyModel, public LaEuScalarTemp
{
private:

    dictionary propsDict_;

    dictionary propsLaEuScalarTempDict_;

    scalar Cp_; 
       
public:

    //- Runtime type information
    TypeName("heatTransfer");

    // Constructors

        //- Construct from components
        heatTransfer
        (
            const dictionary& dict,
            cfdemCloudEnergy& sm
        );


    // Destructor

        virtual ~heatTransfer();


      // Member Functions


        void energyContribution(volScalarField&) const;
	
	void temperatureContribution(volScalarField&) const;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
