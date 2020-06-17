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

#ifndef energyModel_H
#define energyModel_H

#include "fvCFD.H"
#include "cfdemCloudEnergy.H"
#include "probeModel.H"
#include "interpolationCellPoint.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class energyModel Declaration
\*---------------------------------------------------------------------------*/

class energyModel
{

protected:

    // Protected data
        const dictionary& dict_;

        cfdemCloudEnergy& particleCloud_;

        IOdictionary transportProperties_;

        scalar kf0_;                                // fluid thermal conductivity [W/(m*K)]

        scalar Cp_;                                 // specific heat capacity [W*s/(kg*K)]

public:

    //- Runtime type information
    TypeName("energyModel");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            energyModel,
            dictionary,
            (
                const dictionary& dict,
                cfdemCloudEnergy& sm
            ),
            (dict,sm)
        );


    // Constructors

        //- Construct from components
        energyModel
        (
            const dictionary& dict,
            cfdemCloudEnergy& sm
        );


    // Destructor

        virtual ~energyModel();


    // Selector

        static autoPtr<energyModel> New
        (
            const dictionary& dict,
            cfdemCloudEnergy& sm,
            word energyType
        );


    // Member Functions

        virtual void addEnergyContribution(volScalarField&) const = 0;

        virtual void addEnergyCoefficient(volScalarField&) const = 0;

        virtual void calcEnergyContribution() = 0;

        virtual void postFlow() {}

        virtual void solve() {}

        scalar Cp() const;

        virtual scalar aveTpart() const
        {
            return 0.0;
        }


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
