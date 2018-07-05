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
                        M. Efe Kinaci, JKU Linz, Austria

\*---------------------------------------------------------------------------*/

#ifndef chemistryModel_H
#define chemistryModel_H

#include "fvCFD.H"
#include "cfdemCloudEnergy.H"
#include "probeModel.H"
#include "interpolationCellPoint.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class chemistryModel Declaration
\*---------------------------------------------------------------------------*/

class chemistryModel
{

protected:

    // Protected data
        const dictionary& dict_;

        cfdemCloudEnergy& particleCloud_;



public:

    //- Runtime type information
    TypeName("chemistryModel");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            chemistryModel,
            dictionary,
            (
                const dictionary& dict,
                cfdemCloudEnergy& sm
            ),
            (dict,sm)
        );


    // Constructors

        //- Construct from components
        chemistryModel
        (
            const dictionary& dict,
            cfdemCloudEnergy& sm
        );


    // Destructor

        virtual ~chemistryModel();


    // Selector

        static autoPtr<chemistryModel> New
        (
            const dictionary& dict,
            cfdemCloudEnergy& sm,
            word chemistryModelType
        );


    // Member Functions

        virtual void execute() = 0;

        virtual tmp<volScalarField> Smi(const label i) const;

        virtual tmp<volScalarField> Sm() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
