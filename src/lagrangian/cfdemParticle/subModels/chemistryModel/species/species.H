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
    Model to communicate species concentration and changes due to reactions

\*---------------------------------------------------------------------------*/

#ifndef species_H
#define species_H

#include "fvCFD.H"
#include "cfdemCloudEnergy.H"
#include "chemistryModel.H"

#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class species Declaration
\*---------------------------------------------------------------------------*/

// for future use:
// + communicate every N steps

class species
:
    public chemistryModel
{
private:

    dictionary propsDict_;

    bool interpolation_;

    // Species dictionary - To be read where the species are from
    dictionary specDict_;

    const wordList speciesNames_;
    
    PtrList<volScalarField> Y_;

    mutable List<double**> concentrations_;
    
    mutable List<double**> changeOfSpeciesMass_;

    //List<volScalarField> changeOfSpeciesMassFields_;
    
    //volScalarField changeOfGasMassField_;
    
    word tempFieldName_;

    const volScalarField& tempField_;            // ref to gas temperature field

    word partTempName_;

    mutable double **partTemp_;                  // gas temperature at particle positions

    word densityFieldName_;

    const volScalarField& rho_;

    word partRhoName_;
    
    mutable double **partRho_;                   // gas density at particle positions
    
    //word velFieldName_;                          // Can be uncommented when velocity is to be calculated in the species source code

    //const volVectorField& U_;

    //word voidfractionFieldName_;                 // Can be uncommented when voidfraction is to be calculated in the species source code

    //const volScalarField& voidfraction_;


    void allocateMyArrays() const;
       
public:

    //- Runtime type information
    TypeName("species");

    // Constructors

        //- Construct from components
        species
        (
            const dictionary& dict,
            cfdemCloudEnergy& sm
        );


    // Destructor

        virtual ~species();


    // Member Functions

    void execute();
    
//     tmp<fvScalarMatrix> Smi(const label i) const;
//     
//     tmp<fvScalarMatrix> Sm() const;
    
    


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
