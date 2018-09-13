/*---------------------------------------------------------------------------*\
    CFDEMcoupling academic - Open Source CFD-DEM coupling
    
    Contributing authors:
    Thomas Lichtenegger, Gerhard Holzinger
    Copyright (C) 2015- Johannes Kepler University, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling academic.

    CFDEMcoupling academic is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CFDEMcoupling academic is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling academic.  If not, see <http://www.gnu.org/licenses/>.

Description
    Standard Recurrence Model
    
    This class sets limits for the length of a recurrence-sequence. As of now, 
    such a sequence is between one 20-th and 5-th of the total number of data 
    base entries.

\*---------------------------------------------------------------------------*/

#ifndef standardRecModel_H
#define standardRecModel_H

#include "recModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class standardRecModel Declaration
\*---------------------------------------------------------------------------*/

class standardRecModel
:
    public recModel
{
protected:
  
    dictionary propsDict_;
    wordList dataBaseNames_;
    label numDataBases_;
    List<Foam::Time> recTimes_;
    List<instantList> timeDirs_;
    Switch skipZero_;
    labelList numRecFields_;
    labelList cumulativeNumRecFields_;
    label totNumRecFields_;
    
    bool storeAveragedFields_;

    // matrix that contains the recurrence ERROR
    SymmetricSquareMatrix<scalar> recurrenceMatrix_;
    
    // create a data structure for the time indices
    //	constant will not be contained
    //	runTimeIndex -> continuousIndex
    HashTable<label,word> timeIndexList_;
    
    // create a data structure for the time values
    //	constant will not be contained
    //	continuousIndex -> time.value()
    HashTable<label,scalar> timeValueList_;
    
    label contTimeIndex;
    
    label lowerSeqLim_;
    label upperSeqLim_;

    scalar recTimeDilationFactor_;
    
    scalar checkTimeStep();
    
    inline label getVolScalarFieldIndex(word, label) const;
    inline label getVolVectorFieldIndex(word, label) const;
    inline label getSurfaceScalarFieldIndex(word, label) const;

    void init();

    void averageFieldSeries();

    void readFieldSeries();

    void readTimeSeries();

public:

    //- Runtime type information
    TypeName("standardRecModel");


    // Constructors

        //- Construct from components
        standardRecModel
        (
            const dictionary& dict,
            recBase& base
        );

    // Destructor

        ~standardRecModel();

    void exportVolScalarFieldAve(word, volScalarField&, label);
    void exportVolVectorFieldAve(word, volVectorField&, label);
    void exportSurfaceScalarFieldAve(word, surfaceScalarField&, label);
    
    tmp<surfaceScalarField> exportSurfaceScalarFieldAve(word, label = 0);

    void exportVolScalarField(word, volScalarField&);
    void exportVolVectorField(word, volVectorField&);
    void exportSurfaceScalarField(word, surfaceScalarField&);
    
    const volScalarField& exportVolScalarField(word, label);
    const volVectorField& exportVolVectorField(word, label);
    const surfaceScalarField& exportSurfaceScalarField(word, label);
    
  //  tmp<surfaceScalarField> exportAveragedSurfaceScalarField(word, scalar, label index = -1);
    void exportAveragedVolVectorField(volVectorField&, word, scalar, label index = -1) const;
   
    SymmetricSquareMatrix<scalar>& recurrenceMatrix();
    
    const HashTable<label,word>& timeIndexList() const;
    
    label lowerSeqLim() const;
    label upperSeqLim() const;
    
    label numIntervals() const;
    label numRecFields() const;
    label numRecFields(label) const;
    label numDataBaseFields() const;
    
    void updateRecFields();
    
    void writeRecMatrix() const;
	
private:
    
    List<PtrList<volScalarField> > volScalarFieldList_;
    List<PtrList<volVectorField> > volVectorFieldList_;
    List<PtrList<surfaceScalarField> > surfaceScalarFieldList_;
    
    List<PtrList<volScalarField> > aveVolScalarFieldList_;
    List<PtrList<volVectorField> > aveVolVectorFieldList_;
    List<PtrList<surfaceScalarField> > aveSurfaceScalarFieldList_;
    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "standardRecModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
