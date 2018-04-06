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

Class
    Foam::gerhardsRecModel

Description
    A recurrence model for a dataBase with a smaller, or equal number of slots
    than there are snapshots on disk.

    In the case of less slots than snapshots, we follow the Least Recently Used
    page replacement algorithm described in computer science literature.
    A label list is used to keep track of which slot was most recently used,
    the least recently used slot is then vacated if a new snapshot needs to be
    read from disk.

    Reference:
    \verbatim
        "Modern Operating Systems"
        Andrew S. Tannenbaum,
        Prentice Hall, 1992
    \endverbatim

SourceFiles
    gerhardsRecModelI.H
    gerhardsRecModel.C

\*---------------------------------------------------------------------------*/

#ifndef gerhardsRecModel_H
#define gerhardsRecModel_H

#include "recModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class gerhardsRecModel Declaration
\*---------------------------------------------------------------------------*/

class gerhardsRecModel
:
    public recModel
{
protected:

    dictionary propsDict_;
    word dataBaseName_;
    Foam::Time recTime;
    instantList timeDirs;
    Switch skipZero_;
    Switch checkSkipZero_;
    label numRecFields_;

    const Switch verboseVerbose_;
    const Switch verboseVacate_;

    // matrix that contains the recurrence ERROR
    SymmetricSquareMatrix<scalar> recurrenceMatrix_;

    // create a data structure for the time indices
    //  constant will not be contained
    //  runTimeIndex -> continuousIndex
    HashTable<label,word> timeIndexList_;

    HashTable<word, label> revTimeIndexList_;

    // create a data structure for the time values
    //  constant will not be contained
    //  continuousIndex -> time.value()
    HashTable<label,scalar> timeValueList_;

    label contTimeIndex;

    const label lowerSeqLim_;
    const label upperSeqLim_;

    scalar checkTimeStep();

    inline label getVolScalarFieldIndex(word, label) const;
    inline label getVolVectorFieldIndex(word, label) const;
    inline label getSurfaceScalarFieldIndex(word, label) const;

    void readFieldSeries();

    void readTimeSeries();

    void fetchStateForDataBase(label);

    void updateStorageUsageList(label);

    void readNewSnapshot(label, label);

    Switch checkSkipZero();

public:

    //- Runtime type information
    TypeName("gerhardsRecModel");


    // Constructors

        //- Construct from components
        gerhardsRecModel
        (
            const dictionary& dict,
            recBase& base
        );

    // Destructor

        ~gerhardsRecModel();

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

    label numRecFields() const;
    label numDataBaseFields() const;

    void updateRecFields();

    void writeRecMatrix() const;

private:

    /*
        Number of fields kept in the dataBase, this number is smaller, or equal, as the
        number of snapshots.
    */
    const label numDataBaseFields_;

    /*
        As the data base holds less states than there are snapshots, we need to translate
        the index of the snapshot to an appropriate index of the dataBase.
    */
    List<label> storageIndex_;

    /*
        As we repeatedly need to replace fields in the dataBase with new snapshots from
        disk, we need a way to judge, which slot to vacate for the new snapshot.
        Here, we follow the lead of paging algorithms used by OS'es in memory management.

        The storageUsageList_ represents the recentness of a slot's use; with 1 for the
        most recently used element, and higher values for less recently used ones.
    */
    List<label> storageUsageList_;

    // this class does not read all the fields
    List<PtrList<volScalarField> > volScalarFieldList_;
    List<PtrList<volVectorField> > volVectorFieldList_;
    List<PtrList<surfaceScalarField> > surfaceScalarFieldList_;

    // book-keeping
    label nrOfReadsFromDisk_;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "gerhardsRecModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
