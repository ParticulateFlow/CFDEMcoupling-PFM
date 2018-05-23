/*---------------------------------------------------------------------------*\
    CFDEMcoupling academic - Open Source CFD-DEM coupling

    Contributing authors:
    Thomas Lichtenegger
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
\*---------------------------------------------------------------------------*/

#include "error.H"
#include "sqrDiffNorm.H"
#include "recModel.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(sqrDiffNorm, 0);

addToRunTimeSelectionTable
(
    recNorm,
    sqrDiffNorm,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
sqrDiffNorm::sqrDiffNorm
(
    const dictionary& dict,
    recBase& base
)
:
    recNorm(dict, base),
    propsDict_(dict.subDict(typeName + "Props")),
    normConstant_(propsDict_.lookupOrDefault<scalar>("normConstant",-1.0)),
    fieldType_(propsDict_.lookup("fieldType")),
    fieldName_(propsDict_.lookup("fieldName"))
{
    if (propsDict_.found("readRecMat"))
    {
        readRecMat_ = bool(readBool(propsDict_.lookup("readRecMat")));
        if (readRecMat_ && propsDict_.found("recMatName"))
        {
            recMatName_ = word(propsDict_.lookup("recMatName"));
        }
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

sqrDiffNorm::~sqrDiffNorm()
{}

// * * * * * * * * * * * * * protected Member Functions  * * * * * * * * * * * * //

void sqrDiffNorm::computeRecMatrix()
{
    if (readRecMatrix()) return;  

    Info << nl << type() << ": computing recurrence matrix\n" << endl;

    label totNumRecSteps = base_.recM().numRecFields();
    SymmetricSquareMatrix<scalar>& recurrenceMatrix( base_.recM().recurrenceMatrix() );

    scalar normIJ(0.0);
    scalar maxNormIJ(0.0);

    /*
        total number of computed elements: total number of matrix entries
            minus the main diagonal entries, divided by two,
            since the matrix is symmetric
    */
    label size = (totNumRecSteps*(totNumRecSteps-1))/2;
    label counter = 0;
    label percentage = 0;


    label N(this->base_.recM().numRecFields());
    label M(this->base_.recM().numDataBaseFields());

    if (verbose_)
    {
        Info << " N = " << N << ",  M = " << M << endl;
    }

    label ti(0);
    label tj(0);
    label tmp(-1);

    for (int j=0; j<=N/(M-1); j++)
    {
        for (int i=0; i<totNumRecSteps; i++)
        {
            if (verbose_)
            {
                Info << " i = " << i << ",  j = " << j << endl;
            }

            if(counter >= 0.1 * percentage * size)
            {
                Info << "\t" << 10 * percentage << " \% done" << endl;
                percentage++;
            }

            for (int k=i; k<i+(M-1); k++)
            {
                ti = i;
                tj = j*(M-1) + k;

                if (ti > tj)
                {
                    tmp = ti;
                    ti = tj;
                    tj = tmp;
                }

                // skip coordinates outside the recurrence space
                if (ti >= N or tj >= N)
                {
                    continue;
                }

                // start
                // skip main diagonal and upper half
                if (ti >= tj)
                {
                    recurrenceMatrix[ti][tj] = 0;
                    continue;
                }

                if (verbose_)
                {
                    Info << " Doing calculation for element "
                        << ti << " " << tj << endl;
                }

                counter++;

                // compute elements
                if (fieldType_ == "volScalarField")
                {
                    normIJ = normVSF(ti,tj);
                }
                else if (fieldType_ == "volVectorField")
                {
                    normIJ = normVVF(ti,tj);
                }
                else if (fieldType_ == "surfaceScalarField")
                {
                    normIJ = normSSF(ti,tj);
                }
                else
                {
                    FatalError
                        << "sqrDiffNorm: Unknown field type " << fieldType_
                        << abort(FatalError);
                }

                recurrenceMatrix[ti][tj] = normIJ;

                if (normIJ > maxNormIJ)
                {
                    maxNormIJ = normIJ;
                }
                // end
            }
        }
    }


    // normalize matrix and copy lower to upper half
    if(normConstant_ > 0.0) maxNormIJ = normConstant_;

    for(label ti=0;ti<totNumRecSteps;ti++)
    {
        for(label tj=0;tj<totNumRecSteps;tj++)
        {
            if (ti >= tj) continue;

            if (recurrenceMatrix[ti][tj] < 0)
            {
                FatalErrorInFunction << "Error in computation of recurrence matrix!"
                    << nl << "Negative elements encountered. This should not happen!"
                    << abort(FatalError);
            }

            recurrenceMatrix[ti][tj] /= maxNormIJ;
            recurrenceMatrix[tj][ti] = recurrenceMatrix[ti][tj];
        }
    }

    Info << "\nComputing recurrence matrix done\n" << endl;
}

scalar sqrDiffNorm::normVSF(label ti, label tj)
{
    const volScalarField& t1( base_.recM().exportVolScalarField(fieldName_,ti) );
    const volScalarField& t2( base_.recM().exportVolScalarField(fieldName_,tj) );
    dimensionedScalar tNorm( fvc::domainIntegrate( sqr( t1 - t2 ) ) );

    return tNorm.value();
}

scalar sqrDiffNorm::normVVF(label ti, label tj)
{
    const volVectorField& t1( base_.recM().exportVolVectorField(fieldName_,ti) );
    const volVectorField& t2( base_.recM().exportVolVectorField(fieldName_,tj) );
    dimensionedScalar tNorm( fvc::domainIntegrate( magSqr( t1 - t2 ) ) );

    return tNorm.value();
}

scalar sqrDiffNorm::normSSF(label ti, label tj)
{
    const surfaceScalarField& t1( base_.recM().exportSurfaceScalarField(fieldName_,ti) );
    const surfaceScalarField& t2( base_.recM().exportSurfaceScalarField(fieldName_,tj) );
    volVectorField t12 (fvc::reconstruct( t1-t2 ) );
    dimensionedScalar tNorm( fvc::domainIntegrate( magSqr( t12 ) ) );

    return tNorm.value();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
