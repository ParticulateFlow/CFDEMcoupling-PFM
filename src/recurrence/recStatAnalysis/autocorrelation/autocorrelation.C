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
\*---------------------------------------------------------------------------*/

#include "error.H"
#include "Random.H"
#include "autocorrelation.H"
#include "recModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(autocorrelation, 0);

addToRunTimeSelectionTable
(
    recStatAnalysis,
    autocorrelation,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
autocorrelation::autocorrelation
(
    const dictionary& dict,
    recBase& base
)
:
    recStatAnalysis(dict,base),
    propsDict_(dict.subDict(typeName + "Props")),
    fieldname_(propsDict_.lookup("fieldname")),
    fieldtype_(propsDict_.lookup("fieldtype")),
    delaySteps_(propsDict_.lookupOrDefault<label>("delaysteps",0)),
    refCell_(0),
    refPoint_(propsDict_.lookup("refPoint")),
    normalize_(propsDict_.lookupOrDefault<bool>("normalize",true)),
    autoCorrField_
    (   IOobject
        (
            "autoCorrField",
            base.mesh().time().timeName(),
            base.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        base.mesh(),
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0), 0.0)
    ),
    suppressMatrixAndPath_(propsDict_.lookupOrDefault<bool>("suppressMatrixAndPath",false))
{
    if (fieldtype_ != "scalar" && fieldtype_ != "vector")
    {
        FatalError <<"fieldtype needs to be either scalar or vector.\n" << abort(FatalError);
    }

    scalar delayTime = delaySteps_ * base.recM().recTimeStep();
    autoCorrField_.rename("autoCorrField"+name(delayTime));

    refCell_ = base.mesh().findCell(refPoint_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

autocorrelation::~autocorrelation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void autocorrelation::init()
{

}

void autocorrelation::statistics()
{
    autocorr();
}



void autocorrelation::autocorr()
{
    if (fieldtype_ == "scalar") scalarFieldList_.transfer(base_.recM().exportVolScalarFieldList(fieldname_));
    else vectorFieldList_.transfer(base_.recM().exportVolVectorFieldList(fieldname_));

    label tmax = base_.recM().totRecSteps();
    for (label ti = delaySteps_; ti < tmax; ti++)
    {
        forAll(autoCorrField_, cellI)
        {
            autoCorrField_[cellI] += autocorrSummand(ti-delaySteps_,ti,refCell_,cellI);
        }
    }

    autoCorrField_ /= (tmax - delaySteps_);

    autoCorrField_ -= meanProd();


    if (normalize_)
    {
        volScalarField meanProd(autoCorrField_);
        if (fieldtype_ == "scalar")
        {
            volScalarField aveField = base_.recM().exportVolScalarFieldAve(fieldname_);
            forAll(meanProd, cellI)
            {
                meanProd[cellI] = aveField()[cellI] * aveField()[refCell_];
            }
        }
        else
        {
            volVectorField aveField = base_.recM().exportVolVectorFieldAve(fieldname_);
            forAll(meanProd, cellI)
            {
                meanProd[cellI] = aveField()[cellI] & aveField()[refCell_];
            }
        }
        autoCorrField_ /= meanProd;
    }
    else
    {
        dimensionSet fieldDim(0,0,0,0,0);
        if (fieldtype_ == "scalar")
        {
            fieldDim.reset(scalarFieldList_[0].dimensions());
        }
        else
        {
            fieldDim.reset(vectorFieldList_[0].dimensions());
        }

        fieldDim.reset(fieldDim * fieldDim);

        autoCorrField_.dimensions().reset(fieldDim);
    }

    autoCorrField_.write();
}

scalar autocorrelation::autocorrSummand(label t1, label t2, label c1, label c2)
{
    scalar res;
    if (fieldtype_ == "scalar") 
    {
        res = scalarFieldList_[t1][c1] * scalarFieldList_[t2][c2];
    }
    else
    {
        res = vectorFieldList_[t1][c1] & vectorFieldList_[t2][c2];
    }
    return res;
}

volScalarField autocorrelation::meanProd()
{
    volScalarField meanProd(autoCorrField_);
    if (fieldtype_ == "scalar")
    {
        volScalarField aveField = base_.recM().exportVolScalarFieldAve(fieldname_);
        forAll(meanProd, cellI)
        {
            meanProd[cellI] = aveField()[cellI] * aveField()[refCell_];
        }
    }
    else
    {
        volVectorField aveField = base_.recM().exportVolVectorFieldAve(fieldname_);
        forAll(meanProd, cellI)
        {
            meanProd[cellI] = aveField()[cellI] & aveField()[refCell_];
        }
    }
    return meanProd;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
