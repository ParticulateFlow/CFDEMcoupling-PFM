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

#include "error.H"
#include "initMultiLayers.H"
#include "addToRunTimeSelectionTable.H"
#include "dataExchangeModel.H"
#include <sys/stat.h>
#include <unistd.h>
#include <set>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(initMultiLayers, 0);

addToRunTimeSelectionTable
(
        chemistryModel,
        initMultiLayers,
        dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
initMultiLayers::initMultiLayers
(
    const dictionary& dict,
    cfdemCloudEnergy& sm
)
:
    chemistryModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    interpolation_(propsDict_.lookupOrDefault<bool>("interpolation",false)),
    mesh_(sm.mesh()),
    verbose_(propsDict_.lookupOrDefault<bool>("verbose",false)),
    maxNumLayers_(0),
    maxNumParticlesPerType_(propsDict_.lookupOrDefault<label>("maxNumParticlesPerType",1000000)),
    numLayers_(propsDict_.lookupOrDefault<labelList>("numLayers",labelList(1,-1))),
    partTypes_(propsDict_.lookupOrDefault<labelList>("partTypes",labelList(1,-1))),
    searchLayers_(propsDict_.lookupOrDefault<labelList>("searchLayers",labelList(1,-1))),
    relRadiiRegName_(typeName + "relRadii"),
    filepath_(string(propsDict_.lookup("filepath"))),
    initialized_(false)
{
    for (label i=0; i<numLayers_.size(); i++)
    {
        if (numLayers_[i] > maxNumLayers_) maxNumLayers_=numLayers_[i];
    }
    if (maxNumLayers_ > 4)
    {
        FatalError<< "Currently, not more than four layers are supported." << abort(FatalError);
    }
    defaultRelRadii_ = vector(0.998, 0.995, 0.98);
    particleCloud_.checkCG(false);
    particleCloud_.registerParticleProperty<double**>(relRadiiRegName_,maxNumLayers_);

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

initMultiLayers::~initMultiLayers()
{
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //

bool initMultiLayers::init()
{
    double**& relRadii = particleCloud_.getParticlePropertyRef<double**>(relRadiiRegName_);
    // for all types of particles
    for(label types=0; types<partTypes_.size(); types++)
    {
        label type = partTypes_[types];

        labelList indices;
        vectorList positions, relradii;

        std::stringstream ss;
        ss << filepath_ << type;
        std::string filename = ss.str();

        if (access( filename.c_str(), F_OK ) == -1)
        {
            Info << "initMultiLayers: Data file " << filename << "not found. No initialisation of layer radii possible." << endl;
            return false;
        }
        Info << "\nReading" << endl;
        Info << "\t" << filename << endl;
        label numLines = readDump(filename, type, indices, positions, relradii);

        Info << "Binning particle displacements on mesh for type " << type << endl;

        volScalarField particlesInCell
        (
            IOobject
            (
                "particlesInCell",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimensionSet(0,0,0,0,0), 0.0)
        );

        vector position, relrad;
        ss.str("");
        ss << "relRadiiField" << type;
        std::string fieldname = ss.str();
        volVectorField relRadiiField
        (
            IOobject
            (
                fieldname,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("zero", dimensionSet(0,0,0,0,0), vector::zero)
        );

        label cellI = -1;
        for (label lineI = 0; lineI < numLines; lineI++)
        {
            position = positions[lineI];
            relrad = relradii[lineI];
            cellI = mesh_.findCell(position);
            if (cellI < 0)  continue;
            particlesInCell[cellI] += 1;
            relRadiiField[cellI] += relrad;
        }

        forAll (mesh_.C(), cellJ)
        {
            if (particlesInCell[cellJ] > 0.5)
            {
                relRadiiField[cellJ] /= particlesInCell[cellJ];
            }
        }

        // fill field


        // fill particle arrays
        label cellK = -1;
        label cellKoccupied = -1;
        label partType = -1;
        label listIndex = -1;
        vector relRadiiLoopVar = vector::zero;
        for(int index = 0;index < particleCloud_.numberOfParticles(); ++index)
        {
            cellK = particleCloud_.cellIDs()[index][0];
            partType = particleCloud_.particleType(index);
            if(cellK >= 0 && partType == type)
            {
                // look for the nearest occupied cell
                listIndex = getListIndex(partType);
                cellKoccupied = findNearestCellWithValue(cellK, particlesInCell,searchLayers_[listIndex]);
                if (cellKoccupied >= 0)
                {
                    relRadiiLoopVar = relRadiiField[cellKoccupied];
                }
                else
                {
                    relRadiiLoopVar = defaultRelRadii_;
                }
                relRadii[index][0] = 1.0;
                for (label i=1;i<numLayers_[listIndex];i++)
                {
                    relRadii[index][i] = relRadiiLoopVar.component(i-1);
                }
            }
        }
        if (verbose_) relRadiiField.write();
    }
    return true;
}

// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //

void initMultiLayers::execute()
{
    if(!initialized_)
    {
        if (init())
        {
            double**& relRadii = particleCloud_.getParticlePropertyRef<double**>(relRadiiRegName_);
            particleCloud_.dataExchangeM().giveData("relRadii","vector-atom", relRadii);
        }
        initialized_ = true;
    }
}

label initMultiLayers::findNearestCellWithValue(label refCell, volScalarField &particlesInCell, label searchLayers) const
{
    if (particlesInCell[refCell] > 0.5) return refCell;

    std::set<label> neighbors;
    std::set<label> newNeighbors;
    std::set<label> recentNeighbors;

    neighbors.insert(refCell);
    recentNeighbors.insert(refCell);

    label layersSearched = 0;

    while(true)
    {
        for (std::set<label>::iterator it=recentNeighbors.begin(); it!=recentNeighbors.end(); ++it)
        {
            labelList adjacent = mesh_.cellCells()[*it];
            label adj;
            for (label j=0; j<adjacent.size(); j++)
            {
                adj = adjacent[j];
                std::set<label>::iterator it2 = neighbors.find(adj);
                if (it2 == neighbors.end())
                {
                    newNeighbors.insert(adj);
                    neighbors.insert(adj);
                    if (particlesInCell[adj] > 0.5) return adj;
                }
            }
        }

        // if all cells have been searched and none was occupied, return -1; assumes that reasonable default value is available
        if (newNeighbors.size() == 0) return -1;
        recentNeighbors.clear();
        recentNeighbors = newNeighbors;
        newNeighbors.clear();
        layersSearched++;
        if (layersSearched >= searchLayers && searchLayers > 0) return -1;
    }
}

label initMultiLayers::getListIndex(label testElement) const
{
    for(label ind = 0; ind<partTypes_.size(); ind++)
    {
        if (partTypes_[ind] == testElement) return ind;
    }
    // testing
    Pout << "Cannot find list index for element " << testElement << endl;
    return -1;
}

label initMultiLayers::readDump(std::string filename, label type, labelList &indices, vectorList &positions, vectorList &relradii)
{
    #include <fstream>

    const label leadingLines = 9;
    label lineCounter = 0;
    label partIndex, partType;
    scalar x, y, z;
    scalar r0 = 1.0;
    scalar r1 = 1.0;
    scalar r2 = 1.0;
    scalar r3 = 1.0;

    indices.clear();
    positions.clear();
    relradii.clear();

    indices.setSize(maxNumParticlesPerType_);
    positions.setSize(maxNumParticlesPerType_);
    relradii.setSize(maxNumParticlesPerType_);

    std::ifstream file(filename);
    std::string str;
    while (std::getline(file, str))
    {
        if (lineCounter >= leadingLines)
        {
            sscanf(str.c_str(), "%d %d %lf %lf %lf %lf %lf %lf %lf", &partIndex, &partType, &x, &y, &z, &r0, &r1, &r2, &r3);
            if (partType != type)
            {
                FatalError<< "Particle of type " << partType << " detected in " << filename << abort(FatalError);
            }

            indices[lineCounter-leadingLines] = partIndex;
            positions[lineCounter-leadingLines] = vector(x,y,z);
            relradii[lineCounter-leadingLines] = vector(r1,r2,r3);
        }
        lineCounter++;
        if (lineCounter == maxNumParticlesPerType_) break;
    }

    label readLines = lineCounter - leadingLines;
    return readLines;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
