// check which recProperties dicts are present, read them in and construct a PtrList of recBases
// names for dicts can be "recProperties" or "recPropertiesN" where N in {0, 1, ...} 

#include "error.H"

word dictName = "recProperties";
wordList recPropertiesList(0);
PtrList <recBase> recBases(0);
label maxDictNumber = 100;

{
    IOdictionary recPropDict
    (
        IOobject
        (
            dictName,
            mesh.time().constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    if (recPropDict.headerOk())
    {
        recPropertiesList.append(dictName);
    }
}


for (label counter = 0; counter < maxDictNumber; counter++)
{
    word dictNameIter = dictName + Foam::name(counter);
    IOdictionary recPropDict
    (
        IOobject
        (
            dictNameIter,
            mesh.time().constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    if (recPropDict.headerOk())
    {
        recPropertiesList.append(dictNameIter);
    }
}

if (recPropertiesList.size() == 0)
{
    FatalError << "no recProperties dicts found" << endl;
}

else
{
    Info << "found " << recPropertiesList.size() << " dicts with names " << recPropertiesList << endl;
}


for (label counter = 0; counter < recPropertiesList.size(); counter++)
{
    recBases.append( new recBase(mesh, recPropertiesList[counter]));
}
