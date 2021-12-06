//creating the fields according to the recurrence dictionary

IOdictionary recProperties_
(
    IOobject
    (
        "recProperties0",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

List<wordList> fieldsDict_(recProperties_.lookup("fieldsPairs"));

wordList fieldNames(fieldsDict_.size());

for(int i = 0; i < fieldsDict_.size(); i++)
{

    fieldNames[i]= fieldsDict_[i][0];
}

Info<< "\n list of the fields: \n" << fieldNames << endl;

//reading coherent velocity field name
label k = findIndex(fieldNames,"coh_velocity");

if (k < 0)
{
    FatalError <<"\n No field is defiened for the coherent velocity\n" << abort(FatalError);
}
const word Ucoh_pair = fieldsDict_[k][1];

    volVectorField UcohRec
    (
        IOobject
        (
            Ucoh_pair,
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero",dimensionSet(0, 1, -1, 0, 0),vector::zero)
    );

//reading incoherent velocity field name
k = findIndex(fieldNames,"inc_velocity");

if (k < 0)
{
    FatalError <<"\n No field is defiened for the incoherent velocity\n" << abort(FatalError);
}
const word Uinc_pair = fieldsDict_[k][1];

    volVectorField UincRec
    (
        IOobject
        (
            Uinc_pair,
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero",dimensionSet(0, 1, -1, 0, 0),vector::zero)
    );

//reading coherent turb kinetic energy field name
k = findIndex(fieldNames,"kSGS_coh");
if (k < 0)
{
    FatalError <<"\n No field is defiened for the coherent subgrid-scale turbulent kinetic energy\n" << abort(FatalError);
}
const word kSGScoh_pair = fieldsDict_[k][1];

    volScalarField kcohRec
    (
        IOobject
        (
            kSGScoh_pair,
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",dimensionSet(0, 2, -2, 0, 0),0.0)
     );

//reading incoherent turb kinetic energy field name
k = findIndex(fieldNames,"kSGS_inc");
if (k < 0)
{
    FatalError <<"\n No field is defiened for the coherent subgrid-scale turbulent kinetic energy\n" << abort(FatalError);
}
const word kSGSinc_pair = fieldsDict_[k][1];

    volScalarField kincRec
    (
        IOobject
        (
            kSGSinc_pair,
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",dimensionSet(0, 2, -2, 0, 0),0.0)
     );

// calculated fields
Info<< "\nCreating cell volume field\n" << endl;

volScalarField delta
(
    IOobject
    (
        "delta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
   mesh,
   dimensionedScalar("delta", dimLength, 0.0)
);

delta.primitiveFieldRef()=pow(mesh.V(),1.0/3.0);
delta.write();


volVectorField URec
(
    IOobject
    (
        "URec",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
Info<< "\nCreating turb kinetic energy field\n" << endl;

volScalarField kRec
(
    IOobject
    (
        "kRec",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    kcohRec+kincRec
);

// check if there is any negative values
forAll(kRec, cellI)
{
    if (kRec[cellI] < SMALL)
    {
        kRec[cellI] = 0.0;
    }

}

const fvPatchList& patches = mesh.boundary();
forAll(patches, patchI)
{
        kRec.boundaryFieldRef()[patchI] = 0.0;
}


kRec.write();

Info<< "\nCreating turb viscosity field\n" << endl;
volScalarField nutRec
(
    IOobject
    (
        "nutRec",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    sqrt(kRec)*delta*0.094
);

nutRec.write();


Info<< "Calculating face flux field phiRec\n" << endl;
surfaceScalarField phiRec
 (
     IOobject
     (
        "phiRec",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
     ),
     linearInterpolate(URec) & mesh.Sf()
 );

    phiRec.write();

    singlePhaseTransportModel laminarTransport(URec, phiRec);

    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(URec, phiRec, laminarTransport)
    );
  
    dimensionedScalar Sc("Sc", dimless, laminarTransport);
    dimensionedScalar Sct("Sct", dimless, laminarTransport);
  
volScalarField alphat
(
    IOobject
    (
        "alphat",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    nutRec/Sct
);

 // create the scalar field
 Info<< "Creating scalar transport field\n" << endl;

    volScalarField C
    (
        IOobject
        (
            "C",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    fvScalarMatrix CEqn(C, dimless*dimVolume/(dimTime));


 Info<< "reading clockProperties\n" << endl;

    IOdictionary clockProperties
    (
            IOobject
            (
                "clockProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
    );

    autoPtr<clockModel> myClock
    (
            clockModel::New
            (
                clockProperties,
                mesh.time()
            )
    );
