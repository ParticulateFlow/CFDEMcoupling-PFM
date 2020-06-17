# Recurrence-based, multi-phase turbulence modelling

This model implements recurrence-based turbulence models, i.e. the fundamental 
turbulent field quantities are read from the data base and are not solved for.
All derived field quantities are computed just in the same way as the proper 
turbulence models do. By deriving the recurrence-based turbulence models from 
somewhere up the family tree of OpenFOAM's turbulence model class hierarchy, 
the recurrence-based turbulence models are fully compatible with OpenFOAM's 
generic treatment of turbulence modelling, i.e. solvers and libraries interact 
with references to a generic base type of the actual turbulence model. Hence, 
solvers and libraries may remain blissfully ignorant of the actual turbulence 
model in use.

For laminar phases no special treatment is necessary, as the *laminar* 
turbulence model does not compute any fields.


## Development notes

The initial development covers only a small number of turbulence models.


## Notes on usage

The turbulence model in use for the recurrence run must be the recurrence-based 
equivalent of the turbulence model used for generating the data base, i.e. if 
the data base was computed using the *kEpsilon* model, then the recurrence solver 
is to employ the *recurrenceKEpsilon* turbulence model. This model will read 
the relevant model coefficients from the *turbulenceProperties* dictionary, and 
make sure that the turbulent fields `k` and `epsilon` are contained in the data 
base.

Whenever, the solver or a library calls `turbulence->nut()` to access the 
turbulent viscosity, the recurrence-based kEpsilon model will compute `nut` 
according to kEpsilon's relations `nut = Cmu*sqr(k)/epsilon`, with the fields 
`k` and `epsilon` being from the current snapshot provided by the recurrence model.

Thus, the fundamental turbulent field quantities of the employed turbulence model 
have to be added to the *volScalarFields* list in the `recProperties` dictionary 
controlling the recurrence model. This will ensure that the turbulent field 
quantities are read from the data base.


## Notes on the implementation

The base class implements the method `void setRecurrenceBasePtr(recBase*)`, which 
is used to give the recurrence-based turbulence models a reference (technically 
a pointer) to the recurrence model. Thus, after construction of the turbulence 
models and the recurrence model, `setRecurrenceBasePtr()` needs to be called as 
the pointer to the recurrence model is initialized by the constructor with `NULL`.
Trying to access the recurrence model from within the recurrence-based turbulence 
model prior to setting the pointer to the recurrence model with 
`setRecurrenceBasePtr()` will result in a segmentation fault.
In order to be able to call `setRecurrenceBasePtr()`, the generic reference to 
the turbulence model needs to be converted into a reference of the base class' 
type, i.e. `recurrenceTurbulenceModel`.
This unfortunate deviation from good standards, i.e. making full use of C++'s 
polymorphism, should be the only instance of having to use non-pretty hacks.
However, apart from initialisation, i.e. setting the pointer to the recurrence 
model, the recurrence-based turbulence models adhere to the generic interface of 
OpenFOAM's turbulence models, and can be used as any other turbulence model.


The concrete implementations, e.g. *recurrenceKEpsilon*, use the method 
`validate()` to check whether the underlying turbulent quantities are specified 
for use in the data base in the *volScalarFields* list in the `recProperties` 
dictionary. This method is part of the signature of the class `Foam::turbulenceModel`, 
which is the very base class of all turbulence models in OpenFOAM.
In proper turbulence models, this method is used to check whether the internal 
fields are properly initialized and to update all derived quantities.
In the solver, `validate()` must not be called prior to `setRecurrenceBasePtr()`, 
as validate accesses the recurrence model. The wrong order of function calls will 
result in a segmentation fault, as the pointer to the recurrence model is 
initialized by the constructor with `NULL`.


The concrete implementations, e.g. *recurrenceKEpsilon*, use the method 
`correct()` to update the turbulent field quantities from the data base, 
and in turn update the derived quantities, such as `nut`.
This method is part of the signature of the class `Foam::turbulenceModel`, 
which is the very base class of all turbulence models in OpenFOAM.
In proper turbulence models, this method is used to solve for the next time step.


## Compilation

Source OpenFOAM and simply compile with

```bash
./Allwclean
./Allwmake
```

The script `Allwclean` will clear all previous builds. This step is not needed for
first-time compilation. It is, however, recommended for subsequent compilations, as
it completely clears the slate. The script `Allwmake` will run the compilation for
the passive particle model.


## Required software

This model has been tested with the following versions of OpenFOAM:

* OpenFOAM-4.0
* OpenFOAM-5.0


