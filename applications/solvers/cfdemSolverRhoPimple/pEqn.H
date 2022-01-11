rho = thermo.rho();
#if OPENFOAM_VERSION_MAJOR < 5
rho = max(rho, rhoMin);
rho = min(rho, rhoMax);
rho.relax();
rhoeps = rho*voidfraction;
#else
rhoeps = rho*voidfraction;

// Thermodynamic density needs to be updated by psi*d(p) after the
// pressure solution
const volScalarField psip0(psi*p);
#endif

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", (modelType=="A")?fvc::interpolate(voidfraction*rhoeps*rAU):fvc::interpolate(rhoeps*rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

surfaceScalarField phiUs("phiUs", fvc::interpolate(rhoeps*rAU*Ksl*Us)& mesh.Sf());

if (pimple.nCorrPISO() <= 1)
{
    tUEqn.clear();
}

surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::interpolate(rhoeps)*fvc::flux(HbyA)
  + rhorAUf*fvc::ddtCorr(rhoeps, U, phi)
);

if (pimple.transonic())
{
//  transonic version not implemented yet
}
else
{
    // flux without pressure gradient contribution
    phi = phiHbyA + phiUs;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, rhoeps, U, phi, rhorAUf);

#if OPENFOAM_VERSION_MAJOR >= 5
    fvScalarMatrix pDDtEqn
    (
        fvc::ddt(rhoeps)
      + psi*voidfraction*correction(fvm::ddt(p))
      + fvc::div(phi)
     ==
        fvOptions(psi, p, rho.name())
    );
#endif

    while (pimple.correctNonOrthogonal())
    {
        // Pressure corrector
#if OPENFOAM_VERSION_MAJOR < 5
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi*voidfraction, p)
          + fvc::div(phi)
          - fvm::laplacian(rhorAUf, p)
          ==
            fvOptions(psi, p, rho.name())
        );
#else
        fvScalarMatrix pEqn(pDDtEqn - fvm::laplacian(rhorAUf, p));
#endif

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi += pEqn.flux();
        }
    }
}

// Thermodynamic density update
#if OPENFOAM_VERSION_MAJOR >= 5
thermo.correctRho(psi*p - psip0);
#endif

#include "rhoEqn.H"
#include "compressibleContinuityErrsPU.H"

// Explicitly relax pressure for momentum corrector
p.relax();
Info<< "p max/min/ave : " << max(p).value()
    << " " << min(p).value() << " " << average(p).value() << endl;

if (modelType=="A")
{
    U = HbyA - rAU*(voidfraction*fvc::grad(p)-Ksl*Us);
}
else
{
    U = HbyA - rAU*(fvc::grad(p)-Ksl*Us);
}
U.correctBoundaryConditions();
fvOptions.correct(U);
K = 0.5*magSqr(U);

// Recalculate density from the relaxed pressure
#if OPENFOAM_VERSION_MAJOR >= 5
if (pressureControl.limit(p))
{
    p.correctBoundaryConditions();
}
rho = thermo.rho();
#else
rho = thermo.rho();
rho = max(rho, rhoMin);
rho = min(rho, rhoMax);
rho.relax();
#endif
rhoeps = rho*voidfraction;

Info<< "rho max/min/ave : " << max(rho).value()
    << " " << min(rho).value() << " " << average(rho).value() << endl;

if (thermo.dpdt())
{
    dpdt = fvc::ddt(voidfraction,p);
}
