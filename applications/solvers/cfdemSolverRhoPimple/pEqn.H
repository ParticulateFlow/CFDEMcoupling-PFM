rho = thermo.rho();
rho = max(rho, rhoMin);
rho = min(rho, rhoMax);
rho.relax();

volScalarField rAU(1.0/UEqn().A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rhoeps*voidfraction*rAU));

volVectorField HbyA("HbyA", U);
HbyA = rAU*UEqn().H();

surfaceScalarField phiUs("phiUs", fvc::interpolate(rhoeps*rAU*Ksl*Us)& mesh.Sf());

if (pimple.nCorrPISO() <= 1)
{
    UEqn.clear();
}


if (pimple.transonic())
{
//  transonic version not implemented yet
  
/*    surfaceScalarField phid
    (
        "phid",
        fvc::interpolate(psi)
       *(
            (fvc::interpolate(HbyA) & mesh.Sf())
          + rhorAUf*fvc::ddtCorr(rho, U, phi)/fvc::interpolate(rho)
        )
    );

    fvOptions.makeRelative(fvc::interpolate(psi), phid);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          + fvm::div(phid, p)
          - fvm::laplacian(rhorAUf, p)
          ==
            fvOptions(psi, p, rho.name())
        );

        fvOptions.constrain(pEqn);

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi == pEqn.flux();
        }
    }
*/
}
else
{
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (
            (fvc::interpolate(rhoeps*HbyA) & mesh.Sf())
   //       + rhorAUf*fvc::ddtCorr(rho, U, phi)
        )
    );

    fvOptions.makeRelative(fvc::interpolate(rho), phiHbyA);
    
    // shamelessly stolen from interFoam
    setSnGrad<fixedFluxPressureFvPatchScalarField>
    (
        p.boundaryField(),
        (
            phiHbyA.boundaryField()
	  + phiUs.boundaryField()
          - fvOptions.relative((mesh.Sf().boundaryField() & U.boundaryField()))
	    * rhoeps.boundaryField()
        )/(mesh.magSf().boundaryField()*rhorAUf.boundaryField())
    );

    while (pimple.correctNonOrthogonal())
    {
        // Pressure corrector
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi*voidfraction, p)
          + fvc::div(phiHbyA)
          - fvm::laplacian(rhorAUf, p)
	  + fvc::div(phiUs)
          ==
            fvOptions(psi, p, rho.name())
        );

        fvOptions.constrain(pEqn);

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pEqn.flux()+phiUs;

            // Explicitly relax pressure for momentum corrector
            p.relax();

            if (UfromPhi)
	        U = fvc::reconstruct( phi/fvc::interpolate(rhoeps) );
	    else
                U = HbyA + rAU*(voidfraction*fvc::reconstruct((pEqn.flux())/rhorAUf) + Ksl*Us);
	    
            U.correctBoundaryConditions();
            fvOptions.correct(U);
	    
divphi = fvc::div(phi);
divphi2 = fvc::div(linearInterpolate(rho*U*voidfraction) & mesh.Sf());
gP = fvc::grad(p);
gP2 = fvc::reconstruct((pEqn.flux())/rhorAUf);
        }
    }
}



#include "rhoEqn.H"

#include "compressibleContinuityErrsPU.H"

// Explicitly relax pressure for momentum corrector
//p.relax();

// Recalculate density from the relaxed pressure
rho = thermo.rho();
rho = max(rho, rhoMin);
rho = min(rho, rhoMax);
rho.relax();
Info<< "rho max/min : " << max(rho).value()
    << " " << min(rho).value() << endl;

    // old version to get U from p; can have problems at boundaries where grad(p) != snGrad(p)
//U = HbyA - rAU*(voidfraction*fvc::grad(p)-Ksl*Us);
//U.correctBoundaryConditions();
//fvOptions.correct(U);
K = 0.5*magSqr(U);

if (thermo.dpdt())
{
    dpdt = fvc::ddt(voidfraction,p);
}
