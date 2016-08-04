{
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));

    //With included voidfraction
    surfaceScalarField rAUepsf("rAUepsf", fvc::interpolate(rAU*voidfraction));
    volScalarField rAUepsSq("rAUepsSq", rAU * voidfraction * voidfraction);
    surfaceScalarField rAUepsSqf("rAUepsSqf", fvc::interpolate(rAUepsSq));
    volVectorField Ueps("Ueps", U * voidfraction);
    
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA*voidfraction)
      + fvc::interpolate(voidfraction*rho*rAU)*fvc::ddtCorr(U, phi)
    );
    MRF.makeRelative(phiHbyA);    
    adjustPhi(phiHbyA, U, p_rgh);

    surfaceScalarField phig
    (
        (
            mixture.surfaceTensionForce()
          - ghf*fvc::snGrad(rho)
        )*rAUepsSqf*mesh.magSf()
    );
    
    surfaceScalarField phiS
    (
        fvc::flux(voidfraction*Us) * rAUf*fvc::interpolate(Ksl)
        //(fvc::interpolate(voidfraction*Us) & mesh.Sf()) * rAUf*fvc::interpolate(Ksl)
    );

    phiHbyA += phig + phiS;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, Ueps, phiHbyA, rAUepsf, MRF);
    /*setSnGrad<fixedFluxPressureFvPatchScalarField>
    (
        p_rgh.boundaryField(),
        (
            phiHbyA.boundaryField()
          - (mesh.Sf().boundaryField() & U.boundaryField()) * voidfraction.boundaryField()
        )/(mesh.magSf().boundaryField()*rAUepsSqf.boundaryField())
    );*/
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUepsSqf, p_rgh) == fvc::div(phiHbyA) + fvc::ddt(voidfraction)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - p_rghEqn.flux();

            p_rgh.relax();

            U = HbyA + voidfraction*rAU
              * fvc::reconstruct
                (
                    (phig - p_rghEqn.flux())/rAUepsSqf
                )
              + fvc::reconstruct(phiS)/voidfraction;
            
            U.correctBoundaryConditions();
        }
    }

    #include "continuityErrs.H"

    p == p_rgh + rho*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }
}
