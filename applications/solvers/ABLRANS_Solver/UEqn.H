    // Solve the momentum equation

    #include "computeCoriolisForce.H"    //SOWFA implementation

    fvVectorMatrix UEqn
    (
        fvm::ddt(U)                 //time derivative
      + fvm::div(phi, U)            //convection
      + turbulence->divDevReff(U)   //Divergence of shear-rate tensor (div(Tau))
      - fCoriolis                   //Coriolis force 
      - SourceU                   //Driving force (either ideal or from mesoscale)
     ==
        fvOptions(U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                  - ghf*fvc::snGrad(rhok)    // buoyancy force
                  - fvc::snGrad(p_rgh)       // modified pressure gradient
                )*mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
