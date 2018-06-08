{
	// Now the turbulence update is included in turbulenceCorrect.H 
	//alphat = turbulence->nut()/turbulence->Prt();  
    //alphat.correctBoundaryConditions();            

    
    volScalarField alphaEff("alphaEff", turbulence->nu()/Pr + alphat);

    fvScalarMatrix TEqn
    (
        fvm::ddt(T)                  // time derivative
      + fvm::div(phi, T)             // convection
      - fvm::laplacian(alphaEff, T)  // diffusion (molecular + turbulent)
      - fvc::div(qwall)              // temperature flux at boundary
      - SourceT                      // forcing (either ideal of from meso)
     ==
        fvOptions(T)
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve();

    fvOptions.correct(T);

    rhok = 1.0 - ((T - TRef)/TRef);      //Boussinesq bouyancy density

}
