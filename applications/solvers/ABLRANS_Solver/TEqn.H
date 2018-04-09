{
    //alphat = turbulence->nut()/Prt;
	//alphat = turbulence->nut()/turbulence->Prt();        //myyy shit
    //alphat.correctBoundaryConditions();    //commented by MChurch
    
    volScalarField alphaEff("alphaEff", turbulence->nu()/Pr + alphat);

    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
      + fvm::div(phi, T)
      - fvm::laplacian(alphaEff, T)
      - fvc::div(qwall)
      - SourceT                                //M. Church.
     ==
        fvOptions(T)
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve();

    fvOptions.correct(T);

    rhok = 1.0 - ((T - TRef)/TRef);      //Boussinesq bouyancy density

}
