    singlePhaseTransportModel laminarTransport(U, phi);

    // Reference temperature [K]
    dimensionedScalar TRef(laminarTransport.lookup("TRef"));

    // Laminar Prandtl number
    dimensionedScalar Pr(laminarTransport.lookup("Pr"));

    // Turbulent Prandtl number
    //dimensionedScalar Prt(laminarTransport.lookup("Prt"));
    // What I called the "static" turbulent Prandtl number is only a 
    // dummy value as in Sogachev's formulation this is a dynamic value
    // based on atmopheric stability is later updated in turbulence->correct()
    dimensionedScalar Prts(laminarTransport.lookup("Prts"));
