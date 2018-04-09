{
	// update modeled turbulent thermal cond. (as Prt is a dynamic value 
	// in Sogachev model, this step is carried by the turbulence-> correct step
	//alphat = turbulence->nut()/turbulence->Prt();  
	//alphat = turbulence->nut()/Prt;
	
    // update the turbulence model. It also updates the modeled turbulent
    // thermal conductivity
    turbulence->correct();

    // update the boundary conditions on modeled turbulent thermal conductivity.
    alphat.correctBoundaryConditions();

    // update the boundary conditions on surface heat flux.
    qwall.correctBoundaryConditions();

}
