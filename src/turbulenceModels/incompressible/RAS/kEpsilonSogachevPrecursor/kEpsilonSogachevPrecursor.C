/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "kEpsilonSogachevPrecursor.H"
#include "addToRunTimeSelectionTable.H"

#include "backwardsCompatibilityWallFunctions.H"

#define PI 3.14159265359

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(kEpsilonSogachevPrecursor, 0);
addToRunTimeSelectionTable(RASModel, kEpsilonSogachevPrecursor, dictionary);

// * * * * * * * * * * * * Protected Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kEpsilonSogachevPrecursor::kEpsilonSogachevPrecursor
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    RASModel(modelName, U, phi, transport, turbulenceModelName),
    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            coeffDict_,
            0.0256    //default value of Dett and Ettling
        )
    ),
    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            coeffDict_,
            1.13    //default value of Dett and Ettling
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            coeffDict_,
            1.9    //default value of Dett and Ettling
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            coeffDict_,
            1.2987    //default value of Dett and Ettling
        )
    ),
    sigmaK_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaK",
            coeffDict_,
            0.7407    //default value of Dett and Ettling
        )
    ),
    kappa_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kappa",
            coeffDict_,
            0.4    //default value of Dett and Ettling
        )
    ),

    // Construct ABL dictionary and read input data for ABL model
    ABLProperties_
    (
       IOobject
       (
           "ABLProperties",
           runTime_.constant(),
           mesh_,
           IOobject::MUST_READ,
           IOobject::NO_WRITE
       )
    ),

    transportProperties_
    (
       IOobject
       (
           "transportProperties",
           runTime_.constant(),
           mesh_,
           IOobject::MUST_READ,
           IOobject::NO_WRITE
       )
    ),

    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateK("k", mesh_)
    ),
    epsilon_
    (
        IOobject
        (
            "epsilon",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateEpsilon("epsilon", mesh_)
    ),
    nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateNut("nut", mesh_)
    ),
    
     lMY_  // Construct and initialize the maximum mixing length of Mellor-Yamada
     (
        IOobject
        (
            "lMY",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    	mesh_,
    	dimensionedScalar("lMY",dimensionSet(0,1,0,0,0,0,0),40.0)
     ),

    lm_   // Construct and initialize mixing length 
    (
        IOobject
        (
            "lm",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pow(Cmu_,0.75)*pow(k_,1.5)/epsilon_
     ),

     C1ast_   // Construct and Initialize C1 constant according to Apsley & Castro (R. Chavez)
     (
        IOobject
        (
            "C1ast",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        //mesh_,
        //C1_
        C1_ + (C2_-C1_)*lm_/lMY_
     ),

     alphaB_  //Construct alphaB coeff and initialize it as neutral/unstable case
     (
        IOobject
        (
            "alphaB",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        1.0 - lm_/lMY_
     ),

     Prts_
     (
        transportProperties_.lookupOrDefault<dimensionedScalar>("Prts",
        	dimensionedScalar("Prts", dimensionSet(0,0,0,0,0,0,0), 0.74 ))
     ),

     Prt_  //Construct the turb. Prandlt number and initialize it as neutral/unstable case
     (
        IOobject
        (
            "Prt",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    	mesh_,
    	Prts_
     ),

     C3_ 
     (
        IOobject
        (
            "C3",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (C1_-C2_)*alphaB_ + 1.0
     ),

     Rig_
     (
        IOobject
        (
            "Rig",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    	mesh_,
    	dimensionedScalar("Rig",dimensionSet(0,0,0,0,0,0,0),10.0)
     ),

     RiG_
     (
        IOobject
        (
            "RiG",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    	mesh_,
    	dimensionedScalar("RiG",dimensionSet(0,0,0,0,0,0,0),10.0)
     ),

     gradTz_
     (
        IOobject
        (
            "gradTz",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    	mesh_,
    	dimensionedScalar("gradTz",dimensionSet(0,-1,0,1,0,0,0),-1.0)
     ),

     wTz_
     (
        IOobject
        (
            "wTz",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    	mesh_,
    	dimensionedScalar("wTz",dimensionSet(0,1,-1,1,0,0,0),1.0)
     ),

     uStar_
     (
        IOobject
        (
            "uStar",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    	mesh_,
    	dimensionedScalar("uStar",dimensionSet(0,1,-1,0,0,0,0),1.0)
     ),

     L_
     (
        IOobject
        (
            "L",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    	mesh_,
    	dimensionedScalar("L",dimensionSet(0,1,0,0,0,0,0),-1.0E3)
     ),

     g_
     (
        mesh_.lookupObject<uniformDimensionedVectorField>("g")
     ),

     strongLimitForEpsilon_
     (
    	//ABLProperties_.lookupOrDefault<bool>("strongLimitForEpsilon", false)
    	readBool(ABLProperties_.lookup("strongLimitForEpsilon")) // force read to avoid ambiguity
     ),

     limitTopAmbientValues_
     (
    	readBool(ABLProperties_.lookup("limitTopAmbientValues")) // force read to avoid ambiguity
     ),

     minHeightForAmbientValues_
     (
    	readScalar(ABLProperties_.lookup("minHeightForAmbientValues")) // force read to avoid ambiguity
     ),
     
	 lmAmb_
	 (
	 	ABLProperties_.lookupOrDefault<dimensionedScalar>("lmAmb",
	 		dimensionedScalar("lmAmb", lm_.dimensions(), 0.05 ))
     ),

     kAmb_
     (
    	ABLProperties_.lookupOrDefault<dimensionedScalar>("kAmb",
    		dimensionedScalar("kAmb", k_.dimensions(), 1.0E-4 ))
     ),

     epsilonAmb_
     (
    	 "epsilonAmb", pow(Cmu_,0.75)*pow(kAmb_,1.5)/lmAmb_
     ),

     calculatefc_
     (
        readBool(ABLProperties_.lookup("calculatefc")) // force read to avoid ambiguity
     ),

     latitude_
     (
        ABLProperties_.lookupOrDefault<scalar>("latitude",1.0)
     ),

     omegaEarth_
     (
        "omegaEarth", dimensionSet(0,0,-1,0,0,0,0), 7.292115053926e-05
     ),

     fc_
     (
		ABLProperties_.lookupOrDefault<dimensionedScalar>("fc",
				dimensionedScalar("fc", 2.0*omegaEarth_*sin(latitude_*PI/180) ))
     ),

     Ug_
     (
	    ABLProperties_.lookupOrDefault<dimensionedVector>("Ug",
	    	dimensionedVector("Ug", dimensionSet(0,1,-1,0,0,0,0), vector(0,0,0) ))
     ),

     lmax_
     (
    	"lmax", 0.00027*mag(Ug_)/fc_
     ),

     TRef_
     (
        transportProperties_.lookupOrDefault<dimensionedScalar>("TRef",
        	dimensionedScalar("TRef", dimensionSet(0,0,0,1,0,0,0), 300.0  ))
     ),

     lMYreferencePoint_
     (
        ABLProperties_.lookupOrDefault<vector>("lMYreferencePoint",vector(0,0,0))
     ),

     iCellMY
     (
        mesh_.findNearestCell(lMYreferencePoint_)
     ),

     zDir_(0,0,1),

     alphaBcoeff_( "alphaBcoeff", 1 + (C2_-1)/(C2_-C1_) ),

     cellCentres_( "cellCentres" , mesh_.C() ),

     z_( "z" , mesh_.C() & zDir_ )

{

   	 if (calculatefc_)
   	 {
   	 	fc_ = 2.0*omegaEarth_*sin(latitude_*PI/180);
   	 }

   	 Info << " USING kEpsilon MODEL OF SOGACHEV ET AL. 2012 for precursor simulations " << nl
   	 <<	nl << "Coriolis factor fc                    =  " <<    fc_ << nl
   	 << nl << "Graviational constant g               = " << g_ <<nl
   	 << nl << "Reference Temperature Tref            = " << TRef_.value()<<nl
   	 << nl << "Model constant Cmu                    = " << Cmu_.value()<<nl
   	 << nl << "model constant sigmak                 =" << sigmaK_.value()<<nl
   	 << nl << "Turbulent static Prandlt number Prts  = " << Prts_.value() << nl
   	 << nl << "Blackadar value of lmax               =" << lmax_.value()<<nl
   	 << nl << "Coriolis factor being calculated?     = " << calculatefc_<<nl
   	 << nl << "Mellor-Yamada lmax reference location =" << lMYreferencePoint_<<nl
   	 << nl << "Strong limitation of epsilon with lMY?=" << strongLimitForEpsilon_<<nl
   	 << nl << "Limit with ambient values in upper atmosphere? =" << limitTopAmbientValues_<<nl
   	 << endl;
   	 
	 if (limitTopAmbientValues_)
	 {
	     Info << "Strong Clipping so that min mixing length =" << lmAmb_.value()<<"m" << nl
		 << nl << "Upper atmosphere starting at " << minHeightForAmbientValues_ << "m" << nl
   	     << nl << "k ambient =" << kAmb_ << nl
   	     << nl << "epsilon ambient =" << epsilonAmb_ << nl
   	     << endl;
	 }


    bound(k_, kMin_);
    bound(epsilon_, epsilonMin_);

	if (limitTopAmbientValues_)
	 {
	    Info << "Strong Clipping with ambient values" << endl;
		bound(k_, kAmb_);
		bound(epsilon_, epsilonAmb_);
	 }

    nut_ = Cmu_*sqr(k_)/epsilon_;   
    nut_.correctBoundaryConditions();

    printCoeffs();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Reynolds Stress tensor
tmp<volSymmTensorField> kEpsilonSogachevPrecursor::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_)),
            k_.boundaryField().types()
        )
    );
}

tmp<volSymmTensorField> kEpsilonSogachevPrecursor::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -nuEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}

tmp<fvVectorMatrix> kEpsilonSogachevPrecursor::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(T(fvc::grad(U))))
    );
}

tmp<fvVectorMatrix> kEpsilonSogachevPrecursor::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    volScalarField muEff("muEff", rho*nuEff());

    return
    (
      - fvm::laplacian(muEff, U)
      - fvc::div(muEff*dev(T(fvc::grad(U))))
    );
}

bool kEpsilonSogachevPrecursor::read()
{
    if (RASModel::read())
    {
        Cmu_.readIfPresent(coeffDict());
        C1_.readIfPresent(coeffDict());
        C2_.readIfPresent(coeffDict());
        sigmaEps_.readIfPresent(coeffDict());
        sigmaK_.readIfPresent(coeffDict());   //Added RCA
        kappa_.readIfPresent(coeffDict());   //Added RCA

        return true;
    }
    else
    {
        return false;
    }
}


void kEpsilonSogachevPrecursor::correct()
{
    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }

    // calc lMY
    lMY_= 0.075*fvc::domainIntegrate(z_*pow(k_, 0.5)) / fvc::domainIntegrate( pow(k_,0.5));
    //Info << "lMYpointsSum ="<<lMY_.internalField()<<endl;
   
    // compute u*
    volVectorField gradUdotz = T(fvc::grad(U_)) & zDir_;
    uStar_= sqrt(mag(-nut_*gradUdotz));

    //Temperature gradients and L
    const volScalarField& T_ = mesh_.lookupObject<volScalarField>("T");
    
    gradTz_ = fvc::grad(T_) & zDir_;   // dtheta/dz
    wTz_ = -nut_/Prt_ * gradTz_;    //vertical kinematic heat flux [K m s^-1]
    // minimum wTz value
    //volScalarField wTzMin_("wTzMin", dimensionedScalar("wTzMin", wTz_.dimensions(),1.0E-6) );
    dimensionedScalar wTzMin_("wTzMin", wTz_.dimensions(),1.0E-6);

    // Monin obukhov length
    L_ = -pow(uStar_,3) * TRef_/( kappa_*mag(g_)* (sign(wTz_)*max(mag(wTz_),wTzMin_)) );

    //Shear production
    volScalarField G(GName(), nut_*2*magSqr(symm(fvc::grad(U_))));  

    // minimum value for G
    //volScalarField Gmin_("Gmin", dimensionedScalar("Gmin", G.dimensions(),1.0E-12) );
    dimensionedScalar Gmin_("Gmin", G.dimensions(),1.0E-9);

    //Buoyancy production
    volScalarField B("B", (1.0/TRef_) * (nut_ /Prt_) * (g_ & (fvc::grad(T_))) );

    // Update epsilon and G at the wall
    epsilon_.boundaryField().updateCoeffs();  // define appropiately where should I put this term

    // Richardson Number
   	RiG_ = -B/( G + mag( (alphaB_*B/Prt_) ) );
   	
   	Rig_ = -B/ (max(G,Gmin_));

   	scalar gamma = 16.0;

	//models "constants" alphaB and Prt based on stability Rig
    forAll(RiG_, cellI)
    {
        //if (Rig_[cellI] >= 0.0)
        if (RiG_[cellI] >= 0.0)   //stable
        {
        	alphaB_[cellI] = 1.0 - lm_[cellI]/lMY_[cellI];
        	Prt_[cellI]    = Prts_.value();
        }
        else   //unstable
        {
        	alphaB_[cellI] = 1.0 - alphaBcoeff_.value()*lm_[cellI]/lMY_[cellI];

        	Prt_[cellI]    = Prts_.value()*pow( 1.0 - gamma*RiG_[cellI] , -0.25);
        }
    }

    //CHECAR LOS PUTOS GMAX Y GMIN

    // Apsley & Castro 1997 
    C1ast_ = C1_ + (C2_-C1_)* lm_/lMY_;
    
    C3_ = (C1_-C2_)*alphaB_ + 1.0;


    /*
    Info << nl
 	   << " k     : "
       << " min: " << gMin(k_.internalField())
       << " ave: " << k_.weightedAverage(mesh_.V()).value()
       << " max: " << gMax(k_.internalField()) << nl
 	   << " eps   : "
       << " min: " << gMin(epsilon_.internalField())
       << " ave: " << epsilon_.weightedAverage(mesh_.V()).value()
       << " max: " << gMax(epsilon_.internalField()) << nl
       << " T     : "
       << " min: " << gMin(T_.internalField())
       << " ave: " << T_.weightedAverage(mesh_.V()).value()
       << " max: " << gMax(T_.internalField()) << nl
       << " B     : "
       << " min: " << gMin(B.internalField())
       << " ave: " << B.weightedAverage(mesh_.V()).value()
       << " max: " << gMax(B.internalField())  << nl
       << " RiG: "
       << " min: " << gMin(RiG_.internalField())
       << " ave: " << RiG_.weightedAverage(mesh_.V()).value()
       << " max: " << gMax(RiG_.internalField()) << nl
       << " gradTz : "
       << " min: " << gMin(gradTz_.internalField())
       << " ave: " << gradTz_.weightedAverage(mesh_.V()).value()
       << " max: " << gMax(gradTz_.internalField()) << nl
       << endl;
    */

    // Update epsilon and G at the wall
    //epsilon_.boundaryField().updateCoeffs();

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
     ==
        //C1_*G*epsilon_/k_
        (C1ast_*G + C3_*B)*epsilon_/k_        // Shear and buoyancy dissipation terms
      - fvm::Sp(C2_*epsilon_/k_, epsilon_)
    );

    epsEqn().relax();
    epsEqn().boundaryManipulate(epsilon_.boundaryField());

    solve(epsEqn);
    bound(epsilon_, epsilonMin_);       // brute force clipping  (R. Chavez)
    
    // Clipping for ambient values
    if (limitTopAmbientValues_)
    {
        forAll(z_, cellI)
		{
			if (z_[cellI] >= minHeightForAmbientValues_)
			{
				epsilon_[cellI] = max(epsilon_[cellI] , epsilonAmb_.value() );
			}
		}
    }

    // Clipping following Sogachev 2012:"Ensuring epsilon according to lmax"
    if (strongLimitForEpsilon_)
    {
        epsilon_ = max(epsilon_ , pow(Cmu_,0.75)*pow(k_,1.5)/lMY_);
    }

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        G + B                      //shear + buoyancy
      - fvm::Sp(epsilon_/k_, k_)
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, kMin_);
    
    // Clipping for ambient values
    if (limitTopAmbientValues_)
    {   
        forAll(z_, cellI)
		{
			if (z_[cellI] >= minHeightForAmbientValues_)
			{
				k_[cellI] = max(k_[cellI] , kAmb_.value());
			}
		}
    }


    // compute mixing length
    lm_ = pow(Cmu_,0.75)*pow(k_,1.5)/epsilon_;
    
    // Re-calculate viscosity
    nut_ = Cmu_*sqr(k_)/epsilon_;
    nut_.correctBoundaryConditions();
    


    // Update the thermal conductivity.
    volScalarField& alphat_ = const_cast<volScalarField&>(U().db().lookupObject<volScalarField>("alphat"));
    alphat_ = nut_/Prt_;


    /*
    // Update the Turbulent Prandtl Number
    volScalarField& Prt_ = const_cast<volScalarField&>(U().db().lookupObject<volScalarField>("Prt"));
	*/

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
