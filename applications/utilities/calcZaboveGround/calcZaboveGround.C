/******************************************************************************* 
	Application:
		calcZaboveGround.C

	Description:
		Utility to create the Z field with the values of the height above ground

	Author:
		R. Chavez based on the ABLInitialize code of J. Sumner
	
	Date:	
		May 2013

	Last update:
        Ago 2015  Roberto Chavez
        Jan 2018  Roberto Chavez
*******************************************************************************/

#include "fvCFD.H"
//#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"
//#include "RASModel.H"
#include <stdio.h>


#include "wallDist.H"

/************************************ MAIN ************************************/
int main(int argc, char *argv[])
{
	// Default include statements
	#include "setRootCase.H"
	#include "createTime.H"
	#include "createMesh.H"
//    #include "createFields.H"

    volScalarField zAGL
    (
        IOobject
        (
            "zAGL",
            runTime.timeName(),
            mesh,
            //IOobject::NO_READ,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh
	//dimensionedScalar("Z", dimensionSet(0,1,0,0,0,0,0), 0.0)  // by default we assing 0 as ground level at faces
    );

	// Load ABL dictionary
	IOdictionary ABLProperties_
	(
		IOobject
		(
			"ABLProperties",
			runTime.constant(),
			mesh,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		)
	);

	// Read surface-layer parameters
	scalar z0 = readScalar(ABLProperties_.lookup("z0"));

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                           Initialize internal fields
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	const volVectorField& cellCentres = mesh.C();
	const surfaceVectorField& cellFaces = mesh.Cf();

	// Estimate reference ground height
	Info << "Setting internal fields for " << zAGL.size() << " cells with respect to ground patch" << endl;	
	

	// Get reference to ground patch
	label groundPatchID = mesh.boundaryMesh().findPatchID("ground");
	const fvsPatchVectorField& faces = cellFaces.boundaryField()[groundPatchID];
	const faceList& patchFaces = mesh.boundaryMesh()[groundPatchID];


	// Get reference to inlet patch
	//label inletPatchID = mesh.boundaryMesh().findPatchID("inlet");
	//const fvsPatchVectorField& facesInlet = cellFaces.boundaryField()[inletPatchID];
	//const faceList& patchFacesInlet = mesh.boundaryMesh()[inletPatchID];


	//Info << "faceinlet = " << facesInlet << endl;
	//Info << "faceinlet2222 = " << patchFacesInlet << endl;

	//Info << "faceinlet = " << min(facesInlet.z()) << endl;

	/*label face1=0;

	forAll(patchFacesInlet, faceI)
	{
		labelList vertices2 = patchFacesInlet[face1];

		Info << " noooooo " << faceI << " numerco"<<
		vertices2 << endl;
	}*/
	vector up = vector(0,0,1);
	//up.z() = 1.0;

	volScalarField d = mesh.C() & up;
    Info << "Calculating wall distance..." << endl;
    wallDist dWall(mesh);
    d = dWall;



	scalar progress(0.0);	
	int cellCount = 0;
	scalar totalCells = zAGL.size();
	printf ("  %2.1f%%\r",progress);
	fflush(stdout);
	label face = 0;

	forAll(zAGL, celli)
	{
		// Cell position vector (in x-y plane)
		Vector<double> r(cellCentres[celli].x(),cellCentres[celli].y(),0);

		scalar zGround(0.0);		

		// Find patch face with closest centre to current point in x-y plane
		// (This is a very slow but sure method...)
		
		/*scalar currentMin = 1e100;
		forAll(faces, facei)
		{
			// Patch face position vector (in x-y plane)
			Vector<double> p(faces[facei].x(),faces[facei].y(),0);

			if(Foam::mag(r-p) < currentMin)
			{
				zGround = faces[facei].z();
				currentMin = Foam::mag(r-p);
			}
		}*/
		
		//Search for 4 nearest solutions and use linear interpolation for the z value 

		scalar dist0 = 1e100;
		scalar dist1 = 1e100;		
		scalar dist2 = 1e100;
		scalar dist3 = 1e100;
		scalar zGround0 = -999;
		scalar zGround1 = -999;		
		scalar zGround2 = -999;
		scalar zGround3 = -999;				
		scalar temp;
		
		forAll(faces, facei)
		{
			// Patch face position vector (in x-y plane)
			Vector<double> p(faces[facei].x(),faces[facei].y(),0);

			//if p is close enough, then add it to the sorted list
			if(Foam::mag(r-p) < dist0)
			{
				zGround0 = faces[facei].z();
				dist0 = Foam::mag(r-p);
				if(dist0 < dist1)	
				{
					temp = dist1;
					dist1 = dist0;
					dist0 = temp;
					temp = zGround0;
					zGround0 = zGround1;
					zGround1 = temp;
					if(dist1 < dist2)
					{
						temp = dist1;
						dist1 = dist2;
						dist2 = temp;
						temp = zGround2;
						zGround2 = zGround1;
						zGround1 = temp;
						if(dist2 < dist3)
						{
							temp = dist3;
							dist3 = dist2;
							dist2 = temp;
							temp = zGround2;
							zGround2 = zGround3;
							zGround3 = temp;
						}
					}
				}
			}
		}
		//inverse distance weighted, linear interpolation (Shepard interpolation for p=1)
		if(dist3 == 0)
		{
			zGround = zGround3;
		}
		else
		{
			scalar sum1 = zGround0/dist0 + zGround1/dist1 + zGround2/dist2 + zGround3/dist3;
			scalar sum2 = 1/dist0 + 1/dist1 + 1/dist2 + 1/dist3;
			zGround = sum1/sum2;
		}

		scalar deltaZ = cellCentres[celli].z() - zGround;
		if (deltaZ < 0)
		{
	    //Info << "\nCell: " << celli << "--> z < zGround!  z = " << cellCentres[celli].z() << "; zGround = " << zGround << endl;
			cellCount ++;			
			deltaZ = z0;
		}

		zAGL[celli] = deltaZ;

		
		if (celli/totalCells*100.0 - progress > 0.1)
		{
			progress = celli/totalCells*100.0;
			printf ("  %2.1f%%\r",progress);
			fflush(stdout);
		}
		
		//for testing in flat you can use z directly
		//zAGL[celli] = mesh.C()[celli][2];

		Info << "zGround1 = " << zAGL[celli] << "\n dWall  =  " << d[celli]<< endl;

	}
	//Info << endl;
	Info << "\nPercent of cells for which z < zGround = " << cellCount/totalCells*100.0 << "%%\n" << endl;
	
	// Write out Z
	Info << "Writing field\n" << endl;
	zAGL.write();
	Info << "End" << endl;

	return(0);
}
