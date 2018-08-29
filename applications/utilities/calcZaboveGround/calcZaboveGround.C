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
		
		scalar currentMin = 1e100;
		forAll(faces, facei)
		{
			// Patch face position vector (in x-y plane)
			Vector<double> p(faces[facei].x(),faces[facei].y(),0);

			if(Foam::mag(r-p) < currentMin)
			{
				zGround = faces[facei].z();
				currentMin = Foam::mag(r-p);
			}
		}
		
		
		// Find patch face whose x and y extents contain r.x() and r.y()
		// (This is less accurate but should be faster...)
/*
		bool found = false;
		while (!found)
		{
			// Get vertices that define current patch face
			labelList vertices = patchFaces[face];

			scalar minX = 1e100;
			scalar maxX = -1e100;
			scalar minY = minX;
			scalar maxY = maxX;
			
			forAll (vertices, vertexI)
			{
				point vertex = mesh.points()[vertices[vertexI]];
		
				if (vertex.x() >= maxX)
					maxX = vertex.x();
				if (vertex.x() <= minX)
					minX = vertex.x();
				if (vertex.y() >= maxY)
					maxY = vertex.y();
				if (vertex.y() <= minY)
					minY = vertex.y();
			}

//			Info << minX << endl;
//			Info << maxX << endl;
//			Info << minY << endl;
//			Info << maxY << endl;

			// Check if r.x() and r.y() are contained by x and y extents of face
			if (r.x() >= minX && r.x() <= maxX && r.y() >= minY && r.y() <= maxY)
			{
				zGround = faces[face].z();
				found = true;
//				Info << "r = " << r << endl;
//				Info << "patch face centre = " << faces[face] << endl;
//				Info << "zGround = " << zGround << endl;
			}
			else
				face++;

			if (face >= patchFaces.size())
				face = 0;
//			Info << "face = " << face << endl;
		}*/

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
