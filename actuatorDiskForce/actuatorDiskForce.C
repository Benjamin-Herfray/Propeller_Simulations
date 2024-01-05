/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2017 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

	Description
    Process adds volume forces to actuator disk region. Actuator disk parameters are given in fvOptions as follows:

disk
{
    type            actuatorDiskForce;

	fields          (U);    // Names of fields on which to apply source
    startPoint      (0 0 0);
    endPoint		(0.02 0 0);
    interiorRadius	0.0; // Optional; default is 0
    exteriorRadius	0.125;

    Ct      0.155;
    Cp      0.00; // Optional; default is 0
    rpm     5425;

    selectionMode   geometric;
    selection
    {
        rotor
        {
            action      new;
            source      cylinder;
            p1          $startPoint;
            p2          $endPoint;
            radius      $exteriorRadius;
			innerRadius $interiorRadius;
        }
    } 
}

\*--------------------------------------------------------------------------------------------------*/

#include "actuatorDiskForce.H"

#include "faceAreaPairGAMGAgglomeration.H"
#include "fvMesh.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam 
{
namespace fv
{
	defineTypeNameAndDebug(actuatorDiskForce, 0); // Debugging
	addToRunTimeSelectionTable(option, actuatorDiskForce, dictionary);
}
}


Foam::fv::actuatorDiskForce::actuatorDiskForce
(
	const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh

)  // Construct an actuatorDisk object
:
	fv::cellSetOption(name, modelType, dict, mesh), // Initialize fvOptions object with selected cells
	
	// Extract fvOptions parameters
	interiorRadius_(coeffs_.getOrDefault<scalar>("interiorRadius", 0)), // 0 is default if a value is not given
	exteriorRadius_(coeffs_.get<scalar>("exteriorRadius")),
	startPoint_(coeffs_.get<vector>("startPoint")),
	endPoint_(coeffs_.get<vector>("endPoint")),
	Ct_(coeffs_.get<scalar>("Ct")),
	Cp_(coeffs_.getOrDefault<scalar>("Cp", 0)), // 0 is default if a value is not given
	rpm_(coeffs_.get<scalar>("rpm"))
{
	coeffs_.readEntry("fields", fieldNames_); // Read the name of the field to apply source to.

	if (fieldNames_.size() != 1) // The only source should be U
    {
        FatalErrorInFunction
			<< "settings are:" << fieldNames_ << exit(FatalError);
    }

	writeVTK(); // Create a VTK file for visualization in paraView

	fv::option::resetApplied();// Resize/reset applied flag list for all fieldNames_ entries.
}

// ********************************************************************************************************************* //

void Foam::fv::actuatorDiskForce::addSup // Add volume force to equations of motion
(
	fvMatrix<vector>& eqn,
    const label fieldi
)
{
	volVectorField volumeForce // Create vector field
	(
		IOobject // Create readable data object for volume force
		(
			name_ + ":volumeForce",
			mesh_.time().timeName(),
			mesh_
		),
		mesh_,
		dimensionedVector(eqn.dimensions()/dimVolume, Zero)
	);

	calcActuatorDiskVolForce(interiorRadius_, exteriorRadius_, startPoint_, endPoint_, Ct_, Cp_, rpm_, volumeForce); // Calculate the volume force Field

	eqn += volumeForce; // Add vector field to equations of motion

	if (mesh_.time().writeTime()) // Write volume force
    {
    	volumeForce.write();
    }
}

bool Foam::fv::actuatorDiskForce::read(const dictionary& dict) // Not sure, but seems to be important. Best guess is that this ensures routine is not read if class is not implemented in fvOptions
{
	NotImplemented;
 
    return false;
}

scalar Foam::fv::actuatorDiskForce::calcAxialForce // Computes force in axial direction
(
	const scalar& r,
	const scalar& interiorRadius, 
	const scalar& exteriorRadius, 
	const vector& startPoint, 
	const vector& endPoint, 
	const scalar& T
) 
{ 
	scalar PI = constant::mathematical::pi;

	scalar axialForce = 0.0; // Initialize
	scalar x = r / exteriorRadius; // Distance normalized by propeller radius
	scalar xh = interiorRadius / exteriorRadius; // Normalized hub radius
	scalar rs = sqrt((pow(x, 2) - pow(xh, 2)) / (1 - pow(xh, 2))); // Calculate scaled radius: far-wake equivalent distance along from inner radius to outer radius
	scalar Ax = T / (mag(endPoint - startPoint) * PI * (pow(exteriorRadius, 2) - pow(interiorRadius, 2)) * 32 / 105); // Calculate force scaling

	axialForce = Ax * rs * sqrt(1.0 - rs); // Calculate distribution
	
	return axialForce;
}

scalar Foam::fv::actuatorDiskForce::calcTanForce // Tangential force
(
	const scalar& r,
	const scalar& interiorRadius, 
	const scalar& exteriorRadius, 
	const vector& startPoint, 
	const vector& endPoint, 
	const scalar& Q
)
{
	scalar PI = constant::mathematical::pi;

	scalar tangentialForce = 0.0;
	scalar x = r / exteriorRadius; // Distance normalized by propeller radius
	scalar xh = interiorRadius / exteriorRadius; // Normalized hub radius
	scalar rs = sqrt((pow(x, 2) - pow(xh, 2)) / (1 - pow(xh, 2))); // Calculate scaled radius: far-wake equivalent distance along from inner radius to outer radius
	scalar At = Q / (mag(endPoint - startPoint) * PI * (pow(exteriorRadius, 2) - pow(interiorRadius, 2)) * 32 / 105); // Calculate force scaling

	tangentialForce = At * rs * sqrt(1.0 - rs) / r; // Original

	return tangentialForce;
}

void Foam::fv::actuatorDiskForce::calcActuatorDiskVolForce // Calculate the volume force vector field applied to disk region
(
	//const fvMesh& iMesh,
	const scalar& interiorRadius,
	const scalar& exteriorRadius,
    const vector& startPoint,
	const vector& endPoint,
	const scalar& Ct,
	const scalar& Cp,
	const scalar& rpm,
	vectorField& volumeForce
) 
{
	scalar T = Ct * pow(rpm, 2) * pow(exteriorRadius, 4) / 225; // Calculate thrust and torque from propeller parameters
	scalar Q = Cp * pow(rpm, 2) * pow(exteriorRadius, 5) / (225 * constant::mathematical::pi);

	//Loop over all cells in the actuator disk region
	forAll(cells_, i) 
	{
		const label celli = cells_[i]; // New label for the i'th cell
		
		vector vecAxial(endPoint - startPoint); // Create vector along axis of cylinder
		if(mag(vecAxial) != 0.0) { // Normalize
			vecAxial /= mag(vecAxial);
		}
		else {
			FatalErrorInFunction
				<< "Warning: The centerline tangent has zero length." << exit(FatalError);
		}
	
		vector vecStartLineToPoint(mesh_.C()[celli] - startPoint); // Vector from centerline of the starting face to a given point
		scalar radialDist = mag(vecStartLineToPoint - vecAxial * (vecStartLineToPoint & vecAxial)); // Calculate radial component of point by taking the magnitude of the vector component orrthogonal to the axis line found using Gram-Schmidt orthogonalization

		vector tangentialDirection = vecAxial ^ (vecStartLineToPoint - vecAxial * (vecStartLineToPoint & vecAxial)); // Create vector in tangential angular direction by taking the cross-product between the axis unit vector and the radial vector
		if(mag(tangentialDirection) != 0) {
			tangentialDirection /= mag(tangentialDirection); // Normalize if the point is not on the central axis. If it is, we take the tangential force to be 0.		
		}
		
		volumeForce[celli] = vecAxial * calcAxialForce(radialDist, interiorRadius, exteriorRadius, startPoint, endPoint, T); // Add to force vector field by calculating the thrust produced at a given radius and applying it in the axial direction
			
		vector tanForce = tangentialDirection * calcTanForce(radialDist, interiorRadius, exteriorRadius, startPoint, endPoint, Q); // Calculate the density-normalized tangential force vector for the cell.
		volumeForce[celli] += tanForce;
	}
}

void Foam::fv::actuatorDiskForce::writeVTK() // Write the outer surface of the actuator disk to a VTK file so that it can be visualized in Paraview.
{ 
	// Initialization of Variables
	FILE *file; // Create the file
	char fileName[100];

	unsigned int numCells = 100; // The cylindrical surface is visualized as numCells rectangular surfaces
	unsigned int numPoints = 2*numCells; // 2*numCells points are required; numCells points at each end of the cylinder
	unsigned int numInts = 5*numCells; // Number of indices needed in the the VTK file; each surface has 4 corner points, so we need 4 corner indices + the index of the surface = 5 indices per surface
	
	vectorField points(numPoints,vector::zero); // Initialize vector field

	// Create an ON basis for the cylinder
	vector vecAxial(endPoint_ - startPoint_); // vector along the cylinder axis
	if(mag(vecAxial) != 0.0) { // Normalize vector for nontrivial magnitude
		vecAxial /= mag(vecAxial);
	}
	else {
		Info << "Warning: The centerline tangent has zero length.\n"; // Return warning in case axial vector is 0
		return;
	}

	//We need to find a vector in the radial direction. This can be any vector as long as it points in the radial direction.
	//First try with (1 0 0) and see if we can project it onto the normal plane of the actuator disk resulting in a vector in
	//the radial direction.
	vector vecRadial(1.0, 0.0, 0.0);
	vecRadial -= (vecRadial & vecAxial) * vecAxial; // Gram Orthogonalization by removing component in axial direction
	
	if(mag(vecRadial) < SMALL) {
		// If we enter this if statement, our guess (1 0 0) was parallel to the centerline of the actuator disk. Then
		// we try (0 1 0) instead. Since (1 0 0) was parallel to the centerline, (0 1 0) will not be parallel to
		// the centerline.
		vecRadial.x() = 0.0;
		vecRadial.y() = 1.0;
		vecRadial.z() = 0.0;
		vecRadial -= (vecRadial & vecAxial)*vecAxial;			
	}

	vecRadial /= mag(vecRadial); // Normalize vector

	vector vecRadial2 = vecAxial ^ vecRadial; // Second radial vector constructed by cross product
	
	scalar xLocal = 0.0, yLocal = 0.0;
	
	//Compute points on first side of disk region
	double phi = 0.0;
	for(unsigned int i = 0; i < numCells; i++) { // Iterate for numCells number of cells
		xLocal = exteriorRadius_*cos(phi); // Compute x and y components at outer rim of cylinder for a given angle phi
		yLocal = exteriorRadius_*sin(phi);
		
		points[i] = startPoint_ + xLocal*vecRadial + yLocal*vecRadial2; // Add point to list of points
		phi += (1.0/double(numCells))*2*constant::mathematical::pi; // Advance angle by 2pi/numCells each time
	}
	
	//Compute points on second side of disk region
	phi = 0.0;
	for(unsigned int i = 0; i < numCells; i++) {
		xLocal = exteriorRadius_*cos(phi); 
		yLocal = exteriorRadius_*sin(phi);
		
		points[numCells + i] = endPoint_ + xLocal*vecRadial + yLocal*vecRadial2;
		phi += (1.0/double(numCells))*2*constant::mathematical::pi;
	}
	
	// write points to vtk file.
	sprintf(fileName,"actuatorDiskForce.vtk");
	file = fopen(fileName,"w");

	fprintf(file,"# vtk DataFile Version 3.0\n");
	fprintf(file,"Analytical surface of actuator disk. \n");
	fprintf(file,"ASCII\n");

	fprintf(file,"DATASET UNSTRUCTURED_GRID\n");
	fprintf(file,"POINTS %i float\n",numPoints);

	for(int i = 0; i < points.size(); i++) { // Write out locations of all points
		fprintf(file,"%e %e %e\n",points[i].x(),points[i].y(),points[i].z());
	}
	
	fprintf(file,"CELLS %i %i\n",numCells,numInts);

	for(unsigned int i = 0; i < numCells-1; i++) {
		fprintf(file,"%i %i %i %i %i \n",4,i,i+numCells,i+numCells+1,i+1);
	}
	fprintf(file,"%i %i %i %i %i \n",4,numCells-1,2*numCells-1,numCells,0);

	fprintf(file,"CELL_TYPES %i\n",numCells);

	for(unsigned int i = 0; i < numCells; i++) {
		fprintf(file,"%i\n",9);
	}

	fclose(file);
}

// ************************************************************************************************************************************************************************ //