#ifndef ACTUATORDISKFORCEVEL_H_
#define ACTUATORDISKFORCEVEL_H_

#include "fvCFD.H"
#include "fvOption.H"
#include "uniformDimensionedFields.H"
#include "cellSetOption.H"

namespace Foam 
{
namespace fv // Finite Volume OpenFoam class
{
class actuatorDiskForceVel// Create new class for actuatorDiskForce
:
    public fv::cellSetOption // Retrieve class region form cellSetOption operation
{ 

public:
	
    //- Runtime type information
    TypeName("actuatorDiskForceVel"); // Name of class object

    scalar interiorRadius_; // Inputs
    scalar exteriorRadius_;
    vector startPoint_;
    vector endPoint_;
    scalar Ct_;
    scalar Cp_;
    scalar rpm_;
    vector V_;

    // ----------------------------------Constructors----------------------------------- //
    actuatorDiskForceVel
    (
        const word& name,
        const word& modelType,
        const dictionary& dict,
        const fvMesh& mesh
    );

    actuatorDiskForceVel(const actuatorDiskForceVel&) = delete;
 
    // No copy assignment
    void operator=(const actuatorDiskForceVel&) = delete;

    virtual ~actuatorDiskForceVel() = default; // Destructor

    // --------------------------------Member Functions--------------------------------//

    virtual void addSup // Add body force term to equations of motion
    (
        fvMatrix<vector>& eqn,
        const label fieldi
    );

    virtual bool read(const dictionary& dict); // Boolean for whether class has been called upon and should therefore be read? (not sure)

    void calcActuatorDiskVolForce // Calculates volume force at every cell and adds it to the volume force vector field.
    (
        const scalar& interiorRadius, 
        const scalar& exteriorRadius, 
        const vector& startPoint, 
        const vector& endPoint, 
        const scalar& Ct, 
        const scalar& Cp, 
        const scalar& rpm,
        const vector& V,
        const vectorField& U,
        vectorField& VolumeForce
    ); 
    
    void writeVTK(); // writes geometry to a vtk file

    scalar calcAxialForce // Calculates axial force as a function of radius and density
    (
        const scalar& r, 
        const scalar& interiorRadius, 
        const scalar& exteriorRadius, 
        const vector& startPoint, 
        const vector& endPoint, 
        const scalar& T,
        const scalar& V,
        const scalar& U,
        const scalar& W
    ); 
	
    scalar calcTanForce // Calculates tangential force
    (
        const scalar& r, 
        const scalar& interiorRadius, 
        const scalar& exteriorRadius, 
        const vector& startPoint, 
        const vector& endPoint, 
        const scalar& Q,
        const scalar& V,
        const scalar& U
    );
};

}
}
#endif

// ********************************************************************************************* //