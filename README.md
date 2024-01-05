# Propeller_Simulations
A repository containing the propeller class file and simulations using the actuatorDiskSource custom OpenFoam class.

Due to size limitations, I have only included simulation setups, which must be run by users once they have integrated the actuatorDiskForce class.

I have also included Sensitivity and BoundaryTests directories with earlier simulations testing sensitivity to and numarical stability of different loading models and boundary conditions, respectively.

Simulations are classified first by propeller model (e.g. GWS5x4.3) followed by rotation speed in RPM and finally freestream velocity in m/s (written as "mps" is directory names).

Most cases include both a main directory and a preparation directory; the preparation directory (indicated by a "Pre" at the end) run on coarser meshes and are designed to produce approximate field solutions faster. These approximate fields can then be mapped as initial conditions to the corresponding main directory simulation with a finer mesh. Preparation simulations can be run using the Allrun script. Main directories require running the Allprep first (which can be done at any time), followed by the Allrun, which must be launched after the preparation simulation has finished.

Any questions or corrections are welcume at: benjamin (dot) herfray (at) mail (dot) mcgill (dot) ca.
