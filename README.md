# Pipe-Flow-Verification
Verifying a globally explicit workflow of FSI+Erosion using Hagen-Poiseiulle pipe flow

This work implements an FSI scheme which includes OML updates due to erosion at prescribed time increments.
There are three workflows included in this repository. The first uses ABAQUS user subroutines UTRACLOAD and UEXTERNALDB to compute pressure and shear stress loads using Hagen-Poiseiulle pipe flow equations. A UMESHMOTION routine then applies the erosion field, computed as a function of shear stress.
The second workflow is an implementation of cosimulation between ABAQUS and StarCCM+. There is no erosion in this implementation, just strict FSI.
The third and final workflow includes both cosimulation and erosion using UMESHMOTION.

ABAQUS is required for all FEA runs (1,2,3)
StarCCM+ is required for cosimulation (2,3)
Matlab is required for analytical verifications (1,2,3)
