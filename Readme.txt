Developed by Dr.Kularathna, modified by Chihun Sung.

FORMULATION 1: u-p formulation with Euler time integration

Governing equation are
1. Mixture moementum balance
2. Fluid momentum balance

Implement
1. non-incremental projection mehtod
2. incremental projection method
3. rigid contact algorithm
4. contact detection algorithm using edge-to-edge distances among particles
5. GIMP

Stress and Strain are computed by using the B matrix computed at particle level (no enhancement such as BBar is used)

Undrained boundary can be applied by making only the velicity in the direction of the ourward normal to the boundary to zero.

Lumped mass matrix is used in the formulation
