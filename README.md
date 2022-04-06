# MATLAB-Euler-PML
## Perfectly Matched Layer for Linearised Euler Equations

Code developed as part of Master's research project (heavily WIP, code tends to diverge towards end of simulation - could be numerical instability or implementation error)

### Overview

- Linearised Euler Equations in 2D with mean flow in x-direction
- Stable Perfectly Matched Layer formulation for x-layer, y-layer, and corners
- Dispersion Relation Preserving FD scheme for spatial discretisation (TAM, 1992)
- 4th Order Runge Kutta for temporal discretisation (to be extended to LDDRK56 in the future by (HU, 1995))
- 2N Storage scheme
