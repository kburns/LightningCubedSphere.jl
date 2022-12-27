# Lightning conformal cubed sphere

This code solves for the conformal map between a face of the cubed sphere, the gnomic projection of an inscribed or circumscribed cube onto the sphere, and a square in the complex plane.
In the end, this is the same map as computed in [Rancic (1996)](https://doi.org/10.1002/qj.49712253209), but the procedure there is complicated -- the singularities in the map are avoided by performing a stereographic projection through the cube's corners, changing variables to remove the principal singularity at the corner, and computing the remaining mapping terms via a Taylor series.

Here we compute the conformal map using [lightning rational approximations](https://doi.org/10.1007/s40315-020-00325-w) which explicitly include terms with poles clustered near the corners.
The map coefficients are computed simply by least squares.
This much simpler algorithm still gives near machine precision, allowing for independent confirmation of the Rancic coefficients and providing a point of comparison for [discussions of their accuracy](https://github.com/CliMA/CubedSphere.jl/issues/15).
