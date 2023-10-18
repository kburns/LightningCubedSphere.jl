# Lightning conformal cubed sphere

This code solves for the conformal map between a face of the cubed sphere, the gnomic projection of an inscribed or circumscribed cube onto the sphere, and a square in the complex plane.
In the end, this is the same map as computed by [Rančić et al. (1996)](https://doi.org/10.1002/qj.49712253209), but the procedure in that work is a bit complicated -- the singularities in the map are avoided by performing a stereographic projection through the cube's corners, changing variables to remove the principal singularity at the corner, and computing the remaining mapping terms via a Taylor series.

Here, we directly compute the conformal map with a face-centered stereographic projection using [lightning rational approximations](https://doi.org/10.1007/s40315-020-00325-w), which explicitly include terms with poles clustered near the corners.
The map coefficients are computed simply by least squares.
This much simpler algorithm still gives near machine precision, allowing for independent confirmation of the Rančić coefficients and providing a point of comparison for [discussions of their accuracy](https://github.com/CliMA/CubedSphere.jl/issues/15).


## References
- Rančić, M., Purser, R.J. and Mesinger, F. (1996), [A global shallow-water model using an expanded spherical cube: Gnomonic versus conformal coordinates.](https://doi.org/10.1002/qj.49712253209) _Q. J. R. Meteorol. Soc._, **122**, 959-982.
- Trefethen, L.N. (2020), [Numerical conformal mapping with rational functions.](https://doi.org/10.1007/s40315-020-00325-w) _Comput. Methods Funct. Theory_, **20**, 369–387.
