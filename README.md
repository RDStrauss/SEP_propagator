# SEP_propagator

This is a 1D solar energetic particle (SEP) transport model. The focused transport equation for the coordinate along a magnetic field line is solved with a finite difference scheme.

For details see *van den Berg et al. [2020], Space Science Reviews* (https://link.springer.com/article/10.1007/s11214-020-00771-x) (ArXiv link: https://arxiv.org/abs/2012.07570) and the references *[Strauss & Fichnter, 2015; Strauss et al., 2018; Heita, 2019]*.

## Using SEP_propagator

There are a few parameters in the Fortran90 code which can be changed:
- `lambda`: The radial mean free path in AU.
- `energy`: The particle's kinetic energy in MeV.
- `species`: 1 for electrons or 2 for protons.  
- `r_position`: Radial position of observation point in AU.
- `totaltime`: The time the simulation must be computed for in hours.
- `V_sw`: The solar wind speed in km/s.
- `times`: An array with the times in hours when the pitch-angle distribution at the observation point should be printed out.
- `injection_switch`: 1 for a delta or 2 for a Reid-Axford injection in time. For the Reid-Axford profile, there is also an `acceleration_time` and `escape_time` in hours to control its shape.

To compile the Fortran90 code, use either the Intel Fortran compiler with
```
ifort -O2 SEP_propagator.f90
```
*or* the Gfortran compiler with
```
gfortran -O2 SEP_propagator.f90
```
and run the executable `a.out` or `./a.out`, depending on your system. Ater running the Fortran90 code, the Python code, which might not be Python3 compatable, does some basic plotting of the output and can be run with
```
python plot_output.py 
```

The output files are as follows:
- `grid_z.txt`: Values of the spatial grid.
- `grid_mu.txt`: Values of the pitch-angle grid.
- `output.txt`: Output at the observation point. The columns contain the value of the z-grid corresponding to the observation point in AU, the radial distance of the observation point in AU, the time in hours, the omni-directional intensity, and the anisotropy.
- `pitch_angle_distribution.txt`: The pitch-angle distribution at the observation point at different times. The first column contains the pitch-angle grid points, and then the columns alternate between the time and the pitch-angle distribution.
- `model_setup.txt`: The model setup for future reference. The parameters included, in order, are `lambda`, `energy`, `species`, `injection_switch`, `V_sw`, `acceleration_time`, `escape_time`, `N`, `M`, `Z_max`, `mulimiter`, and `zlimiter`. The last 5 parameters are technical parameters needed in the numerical scheme, *which should only be changed if you have experience with the numerical scheme*. `N` is the number of position grid points, `M` is the number of pitch-angle grid points, `Z_max` is the outer simulation boundary, and `mulimiter` and `zlimiter` are flux limiters.

## Examples

The examples include two figures presented in *van den Berg et al. [2020]* to show case how the model works. The first example (Figure 12) is a fitting of the model to an SEP event and also includes an example of how the code can be used to calculate the ratio of particles propagating away from and towards the Sun at the observation point. The second example (Figure 21) illustrates how the solution behaves if the mean free path and acceleration or injection time is changed.

## Disclosure and Notice

The  code  is published under the Creative Commons license, but is not intended to be used for commercial applications. We ask anyone using this model to reference *van den Berg et al. [2020]* in all research outputs and to contact the authors when used extensively.

## Changes to be implemented

-- Integral over pitch-angle dependence of Duu should be re-done with a finer u-grid!

## References

[Heita, P.K.N. 2019. Numerical investigation of solar energetic particle transport between the Sun, Earth, and Mars. *MSc thesis*. North-West University, South Africa.](https://dspace.nwu.ac.za/handle/10394/33865)

[Strauss, R.D. and Fichnter, H. 2015. On aspects pertaining to the perpendicular diffusion of solar energetic particles. *The Astrophysical Journal*, 801: 29.](https://ui.adsabs.harvard.edu/abs/2015ApJ...801...29S/abstract)

[Strauss, R.D., Ogunjobi, O., Moraal, H., McCracken, K.G., and Caballero-Lopez, R.A. 2018. On the pulse shape of ground-level enhancements. *Solar Physics*, 292(4): 51.](https://ui.adsabs.harvard.edu/abs/2017SoPh..292...51S/abstract)

van den Berg, J.P., Strauss, R.D., and Effenberger, F. 2020. A primer on focused solar energetic particle transport. *Space Science Reviews*, https://link.springer.com/article/10.1007/s11214-020-00771-x
