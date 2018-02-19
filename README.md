# plot_lightcone

This code receives an input catalogue with positions and velocities at a given redshift and it displays galaxies in a lightcone. You can choose to display galaxies in real and redshift space. Also, redshift uncertainties can be specified with a value of sigma_z, so that

![](https://latex.codecogs.com/gif.latex?%5Cdelta%20z%20%3D%20%5Csigma_z%281&plus;z%29)

The redshift of each galaxy z_gal = z + G(z), where z is given by the line-of-sight position, and G(z) is a Gaussian centered in zero with variance dz^2

Here are some examples of the output of the code:

![Example of output](https://github.com/aaorsi/plot_lightcone/blob/master/zspace.png)
![Example of output](https://github.com/aaorsi/plot_lightcone/blob/master/zspace_2.png)
