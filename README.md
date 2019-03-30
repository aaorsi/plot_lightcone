# plot_lightcone

This code receives an input catalogue with positions and velocities at a given redshift and it displays galaxies in a lightcone. You can choose to display galaxies in real and redshift space. Also, redshift uncertainties can be specified with a value of sigma_z, so that

![](https://latex.codecogs.com/gif.latex?%5Cdelta%20z%20%3D%20%5Csigma_z%281&plus;z%29)

The redshift of each galaxy z_gal = z + G(z), where z is given by the line-of-sight position, and G(z) is a Gaussian centered in zero with variance dz^2.

The code makes use of numpy, matplotlib, scipy and astropy

The following examples of usage of this code have been published in Figure 1 of Dickinson et al. 2019, [Observing Galaxy Evolution in the Context of Large-Scale Structure](https://arxiv.org/abs/1903.07409). The figures show the impact of adding progresively larger uncertainties in the redshift of galaxies, thus, eventually smearing out the spatial structure of the galaxy distribution.

![Example of output](https://github.com/aaorsi/plot_lightcone/blob/master/zspace_sigmaz0.0001_dickinson.png=250x)
![Example of output](https://github.com/aaorsi/plot_lightcone/blob/master/zspace_sigmaz0.0010_dickinson.png=250x)
![Example of output](https://github.com/aaorsi/plot_lightcone/blob/master/zspace_sigmaz0.0100_dickinson.png=250x)
