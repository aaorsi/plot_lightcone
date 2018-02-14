# Plots to show z-errors using galaxy catalogues directly


from astropy.cosmology import FlatLambdaCDM
from scipy.interpolate import interp1d
cosmo = FlatLambdaCDM(H0 = 73.0,Om0=0.25)

from matplotlib.transforms import Affine2D
import mpl_toolkits.axisartist.floating_axes as floating_axes
import numpy as np
import mpl_toolkits.axisartist.angle_helper as angle_helper
from matplotlib.projections import PolarAxes
from mpl_toolkits.axisartist.grid_finder import (FixedLocator, MaxNLocator,
                                                 DictFormatter)
import matplotlib.pyplot as plt
from matplotlib import pyplot, transforms


redshift = 2.0
galfile = 'gal_iz30_atlascat_vel.cat'
x, y, z, vy, loglum = np.loadtxt(galfile, unpack = True)

zc = z < 20

aexp = 1./(1 + redshift)
Hz = cosmo.H(redshift).value

sz = 1e-2
PlotName = 'sigmax' # zspace # sigmaX


y_ = y if PlotName == 'real' else y + vy / (aexp*Hz) # zspace

pname = PlotName

if PlotName == 'sigmax':
  pname = 'zspace_sigmaz%.4f_thick' % sz


ym = np.mean(y_[zc])
xm = np.mean(x[zc])

y1 = y_[zc] - ym
dc = cosmo.comoving_distance(redshift) * cosmo.h
rad = dc.value + y1  # Radial distance along the y-axis w/r to observer at z=0

zarr = np.linspace(1.0,3.0,num=1000)
dcarr = cosmo.comoving_distance(zarr)*cosmo.h

zmap = interp1d(dcarr,zarr)

zgal = zmap(rad)                   # redshift of galaxies


if PlotName == 'sigmax':
  for i in range(len(zgal)):
    zgal[i] += np.random.normal(0.0,sz * (1 + zgal[i]))

tgal = np.arcsin( (x[zc] - xm) / rad)  # angular position


tmax = tgal.max()*.75
tmin = tgal.min()*.75

zmax = 2.13
zmin = 1.78




def setup_axes2(fig, rect,tmin, tmax,zmin,zmax):
  """
  With custom locator and formatter.
  Note that the extreme values are swapped.
  """
    # rotate a bit for better orientation
  #tr_rotate = Affine2D().scale(2, 1).rotate_deg(90)
#  tr_rotate = Affine2D().translate(10, 1)
  
  # scale degree to radians
 # tr_scale = Affine2D().scale(np.pi/180., 1.)

  tr =PolarAxes.PolarTransform()


  pi = np.pi
#  angle_ticks = [(0, r"$0$"),
#                 (.25*pi, r"$\frac{1}{4}\pi$"),
#                 (.5*pi, r"$\frac{1}{2}\pi$")]

  angle_ticks = [(tmin, '%.2f' % tmin), (0,r'$0$'), (tmax, '%.2f' % tmax)]

  grid_locator1 = FixedLocator([v for v, s in angle_ticks])
  tick_formatter1 = DictFormatter(dict(angle_ticks))

  grid_locator2 = MaxNLocator(4)

  grid_helper = floating_axes.GridHelperCurveLinear(
      tr, extremes=(tmax, tmin, zmax, zmin),
      grid_locator1=grid_locator1,
      grid_locator2=grid_locator2,
      tick_formatter1=tick_formatter1,
      tick_formatter2=None)

  ax1 = floating_axes.FloatingSubplot(fig, rect, grid_helper=grid_helper)
  fig.add_subplot(ax1)

  # create a parasite axes whose transData in RA, cz
  aux_ax = ax1.get_aux_axes(tr)

  aux_ax.patch = ax1.patch  # for aux_ax to have a clip path as in ax
  ax1.patch.zorder = 0.95  # but this has a side effect that the patch is
  # drawn twice, and possibly over some other
  # artists. So, we decrease the zorder a bit to
  # prevent this.

  return ax1, aux_ax


def setup_axes3(fig, rect):
  """
  Sometimes, things like axis_direction need to be adjusted.
  """

  # rotate a bit for better orientation
  tr_rotate = Affine2D().translate(-95, 0)

  # scale degree to radians
  tr_scale = Affine2D().scale(np.pi/180., 1.)

  tr = tr_rotate + tr_scale + PolarAxes.PolarTransform()

  grid_locator1 = angle_helper.LocatorHMS(4)
  tick_formatter1 = angle_helper.FormatterHMS()

  grid_locator2 = MaxNLocator(3)

  ra0, ra1 = 8.*15, 14.*15
  cz0, cz1 = 0, 14000
  grid_helper = floating_axes.GridHelperCurveLinear(
      tr, extremes=(ra0, ra1, cz0, cz1),
      grid_locator1=grid_locator1,
      grid_locator2=grid_locator2,
      tick_formatter1=tick_formatter1,
      tick_formatter2=None)

  ax1 = floating_axes.FloatingSubplot(fig, rect, grid_helper=grid_helper)
  fig.add_subplot(ax1)

  # adjust axis
  ax1.axis["left"].set_axis_direction("bottom")
  ax1.axis["right"].set_axis_direction("top")

  ax1.axis["bottom"].set_visible(False)
  ax1.axis["top"].set_axis_direction("bottom")
  ax1.axis["top"].toggle(ticklabels=True, label=True)
  ax1.axis["top"].major_ticklabels.set_axis_direction("top")
  ax1.axis["top"].label.set_axis_direction("top")

  ax1.axis["left"].label.set_text(r"cz [km$^{-1}$]")
  ax1.axis["top"].label.set_text(r"$\alpha_{1950}$")

  # create a parasite axes whose transData in RA, cz
  aux_ax = ax1.get_aux_axes(tr)

  aux_ax.patch = ax1.patch  # for aux_ax to have a clip path as in ax
  ax1.patch.zorder = 0.9  # but this has a side effect that the patch is
  # drawn twice, and possibly over some other
  # artists. So, we decrease the zorder a bit to
  # prevent this.

  return ax1, aux_ax



fig = plt.figure(1,figsize=(15,15))

ax, aux_ax = setup_axes2(fig, 111,tmin,tmax,zmin,zmax)

#import ipdb ; ipdb.set_trace()

#base = pyplot.gca().transData
#rot = transforms.Affine2D().rotate_deg(90)

aux_ax.plot(tgal,zgal,'.',color='k',markersize=3)
#aux_ax.scatter(tgal,zgal)
plt.ylabel('redshift')



plt.savefig(pname+'.png',bbox_inches='tight',dpi=300)



