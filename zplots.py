# Plots to show z-errors using galaxy catalogues directly


from astropy.cosmology import FlatLambdaCDM
from scipy.interpolate import interp1d
cosmo = FlatLambdaCDM(H0 = 73.0,Om0=0.25)

import matplotlib
from matplotlib.transforms import Affine2D
import mpl_toolkits.axisartist.floating_axes as floating_axes
import numpy as np
import mpl_toolkits.axisartist.angle_helper as angle_helper
from matplotlib.projections import PolarAxes
from mpl_toolkits.axisartist.grid_finder import (FixedLocator, MaxNLocator,
                                                 DictFormatter)
import matplotlib.pyplot as plt
from matplotlib import pyplot, transforms


font = {'family' : 'STIXGeneral',
        'size'   : 22}

matplotlib.rc('font', **font)


redshift = 2.0
galfile = 'gal_iz30_atlascat_vel.cat'
x, y, z, vy, loglum = np.loadtxt(galfile, unpack = True)

zc = z < 30

aexp = 1./(1 + redshift)
Hz = cosmo.H(redshift).value

sz = 1e-4
PlotName = 'sigmax' # zspace # sigmax # real


y_ = y if PlotName == 'real' else y + vy / (aexp*Hz) # zspace

pname = PlotName

if PlotName == 'sigmax':
  pname = 'zspace_sigmaz%.4f_andrea_2' % sz


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


tmax = tgal.max()*.25
tmin = tgal.min()*.25

zmax = 2.15
zmin = 1.80


def setup_axes2(fig, rect,tmin, tmax,zmin,zmax):
  """
  With custom locator and formatter.
  Note that the extreme values are swapped.
  """

  tr =PolarAxes.PolarTransform()
  pi = np.pi

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

rmin = cosmo.comoving_distance(zmin)*cosmo.h
rmax = cosmo.comoving_distance(zmax)*cosmo.h
#ax.set_ylabel('redshift',fontsize=20)

#plt.title('hello')

#ax.text(0.5,-0.01,'hello',transform=ax.transAxes, rotation=+155)


ax.axis["left"].set_axis_direction("top")
ax.axis["left"].toggle(ticklabels=True, label=True)
ax.axis["left"].major_ticklabels.set_axis_direction("top")
ax.axis["left"].label.set_axis_direction("top")

ax.axis["right"].major_ticklabels.set_color("white")
ax.axis["right"].set_axis_direction("bottom")
ax.axis["right"].toggle(ticklabels=True, label=True)


ax.axis['right'].label.set_visible(True)
ax.axis['right'].label.set_text(r'${\rm Mpc/}h$')

ax.axis['left'].label.set_visible(True)
ax.axis['left'].label.set_text('redshift')
#ax2 = fig.add_subplot(111, sharex=ax, frameon=False)
#ax2.yaxis.tick_right()
#ax2.yaxis.set_label_position("right")
#ax2.set_ylabel(r'${\rm Mpc}/h$',fontsize=20)
#ax2.set_ylim([rmin.value, rmax.value])



zt_arr = np.linspace(zmin,zmax*.99,num=3)
rt_arr = cosmo.comoving_distance(zt_arr)*cosmo.h
#ax.axis.set_xticks(zt_arr, rt_arr)

for zz, rr in zip(zt_arr, rt_arr):
  aux_ax.text(tmin*1.25,zz,'%ld' % rr.value, rotation=-1.5)

print tmin
#ax.text(0.45, 0.005, '%ld' % rt_arr[1].value, transform=ax.transAxes,rotation=-3)
#ax.text(0.73, -0.003, '%ld' % rt_arr[2].value, transform=ax.transAxes,rotation=-3)




plt.savefig(pname+'.eps',bbox_inches='tight')



