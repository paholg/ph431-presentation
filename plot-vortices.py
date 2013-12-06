#!/usr/bin/python
from __future__ import division
import matplotlib, sys
if not 'show' in sys.argv:
  matplotlib.use('Agg')
from pylab import *
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.tri as mtri

import optical_vortex as vort

nframes = 30

diff = .05
xlo = .5 - diff
xhi = .5 + diff

cdict = { 'red':  [(0.0,  0.0, 0.0),
                   (xlo,  1.0, 1.0),
                   (xhi,  1.0, 1.0),
                   (1.0,  1.0, 1.0)],

          'green': [(0.0, 0.0, 0.0),
                    (xlo, 1.0, 1.0),
                    (xhi, 1.0, 1.0),
                    (1.0, 0.0, 0.0)],

          'blue':   [(0.0,  1.0, 1.0),
                     (xlo,  1.0, 1.0),
                     (xhi,  1.0, 1.0),
                     (1.0,  0.0, 0.0)]}

cmap = matplotlib.colors.LinearSegmentedColormap('mine', cdict)

cdictpos = { 'red':  [(0.0,  1.0, 1.0),
                      (diff, 1.0, 1.0),
                      (1.0,  1.0, 1.0)],

             'green': [(0.0, 1.0, 1.0),
                       (diff, 1.0, 1.0),
                       (1.0, 0.0, 0.0)],

             'blue':   [(0.0,  1.0, 1.0),
                        (diff, 1.0, 1.0),
                        (1.0,  0.0, 0.0)]}

cmappos = matplotlib.colors.LinearSegmentedColormap('mine', cdictpos)


def plot_all(x, y, z, t, plane='z'):
  for F, tit in [(vort.E(x, y, z, t), 'E'), (vort.B(x, y, z, t), "B"),
                 (vort.S(x, y, z, t), 'S'), (vort.u(x, y, z, t), 'u'), (vort.f(x, y, z, t), 'f')]:
    for i in xrange(4):
      if i < 3:
        if 'u' in tit:
          continue
        f = F[i]
      else:
        if 'u' in tit:
          f = F
        else:
          f = sqrt(vort.sqNorm(F))
      vmax = max(f.flat)
      if vmax < dx:
        continue
      rcParams.update({'font.size':14, 'legend.fontsize':8})
      fig = figure()
      ax = fig.add_subplot(111)
      if plane == 'z':
        planeval = z
        ylab = '$y/w_0$'
        xlab = '$x/w_0$'
        xlim(-xmax, xmax)
        ylim(-xmax, xmax)
        yax = y
        xax = x
      elif plane == 'y':
        planeval = y
        ylab = '$z/w_0$'
        xlab = '$x/w_0$'
        xlim(-xmax, xmax)
        ylim(0, zmax)
        yax = z
        xax = x
      elif plane == 'x':
        planeval = x
        ylab = '$y/w_0$'
        xlab = '$z/w_0$'
        xlim(0, zmax)
        ylim(-xmax, xmax)
        yax = y
        xax = z

      xlabel(xlab)
      ylabel(ylab)

      top = vmax
      bot  = -vmax
      if i == 0:
        coord = 'x'
      elif i == 1:
        coord = 'y'
      elif i == 2:
        coord = 'z'
      else:
        if 'u' in tit:
          coord = ''
        else:
          coord = 'mag'
      if 'u' in tit:
        titstring = '$' + tit + '\quad' + plane + '= %g\quad t = %g$' %(planeval, t)
      else:
        titstring = '$' + tit + '_{%s} '%coord + '\quad' + plane + '= %g\quad t = %g$' %(planeval, t)
      title(titstring)
      print tit, coord, max(f.flat)
      pcolormesh(xax, yax, f, cmap=cmap, vmin = bot, vmax = top)
      ax.set_aspect('equal')
      #colorbar()
      fname =  tit+coord+'-'+plane+'-'+str(planeval)+'-t-'+str(t)+'.png'
      tight_layout()
      if 'save' in sys.argv:
        savefig('figs/'+fname)

xmax = 5
zmax = 45
dx = .05
x = arange(-xmax, xmax+dx/2, dx)
y = arange(-xmax, xmax+dx/2, dx)
x, y = meshgrid(x, y)
z = 0

A = vort.A(x, y, 0, 0)
E = vort.E(x, y, 0, 0)
vmax = max(E[0].flat)
vmin = -vmax

if 'paper' in sys.argv:
  x = arange(-xmax, xmax+dx/2, dx)
  y = arange(-xmax, xmax+dx/2, dx)
  x, y = meshgrid(x, y)
  t = 0
  for z in xrange(0,1):
    plot_all(x, y, z, t)


  # y = arange(-xmax, xmax+dx/2, dx)
  # z = arange(0, zmax+dx/2, dx)
  # z, y = meshgrid(z, y)
  # for x in [0]:
  #   plot_all(x, y, z, t, 'x')
  show()

if 'B' in sys.argv:
  name = 'B'
elif 'E' in sys.argv:
  name = 'E'
else:
  name = 'S'

if 'x' in sys.argv:
  name += 'x'
elif 'y' in sys.argv:
  name += 'y'
elif 'z' in sys.argv:
  name += 'z'

xz_axis = False
if 'Z' in sys.argv:
  xz_axis = True
  z = arange(0, zmax+dx/2, dx)
  x = arange(-xmax, xmax+dx/2, dx)
  z, x = meshgrid(z, x)
  y = .75


if 'slice' in sys.argv:
  def update(x, y, z, t):
    if 'B' in name:
      F = vort.B(x,y,z,t)
    elif 'E' in name:
      F = vort.E(x,y,z,t)
    else:
      F = vort.S(x,y,z,t)

    if 'x' in name:
      return F[0]
    elif 'y' in name:
      return F[1]
    return F[2]

  def update_contour_plot(t, i):
    #cla()
    size = (12, 6) if xz_axis else (4, 4)
    fig = figure(figsize=size)
    ax = fig.add_subplot(111)
    F = update(x, y, z, t)
    # vmax = max(F.flat)
    if xz_axis:
      xax = z
      yax = x
    else:
      xax = x
      yax = y
    im = contourf(xax, yax, F, 20, cmap=cmap, vmax=vmax, vmin=-vmax)

    rcParams.update({'font.size':16, 'legend.fontsize':8})
    tit = '$' + name[0] + '_' + name[1] + '$'
    if xz_axis:
      tit = '$' + name[0] + '_' + name[1] + '\quad y/w_0 = %0.2f$' %y
    title(tit)
    #ax.text(0, 5.6, tit, horizontalalignment='center', fontsize='30')
    axes().set_aspect('equal', 'datalim')
    tight_layout()

    if xz_axis:
      ax.set_xlabel('$z/w_0$')
      ax.set_ylabel('$x/w_0$')
      xlim(0, zmax)
      ylim(-5, 5)
      fname = 'anim/sliceZ-%s-%02i.png' %(name, i)
    else:
      ax.set_xlabel('$x/w_0$')
      ax.set_ylabel('$y/w_0$')
      ylim(-1, 1)
      fname = 'anim/slice-%s-%02i.png' %(name, i)
    print 'saving', fname
    savefig(fname)
    if not 'show' in sys.argv:
      close()

  t = 0


  F = update(x, y, z, t)
  vmax = max(F.flat)

  tmax = 2*pi/vort.omega
  i = 0

  update_contour_plot(0, 0)
  show()

  for t in linspace(0, tmax, nframes):
    update_contour_plot(t, i)

    i += 1


if '3d' in sys.argv:
  def field(x, y, z, t):
    if 'B' in name:
      return vort.B(x, y, z, t)
    elif 'E' in name:
      return vort.E(x, y, z, t)
    return vort.S(x, y, z, t)

  frame = 0
  t = 0

  dt = 2*pi/vort.omega/nframes
  import mayavi.mlab as mlab
  from tvtk.api import tvtk

  figure = mlab.figure(bgcolor=(1,1,1), fgcolor=(0,0,0), size=(500,768))

  xmax = 10
  dx = 1
  zmax = 25
  x,y,z = mgrid[-xmax:xmax:dx, -xmax:xmax:dx, 0:zmax:dx]

  F = field(x, y, z, t)

  size = F[0].shape[0]*F[0].shape[1]*F[0].shape[2]

  q = mlab.quiver3d(x, y, z, *F)
  s = q.mlab_source

  az = 15
  mlab.outline(extent=[-xmax/2, xmax/2, -xmax/2, xmax/2, 0, zmax])
  mlab.view(azimuth=az, elevation=70, roll=90, distance=75, focalpoint=(0, 0, zmax/2))

  spin_dir = 1
  @mlab.animate(delay=100, ui=False)
  def anim():
    global t, frame, F, az, spin_dir
    while True:
      print mlab.view(), mlab.roll()
      print frame, t
      t += dt
      az += spin_dir*.5

      mlab.view(azimuth=az, distance=60, focalpoint=(0, 0, zmax/2))

      F = field(x, y, z, t)

      s.set(u=F[0].reshape(size), v=F[1].reshape(size), w=F[2].reshape(size))
      yield
      fname = "anim/3d-%s-%02i.png" %(name, frame)
      print 'saving', fname
      figure.scene.save(fname)
      frame += 1
      if frame > nframes/2:
        spin_dir = -1
      if frame > nframes:
        exit(0)


  a = anim()

  mlab.show()

t = 7
if 'test' in sys.argv:
  import mayavi.mlab as mlab
  from tvtk.api import tvtk

  figure = mlab.figure()
  x,y,z = mgrid[-10:10, -10:10, -10:10]
  E = E_x(x, y, z, t)
  E = (E, 0*E, 0*E)

  B = B_y(x, y, z, t)
  B = (0*B, B, 0*B)

  S = S_z(x, y, z, t)
  S = (0*S, 0*S, S)
  print t

  #mlab.contour3d(x, y, z, E[0])
  s = E[0]
  #mlab.pipeline.volume(mlab.pipeline.scalar_field(E[0]), vmin=-.1, vmax=.1)


  mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(s),
                                              plane_orientation='x_axes',
                                              slice_index=10,
                                            )
  mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(s),
                                   plane_orientation='y_axes',
                                   slice_index=10,
                                 )
  mlab.outline()
  mlab.show()
