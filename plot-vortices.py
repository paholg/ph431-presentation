#!/usr/bin/python
from __future__ import division
import matplotlib, sys
if 'hide' in sys.argv:
  matplotlib.use('Agg')
from pylab import *
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.tri as mtri

import optical_vortex as vort
# X = 0 # coordinate indices
# Y = 1
# Z = 2

# delta = 1e-5 # finite delta

# # physical constants
# c = 1
# epsilon = 1
# mu = 1/(epsilon*c*c)

# # beam properties
# omega = 100 # angular frequency
# k = omega/c # wave vector
# l = 3 # topological charge
# p = 1 # radial index
# A_0 = 1 # vector potential seed value
# w_0 = 1 # minimum beam waist
# z_R = k*w_0**2/2 # Rayleigh length

# # actual beam waist
# def w(z):
#   return w_0*sqrt(1+(z/z_R)**2)

# i = complex(0,1)

# def factorial(n):
#   return 1 if n == 0 else n * factorial(n-1)

# # Associated Laguerre polynomial
# def L(l,p,x):
#   m = 0
#   value = 0
#   while m <= p:
#     value += x**m * (-1)**m * factorial(p+l) \
#               / ( factorial(p-m)*factorial(l+m)*factorial(m) )
#     m += 1
#   return value

# # complex vector potential
# # assumes that beam of light travels in the z direction
# # vector potential is polarized in x direction
# def A_x(l,p,x,y,z,t):
#   r = sqrt(x**2 + y**2)
#   phi = arctan(y/x)# if x != 0 else 0
#   wz = w(z)
#   return A_0 * w_0/wz * (r*sqrt(2)/wz)**l * L(l,p,2*(r/wz)**2) \
#     * exp(-(r/wz)**2) * exp(-i*k*r**2*z/(2*(z**2+z_R**2))) \
#     * exp(-i*l*phi) * exp(i*(2*p+l+1)*arctan(z/z_R)) * exp(-i*omega*t)

# def A_I(x,y,z,t): # imaginary component of vector potential
#   A = imag(A_x(l,p,x,y,z,t))
#   A[x<0] *= -1
#   return A

# def E_x(x,y,z,t): # electric field, in x direction
#   return -omega * A_I(x,y,z,t)

# def B_y(x,y,z,t): # magnetic field, in y direction
#   return -k * A_I(x,y,z,t)

# def u(x,y,z,t): # energy density
#   return epsilon * omega**2 * A_I(x,y,z,t)**2

# def I(x,y,z,t): # intensity
#   return c / 2 * u(x,y,z,t)

# def S_z(x,y,z,t): # poynting vector, in z direction
#   return c * u(x,y,z,t)

# def sigma(i,j,x,y,z,t): # Maxwell stress tensor
#   return u(x,y,z,t) if i == j else 0

# def f(x,y,z,t): # force per unit volume
#   f = zeros(3)
#   u_0 = u(x,y,z,t)
#   f[X] = (u(x+delta,y,z,t)-u_0)/delta
#   f[Y] = (u(x,y+delta,z,t)-u_0)/delta
#   f[Z] = (u(x,y,z+delta,t)-u_0)/delta - (u(x,y,z,t+delta)-u_0)/delta/c
#   return f,sqrt(f[X]**2 + f[Y]**2 + f[Z]**2)

diff = .001
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
zmax = 10
dx = .05
x = arange(-xmax, xmax+dx/2, dx)
y = arange(-xmax, xmax+dx/2, dx)
x, y = meshgrid(x, y)

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



if 'anim' in sys.argv:
  fig = figure()
  ax = fig.add_subplot(111)
  z = 0
  t = 0
  im = pcolormesh(x, y, vort.E(x, y, z, t)[0], cmap=cmap)

  def update(x, y, z, t):
    return vort.S(x, y, z, t)[1]

  F = update(x, y, z, t)
  vmax = max(F.flat)

  def update_contour_plot(t):
    clf()
    F = update(x, y, z, t)
    # vmax = max(F.flat)
    im = pcolormesh(x, y, F, cmap=cmap, vmax=vmax, vmin=-vmax)
    #colorbar(im)
    ax.set_aspect('equal')
    return im,


  ani = animation.FuncAnimation(
    fig, update_contour_plot, frames=xrange(10000), fargs=(t), interval=1)

  show()

if 'slices' in sys.argv:
  zmin = 0
  zmax = 9
  zstep = 3

  fig = figure()
  ax = fig.add_subplot(111, projection='3d')
  x = arange(-xmax, xmax, .1)
  y = arange(-xmax, xmax, .1)
  x, y = meshgrid(x, y)

  for z in xrange(zmin, zmax, zstep):
    E = E_x(x, y, z, t)
    E[abs(E) < .1] = nan
    ax.contourf(x, y, E, cmap=cmap, offset = z)
  show()

if '3d' in sys.argv:
  frame = 0
  t = 0
  dt = 0.001
  import mayavi.mlab as mlab
  from tvtk.api import tvtk

  figure = mlab.figure(bgcolor=(1,1,1), fgcolor=(0,0,0), size=(1024,768))
  mlab.view(azimuth=-153, elevation=-97, roll=-95, distance=61, focalpoint=(0, 0, 0))

  xmax = 10
  dx = 1
  zmax = 30
  x,y,z = mgrid[-xmax:xmax:dx, -xmax:xmax:dx, 0:zmax:dx]
  E = vort.E(x, y, z, t)
  B = vort.B(x, y, z, t)
  S = vort.S(x, y, z, t)
  u = vort.u(x, y, z, t)

  size = E[0].shape[0]*E[0].shape[1]*E[0].shape[2]

  # qe = mlab.quiver3d(x, y, z, *E)
  # se = qe.mlab_source

  # qb = mlab.quiver3d(x, y, z, *B)
  # sb = qb.mlab_source

  qs = mlab.quiver3d(x, y, z, *S)
  ss = qs.mlab_source


  @mlab.animate(delay=100, ui=False)
  def anim():
    global t, frame, E, B, S
    while True:
      print mlab.view(), mlab.roll()
      print frame, t
      t += dt
      E = vort.E(x, y, z, t)
      B = vort.B(x, y, z, t)
      S = vort.S(x, y, z, t)

      # se.set(u=E[0].reshape(size), v=E[1].reshape(size), w=E[2].reshape(size))
      # sb.set(u=B[0].reshape(size), v=B[1].reshape(size), w=B[2].reshape(size))
      ss.set(u=S[0].reshape(size), v=S[1].reshape(size), w=S[2].reshape(size))
      yield
      # figure.scene.save("anim/test-%02i.png" %frame)
      frame += 1


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
