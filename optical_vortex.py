from __future__ import division
from pylab import *

X = 0 # coordinate indices
Y = 1
Z = 2
T = 4

XYZ = 3 # spacial dimensionality

delta = 1e-5 # finite delta

# physical constants
c = 1
epsilon = 1
mu = 1/(epsilon*c*c)

# beam properties

omega = 100 # angular frequency
k = omega/c # wave vector

l = 3 # topological charge
p = 1 # radial index
A_0 = 1 # vector potential seed value
w_0 = 1 # minimum beam waist
z_R = k*w_0**2/2 # Rayleigh length

# actual beam waist
def w(z):
  return w_0*sqrt(1+(z/z_R)**2)

i = complex(0,1)

def factorial(n):
  return 1 if n == 0 else n * factorial(n-1)

# Associated Laguerre polynomial
def L(l,p,x):
  m = 0
  value = 0
  while m <= p:
    value += x**m * (-1)**m * factorial(p+l) \
              / ( factorial(p-m)*factorial(l+m)*factorial(m) )
    m += 1
  return value

def sqNorm(F):
  return F[X]*F[X] + F[Y]*F[Y] + F[Z]*F[Z]

# complex vector potential
# assumes that beam of light travels in the z direction
# vector potential is polarized in x direction
def A(x,y,z,t):
  r = sqrt(x**2 + y**2)
  phi = arctan(y/x)
  wz = w(z)
  A_x = A_0 * w_0/wz * (r*sqrt(2)/wz)**l * L(l,p,2*(r/wz)**2) \
      * exp(-(r/wz)**2) * exp(-i*k*r**2*z/(2*(z**2+z_R**2))) \
      * exp(-i*l*phi) * exp(i*(2*p+l+1)*arctan(z/z_R)) * exp(-i*omega*t)
  A_x[x<0] *= -1
  return array([A_x, 0*A_x, 0*A_x])

def E(x,y,z,t):
  Ap = A(x,y,z,t)
  E_x = -omega*imag(Ap[0])
  E_y = -omega*imag(Ap[1])
  E_z = -omega*imag(Ap[2])
  return array([E_x, E_y, E_z])

def B(x,y,z,t): # magnetic field, in y direction
  Ap = A(x,y,z,t)
  A_d = [None]*XYZ
  A_d[X] = A(x+delta,y,z,t)
  A_d[Y] = A(x,y+delta,z,t)
  A_d[Z] = A(x,y,z+delta,t)

  BX = ( (A_d[Y][Z] - Ap[Z]) - (A_d[Z][Y] - Ap[Y]) ) / delta
  BY = ( (A_d[Z][X] - Ap[X]) - (A_d[X][Z] - Ap[Z]) ) / delta
  BZ = ( (A_d[X][Y] - Ap[Y]) - (A_d[Y][X] - Ap[X]) ) / delta

  return array([real(BX), real(BY), real(BZ)])

def S(x,y,z,t): # Poynting vector
  S = array(zeros(XYZ))
  Ep = E(x,y,z,t)
  Bp = B(x,y,z,t)
  Sx = (Ep[Y]*Bp[Z] - Ep[Z]*Bp[Y])/mu
  Sy = (Ep[Z]*Bp[X] - Ep[X]*Bp[Z])/mu
  Sz = (Ep[X]*Bp[Y] - Ep[Y]*Bp[X])/mu
  return array([Sx, Sy, Sz])

def u(x,y,z,t): # energy density
  return 1/2 * ( epsilon*sqNorm(E(x,y,z,t)) + sqNorm(B(x,y,z,t))/mu )

def I(x,y,z,t): # intensity
  return c / 2 * u(x,y,z,t)

def sigma(i,j,x,y,z,t): # Maxwell stress tensor
  return u(x,y,z,t) if i == j else 0

def f(x,y,z,t): # force per unit volume
  f = zeros(XYZ)
  u_0 = u(x,y,z,t)
  fx = (u(x+delta,y,z,t)-u_0)/delta
  fy = (u(x,y+delta,z,t)-u_0)/delta
  fz = (u(x,y,z+delta,t)-u_0)/delta - (u(x,y,z,t+delta)-u_0)/delta/c
  return array([fx, fy, fz])
