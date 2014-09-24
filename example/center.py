import sys
sys.path.append('../')

from scipy.integrate._ode import *
from src.ode import dirField_2, Euler
from pylab import *

def f1(x1, x2, *fparams):
  """
  INPUT:
    x1 - position
    x2 - position
  OUTPUT:
    F - dx2/dx1
  """
  F = 2*x1 - 5*x2 
  return array(F)

def f2(x1, x2, *fparams):
  """
  INPUT:
    x1 - position
    x2 - position
  OUTPUT:
    F - dx2/dx1
  """
  F = 1*x1 - 2*x2 
  return array(F)


def plot_solutions(cvec, x):
  '''
  plot the solutions to the matrix x for initial conditions c.
  '''
  for c in cvec:
    x1 = c*x[0,0] + c*x[0,1]
    x2 = c*x[1,0] + c*x[1,1]
    plot(x1, x2, 'r', lw=1.5)


# times to evaluate a solution :
x1min = -3
x1max = 3
x2min = -3
x2max = 3
dx    = 0.1

x1 = arange(x1min, x1max, dx)
x2 = arange(x2min, x2max, dx)

# time range to plot the solution curves :
t = arange(-200, 200, 0.01)

# set up the plot window :
fig = figure()
ax = fig.add_subplot(111) 
ax.set_ylim(x2min, x2max)
ax.set_xlim(x1min, x1max)

# plot the direction field for the problem :
dirField_2(f1, f2, ax)

# computed eigenvalues with numpy :
X = array([[2, -5],[1, -2]])
val, vec = eig(X)
a   = array([vec[0,0].real, vec[1,0].real])
b   = array([vec[0,0].imag, vec[1,0].imag])
lam = val[0].real
mu  = val[0].imag

u0 = exp(lam*t) * (a[0]*cos(mu*t) - b[0]*sin(mu*t))
u1 = exp(lam*t) * (a[1]*cos(mu*t) - b[1]*sin(mu*t))
v0 = exp(lam*t) * (a[0]*sin(mu*t) + b[0]*cos(mu*t))
v1 = exp(lam*t) * (a[1]*sin(mu*t) + b[1]*cos(mu*t))
x = array([[u0, v0], [u1, v1]])

# plot the solutions :
#c = [-3, -2, -1, -.5, -.1, .1, .5, 1, 2, 3]  # initial values
c = [1,2,3]  # initial values
plot_solutions(c, x)

legend(loc='lower right')
xlabel(r'$x_1$')
ylabel(r'$x_2$')
title(r"Phase plane")
#savefig('image.png', dpi=150)
show()


