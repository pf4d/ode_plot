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
  F = -1*x1 + 0*x2 
  return array(F)

def f2(x1, x2, *fparams):
  """
  INPUT:
    x1 - position
    x2 - position
  OUTPUT:
    F - dx2/dx1
  """
  F = 0*x1 + -1*x2 
  return array(F)


def plot_solutions(cvec, x):
  '''
  plot the solutions to the matrix x for initial conditions c.
  '''
  for c1 in cvec:
    for c2 in cvec:
      x1 = c1*x[0,0] + c2*x[0,1]
      x2 = c1*x[1,0] + c2*x[1,1]
      plot(x1, x2, 'r', lw=1.5)


# times to evaluate a solution :
x1min = -1.0001
x1max = 1
x2min = -1
x2max = 1
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
X = array([[-1, 0],[0, -1]])
val, vec = eig(X)

x11 = vec[0,0]*exp(val[0]*t)
x12 = vec[0,1]*exp(val[1]*t)
x21 = vec[1,0]*exp(val[0]*t)
x22 = vec[1,1]*exp(val[1]*t)
x = array([[x11, x12], [x21, x22]])

# plot the solutions :
c = [-3,-1,1,3]  # initial values
plot_solutions(c, x)

legend(loc='lower right')
xlabel(r'$x_1$')
ylabel(r'$x_2$')
title(r"Phase plane")
#savefig('image.png', dpi=150)
show()


