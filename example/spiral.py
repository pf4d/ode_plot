import sys
sys.path.append('../')

from scipy.integrate.ode import *
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
  F = 2*x1 - 5/2.*x2 
  return array(F)

def f2(x1, x2, *fparams):
  """
  INPUT:
    x1 - position
    x2 - position
  OUTPUT:
    F - dx2/dx1
  """
  F = 9/5.*x1 - 1*x2 
  return array(F)


# times to evaluate a solution :
x1min = -3
x1max = 3
x2min = -3
x2max = 3
dx    = 0.1

x1 = arange(x1min, x1max, dx)
x2 = arange(x2min, x2max, dx)

# time range to plot the solution curves :
t = arange(-20, 10, 0.01)

# set up the plot window :
fig = figure()
ax = fig.add_subplot(111) 
ax.set_ylim(x2min, x2max)
ax.set_xlim(x1min, x1max)

# plot the direction field for the problem :
dirField_2(f1, f2, ax)

# coefficient choices are made log to plot better :
c = [.5, 1, 2, 3]

# formation of solution matrix X :
x11 = exp(t/2)*5*cos(3/2.*t)
x12 = exp(t/2)*5*sin(3/2.*t)
x21 = exp(t/2)*3*(cos(3/2.*t) + sin(3/2.*t))
x22 = exp(t/2)*3*(-cos(3/2.*t) + sin(3/2.*t))
x = array([[x11, x12], [x21, x22]])

def plot_solutions(cvec, x):
  '''
  plot the solutions to the matrix x for initial conditions c.
  '''
  for c in cvec:
    x1 = c*x[0,0] + c*x[0,1]
    x2 = c*x[1,0] + c*x[1,1]
    plot(x1, x2, 'r', lw=1.5)

plot_solutions(c, x)

## can we do this with a contour function of some Z?
#X1, X2 = np.meshgrid(x1,x2)
# Z = X1**2 + X2**2
#contour(X1, X2, Z)

legend(loc='lower right')
xlabel(r'$x_1$')
ylabel(r'$x_2$')
title(r"Phase Plane")
savefig('image.png', dpi=150)
show()


