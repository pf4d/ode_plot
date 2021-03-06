import sys
sys.path.append('../')

from scipy.integrate._ode import *
from src.ode import dirField_2, Euler
from pylab import *
import sys

def f1(x1, x2, *fparams):
  """
  INPUT:
    x1 - position
    x2 - position
  OUTPUT:
    F - dx2/dx1
  """
  F = 2*x1 - 1*x2 
  return array(F)

def f2(x1, x2, *fparams):
  """
  INPUT:
    x1 - position
    x2 - position
  OUTPUT:
    F - dx2/dx1
  """
  F = 3*x1 - 2*x2 
  return array(F)


# times to evaluate a solution :
x1min = -3
x1max = 3
x2min = -3
x2max = 3
dx    = 0.1

x1 = arange(x1min, x1max, dx)
x2 = arange(x2min, x2max, dx)

# set up the plot window :
fig = figure()
ax = fig.add_subplot(111) 
ax.set_ylim(x2min, x2max)
ax.set_xlim(x1min, x1max)

# plot equilibrium solutions :
plot(x2, x1,  'k--', lw=2, label=r'$x^1(t)$')
plot(x2, 3*x1, 'k',   lw=2, label=r'$x^2(t)$')

# plot the direction field for the problem :
dirField_2(f1, f2, ax)

# time range to plot the solution curves :
t = arange(-10, 10, 0.01)

# coefficient choices are made log to plot better :
c = log(arange(1.3, 10, 5))
c = append(-c[::-1], c)

# formation of solution matrix X :
x11 = exp(t)
x12 = exp(-t)
x21 = exp(t)
x22 = 3*exp(-t)
x = array([[x11, x12], [x21, x22]])

def plot_solutions(cvec, x):
  '''
  plot the solutions to the matrix x for initial conditions c.
  '''
  for c1 in cvec:
    for c2 in cvec:
      x1 = c1*x[0,0] + c2*x[0,1]
      x2 = c1*x[1,0] + c2*x[1,1]
      plot(x1, x2, 'r', lw=1.5)

plot_solutions(c, x)

## can we do this with a contour function of some Z?
#X1, X2 = np.meshgrid(x1,x2)
# Z = X1**2 + X2**2
#contour(X1, X2, Z)

legend(loc='lower right')
xlabel(r'$x_1$')
ylabel(r'$x_2$')
#title(r"Euler's Method")
savefig('image.png', dpi=150)
show()


