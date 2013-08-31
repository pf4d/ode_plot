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
  F = 1*x1 + -1*x2 
  return array(F)

def f2(x1, x2, *fparams):
  """
  INPUT:
    x1 - position
    x2 - position
  OUTPUT:
    F - dx2/dx1
  """
  F = 1*x1 + 3*x2 
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
x1min = -3.00000001
x1max = 3
x2min = -3
x2max = 3
dx    = 0.1

x1 = arange(x1min, x1max, dx)
x2 = arange(x2min, x2max, dx)

# time range to plot the solution curves :
t = arange(-10, 10, 0.01)

# set up the plot window :
fig = figure()
ax = fig.add_subplot(111) 
ax.set_ylim(x2min, x2max)
ax.set_xlim(x1min, x1max)

# plot the direction field for the problem :
dirField_2(f1, f2, ax)

# computed eigenvalues with numpy :
X = array([[1, -1],[1, 3]])
val, vec = eig(X)
r   = val[0]
xi  = array([1,-1])    # eigenvector
eta = array([0,-1])    # generalized eigenvector TODO make determination 
                       #                              of this automated

x11 = xi[0]*exp(r*t)
x12 = xi[0]*t*exp(r*t) + eta[0]*exp(r*t)
x21 = xi[1]*exp(r*t)
x22 = xi[1]*t*exp(r*t) + eta[1]*exp(r*t)
x = array([[x11, x12], [x21, x22]])

# plot the solutions :
c = [-3, -2, -1, -.5, -.1, .1, .5, 1, 2, 3]  # initial values
#c = [5]  # initial values
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


