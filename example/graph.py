import sys
sys.path.append('../')

from scipy.integrate.ode import *
from src.ode import dirField_2, Euler
from pylab import *

def f1(x1, x2, *fparams):
  """
  INPUT:
    y - position
    t - time
  OUTPUT:
    F - dy/dt
  """
  F = -3*x1 + sqrt(2)*x2 
  return array(F)

def f2(x1, x2, *fparams):
  """
  INPUT:
    y - position
    t - time
  OUTPUT:
    F - dy/dt
  """
  F = sqrt(2)*x1 - 2*x2 
  return array(F)

## Initial Conditions
#y0  = 0.0
#t0  = 1.0
#idt = float(sys.argv[1])
#dt  = idt#0.1
#tf  = 1.8
#tf += dt
#
## Position arrays
#yf  = [y0]
#
## Times to evaluate a solution. 
#t   = arange(t0,tf,dt)
#
## CREATE ODE OBJECTS
#i   = ode(f1)
#
## Main loops for the integration
## Euler Method:
#i.set_integrator('Euler', dt=idt)
#i.set_initial_value(y0,t0)
#for time in t[:-1]:
#  i.integrate(i.t+dt)
#  yf.append(i.y[0])
#yf = array(yf)

# analytic solution:
#ya = 3/5.*(sin(t) + 2*cos(t) - 2*exp(-2*t))

# Plot the results
fig = figure()
ax = fig.add_subplot(111) 
ax.set_xlim(-3,3)
ax.set_ylim(-3,3)

#plot(t, yf, label='Euler')
#plot(t, ya, 'r', label='Analytic')

# plot the direction field for the problem
dirField_2(f1, f2, ax)

legend(loc='lower right')
xlabel(r'$x_1$')
ylabel(r'$x_2$')
#title(r"Euler's Method")
#savefig('image.png', dpi=150)
show()


