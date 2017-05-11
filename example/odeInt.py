from scipy.integrate._ode import *
from ode_plot import dirField, Euler
from pylab import *

def func(t, y, *fparams):
  """
  INPUT:
    y - position
    t - time
  OUTPUT:
    F - dy/dt
  """
  F = y*(3 - t*y) 
  return array(F)


# Initial Conditions
y0  = 0.5
t0  = 0.0
dt  =  0.5
idt = float(sys.argv[1])
tf  = 3.0
tf += dt

# Position arrays
yf  = [y0]

# Times to evaluate a solution. 
t   = arange(t0,tf,dt)

# CREATE ODE OBJECTS
i   = ode(func)

# Main loops for the integration
# Euler Method:
i.set_integrator('Euler', dt=idt)
i.set_initial_value(y0,t0)
for time in t[:-1]:
  i.integrate(i.t+dt)
  yf.append(i.y[0])
yf = array(yf)

# analytic solution:
#ya = 3/5.*(sin(t) + 2*cos(t) - 2*exp(-2*t))

# Plot the results
fig = figure()
ax = fig.add_subplot(111) 

plot(t, yf, label='Euler')
#plot(t, ya, 'r', label='Analytic')

# plot the direction field for the problem
dirField(func, ax)

legend(loc='lower right')
xlabel(r'$t$')
ylabel(r'$y$')
title(r"Euler's Method")
grid()
show()


