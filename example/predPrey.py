import sys
sys.path.append('../')

from scipy.integrate._ode import ode
from src.ode import dirField, dirField_2, Euler
from pylab import *

#---------------------------------------------------------------------
# ODE function to be integrated
def dudt(u, v, params):
  """
  INPUT:
    u - population of species u.
    v - population of species v.
  OUTPUT:
    du/dt - time derivative of population u as a function of u and v.
  """
  alpha = params[0]
  beta  = params[1]
  dudt = u - alpha*u**2 - beta*u*v
  return array(dudt)

def dvdt(u, v, params):
  """
  INPUT:
    u - population of species u.
    v - population of species v.
  OUTPUT:
    dv/dt - time derivative of population v as a function of u and v.
  """
  gamma = params[0]
  delta = params[1]
  dvdt = 2*v - gamma*v**2 - delta*u*v
  return array(dvdt)

def f(t, y, dudt, dvdt, alpha, beta, gamma, delta):
  """
  INPUT: 
    t    - time array
    dudt - function 
    uvdt - function
    alpha, beta, gamma, delta - constants 
  OUTPUT:
    ydot[0] = time derivative of y[0],
    ydot[1] = time derivative of y[1].
  """
  u     = y[0]
  v     = y[1]
  # For small problems the following syntax may be used:
  rhs1 = dudt(u,v, (alpha, beta))    # right hand side 1st eqn
  rhs2 = dvdt(u,v, (gamma, delta))   # right hand side 2nd eqn
  ydot = array([rhs1, rhs2])
  return ydot
   

# Initial conditions
y0 = [20.0, 10.0]
      
# Additional parameters being passed to the ODE function
alpha   = 0.01
beta    = 0.02
gamma   = 0.02
delta   = 0.04

t0  = 0                            # Starting time for simulation
tf  = 20                           # Time to end the simulation
dt  = 0.01                         # time step for integrator
ta  = arange(t0, tf+dt, dt)        # Time span
ys  = linspace(0, max(y0), 1000)   # ys for plotting dy/dt

# Call function that integrates the ODE:
r = ode(f)
r.set_integrator('dopri5', atol=1e-6, rtol=1e-3)
r.set_initial_value(y0, t0)
r.set_f_params(dudt, dvdt, alpha, beta, gamma, delta)

sol = []
sol.append(y0)
for t in ta[:-1]:
  r.integrate(r.t + dt)
  sol.append(r.y)
sol = array(sol).T

#Plot the solution:
fig = figure(figsize=(12,5))
ax1 = fig.add_subplot(121) 
ax1.plot(ta, sol[0], '-', lw=2.0, label=r'$u$')
ax1.plot(ta, sol[1], '-', lw=2.0, label=r'$v$')

ax1.set_xlabel(r'$u,v$')
ax1.set_ylabel(r'$f(u),f(v)$')
ax1.set_title('Solution')
leg = ax1.legend(loc='right center')
leg.get_frame().set_alpha(0.5)
ax1.grid()


# Plot the results
ax2 = fig.add_subplot(122)
xmin = sol[0].min() 
xmax = sol[0].max() 
ymin = sol[1].min() 
ymax = sol[1].max() 
ax2.set_xlim([xmin, xmax])
ax2.set_ylim([ymin, ymax])

# plot the direction field for the problem
dirField_2(dudt, dvdt, (alpha,beta), (gamma, delta), ax2)
ax2.plot(sol[0], sol[1])

ax2.set_title(r'$(u,v)$ phase plane')
ax2.set_xlabel(r'$u$')
ax2.set_ylabel(r'$v$')
tight_layout()
show()



