import sys
sys.path.append('../')

from scipy.integrate._ode import ode
from src.ode import dirField, dirField_2, Euler
from pylab import *

#---------------------------------------------------------------------
# ODE function to be integrated
def dSdt(S, I, R, params):
  """
  INPUT:
    S - population of susceptible
    I - population of infected.
    R - population of recovered.
  OUTPUT:
    dS/dt - time derivative of susceptible as a function of S and I.
  """
  alpha = params[0]
  gamma = params[1]
  dSdt  = -alpha*I*S + gamma*R
  return array(dSdt)

def dIdt(S, I, R, params):
  """
  INPUT:
    S - population of susceptible
    I - population of infected.
    R - population of recovered.
  OUTPUT:
    dI/dt - time derivative of infected as a function of S and I.
  """
  alpha = params[0]
  beta  = params[1]
  dIdt  = alpha*I*S - beta*I
  return array(dIdt)

def dRdt(S, I, R, params):
  """
  INPUT:
    S - population of susceptible
    I - population of infected.
    R - population of recovered.
  OUTPUT:
    dR/dt - time derivative of infected as a function of S and I.
  """
  beta  = params[0]
  gamma = params[1]
  dRdt  = beta*I - gamma*R
  return array(dRdt)

def f(t, y, dSdt, dIdt, dRdt, S_params, I_params, R_params):
  """
  INPUT: 
    t        - time array
    dSdt     - function 
    uIdt     - function
    dRdt     - function
    S_params - parameters for dSdt
    I_params - parameters for dIdt
    R_params - parameters for dRdt
  OUTPUT:
    ydot[0] = time derivative of y[0],
    ydot[1] = time derivative of y[1],
    ydot[2] = time derivative of y[2].
  """
  S     = y[0]
  I     = y[1]
  R     = y[2]
  # For small problems the following syntax may be used:
  rhs1 = dSdt(S, I, R, S_params)    # right hand side 1st eqn
  rhs2 = dIdt(S, I, R, I_params)    # right hand side 2nd eqn
  rhs3 = dRdt(S, I, R, R_params)    # right hand side 3rd eqn
  return array([rhs1, rhs2, rhs3])
   

# Initial conditions
y0 = [490.0, 10.0, 0.0]   # no recovered or immunized
#y0 = [290.0, 10.0, 200.0]  # 200 recovered or immunized
      
# Additional parameters being passed to the ODE function
N       = 500.0
alpha   = 0.001
gamma   = 0.01
beta    = 0.01

S_params = [alpha, gamma]
I_params = [alpha, beta]
R_params = [beta,  gamma]

t0  = 0                            # Starting time for simulation
tf  = 200                          # Time to end the simulation
dt  = 0.01                         # time step for integrator
ta  = arange(t0, tf+dt, dt)        # Time span
ys  = linspace(0, max(y0), 1000)   # ys for plotting dy/dt

# Call function that integrates the ODE:
r = ode(f)
r.set_integrator('dopri5', atol=1e-6, rtol=1e-3)
r.set_initial_value(y0, t0)
r.set_f_params(dSdt, dIdt, dRdt, S_params, I_params, R_params)

sol = []
sol.append(y0)
for t in ta[:-1]:
  r.integrate(r.t + dt)
  sol.append(r.y)
sol = array(sol).T

#Plot the solution:
fig = figure(figsize=(12,5))
ax1 = fig.add_subplot(121) 
ax1.plot(ta, sol[0], '-', lw=2.0, label=r'$S$')
ax1.plot(ta, sol[1], '-', lw=2.0, label=r'$I$')
ax1.plot(ta, sol[2], '-', lw=2.0, label=r'$R$')

ax1.set_xlabel(r'$S,I$')
ax1.set_ylabel(r'$f(S),f(I)$')
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
dirField_2(dSdt, dIdt, ax2, S_params, I_params)
ax2.plot(sol[0], sol[1])

ax2.set_title(r'$(S,I)$ phase plane')
ax2.set_xlabel(r'$S$')
ax2.set_ylabel(r'$I$')
tight_layout()
show()




