import sys
sys.path.append('../')

from scipy.integrate._ode import ode
from src.ode              import dirField, dirField_2, Euler
from pylab                import *
from mpl_toolkits.mplot3d import Axes3D
from scipy.io             import loadmat
from scipy.optimize       import fmin

#---------------------------------------------------------------------
# ODE function to be integrated :
def fungenlogistic(t, y, p):
  ydot = p[2]*p[1]*y * (p[1]*t)**(p[2]-1) * (1-y/p[0])
  return ydot

def bootOptimize(params, f, t, u, dt):
  """
  FUNCTION BOOTOPTIMIZE: This function calls an ODE solver to fit a 
    specific ODE function "odemodel" to a set of data given by (ydat,
    xdat).
  """
  t0   = t[0]
  tf   = t[-1]
  u0   = params[0]
  pars = parameters[1:]
  
  ta  = arange(t0, tf+dt, dt)        # Time span
  ys  = linspace(0, max(y0), 1000)   # ys for plotting dy/dt
  
  # Call function that integrates the ODE:
  r = ode(f)
  r.set_integrator('dopri5', atol=1e-8, rtol=1e-5)
  r.set_initial_value(u0, t0)
  r.set_f_params(pars)
  
  sol = []
  sol.append(u0)
  for t in ta[:-1]:
    r.integrate(r.t + dt)
    sol.append(r.y)
  sol = array(sol).T
  
  # Sum of Square Errors :
  SSE = sum((u - sol)**2)
  return SSE 


# data:
data = loadmat('data/blood.mat')   
t    = data['blood'][0][0][0].T[0]
u    = data['blood'][0][0][1].T[0]

# Initial conditions
u0 = [27, 221, 0.33, 1]
dt = 0.01


# ================================================================ #
# Fits the ODE to the data using the simplex method with objective #
# function "bootoptimize" (usually least squares), starting values #
# "startvals" and data in "y,x" to fit the data to the model       #
# ================================================================ #
ftn  = lambda p: bootOptimize(p, fungenlogistic, t, u, dt)

phat = fmin(func=ftn, x0=u0)

#Plot the solution:
fig = figure(figsize=(12,5))
ax1 = fig.add_subplot(121)
ax1.plot(ta, sol[0], '-', lw=2.0, label=r'$S$')
ax1.plot(ta, sol[1], '-', lw=2.0, label=r'$I$')
ax1.plot(ta, sol[2], '-', lw=2.0, label=r'$R$')

ax1.set_xlabel(r'$S,I$')
ax1.set_ylabel(r'$f(S),f(I)$')
ax1.set_title('Solution')
leg = ax1.legend(loc='upper right')
leg.get_frame().set_alpha(0.5)
ax1.grid()


# Plot the results
ax2 = fig.add_subplot(122)
xmin = sol[0].min()
xmax = sol[0].max()
ymin = sol[1].min()
ymax = sol[1].max()
zmin = sol[2].min()
zmax = sol[2].max()
ax2.set_xlim([xmin, xmax])
ax2.set_ylim([ymin, ymax])

dirField_2(dSdt, dIdt, ax2, S_params, I_params)
ax2.plot(sol[0], sol[1])

ax2.set_title(r'$(S,I)$ phase plane')
ax2.set_xlabel(r'$S$')
ax2.set_ylabel(r'$I$')
tight_layout()
savefig('prb4b.png', dpi=300)
show()




