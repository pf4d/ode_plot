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
def f(t, y, p):
  theta1 = p[0]
  theta2 = p[1]
  ydot   = theta1*exp(-theta1*t) - theta2*y
  return ydot

def integrate(params, f, ta, u, u0, dt):
  t0 = ta[0]
  
  # Call function that integrates the ODE:
  r = ode(f)
  r.set_integrator('dopri5', atol=1e-8, rtol=1e-5)
  r.set_initial_value(u0, t0)
  r.set_f_params(params)
  
  sol = []
  sol.append(u0)
  for t in ta[:-1]:
    r.integrate(r.t + dt)
    sol.append(r.y)
  sol = array(sol).T
  return sol

# function to be optimized :
def bootOptimize(params, f, ta, u, u0):
  """
  FUNCTION BOOTOPTIMIZE: This function calls an ODE solver to fit a 
    specific ODE function "odemodel" to a set of data given by (ydat,
    xdat).
  """
  dt  = ta[1] - ta[0]
  sol = integrate(params, f, ta, u, u0, dt)
  
  # Sum of Square Errors :
  SSE = sum((u - sol)**2)
  return SSE 

# non-linear model bootstrap :
def bootNLM(t, u, u0, p0, f, B, dt):
  n   = len(t)
  np  = len(p0)
  t0  = t[0]                         # Starting time for simulation
  tf  = t[-1]                        # Time to end the simulation
  ta  = arange(t0, tf+dt, dt)        # Time span
  
  yhat = integrate(p0, f, ta, u, u0, dt)
  
  Ax_min = array([])
  ys     = array([])
  Fhat   = array([])
  for i in range(np):
    Ax_min      = append(Ax_min, p0)
    Ax_min[i,i] = Ax_min[i,i] + dt
    p0_new      = Axmin[i,:]
    yhat_new    = integrate(p0_new, f, ta, u, u0, dt)
    
    ys   = append(ys, yhat_new)
    Jhat = (yhat_new - yhat) / dt
    Fhat = append(Fhat, Jhat) 




# data:
data = loadmat('data/blood.mat')   
t    = data['blood'][0][0][0].T[0]
u    = data['blood'][0][0][1].T[0]

# Initial conditions
p0 = [1, 1]
u0 = 0.0


# ================================================================ #
# Fits the ODE to the data using the simplex method with objective #
# function "bootoptimize" (usually least squares), starting values #
# "startvals" and data in "y,x" to fit the data to the model       #
# ================================================================ #
ftn  = lambda p: bootOptimize(p, f, t, u, u0)

phat = fmin(func=ftn, x0=p0)

dt = 0.0001


##Plot the solution:
#fig = figure(figsize=(12,5))
#ax1 = fig.add_subplot(121)
#ax1.plot(ta, sol[0], '-', lw=2.0, label=r'$S$')
#ax1.plot(ta, sol[1], '-', lw=2.0, label=r'$I$')
#ax1.plot(ta, sol[2], '-', lw=2.0, label=r'$R$')
#
#ax1.set_xlabel(r'$S,I$')
#ax1.set_ylabel(r'$f(S),f(I)$')
#ax1.set_title('Solution')
#leg = ax1.legend(loc='upper right')
#leg.get_frame().set_alpha(0.5)
#ax1.grid()
#
#
## Plot the results
#ax2 = fig.add_subplot(122)
#xmin = sol[0].min()
#xmax = sol[0].max()
#ymin = sol[1].min()
#ymax = sol[1].max()
#zmin = sol[2].min()
#zmax = sol[2].max()
#ax2.set_xlim([xmin, xmax])
#ax2.set_ylim([ymin, ymax])
#
#dirField_2(dSdt, dIdt, ax2, S_params, I_params)
#ax2.plot(sol[0], sol[1])
#
#ax2.set_title(r'$(S,I)$ phase plane')
#ax2.set_xlabel(r'$S$')
#ax2.set_ylabel(r'$I$')
#tight_layout()
#savefig('prb4b.png', dpi=300)
#show()




