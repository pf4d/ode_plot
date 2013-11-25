import sys
sys.path.append('../')

from scipy.integrate._ode import ode
from src.ode              import dirField, dirField_2, Euler
from pylab                import *
from mpl_toolkits.mplot3d import Axes3D
from scipy.io             import loadmat
from scipy.optimize       import fmin
from numpy.random         import choice

#---------------------------------------------------------------------
# ODE function to be integrated :
def f(t, y, p):
  theta1 = p[0]
  theta2 = p[1]
  ydot   = theta1*exp(-theta1*t) - theta2*y
  return ydot

def integrate(params, f, ta, u, dt):
  t0 = ta[0]
  u0 = params[0]
  p0 = params[1:]
    
  # Call function that integrates the ODE:
  r = ode(f)
  r.set_integrator('dopri5', atol=1e-8, rtol=1e-5)
  r.set_initial_value(u0, t0)
  r.set_f_params(p0)
  
  sol = []
  sol.append(u0)
  for t in ta[:-1]:
    r.integrate(r.t + dt)
    sol.append(r.y)
  sol = array(sol).T
  return sol

# function to be optimized :
def bootOptimize(params, f, ta, u):
  """
  FUNCTION BOOTOPTIMIZE: This function calls an ODE solver to fit a 
    specific ODE function "odemodel" to a set of data given by (ydat,
    xdat).
  """
  dt  = ta[1] - ta[0]
  sol = integrate(params, f, ta, u, dt)
  
  # Sum of Square Errors :
  SSE = sum((u - sol)**2)
  return SSE 

# non-linear model bootstrap :
def bootNLM(t, u, params, f, B, dtmin):
  """
  FUNCTION BOOTNLM: This function reads in the data (x,y), the parameter 
    starting values (generally the parameter estimates), and the number 
    of bootstrap samples desired (B).  It returns the bootstrap estimates 
    of the parameters in a B x p vector (bootest) and the Delta method 
    standard errors (sedelta).
  """
  n   = len(t)
  np  = len(p0)
  dt  = t[1] - t[0]
  u0  = params[0]
  p0  = params[1:]

  uhat = integrate(p0, f, ta, u, dt)
  
  Ax_min = array([])
  us     = array([])
  Fhat   = array([])
  for i in range(np):
    Ax_min      = append(Ax_min, params)
    Ax_min[i,i] = Ax_min[i,i] + dtmin
    p0_new      = Ax_min[i,:]
    uhat_new    = integrate(p0_new, f, t, u, dt)
    
    us   = append(us, uhat_new)
    Jhat = (uhat_new - uhat) / dtmin
    Fhat = append(Fhat, Jhat)

  Fhat = Fhat.T
  us   = us.T
  
  resid  = u - uhat
  MSE    = sum(resid**2) / (n-np)
  FTF    = dot(Fhat.T, Fhat)
  varMat = inv(FTF * MSE)
  sedelt = sqrt(diag(varMat))

  leverage = diag(dot(Fhat, dot(inv(FTF), Fhat.T)))
  uhatmat  = dot(ones(B,1), yhat)

  modres   = resid / sqrt(1 - leverage)
  modres  -= mean(modres)
  residmat = choice(modres, B*n, replace=True)
  residmat = reshape(residmat, (B,n))

  bsamp    = uhatmat + residmat
  
  for i in range(B):
    btdat = t
    budat = bsamp[i,:]

    ftn = lambda p: bootOptimize(p, f, btdat, budat, u0) 
    bootest[i,np] = fmin()

  



# data:
data = loadmat('data/blood.mat')   
t    = data['blood'][0][0][0].T[0]
u    = data['blood'][0][0][1].T[0]

# Initial conditions
p0 = [0, 1, 1]


# ================================================================ #
# Fits the ODE to the data using the simplex method with objective #
# function "bootoptimize" (usually least squares), starting values #
# "startvals" and data in "y,x" to fit the data to the model       #
# ================================================================ #
ftn  = lambda p: bootOptimize(p, f, t, u)

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




