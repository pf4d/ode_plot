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

def integrate(params, f, ta, dt):
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
def bootOptimize(params, f, t, u, dt):
  """
  FUNCTION BOOTOPTIMIZE: This function calls an ODE solver to fit a 
    specific ODE function "odemodel" to a set of data given by (ydat,
    xdat).
  """
  t0  = t[0]
  tf  = t[-1]
  ta  = arange(t0, tf+dt, dt)
  sol = integrate(params, f, ta, dt)
  
  t_indx = []
  for ti in t:
    t_indx.append(where(ta == ti)[0][0])

  # Sum of Square Errors :
  SSE = sum((u - sol[t_indx])**2)
  return SSE 

# non-linear model bootstrap :
def bootNLM(t, u, params, f, B, dt):
  """
  FUNCTION BOOTNLM: This function reads in the data (x,y), the parameter 
    starting values (generally the parameter estimates), and the number 
    of bootstrap samples desired (B).  It returns the bootstrap estimates 
    of the parameters in a B x p vector (bootest) and the Delta method 
    standard errors (sedelta).
  """
  n   = len(t)
  np  = len(params)
  u0  = params[0]
  p0  = params[1:]
  print u0, p0

  uhat = integrate(params, f, t, dt)
  
  Ax_min = tile(params, (np, 1))
  us     = array([])
  Fhat   = array([])
  for i in range(np):
    Ax_min[i,i] = Ax_min[i,i] + dt
    p0_new      = Ax_min[i,:]
    uhat_new    = integrate(p0_new, f, t, dt)
    Jhat        = (uhat_new - uhat) / dt
    if i == 0:
      us   = uhat_new
      Fhat = Jhat
    else:
      us   = vstack((us, uhat_new))
      Fhat = vstack((Fhat, Jhat))

  Fhat = Fhat.T
  us   = us.T
  
  resid  = u - uhat
  MSE    = sum(resid**2) / (n-np)
  FTF    = dot(Fhat.T, Fhat)
  varMat = inv(FTF * MSE)
  sedelt = sqrt(diag(varMat))

  leverage = diag(dot(Fhat, dot(inv(FTF), Fhat.T)))
  uhatmat  = tensordot(ones(B), uhat, 0)

  modres   = resid / sqrt(1 - leverage)
  modres  -= mean(modres)
  residmat = choice(modres, B*n, replace=True)
  residmat = reshape(residmat, (B,n))

  bsamp    = uhatmat + residmat
  bootest  = array([]) 
  for i in range(B):
    btdat = t
    budat = bsamp[i,:]

    ftn = lambda p: bootOptimize(p, f, btdat, budat, dt) 
    bootnew = fmin(func=ftn, x0=params)
    bootest = append(bootest, bootnew)
    print 'Bootstrap iteration %i' % i
  
  return bootest, sedelt
  



# data:
data = loadmat('data/blood.mat')   
t    = data['blood'][0][0][0].T[0]
u    = data['blood'][0][0][1].T[0]

# initial conditions :
p0 = [0, 1, 1]

# time parameters :
t0  = t[0]
tf  = t[-1]
dt  = 0.25
pdt = 0.001
ta  = arange(t0, tf, pdt)
B   = 2


# ================================================================ #
# Fits the ODE to the data using the simplex method with objective #
# function "bootoptimize" (usually least squares), starting values #
# "startvals" and data in "y,x" to fit the data to the model       #
# ================================================================ #
ftn  = lambda p: bootOptimize(p, f, t, u, dt)
phat = fmin(func=ftn, x0=p0)
uhat = integrate(phat, f, ta, pdt)

betaBoot, seDelta = bootNLM(t, u, phat, f, B, dt)

fig = figure()
plot(t,  u, 'ro')
plot(ta, uhat, 'k-', lw=2.0)
grid()
xlabel(r'$t$')
ylabel(r'concentration')
show()

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




