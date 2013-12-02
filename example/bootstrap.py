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

  uhat = integrate(params, f, t, dt)
  
  Ax_min = tile(params, (np, 1))
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
  for i in range(B):
    btdat = t
    budat = bsamp[i,:]

    ftn = lambda p: bootOptimize(p, f, btdat, budat, dt) 
    bootnew = fmin(func=ftn, x0=params)
    if i == 0:
      bootest = bootnew
    else:
      bootest = vstack((bootest, bootnew))
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

# Fits the ODE to the data :
ftn  = lambda p: bootOptimize(p, f, t, u, dt)
phat = fmin(func=ftn, x0=p0)
uhat = integrate(phat, f, ta, pdt)

# perform bootstrap method 2 :
B = 5000
betaBoot, seDelta = bootNLM(t, u, phat, f, B, dt)

fig = figure()
plot(t,  u, 'ro')
plot(ta, uhat, 'k-', lw=2.0)
grid()
xlabel(r'$t$')
ylabel(r'concentration')
tight_layout()
show()

bm    = mean(betaBoot, axis=0)
bias  = bm - phat
bse   = std(betaBoot, axis=0)
p     = len(phat)
nbins = max(10, round(B/50.0))

fig = figure(figsize=(12,5))
tit = [r'$u_0$', r'$\theta_1$', r'$\theta_2$']
for i in range(p):
  ax = fig.add_subplot(131 + i)
  ax.hist(betaBoot[:,i], nbins)
  ax.set_xlabel('Bootstrap Estimates')
  ax.set_ylabel('Frequency')
  ax.set_title('Bootstrap ' + tit[i] + 'Estimates')
  ax.grid()
tight_layout()
show()

sbeta = betaBoot.copy()
for i in range(p):
  sbeta[:,i] = sort(sbeta[:,i])

alpha = 0.05
cilow = sbeta[round(alpha*(B-1)/2.0), :]
cihi  = sbeta[round(B-alpha*(B-1)/2.0 - 1), :]

out   = {'true'  : phat,
         'mean'  : bm,
         'se'    : bse,
         'bias'  : bias,
         'cilow' : cilow,
         'cihi'  : cihi}



