import sys
sys.path.append('../')

from scipy.integrate._ode import ode
from src.ode import dirField, dirField_2, Euler
from pylab import *
from mpl_toolkits.mplot3d import Axes3D

#---------------------------------------------------------------------
# ODE function to be integrated
def dLdt(L, S, C, params):
  """
  INPUT:
  OUTPUT:
  """
  K_L   = params[0]
  F     = params[1]
  V     = params[2]
  dLdt  = K_L*C*L - F*L/V
  return array(dLdt)

def dSdt(L, S, C, params):
  """
  INPUT:
  OUTPUT:
  """
  K_S   = params[0]
  F     = params[1]
  V     = params[2]
  beta  = params[3]
  dSdt  = K_S*C*S - F*S/V - beta*S*L
  return array(dSdt)

def dCdt(L, S, C, params):
  """
  INPUT:
  OUTPUT:
  """
  a_L  = params[0]
  a_S  = params[1]
  K_L  = params[2]
  K_S  = params[3]
  F    = params[4]
  V    = params[5]
  C_0  = params[6]
  dCdt = -a_L*K_L*C*L - a_S*K_S*C*S - F*C/V + F*C_0/V
  return array(dCdt)

def f(t, y, dLdt, dSdt, dCdt, L_params, S_params, C_params):
  """
  INPUT: 
    t        - time array
    dLdt     - function 
    uSdt     - function
    dCdt     - function
    L_params - parameters for dLdt
    S_params - parameters for dSdt
    C_params - parameters for dCdt
  OUTPUT:
    ydot[0] = time derivative of y[0],
    ydot[1] = time derivative of y[1],
    ydot[2] = time derivative of y[2].
  """
  L = y[0]
  S = y[1]
  C = y[2]
  rhs1 = dLdt(L, S, C, L_params)    # right hand side 1st eqn
  rhs2 = dSdt(L, S, C, S_params)    # right hand side 2nd eqn
  rhs3 = dCdt(L, S, C, C_params)    # right hand side 3rd eqn
  return array([rhs1, rhs2, rhs3])
 

# Initial conditions
y0 = [0.3, 0.3, 0.001]  

def model(F, beta, tf, dt):
  """
  Run model for given volume flow rate <F> and toxin coef <beta> for total
  time <tf> in hours at timestep <dt>, also in hours.
  """ 
  # Additional parameters being passed to the ODE function
  a_L  = 0.1124
  a_S  = 0.0325
  K_L  = 19.0288
  K_S  = 20.1818
  V    = 1.0
  C_0  = 0.02
  
  L_params = [K_L, F, V]
  S_params = [K_S, F, V, beta]
  C_params = [a_L, a_S, K_L, K_S, F, V, C_0]
 
  t0  = 0.0                          # initial time 
  ta  = arange(t0, tf+dt, dt)        # Time span
  ys  = linspace(0, max(y0), 1000)   # ys for plotting dy/dt
  
  # Call function that integrates the ODE:
  r = ode(f)
  r.set_integrator('dopri5', atol=1e-6, rtol=1e-5)
  r.set_initial_value(y0, t0)
  r.set_f_params(dLdt, dSdt, dCdt, L_params, S_params, C_params)
  
  sol = []
  sol.append(y0)
  for t in ta[:-1]:
    r.integrate(r.t + dt)
    sol.append(r.y)
  sol = array(sol).T
  
  return ta, sol

ta, sol = model(F=0.02, beta=0.009134, tf=40000, dt=2)

#Plot the solution:
fig = figure()
ax1 = fig.add_subplot(111)
ax1.plot(ta, sol[0], 'r-',  lw=2.0, label=r'$L$')
ax1.plot(ta, sol[1], 'k-',  lw=2.0, label=r'$S$')
ax1.plot(ta, sol[2], 'k--', lw=2.0, label=r'$C$')

ax1.set_xlabel(r'$t$')
ax1.set_ylabel(r'Optical Density')
ax1.set_title('Solution')
leg = ax1.legend(loc='upper right')
leg.get_frame().set_alpha(0.5)
ax1.grid()
show()


