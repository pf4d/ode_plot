from scipy.integrate._ode    import ode
from scipy.io                import savemat
from pylab                   import *
from mpl_toolkits.mplot3d    import Axes3D
from multiprocessing         import Queue, cpu_count, Process
from mpl_toolkits.axes_grid1 import make_axes_locatable

#---------------------------------------------------------------------
# ODE function to be integrated
def dLdt(L, S, C, params):
  """
  INPUT:
    L - population of killer yeast.
    S - population of sensitive yeast.
    C - population of nutrient.
    params - dLdt equation parameters.
  OUTPUT:
    dLdt - time derivative of L.
  """
  K_L   = params[0]
  F     = params[1]
  V     = params[2]
  dLdt  = K_L*C*L - F*L/V
  return array(dLdt)

def dSdt(L, S, C, params):
  """
  INPUT:
    L - population of killer yeast.
    S - population of sensitive yeast.
    C - population of nutrient.
    params - dSdt equation parameters.
  OUTPUT:
    dSdt - time derivative of S.
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
    L - population of killer yeast.
    S - population of sensitive yeast.
    C - population of nutrient.
    params - dCdt equation parameters.
  OUTPUT:
    dCdt - time derivative of C.
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
 
def model(F, beta, y0, ta, dt):
  """
  Run model for given volume flow rate <F> and toxin coef <beta> for total
  time array <ta> in hours at timestep <dt>, also in hours.  Returns the 
  last solution 3-tuple for L, S, and C.
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
  
  # Call function that integrates the ODE:
  r = ode(f)
  r.set_integrator('dopri5', atol=1e-6, rtol=1e-5)
  r.set_initial_value(y0, ta)
  r.set_f_params(dLdt, dSdt, dCdt, L_params, S_params, C_params)
  
  sol = []
  sol.append(y0)
  for t in ta[:-1]:
    r.integrate(r.t + dt)
    sol.append(r.y)
  sol = array(sol).T
  
  return sol[:,-1]


class solveProcess(Process):
  """
  Process to solve the model function.
  """
  def __init__(self, i, queue, beta_a, F_a, y0, ta, dt, p):
    """
    Initialize the Process with ID <i>, processing queue <queue>, beta array
    <beta_a>, flow array <F_a>, time array <ta>, timestep <dt>, and number of
    parameters <p>.
    """
    Process.__init__(self)
    self.i      = i
    self.q      = queue
    self.beta_a = beta_a
    self.F_a    = F_a
    self.y0     = y0
    self.ta     = ta
    self.dt     = dt
    self.m      = len(beta_a)
    self.n      = len(F_a)
    self.p      = p
  
  def run(self):
    """
    solve the differential equations for all beta_a and F_a.
    """
    p = self.p
    m = self.m
    n = self.n

    SS_L = zeros((m,n))
    SS_S = zeros((m,n))
    SS_C = zeros((m,n))

    for i, beta in enumerate(self.beta_a):
      for j, F in enumerate(self.F_a):
        print 'Process %i solving: beta=%f, F=%f' % (self.i, beta, F)
        sol = model(F=F, beta=beta, y0=self.y0, ta=self.ta, dt=self.dt)
        SS_L[i,j] = sol[0]
        SS_S[i,j] = sol[1]
        SS_C[i,j] = sol[2]

    self.q.put(array([SS_L, SS_S, SS_C]))  # add the result to the queue.

def plot_sol(ax, f, extent, tit, cmap='Greys'):
  """
  plot the 2D solution <f> to axes <ax>.
  """
  im      = ax.imshow(f[::-1,:], extent=extent, cmap=cmap)
  divider = make_axes_locatable(ax)
  cax     = divider.append_axes("right", size="5%", pad=0.05)
  ax.set_title(tit)
  ax.set_ylabel(r'$\beta$')
  ax.set_xlabel(r'$F$')
  colorbar(im, cax=cax)

# parameters :
m  = 250     # number of beta discretizations.
n  = 250     # number of F discretizations.
p  = 3       # number of parameters
t0 = 0.0     # initial time
tf = 40000   # final time
dt = 500     # time step

# Initial conditions
y0 = [0.3, 0.3, 0.001]  

# range of beta, flow, and time to model :
betaMin = 0.0
betaMax = 1.0
Fmin    = 0.0
Fmax    = 0.5

beta_a  = linspace(betaMin, betaMax, m)
F_a     = linspace(Fmin,    Fmax,    n)
ta      = arange(t0, tf+dt, dt)

# multiprocessing data structures :
solvers = []
queue   = []
numCpus = cpu_count()
Fs      = array_split(F_a,    numCpus)

# create a solver for each processor and begin solving each :
for i in range(numCpus):
  q = Queue()
  queue.append(q)
  solver = solveProcess(i, q, beta_a, Fs[i], y0, ta, dt, p)
  solvers.append(solver)
  solver.start()

# wait until solver (started above) finishes :
for s in solvers:
  s.join()

# retrive the results :
sols = []
for q in queue:
  while q.empty() == False:
    sols.append(q.get())

# put the results from the individual cores back together :
for i, s in enumerate(sols):
  if i == 0:
    L_sol = s[0]
    S_sol = s[1]
    C_sol = s[2]
  else:
    L_sol = hstack((L_sol, s[0]))
    S_sol = hstack((S_sol, s[1]))
    C_sol = hstack((C_sol, s[2]))

Beta, Flow = meshgrid(beta_a, F_a)

data = {'Beta'  : Beta,
        'Flow'  : Flow,
        'L_sol' : L_sol,
        'S_sol' : S_sol,
        'C_sol' : C_sol}
savemat('../../killer_yeast/data/results.mat', data)

# plot the results :
fig = plt.figure()
ax  = fig.add_subplot(111, projection='3d')

ax.plot_wireframe(Beta, Flow, L_sol, color='r', lw=2.0, rstride=5, cstride=5)
ax.plot_wireframe(Beta, Flow, S_sol, color='k', lw=2.0, rstride=5, cstride=5)
ax.set_ylabel(r'$\beta$')
ax.set_xlabel(r'$F$')
show()


fig = plt.figure(figsize=(15,5))
ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)

extent = [betaMin, betaMax, Fmin, Fmax]
plot_sol(ax1, L_sol, extent, r'Killer',    cmap='Greys')
plot_sol(ax2, S_sol, extent, r'Sensitive', cmap='Greys')
plot_sol(ax3, C_sol, extent, r'Nutrient',  cmap='Greys')
tight_layout()
savefig('../../killer_yeast/doc/images/sols.png', dpi=300)
show()



