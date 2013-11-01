from pylab import *
from scipy.integrate._ode import IntegratorBase
from numpy import array, isfinite

mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['legend.fontsize'] = 'medium'

class Euler(IntegratorBase):
    runner = True

    def __init__(self,dt=.01):
      self.dt = dt

    def reset(self,n,has_jac):
      pass

    def run(self,f,jac,y0,t0,t1,f_params,jac_params):
      # this method is called to integrate from t=t0 to t=t1
      # with initial condition y0. f and jac are user-supplied functions
      # that define the problem. f_params,jac_params are additional
      # arguments to these functions.

      yo = array(y0) # Initial condition
      t = t0
      dt = self.dt

      while t < t1:
   
        # For last value of dt:
        if t + dt > t1:
          dt = t1 - t
        
        yn = yo + f(t,yo,*f_params) * dt
        yo = yn.copy()
        t += dt

      if isfinite(yn[-1]): self.success = True # Check for success
      return yn,t

if Euler.runner:
    IntegratorBase.integrator_classes.append(Euler)


def dirField_2(f1, f2, f1_params, f2_params, ax):
  """
  Plot the direction field of dx2/dx1 = f2/f1 on an axes object ax.
  """
  xmin = ax.get_xlim()[0]
  xmax = ax.get_xlim()[1]
  ymin = ax.get_ylim()[0]
  ymax = ax.get_ylim()[1]

  xstep  = (xmax - xmin)/30.
  ystep  = (ymax - ymin)/30.
  
  # coordinate initialization :
  x1     = arange(xmin, xmax, xstep)
  x2     = arange(xmin, xmax, xstep/4)    # higher-res. plot x coord.
  y1     = arange(ymin, ymax, ystep)
  y2     = arange(ymin, ymax, ystep/4)
  coord  = meshgrid(x1, y1)
  x      = coord[0]
  y      = coord[1]
  r      = 0.5                            # length of arrows - irrelevant
  
  #=============================================================================
  # differential function :
  f1  = f1(x, y, f1_params)
  f2  = f2(x, y, f2_params)
  f   = f2/f1

  v  = sqrt((r**2) / (1 + 1/f**2))        # length of arrow in y-dir
  u  = v/f                                # length of arrow in x-dir
  
  #=============================================================================
  # plotting :
  gray = '#5f5f5f'                           # color for axes
  ax.axhline(lw=1, c=gray)                   # x-axis
  ax.axvline(lw=1, c=gray)                   # y-axis
  ax.quiver(x, y, u, v, pivot='middle',
                        scale=None,
                        angles='xy',
                        headwidth=0.0, 
                        headlength=0.0, 
                        headaxislength=0.0)  # plot the dir. field
  ylim([ymin, ymax])                      # plotting y-axis limits
  xlim([xmin, xmax])                      # plotting x-axis limits


def dirField(f, ax):
  """
  Plot the direction field of f on an axes object ax.
  """
  xmin = ax.get_xlim()[0]
  xmax = ax.get_xlim()[1]
  ymin = ax.get_ylim()[0]
  ymax = ax.get_ylim()[1]

  xstep  = (xmax - xmin)/30.
  ystep  = (ymax - ymin)/30.
  
  # coordinate initialization :
  x1     = arange(xmin, xmax, xstep)
  x2     = arange(xmin, xmax, xstep/4)    # higher-res. plot x coord.
  y1     = arange(ymin, ymax, ystep)
  y2     = arange(ymin, ymax, ystep/4)
  coord  = meshgrid(x1, y1)
  x      = coord[0]
  y      = coord[1]
  r      = 0.5                            # length of arrows
  
  #=============================================================================
  # differential function :
  f  = f(x, y)
  
  v  = sqrt((r**2) / (1 + 1/f**2))        # length of arrow in y-dir
  u  = v/f                                # length of arrow in x-dir
  
  #=============================================================================
  # plotting :
  gray = '#5f5f5f'                           # color for axes
  ax.axhline(lw=1, c=gray)                   # x-axis
  ax.axvline(lw=1, c=gray)                   # y-axis
  ax.quiver(x, y, u, v, pivot='middle',
                        scale=None,
                        angles='xy',
                        headwidth=0.0, 
                        headlength=0.0, 
                        headaxislength=0.0)  # plot the dir. field
  ylim([ymin, ymax])                      # plotting y-axis limits
  xlim([xmin, xmax])                      # plotting x-axis limits

