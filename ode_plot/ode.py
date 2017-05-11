#
#    Copyright (C) <2012>  <cummings.evan@gmail.com>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
from pylab                import *
from scipy.integrate._ode import IntegratorBase
from numpy                import array, isfinite
from matplotlib           import colors
from matplotlib.ticker    import LogFormatter

mpl.rcParams['font.family']     = 'serif'
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


def dirField_2(dxdt, dydt, fig, ax, nx=30, norm_bkg=False, Umin=0.05,
               cmap='jet', scale='lin', x_params=None, y_params=None):
  """
  Plot the direction field of dx2/dx1 = f2/f1 on an axes object ax.
  """
  cmap=get_cmap(cmap)

  print ax

  xmin = ax.get_xlim()[0]
  xmax = ax.get_xlim()[1]
  ymin = ax.get_ylim()[0]
  ymax = ax.get_ylim()[1]

  xstep  = (xmax - xmin)/float(nx)
  ystep  = (ymax - ymin)/float(nx)
  
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
  u   = dxdt(x, y, x_params) + 1e-16
  v   = dydt(x, y, y_params) + 1e-16
  
  Unorm = sqrt(u**2 + v**2) + 1e-16
  u /= Unorm
  v /= Unorm

  #=============================================================================
  # plotting :
  gray = '#5f5f5f'                           # color for axes
  ax.axhline(lw=1, c=gray)                   # x-axis
  ax.axvline(lw=1, c=gray)                   # y-axis
 
  if scale == 'log':
    norm = colors.LogNorm()
    Unorm[Unorm < Umin] = Umin
  else:
    norm = None
  
  if norm_bkg:
    axn = ax.imshow(Unorm, aspect='auto', 
                           cmap=cmap,
                           extent=(xmin,xmax,ymin,ymax),
                           norm=norm)
    axq = ax.quiver(x, y, u, v, pivot='middle',
                                cmap=cmap,
                                headwidth=3.5, 
                                headlength=2.0, 
                                headaxislength=2.0)  # plot the dir. field
  else:
    c   = Unorm
    axn = ax.quiver(x, y, u, v, c, pivot='middle',
                                   cmap=cmap,
                                   norm=norm,
                                   headwidth=3.5, 
                                   headlength=2.0, 
                                   headaxislength=2.0)  # plot the dir. field
  fig.colorbar(axn)


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
  f  = f(x, y) + 1e-16
  
  v  = sqrt((r**2) / (1 + 1/f**2))           # length of arrow in y-dir
  u  = v/f                                   # length of arrow in x-dir
  
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

