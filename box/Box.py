"""
Main box module
"""

import numpy as np
import ctypes as C
Cbox = C.CDLL('box/Cbox_fixed.so')

class Box(object):
  """
  Box class
  """
  def __init__(self, x0, xf, t='Periodic', vecinos=0 ):
    """
    Parameters
    ----------

    x0 : NumPy array
        Initial vertex of the box

    xf : NumPy array
        Final vertex of the box

    t : {'Periodic', 'Fixed', 'Fixed Bounce'}
        Type of boundary
    """
    self.x0 = np.array(x0)
    self.xf = np.array(xf)
    self.xl = np.array(xf) -  np.array(x0)
    self.t = t
    self.fltp = C.POINTER(C.c_float)

  def wrap_boundary(self, x, v):
    """
    Apply boundary conditions

    Parameters
    ----------

    x, v : NumPy array
        Positions and velocities of the particles

    Returns
    -------

    x, v : NumPy array
        Positions and velocities updated
    """
    x_b = np.copy(x)
    v_b = np.copy(v)
    if self.t == 'Periodic':
      l = self.xf - self.x0
      for i, pos in enumerate(x):
        for j, p in enumerate(pos):
          if p > self.xf[j]:  p -= l[j]
          if p < self.x0[j]:  p += l[j]
          x_b[i, j] = p

    elif self.t == 'Fixed':
      l = self.xf - self.x0
      for i, pos in enumerate(x):
        for j, p in enumerate(pos):
          m = 1
          while (p > self.xf[j]) or (p < self.x0[j]):
            if p > self.xf[j]:
              m *= -1
              p = 2*self.xf[j] - p
            if p < self.x0[j]:
              p = 2*self.x0[j] - p
              m *= -1
          x_b[i, j] = p
          v_b[i, j] *= m
    return x_b, v_b

  def C_wrap_boundary(self, x, v):
    """
    Apply boundary conditions

    Parameters
    ----------

    x, v : NumPy array
        Positions and velocities of the particles

    Returns
    -------

    x, v : NumPy array
        Positions and velocities updated
    """

    if self.t == 'Periodic':
      x_b = np.copy(x)
      v_b = np.copy(v)
      l = self.xf - self.x0
      for i, pos in enumerate(x):
        for j, p in enumerate(pos):
          if p > self.xf[j]:  p -= l[j]
          if p < self.x0[j]:  p += l[j]
          x_b[i, j] = p
      x=x_b

    elif self.t == 'Fixed':
      x = x.astype( C.c_float)
      v = v.astype( C.c_float)
      Cbox.box_fixed( 
               x.ctypes.data_as(self.fltp), 
               v.ctypes.data_as(self.fltp),
               self.x0.astype( C.c_float).ctypes.data_as(self.fltp),
               self.xf.astype( C.c_float).ctypes.data_as(self.fltp),
               C.c_int( x.shape[0] ) )

    elif self.t == 'Fixed Bounce':
      x = x.astype( C.c_float)
      v = v.astype( C.c_float)
      Cbox.box_fixed( 
               x.ctypes.data_as(self.fltp), 
               v.ctypes.data_as(self.fltp),
               self.x0.astype( C.c_float).ctypes.data_as(self.fltp),
               self.xf.astype( C.c_float).ctypes.data_as(self.fltp),
               C.c_int( x.shape[0] ) )
    return x, v











