"""
Main Integrator module
"""

import numpy as np
import ctypes as C
velverlet = C.CDLL('integrator/c_first_step.so')

class Integrator(object):
  """
  Integrator base class
  """

  def __init__(self, dt):
    self.dt = dt

  def first_step(self, x, v, a):
    return x, v

  def last_step(self, x, v, a):
    return x, v


class NVE(Integrator):
  """
  NVE base class
  """

class VelVerlet(NVE):
  """
  Velocity Verlet Integrator
  """
  def first_step(self, x, v, a, constrain):
    x[constrain==1] = x[constrain==1] + v[constrain==1]*self.dt + 0.5*a[constrain==1]*self.dt**2
    v[constrain==1] = v[constrain==1] + 0.5*a[constrain==1]*self.dt
    return x, v

  def C_first_step(self, x ,v, a):
    in1 = x.astype( C.c_float )
    in2 = v.astype( C.c_float )
    in3 = a.astype( C.c_float )
    in4 = C.c_float( self.dt )
    out1 = np.zeros(x.shape, dtype=np.float32)
    out2 = np.zeros(x.shape, dtype=np.float32)


    fltp = C.POINTER(C.c_float)
    velverlet.first_step( in1.ctypes.data_as(fltp),
               in2.ctypes.data_as(fltp),
               in3.ctypes.data_as(fltp),
              
               out1.ctypes.data_as(fltp),
               out2.ctypes.data_as(fltp),
               C.byref( in4 ),

               C.c_int(x.shape[0]),         
           )
    return out1, out2  


  def last_step(self, x, v, a, constrain):
    v[constrain==1] = v[constrain==1] + 0.5*a[constrain==1]*self.dt
    return x, v

class NVT(Integrator):
  """
  NVT base class
  """
  def __init__(self, dt, temperature):
    self.temperature = temperature
    Integrator.__init()

class Andersen(NVT):
      """
      Andersen Thermostat
      """
      def __init__(self, dt, temperature, freq):
         self.freq = freq
         NVT.__init__(dt, temperature)

      def first_step(self, x, v, a):
         x = x + v*self.dt + 0.5*a*self.dt**2
         v = v + 0.5*a*self.dt
         return x, v

      def last_step(self, x, v, a):
         v = v + 0.5*a*self.dt
         p_i = np.random.rand()
         nparticles = len(x)
         if p_i < self.freq * self.dt:
              i_ran = np.random.randint(0,nparticles)
              v [i_ran, :] = np.random.normal(0.,np.sqrt(self.temperature),size = 3)
         return x, v




