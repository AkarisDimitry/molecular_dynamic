"""
Main box module
"""

import numpy as np

class Neiboarhood(object):
  """
  Neiboarhood class
  """
  def __init__(self, x0, xf, p_cell, rcut):
    """
    Parameters
    ----------

    x0 : NumPy array
        Initial vertex of the box

    xf : NumPy array
        Final vertex of the box

    p_cell :  Numpy array
        Cell parameters

    rcut : float
            The cut radius parameter

    """

    self.x0 = np.array(x0)
    self.xf = np.array(xf)
    self.p_cell = p_cell
    self.l = (self.xf - self.x0)/self.p_cell
    
    self.rcut = rcut
    self.neiboarhoods = []
    
    for n in range(int(self.l[0])):
        vec02 = []
        for m in range(int(self.l[1])):
            vec01 = []
            for o in range(int(self.l[2])):
                vec01.append([])
            vec02.append(vec01)
        self.neiboarhoods.append(vec02)



  def new_neiboarhood(self, particles ):
      nh = []
      for n in range(int(self.l[0])+1):
        vec02 = []
        for m in range(int(self.l[1])+1):
           vec01 = []
           for o in range(int(self.l[2])+1):
               vec01.append([])
           vec02.append(vec01)
        nh.append(vec02)

      for i, n in enumerate(particles/self.rcut):
          nh[int(n[0])][int(n[1])][int(n[2])].append( i )
      self.neiboarhood = nh
      return nh




