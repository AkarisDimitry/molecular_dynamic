"""
Main Particles module.
"""
import copy 
import numpy as np

class Base(object):
  """
  Base Particles class. It is abstract and we should specify which
  type of particle we actually want in order to fill it
  """
  def __init__(self, N=1, x=None, v=None, t=None, f=None, mass=None, idn=None):

    # == alocate memory == #
    self.N = N
    self._x = x
    self._v = v
    self._t = t
    self._f = f

    self._mass = mass
    self.idn = idn
    self.idx = np.arange(N)

    self.store_data = { 'x':  None,
                        'v':  None,
                        'f':  None,
                        'Ep': None,
                        'Ec': None,}
    
    self.id_dict = { i:n for i, n in enumerate(['H',  'He', 
                                                'Li', 'Be', 'B',  'N',  'C', 'O',  'F',  'Ne',
                                                'Na', 'Mg', 'Al', 'Si', 'P', 'S',  'Cl', 'Ar', 
                                                'K',  'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
                                                'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I',  'Xe', 
                                                'Cs', 'Ba', ]) }   

  @property
  def x(self):
    return self._x

  @x.setter
  def x(self, value):
    """
    Set positions of particles.

    Parameters
    ----------

    value : 2D NumPy array
        The new positons of particles in an Nx3 array
    """
    number = np.shape(value)[0]
    if self.N == number:
      self._x = value
    else:
      msg = "Trying to set {0} positions for a system with {1} particles"
      raise ValueError(msg.format(number, self.N))

  @property
  def v(self):
    return self._v

  @v.setter
  def v(self, value):
    """
    Set velocities of particles.

    Parameters
    ----------

    value : 2D NumPy array
        The new velocities of particles in an Nx3 array
    """
    number = np.shape(value)[0]
    if self.N == number:
      self._v = value
    else:
      msg = "Trying to set {0} velocities for a system with {1} particles"
      raise ValueError(msg.format(number, self.N))

  @property
  def f(self):
    return self._f

  @f.setter
  def f(self, value):
    """
    Set forces of particles.

    Parameters
    ----------

    value : 2D NumPy array
        The new forces of particles in an Nx3 array
    """
    number = np.shape(value)[0]
    if self.N == number:
      self._f = value
    else:
      msg = "Trying to set {0} forces for a system with {1} particles"
      raise ValueError(msg.format(number, self.N))

  @property
  def a(self):
    return self._f/self._mass[:, np.newaxis]

  @property
  def t(self):
    return self._t

  @t.setter
  def t(self, value):
    """
    Set types of particles.

    Parameters
    ----------

    value : 1D NumPy array or integer
        The new types of particles in an Nx3 array
    """
    if type(value) == int:
      self._t = np.array([value]*self.N)
    else:
      number = np.shape(value)[0]
      if self.N == number:
        self._f = value
      else:
        msg = "Trying to set {0} types for a system with {1} particles"
        raise ValueError(msg.format(number, self.N))

  @property
  def mass(self):
    return self._mass

  @mass.setter
  def mass(self, value):
    """
    Set types of particles.

    Parameters
    ----------

    value : 1D NumPy array or integer
        The new types of particles in an Nx3 array
    """
    if type(value) == float:
      self._mass = np.array([value]*self.N)
    else:
      number = np.shape(value)[0]
      if self.N == number:
        self._mass = value
      else:
        msg = "Trying to set {0} masses for a system with {1} particles"
        raise ValueError(msg.format(number, self.N))

  @property
  def idn(self):
    return self._idn

  @idn.setter
  def idn(self, value):
    """
    Set types of particles.

    Parameters
    ----------

    value : 1D NumPy array or integer
        The new types of particles in an Nx3 array
    """
    try:
      if type(value) == float:
        self._idn = np.array([value]*self.N)
      else:
        number = np.shape(value)[0]
        if self.N == number:
          self._idn = value
        else:
          msg = "Trying to set {0} idn for a system with {1} particles"
          raise ValueError(msg.format(number, self.N))
    except:
      pass

  def alocate_memory(self, steps=None, N=None):
    N = self.N if type(N) == type(None) else N
    steps = self.steps if type(steps) == type(None) else steps

    if type(self.x) == None:         self.x = np.zeros((N, 3), dtype=np.float32)
    if type(self.v) == None:         self.v = np.zeros((N, 3), dtype=np.float32)
    if type(self.f) == None:         self.f = np.zeros((N, 3), dtype=np.float32)
    if type(self.t) == None:         self.t    = np.zeros(N, dtype=np.int32)
    if type(self.mass) == None:      self.mass = np.ones(N, dtype=np.float32)
    if type(self.idn) == None:       self.idn  = np.ones(N, dtype=np.float32)
    if type(self.idx) == None:       self.idx  = np.ones(N, dtype=np.float32)

    self.store_data = {
                        'x':  np.zeros((steps, N, 3), dtype=np.float32),
                        'v':  np.zeros((steps, N, 3), dtype=np.float32),
                        'f':  np.zeros((steps, N, 3), dtype=np.float32),
                        'Ep': np.zeros((steps, ), dtype=np.float32),
                        'Ec': np.zeros((steps, ), dtype=np.float32),}
    return None

  def store_step(self, Ec=None, Ep=None, step=None):

    self.store_data['x'][step, :, :] = self.x
    self.store_data['v'][step, :, :] = self.v
    self.store_data['f'][step, :, :] = self.f

    self.store_data['Ec'][step] = np.sum(self.v**2)/2
    self.store_data['Ep'][step] = self.E

    return None

  def load(self, particles):

   if    particles['init'] == 'random':   self.load_random( N=              particles['n'], 
                                                            dimentionality= particles['dimentionality'], 
                                                            velocities=     'zeros',
                                                            position_radius=particles['position_radius'])

   elif  particles['init'] == 'load':     self.load_POSFILE(particles['filename'])

   elif  particles['init'] == 'material': self.load_material( RD=particles['RD'], rep=particles['rep'], )

   self.N = self.x.shape[0]

   return None

  def load_VASP(self, file_name, direct=1):
    """
    Load from VASP.
    POSCAR or CONTCAR files

    Parameters
    ----------

    value : 1D NumPy array or integer
        The new types of particles in an Nx3 array

    """
    file = open(file_name, 'r')
    cell = np.zeros((3,3)) ; atom_name = [] ; atom_numbers = [] ; N = 0 ; atom_position = [] ; atom_allposition = []; direct = direct

    for i, n in enumerate(file):
        vec = [m for m in n.split(' ') if m != '']
        vec[-1] = vec[-1][:-1]
             
        for m in range(3):
            if i == 2+m:  cell[m,0]=float(vec[0]);  cell[m,1]=float(vec[1]); cell[m,2]=float(vec[2])
       
        # add atoms names categories
        if i == 5: 
            for m in vec: 
              if m != '': atom_name.append(m)
       
        # atom numbers
        if i == 6:
            for m in vec:  N += int(m); atom_numbers.append(int(m)); atom_position.append( np.zeros((int(m), 3)) )
            atom_allposition = np.zeros((N, 3))

        if i == 8 and vec[0]=='Direct': direct = 1

        if i > 8:
            var = 9
            for j, m in enumerate(atom_numbers):
                if i >= var and i < var+m:  
                    atom_position[j][i-var,0]=float(vec[0]) ;  atom_position[j][i-var,1]=float(vec [1]) ;   atom_position[j][i-var,2]=float(vec[2])
                    atom_allposition[i-9,0]=float(vec[0])   ;  atom_allposition[i-9,1]=float(vec [1])   ;   atom_allposition[i-9,2]=float(vec[2])
                var += m
    if direct==1: atom_allposition = np.dot(atom_allposition, cell)

    return atom_name, atom_numbers, atom_position, atom_allposition, cell, N

  def load_material(self, velocities='rand', 
                          RD={'a1':np.array([1,0,0]), 'a2':[0,1,0], 'a3':[0,0,1]}, 
                          rep={'x':5, 'y':5, 'z':1}, ):

    RD = { key:np.array(vector) for key, vector in RD.items() }
    a1, a2, a3 = RD['a1'], RD['a2'], RD['a3']
    N = rep['x']*rep['y']*rep['z']

    self.N = N
    self.x = np.array( [ nx*a1+ny*a2+nz*a3 for nx in range(rep['x']) for ny in range(rep['y']) for nz in range(rep['z'])] )
    self.v = np.zeros((N, 3))

    # ---- prop ---- #
    self.t = np.zeros(N, dtype=np.int32)
    self.f = np.zeros((N, 3), dtype=np.float32)
    self.mass = np.ones(N, dtype=np.float32) 
    self.charge = np.ones(N, dtype=np.float32) * 10
    self.constrains = np.ones((N, 3), dtype=np.float32)
    self.idx = np.ones(N, dtype=np.int32)
    print(self.id_dict)
    self.idn = np.array([self.id_dict[x] for x in self.idx])

  def load_random(self, N=0, dimentionality=3, 
                        position='rand', velocities='rand', position_radius=1.0):
    if N == 0: print('Warnning :: Partibles.load_random :: particles number == 0')
    
    self.N = N

    # ---- position ---- #
    if position == 'non':     pass
    if position == 'rand':    self.x = np.ndarray.astype(np.random.rand(N, 3)*position_radius, dtype = np.float32)

    if dimentionality <= 1: self.x[:,1] = 0
    if dimentionality <= 2: self.x[:,2] = 0

    # ---- velocities ---- #
    if velocities == 'non':       pass
    elif velocities == 'zeros':   self.v = np.zeros((N, 3))
    elif velocities == 'rand':    self.v = np.random.rand(N, 3) 

    if dimentionality <= 1: self.v[:,1] = 0
    if dimentionality <= 2: self.v[:,2] = 0

    # ---- prop ---- #
    self.t = np.zeros(N, dtype=np.int32)
    self.f = np.zeros((N, 3), dtype=np.float32)
    self.mass = np.ones(N, dtype=np.float32) 
    self.charge = np.ones(N, dtype=np.float32) * 10
    self.constrains = np.ones((N, 3), dtype=np.float32)
    self.idx = np.ones(N, dtype=np.int32)
    self.idn = np.array([self.id_dict[x] for x in self.idx])

  def load_POSFILE(self, file_name):
    file = open(file_name, 'r')
    x, v, charge, constrains = [], [], [], []
    for i, n in enumerate(file):
      # primogenial embeding
      vec = n.replace('\t', ' ').split(' ')
      # eliminate '\n'
      if vec[-1][-1:] == '\n':  vec[-1] = vec[-1][:-1] 

      if i > 1 and len(vec) > 6:
        x.append( [float(vec[1]),float(vec[2]),float(vec[3])] )
        v.append( [float(vec[4]),float(vec[5]),float(vec[6])] ) 
        charge.append( float(vec[7]) )
        constrains.append( [float(vec[8]), float(vec[9]), float(vec[10])] )

    self.N = len(x)
    self.x = np.array(x)
    self.x.astype( dtype = np.float32)
    self.v = np.array(v)
    self.v.astype( dtype = np.float32)


    self.t = np.zeros(self.N, dtype=np.int32)
    self.f = np.zeros((self.N, 3), dtype=np.float32)
    self.mass = np.ones(self.N, dtype=np.float32) 

    self.charge = np.array(charge)
    self.charge.astype( dtype = np.float32)

    self.constrains = np.array(constrains)
    self.charge.astype( dtype = np.float32)

    self.idx = np.ones(N, dtype=np.int32)
    self.idn = np.array([self.id_dict[x] for x in self.idx])

  def save_OUTFILE(self, file_name):
    file = open(file_name, 'w')
    for i, step in enumerate(self.store_data['x']):
      file.write( '{} \n'.format(int(self.N)) )
      file.write( ' TIme: {} fs -- step count: {} \n'.format(i, i) )
      for j, atom in enumerate(step):
        file.write( '{} \t {} \t {} \t {} 10 0 0 \n'.format(self.idn[j], float(atom[0]), float(atom[1]), float(atom[2])) )
           
class PointParticles(Base):
  """
  PointParticles class.
  Typical point particles used, for example, in LJ potential.
  """
  pass

