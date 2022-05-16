"""
Main Interaction module.
"""
import numpy as np
import itertools
#import sklearn.neighbors
import ctypes as C

LJforces = C.CDLL('interaction/LennardJones.so')
COforces = C.CDLL('interaction/Coulomb.so')
Fforces = C.CDLL('interaction/Flocking.so')

class Interaction(object):
  """
  Base Interaction class.
  """
  def __init__(self, types):
    self.types = types

  def forces(self, x, v, t):
    """
    Main loop calculation.

    NOTE: This is highly experimental and slow.
    It is just meant to be a proof of concept for the main loop, but
    has to change *dramatically*, even in arity, when we want to add
    lists of neighbors, parallelization and so on.
    """
    return np.zeros_like(x), 0.0

class ShortRange(Interaction):
  """
  Base short-range class
  """
  def __init__(self, types, rcut, shift_style='None'):
    """
    Base short-range class

    Parameters
    ----------

    rcut : float
        The cut radius parameter

    shift_style: {'None', 'Displace', 'Splines'}
        Shift style when approaching rcut

    .. note:: 'Splines' not implemented yet
    """
    self.rcut = np.array(rcut)
    self.shift_style = shift_style
    Interaction.__init__(self, types)

class FlockingSimulation(ShortRange):
  """
  Flocking Simulation potential
  """
  def __init__(self, types, rcut, coef01, coef02, coef03, shift_style='None'):
    self.coef01 = coef01
    self.coef02 = coef02
    self.coef03 = coef03
    self.rcut = rcut
    ShortRange.__init__(self, types, rcut, shift_style)

  def forces(self, x, v, charge, t, box):
    """
    Calculos de las fuerzas segun el potencial de  Flocking Simulation utilizando python
    """
    forces = np.zeros_like(x)
    energ = 0
    coef_separation = self.coef01
    coef_alignment =  self.coef02
    coef_cohesion =   self.coef03

    for i, s1 in enumerate(x):
      for j, s2 in enumerate(x[i+1:]):
          d = np.linalg.norm(s1-s2)
          if d < self.rcut:
            r = (s1-s2)/d
            dv = np.nan_to_num(np.linalg.norm(v[i]-v[i+1+j]))
            rv = np.nan_to_num((v[i]-v[i+1+j])/dv)

            f_sep = self.separation( d, r, coef_separation)
            f_alig = self.alignment( dv, rv, coef_alignment)
            f_cohe = self.cohesion(  d, r, coef_cohesion)

            forces[i]     += f_sep + f_alig + f_cohe
            forces[i+1+j] -= f_sep + f_alig + f_cohe

    forces = np.nan_to_num(forces / np.linalg.norm(forces, axis=1)[:,np.newaxis])
    return forces, energ

  def separation(self, d, r, coef):
    force = coef[0] + coef[1]*r
    return force

  def alignment(self, dv, rv, coef):
    force = coef[0] - coef[1]*rv
    return force

  def cohesion(self, d, r, coef):
    force = coef[0] - coef[1]*r 
    return force

  def Cforces(self, x, v, charge, t, box):
    """
    Calculos de las fuerzas segun el potencial de LennardJones utilizando la libreia dinamica floking.so 
    """
    in1 = x.astype( C.c_float)
    in2 = v.astype( C.c_float)
    in_charge = charge.astype( C.c_float)

    out1 = np.zeros(x.shape, dtype=np.float32)
    out2 = np.zeros(1, dtype=np.float32)

    in3 = self.rcut.astype( C.c_float)
    in4 = box.astype( C.c_float)

    coef_separation = np.array([0,1,0]).astype( C.c_float)
    coef_alignment  = np.array([0,-100,0]).astype( C.c_float)
    coef_cohesion   = np.array([0,-2,0]).astype( C.c_float)

    fltp = C.POINTER(C.c_float)
    Fforces.forces( 
               in1.ctypes.data_as(fltp),
               in2.ctypes.data_as(fltp),
               in_charge.ctypes.data_as(fltp),

               out1.ctypes.data_as(fltp),
               out2.ctypes.data_as(fltp),
               
               C.c_int(x.shape[0]),
               in3.ctypes.data_as(fltp),
               in4.ctypes.data_as(fltp),

               coef_separation.ctypes.data_as(fltp),
               coef_alignment.ctypes.data_as(fltp),
               coef_cohesion.ctypes.data_as(fltp),
           )

    return out1, out2[0]

class Coulomb(ShortRange):
  """
  Coulomb potential
  """
  def __init__(self, types, rcut, eps, sigma, shift_style='None'):
    self.K = eps
    ShortRange.__init__(self, types, rcut, shift_style)

  def Cforces(self, x, v, charge, t, box):
    """
    Calculos de las fuerzas segun el potencial de LennardJones utilizando la libreia dinamica Lennajones.so 
    """
    in1 = x.astype( C.c_float)
    in2 = v.astype( C.c_float)
    in_charge = charge.astype( C.c_float)
    out1 = np.zeros(x.shape, dtype=np.float32)
    out2 = np.zeros(1, dtype=np.float32)
    in3 = self.rcut.astype( C.c_float)
    in4 = box.astype( C.c_float)

    fltp = C.POINTER(C.c_float)
    COforces.forces( 
               in1.ctypes.data_as(fltp),
               in2.ctypes.data_as(fltp),
               in_charge.ctypes.data_as(fltp),

               out1.ctypes.data_as(fltp),
               out2.ctypes.data_as(fltp),
               
               C.c_int(x.shape[0]),
               in3.ctypes.data_as(fltp),
               in4.ctypes.data_as(fltp),
           )
    
    return out1, out2[0]

  def Cudaforces(self, x, v, charge, t, box):
    @vectorize(['float32(float32, float32)'], target='gpu')
    def dist(x ):
      return 
    """
    Calculos de las fuerzas segun el potencial de LennardJones utilizando la libreia dinamica Lennajones.so 
    """
    in1 = x.astype( C.c_float)
    in2 = v.astype( C.c_float)
    in_charge = charge.astype( C.c_float)
    out1 = np.zeros(x.shape, dtype=np.float32)
    out2 = np.zeros(1, dtype=np.float32)
    in3 = self.rcut.astype( C.c_float)
    in4 = box.astype( C.c_float)

    fltp = C.POINTER(C.c_float)
    COforces.forces( 
               in1.ctypes.data_as(fltp),
               in2.ctypes.data_as(fltp),
               in_charge.ctypes.data_as(fltp),

               out1.ctypes.data_as(fltp),
               out2.ctypes.data_as(fltp),
               
               C.c_int(x.shape[0]),
               in3.ctypes.data_as(fltp),
               in4.ctypes.data_as(fltp),
           )
    
    return out1, out2[0]

class LennardJones(ShortRange):
  """
  Lennard-Jones potential
  """
  def __init__(self, types, rcut, eps, sigma, shift_style='None'):
    self.eps = eps
    self.sigma = sigma
    ShortRange.__init__(self, types, rcut, shift_style)

    self.fltp = C.POINTER(C.c_float)
    
    
    self.E = np.zeros(1, dtype=np.float32)
    self.in_rcut = self.rcut.astype( C.c_float).ctypes.data_as(self.fltp)

  def Cforces(self, x, v, charge, t, box):
    """
    Calculos de las fuerzas segun el potencial de LennardJones utilizando la libreia dinamica Lennajones.so 
    """
    in1 = x.astype( C.c_float)
    in2 = v.astype( C.c_float)
    #out1 = np.zeros(x.shape, dtype=np.float32)
    #out2 = np.zeros(1, dtype=np.float32)
    
    in4 = box.astype( C.c_float)
    self.forces = np.zeros(x.shape, dtype=np.float32)
    
    LJforces.forces( in1.ctypes.data_as(self.fltp),
                    in2.ctypes.data_as(self.fltp),
                    self.forces.ctypes.data_as(self.fltp),
                    self.E.ctypes.data_as(self.fltp),
               
                    C.c_int(x.shape[0]),

                    self.in_rcut,
                    in4.ctypes.data_as(self.fltp),
           )
    
    return self.forces , self.E[0]

  def forces(self, x, v, t, box):
    """
    Calculate Lennard-Jones force
  
    """
    x1 = x[t == self.types[0]]
    x2 = x[t == self.types[1]]
    i1 = np.arange(len(x))[t == self.types[0]]
    i2 = np.arange(len(x))[t == self.types[1]]
    forces = np.zeros_like(x)
    energ = 0
    # I have to split it to avoid double-counting. Don't want to get
    # too fancy since it will change when creating neighbor lists
    if self.types[0] == self.types[1]:
      for i, s1 in enumerate(x1):
        for j, s2 in enumerate(x2[i+1:]):
          f = self.pair_force(s1, s2)
          ii = i1[i]
          #print s1 
          #if np.abs(np.sum(f)) > 0.001:
               #print i, j+i+1, f
          jj = i2[j+i+1]
          forces[ii] += f
          forces[jj] -= f 
          energ += self.pair_energ(s1, s2)
    else:
      for i, s1 in enumerate(x1):
        for j, s2 in enumerate(x2):
          f = self.pair_force(s1, s2)
          ii = i1[i]
          jj = i2[j]
          forces[ii] += f
          forces[jj] -= f
          energ += self.pair_energ(s1, s2)
    return forces, energ

  def forces_neiboarhood(self, x, v, t, nh, l, rcut):
      """
      Calculate Lennard-Jones force
      Neiboarhood list implementation

      """
      x1 = x[t == self.types[0]]
      x2 = x[t == self.types[1]]
      i1 = np.arange(len(x))[t == self.types[0]]
      i2 = np.arange(len(x))[t == self.types[1]]
      forces = np.zeros_like(x)
      energ = 0
      pairs = [] 

      if self.types[0] == self.types[1]:
          for i, s1 in enumerate(x):
              for n in itertools.product( [0,1], repeat=3):
                   
                    for m in nh[int(s1[0]/rcut)+n[0]][int(s1[1]/rcut)+n[1]][int(s1[2]/rcut)+n[2]]:
                        if i != m  and not [i, m] in pairs and not [m, i] in pairs :
                          s2 = x[m,:]
                          #print m, i                            
                          f = self.pair_force(s1, s2)
                          forces[i] += f
                          forces[m] -= f
                    
                          energ += self.pair_energ(s1, s2)
                          pairs.append([i, m])

      return forces, energ

  def force_tree(self, x, v, t, Mtree, rcut):
      """
      sklearn neighbors KDTree implementation
      Mtree : is a KDTree (from library SKlearn)
      """

      forces = np.zeros_like(x)
      energ = 0
      for i, s1 in enumerate(x):
          for m in Mtree.query_radius(s1.reshape(1,-1), r=rcut)[0]:
              if m != i and m>i:
                 s2 = x[m,:]
                 f = self.pair_force(s1, s2)
                 forces[i] += f
                 forces[m] -= f
                 energ += self.pair_energ(s1, s2)
      return forces, energ

  def pair_force(self, s1, s2):
    d = np.linalg.norm(s1-s2)
    if d > self.rcut:
      return np.zeros_like(s1)
    ljf = 24*self.eps*(2*self.sigma**12/d**14 - self.sigma**6/d**8)*(s1-s2)
    if self.shift_style == 'None':
      return ljf
    elif self.shift_style == 'Displace':
      return ljf

  def pair_energ(self, s1, s2):
    vcut = 4*self.eps*(self.sigma**12/self.rcut**12 - self.sigma**6/self.rcut**6)
    d = np.linalg.norm(s1-s2)
    if d >= self.rcut:
      return 0
    ljf = 4*self.eps*(self.sigma**12/d**12 - self.sigma**6/d**6)
    if self.shift_style == 'None':
      return ljf
    elif self.shift_style == 'Displace':
      return ljf - vcut



class Morse(ShortRange):
  """
  Morse  potential
  """
  def __init__(self, types, rcut, bener, blem, beta, shift_style='None'):
     self.bener = bener
     self.beta = beta
     self.blem = blem
     ShortRange.__init__(self, types, rcut, shift_style)

  def forces(self, x, v, t):
     """
     Calculate Lennard-Jones force
     """
     x1 = x[t == self.types[0]]
     x2 = x[t == self.types[1]]
     i1 = np.arange(len(x))[t == self.types[0]]
     i2 = np.arange(len(x))[t == self.types[1]]
     forces = np.zeros_like(x)
     energ = 0
            # I have to split it to avoid double-counting. Don't want to get
            # too fancy since it will change when creating neighbor lists
     if self.types[0] == self.types[1]:
          for i, s1 in enumerate(x1):
               for j, s2 in enumerate(x2[i+1:]):
                   f = self.pair_force(s1, s2)
                   ii = i1[i]
                   jj = i2[j+i+1]
                   forces[ii] += f
                   forces[jj] -= f
                   energ += self.pair_energ(s1, s2)
     else:
         for i, s1 in enumerate(x1):
              for j, s2 in enumerate(x2):
                  f = self.pair_force(s1, s2)
                  ii = i1[i]
                  jj = i2[j]
                  forces[ii] += f
                  forces[jj] -= f
                  energ += self.pair_energ(s1, s2)
     return forces, energ



  def pair_force(self, s1, s2):
      d = np.linalg.norm(s1-s2)
      if d > self.rcut:
         return np.zeros_like(s1)
      morf = 2.0*self.bener*(1.0-np.exp(-self.beta*(d-self.blem)))*np.exp(-self.beta*(d-self.blem))*self.beta*(s1-s2)/d
      if self.shift_style == 'None':
         return morf
      elif self.shift_style == 'Displace':
         return morf
  
  def pair_energ(self, s1, s2):
      vcut = -self.bener*(1.0-np.exp(-self.beta*(self.rcut-self.blem)))**2
      d = np.linalg.norm(s1-s2)
      if d >= self.rcut:
         return 0
      morf = -self.bener*(1.0-np.exp(-self.beta*(d-self.blem)))**2
      if self.shift_style == 'None':
         return morf
      elif self.shift_style == 'Displace':
         return morf - vcut

class NeuralNetwork(ShortRange):

    """Class that implements a basic feed-forward neural network.

    Parameters
    ----------
    hiddenlayers : dict
        Dictionary of chemical element symbols and architectures of their
        corresponding hidden layers of the conventional neural network. Number
        of nodes of last layer is always one corresponding to energy.  However,
        number of nodes of first layer is equal to three times number of atoms
        in the system in the case of no descriptor, and is equal to length of
        symmetry functions of the descriptor. Can be fed using tuples as:

        >>> hiddenlayers = (3, 2,)

        for example, in which a neural network with two hidden layers, the
        first one having three nodes and the second one having two nodes is
        assigned (to the whole atomic system in the no descriptor case, and to
        each chemical element in the atom-centered mode). When setting only one
        hidden layer, the dictionary can be fed as:

        >>> hiddenlayers = (3,)

        In the atom-centered mode, neural network for each species can be
        assigned seperately, as:

        >>> hiddenlayers = {"O":(3,5), "Au":(5,6)}

        for example.

    activation : str
        Assigns the type of activation funtion. "linear" refers to linear
        function, "tanh" refers to tanh function, and "sigmoid" refers to
        sigmoid function.
    weights : dict
        In the case of no descriptor, keys correspond to layers and values are
        two dimensional arrays of network weight.  In the atom-centered mode,
        keys correspond to chemical elements and values are dictionaries with
        layer keys and network weight two dimensional arrays as values. Arrays
        are set up to connect node i in the previous layer with node j in the
        current layer with indices w[i,j]. The last value for index i
        corresponds to bias. If weights is not given, arrays will be randomly
        generated.
    scalings : dict
        In the case of no descriptor, keys are "intercept" and "slope" and
        values are real numbers. In the fingerprinting scheme, keys correspond
        to chemical elements and values are dictionaries with "intercept" and
        "slope" keys and real number values. If scalings is not given, it will
        be randomly generated.
    fprange : dict
        Range of fingerprints of each chemical species.  Should be fed as
        a dictionary of chemical species and a list of minimum and maximun,
        e.g.:

        >>> fprange={"Pd": [0.31, 0.59], "O":[0.56, 0.72]}

    regressor : object
        Regressor object for finding best fit model parameters, e.g. by loss
        function optimization via amp.regression.Regressor.
    mode : str
        Can be either 'atom-centered' or 'image-centered'.
    lossfunction : object
        Loss function object, if at all desired by the user.
    checkpoints : int
        Frequency with which to save parameter checkpoints upon training. E.g.,
        100 saves a checpoint on each 100th training setp.  Specify None for no
        checkpoints.

    RuntimeError, NotImplementedError
    """

    def __init__(self, hiddenlayers=(5, 5), activation='tanh', weights=None, descriptor=None,
                 scalings=None, fprange=None, mode=None, checkpoints=100):

        # The parameters dictionary contains the minimum information
        # to produce a compatible model; e.g., one that gives
        # the identical energy (and/or forces) when fed a fingerprint.
        self.importname = '.model.neuralnetwork.NeuralNetwork'
        self.hiddenlayers = hiddenlayers
        self.weights = weights
        self.scalings = scalings
        self.fprange = fprange
        self.activation = activation
        self.mode = mode

        # Checking that the activation function is given correctly:
        if activation not in ['linear', 'tanh', 'sigmoid']:
            _ = ('Unknown activation function %s; must be one of '
                 '"linear", "tanh", or "sigmoid".' % activation)
            raise NotImplementedError(_)

        self.parent = None  # Can hold a reference to main Amp instance.
        self.checkpoints = checkpoints

    def calculate_ALL(self, x, t, f, E, G, box):
      in1 = x.astype( C.c_float) # positions
      in2 = t.astype( C.c_float)  # atom tipe
      in3 = f.astype( C.c_float)  # forces 
      in4 = E.astype( C.c_float)   # energy
      
      in5 = G.G2.astype( C.c_float)   # G2
      in6 = G.nG2.astype( C.c_int)

      in7 = G.G4_eta.astype( C.c_float)  # G4
      in7 = G.nG4_eta.astype( C.c_float)
      in7 = G.G4_gamma.astype( C.c_float)
      in7 = G.nG4_gamma.astype( C.c_float)
      in7 = G.G4_zeta.astype( C.c_float)
      in7 = G.nG4_zeta.astype( C.c_float)
      in8 = G.nG4.astype( C.c_int)

      fltp = C.POINTER(C.c_float)


      pass

    def calculate_fingerprints(self,):
      pass

    def calculate_force(self,):
      pass

    def calculate_energy(self,):
      pass

def Cforces(self, x, v, t, box):
    """
    Calculos de las fuerzas segun el potencial de LennardJones utilizando la libreia dinamica Lennajones.so 
    """
    in1 = x.astype( C.c_float)
    in2 = v.astype( C.c_float)
    out1 = np.zeros(x.shape, dtype=np.float32)
    out2 = np.zeros(1, dtype=np.float32)
    in3 = self.rcut.astype( C.c_float)
    in4 = box.astype( C.c_float)

    fltp = C.POINTER(C.c_float)
    LJforces.forces( in1.ctypes.data_as(fltp),
               in2.ctypes.data_as(fltp),
               out1.ctypes.data_as(fltp),
               out2.ctypes.data_as(fltp),
               
               C.c_int(x.shape[0]),
               in3.ctypes.data_as(fltp),
               in4.ctypes.data_as(fltp),
           )
    
    return out1, out2[0]


def interaction_config( parameters ):
  # make the correct interaction object with the required parameters #
  if parameters['type'] == 'LJ':
      ''' eg.
                           {  'type':  'LJ', 
                              'rcut':  3.5, 
                              'eps':   1.0, 
                              'sigma': 1.0, 
                              'shift_style':'Displace' }, ;
      ''' 
      interaction = LennardJones(   (0,0), 
                                    rcut=       parameters['rcut'] , 
                                    eps=        parameters['eps'], 
                                    sigma=      parameters['sigma'], 
                                    shift_style= parameters['shift_style'])
   
  elif parameters['type'] == 'Coulomb':
      ''' eg.
                           {  'type':  'Coulomb', 
                              'rcut':  5.5, 
                              'eps':   1.0, 
                              'sigma': 1.0, 
                              'shift_style':'Displace' },  
      '''
      interaction = Coulomb(  (0,0), 
                              rcut=       parameters['rcut'], 
                              eps=        parameters['eps'], 
                              sigma=      parameters['sigma'], 
                              shift_style= parameters['shift_style'])

  elif parameters['type'] == 'flocking':
      ''' eg.
                           {  'type':  'flocking', 
                              'rcut':  3.5, 
                              'coef01':  [0,1,0],
                              'coef02':  [0,-100,0],
                              'coef03':  [0,-3,0], 
                              'correction':'Displace' },  
      '''
      interaction = FlockingSimulation(   (0,0), 
                                          rcut=       parameters['rcut'], 
                                          coef01=     parameters['coef01'] , 
                                          coef02=     parameters['coef02'], 
                                          coef03=     parameters['coef03'], 
                                          shift_style= parameters['shift_style'])
  else: interaction = None

  if interaction == None: print('Error :: Can NOT load interaction')
  
  return interaction

