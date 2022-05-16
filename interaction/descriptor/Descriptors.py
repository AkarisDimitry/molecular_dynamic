"""
Main Descriptor module
"""
from cutoffs import Cosine, dict2cutoff
import numpy as np
import ctypes as C

class Descriptor(object):
  """
  Descriptor class
  """
  def __init__(self, types):
    """
    Parameters
    ----------
    """
    self.types = types

class Gaussian(Descriptor):
    """
    Base short-range class
    """
    """Class that calculates Gaussian fingerprints (i.e., Behler-style).

    Parameters
    ----------
    cutoff : object or float
        Cutoff function, typically from amp.descriptor.cutoffs.  Can be also
        fed as a float representing the radius above which neighbor
        interactions are ignored; in this case a cosine cutoff function will be
        employed.  Default is a 6.5-Angstrom cosine cutoff.
    Gs : dict
        Dictionary of symbols and lists of dictionaries for making symmetry
        functions. Either auto-genetrated, or given in the following form, for
        example:

               >>> Gs = {"O": [{"type":"G2", "element":"O", "eta":10.},
               ...             {"type":"G4", "elements":["O", "Au"],
               ...              "eta":5., "gamma":1., "zeta":1.0}],
               ...       "Au": [{"type":"G2", "element":"O", "eta":2.},
               ...              {"type":"G4", "elements":["O", "Au"],
               ...               "eta":2., "gamma":1., "zeta":5.0}]}

    dblabel : str
        Optional separate prefix/location for database files, including
        fingerprints, fingerprint derivatives, and neighborlists. This file
        location can be shared between calculator instances to avoid
        re-calculating redundant information. If not supplied, just uses the
        value from label.
    elements : list
        List of allowed elements present in the system. If not provided, will
        be found automatically.
    mode : str
        Can be either 'atom-centered' or 'image-centered'.

    Raises
    ------
        RuntimeError
    """

    def __init__(self, cutoff=Cosine(6.5), G2=None, G4=None, 
                        elements=None, mode='atom-centered'):

        # Check that the mode is atom-centered.
        if mode != 'atom-centered':
            raise RuntimeError('Gaussian scheme only works '
                               'in atom-centered mode. %s '
                               'specified.' % mode)

        # If the cutoff is provided as a number, Cosine function will be used
        # by default.
        if isinstance(cutoff, int) or isinstance(cutoff, float):
            cutoff = Cosine(cutoff)


        # The parameters dictionary contains the minimum information
        # to produce a compatible descriptor; that is, one that gives
        # an identical fingerprint when fed an ASE image.

        self.cutoff = cutoff

        self.descriptor = None

        self.G2 = None
        self.nG2 = None
        self.G4_eta = None
        self.G4_gamma = None
        self.G4_zeta = None
        self.nG4_eta = None
        self.nG4_gamma = None
        self.nG4_zeta = None 
        self.elements = None

    def load_fromAMP(self, file_name):
        if hasattr(file_name, 'read'):
          text = file_name.read()
        else:
          with open(file_name) as f:
            text = f.read()

        p = self.string2dict(text)
        for key in ['descriptor', 'model']:
          p[key] = self.string2dict(p[key])

        self.descriptor = p['descriptor']
        if p['descriptor']['cutoff']['name'] == 'Cosine':
          self.cutoff = p['descriptor']['cutoff']['kwargs']['Rc']

        self.elements = p['descriptor']['elements']
        self.G2_eta = []; self.G4_eta = []; self.G4_gamma = [];  self.G4_zeta = []
        for n in self.elements:
          for m in p['descriptor']['Gs'][n]:
            if m['type'] == 'G2':  self.G2_eta.append( m['eta'] )
            if m['type'] == 'G4':  self.G4_eta.append( m['eta'] ); self.G4_gamma.append( m['gamma'] ); self.G4_zeta.append( m['zeta'] )

        self.G2_eta = np.array(list(set(self.G2_eta)))
        self.G4_eta = np.array(list(set(self.G4_eta)))
        self.G4_gamma = np.array(list(set(self.G4_gamma)))
        self.G4_zeta = np.array(list(set(self.G4_zeta)))
    
        self.nG2_eta = self.G2_eta.shape
        self.nG4_eta, self.nG4_gamma, self.nG4_zeta = self.G4_eta.shape, self.G4_gamma.shape, self.G4_zeta.shape 

        return p['model'], p['descriptor']

    def string2dict(self, text):
      """ converts a string into a dicctionary """
      try:
        dicctionary = eval(text)
      except:
        from collections import OrderedDict
        from numpy import array, matrix
        dicctionary = eval(text)

      return dicctionary

    def get_fingerprint(self, x, box, ):
      pass
















