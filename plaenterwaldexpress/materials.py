from dolfin import *

class ElementaryBands:
  """Band edge description of a material consisting of a single atomic type."""
  
  def __init__(self, param):
    self.param = dict()

    # Band energy parameters
    self.param['EG'] = param['EG']
    self.param['EL'] = param['EL']
    self.param['EX'] = param['EX']
    self.param['Eav'] = param['Eav']
    self.param['delta_so'] = param['delta_so']
  
    # Lattice parameters
    self.param['a0'] = param['a0']
    self.param['C11'] = param['C11']
    self.param['C12'] = param['C12']
    self.param['C44'] = param['C44']
  
    # Deformation potentials
    self.param['av'] = param['av'] 
    self.param['bv'] = param['bv']
    self.param['aG'] = param['aG']
    self.param['aL'] = param['aL']
    self.param['aX'] = param['aX'] 
    self.param['xiX'] = param['xiX']

  def get_value(self, property_id):
    """Returns specified material property."""
    return self.param[property_id]

class ElementaryEMA(ElementaryBands):
  """EMA description of a material consisting of a single atomic type."""
  def __init__(self, param):
    ElementaryBands.__init__(self, param)

    # Effective masses
    self.param['meffG'] = param['meffG']
    self.param['meffL'] = param['meffL']
    self.param['meffX'] = param['meffX']
    self.param['meffhh'] = param['meffhh']
    self.param['mefflh'] = param['mefflh']
    self.param['meffso'] = param['meffso']

class ElementaryKP6(ElementaryBands):
  """k.p 6-band description of a material consisting of a single atomic type."""
  def __init__(self, param):
    ElementaryBands.__init__(self, param)

    # Kane k.p 6-band parameters
    self.param['L'] = param['L']
    self.param['M'] = param['M']
    self.param['N'] = param['N']
    self.param['Np'] = param['Np']
    self.param['Nm'] = param['Nm']

class ElementaryKP8(ElementaryBands):
  """k.p 8-band description of a material consisting of a single atomic type."""
  def __init__(self, param):
    ElementaryBands.__init__(self, param)

    # Kane k.p 8-band parameters
    self.param['Ac'] = param['Ac']
    self.param['P'] = param['P']
    self.param['Lp'] = param['Lp']
    self.param['M'] = param['M']
    self.param['Np'] = param['Np']
    self.param['Npp'] = param['Npp']
    self.param['Nm'] = param['Nm']

class Binary:
  """Material consisting of two different types of atoms."""

  def __init__(self, elem_a, elem_b, bowing = None):
    """Alloy consisting of two elements with variable composition."""
    self.elem_a = elem_a.param
    self.elem_b = elem_b.param
    self.bowing = bowing
    self.parent_class = elem_a.__class__
  
  def make_elementary(self, x):
    """Converts binary alloy into a single atomic material."""
    elem_param = dict()
    
    for key in self.elem_a:
      elem_param[key] = self.get_value(key, x)
    

    return self.parent_class(elem_param)

  def get_value(self, property_id, x):
    """Computes specific alloy property based on material composition."""
    if isinstance(self.elem_a[property_id], tuple):
      # Calculate bowing term
      if self.bowing and property_id in self.bowing:
        bowing = tuple(x * (1 - x) * b for b in self.bowing[property_id])
      else:
        bowing = (0,) * len(self.elem_a[property_id]) 
      
      # Vegard's law plus bowing term
      value = tuple(x * va + (1 - x) * vb - b \
                for va, vb, b in \
 	        zip(self.elem_a[property_id], self.elem_b[property_id], bowing))

    else:
      # Calculate bowing term
      if self.bowing and property_id in self.bowing:
        bowing = self.bowing[property_id] * x * (1 - x)
      else:
        bowing = 0
      
      # Vegard's law plus bowing term
      value = self.elem_a[property_id] * x + \
        self.elem_b[property_id] * (1 - x) - bowing

    return value

