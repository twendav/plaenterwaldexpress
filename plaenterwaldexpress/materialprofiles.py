from dolfin import *
from utilities import *

class MatParamBasics():
  """Parameter class with basic functionality."""
  def __init__(self, material, gradient):
    self.material = material
    self.gradient = gradient
  
  def get_value(self, coord, argument):
    if self.gradient:
      x = self.gradient(coord)
      param = self.material.get_value(argument, x)
    else:
      param = self.material.get_value(argument)

    return param

def materialKP8_profile_1d(material, gradient = None, asub = None):
  
  # Create and initiate class for valence band (potential)
  class MaterialPotential1D(Expression):
    def __init__(self, material, band, gradient):
      self.material = MatParamBasics(material, gradient)
      self.band = band

    def eval(self, value, coord):
      # Return band parameter
      if self.band == 'vb':
        value[0] = self.material.get_value(coord, 'Eav') + \
                   self.material.get_value(coord, 'delta_so') / 3
      else:
        value[0] = self.material.get_value(coord, 'EG')
  
  cpot = MaterialPotential1D(material, 'cb', gradient)
  vpot = MaterialPotential1D(material, 'vb', gradient)

  # Create and initiate class for spin-orbit splitting
  class MaterialDeltaSo1D(Expression):
    def __init__(self, material, gradient):
      self.material = MatParamBasics(material, gradient)

    def eval(self, value, coord):
      # Return band parameter
      value[0] = self.material.get_value(coord, 'delta_so')
  
  delta_so = MaterialDeltaSo1D(material, gradient)

  # Create and initiate class for band parameters
  class MaterialBand1D(Expression):
    def __init__(self, material, parameter, gradient):
      self.material = MatParamBasics(material, gradient)
      self.parameter = parameter

    def eval(self, value, coord):
      # Return band parameter
      value[0] = self.material.get_value(coord, self.parameter)
  
  Ac = MaterialBand1D(material, 'Ac', gradient) 
  P = MaterialBand1D(material, 'P', gradient)
  Lp = MaterialBand1D(material, 'Lp', gradient)
  M = MaterialBand1D(material, 'M', gradient)
  Np = MaterialBand1D(material, 'Np', gradient)
  Npp = MaterialBand1D(material, 'Npp', gradient)
  Nm = MaterialBand1D(material, 'Nm', gradient)
  kane = [Ac, P, Lp, M, Np, Npp, Nm]

  # Create and initiate class for strain tensor
  class MaterialEpsilon1D(Expression):
    def __init__(self, material, asub, direction, gradient):
      self.material = MatParamBasics(material, gradient)
      self.asub = asub
      self.direction = direction

    def eval(self, value, coord):
      eps_para = self.asub / self.material.get_value(coord, 'a0') - 1
        
      if self.direction == 'para':
        value[0] = eps_para
      else:
        cratio = self.material.get_value(coord, 'C12') / \
                   self.material.get_value(coord, 'C11')
        value[0] = -2 * cratio * eps_para
  
  if asub:
    eps_para = MaterialEpsilon1D(material, asub, 'para', gradient)
    eps_perp = MaterialEpsilon1D(material, asub, 'perp', gradient)
    epsilon = [eps_para, eps_perp]
  else:
    epsilon = None

  # Create and initiate class for deformation potentials
  class MaterialDeformation1D(Expression):
    def __init__(self, material, parameter, gradient):
      self.material = MatParamBasics(material, gradient)
      self.parameter = parameter

    def eval(self, value, coord):
      aG = self.material.get_value(coord, 'aG')
      av = self.material.get_value(coord, 'av')
      bv = self.material.get_value(coord, 'bv')

      if self.parameter == 'aG':
	value[0] = aG
      elif self.parameter == 'l':
        value[0] = av + 2 * bv
      elif self.parameter == 'm':
        value[0] = av - bv
      elif self.parameter == 'n':
        value[0] = 0
  
  def_aG = MaterialDeformation1D(material, 'aG', gradient)
  def_l = MaterialDeformation1D(material, 'l', gradient)
  def_m = MaterialDeformation1D(material, 'm', gradient)
  def_n = MaterialDeformation1D(material, 'n', gradient)
  defpot = [def_aG, def_l, def_m, def_n]
  
  return [cpot, vpot, delta_so, kane, epsilon, defpot]


def materialKP6_profile_1d(elementary, gradient = None, asub = None):

  # Create and initiate class for valence band (potential)
  class MaterialPotential1D(Expression):
    def __init__(self, material, gradient):
      self.material = MatParamBasics(material, gradient)

    def eval(self, value, coord):
      # Return band parameter
      value[0] = self.material.get_value(coord, 'Eav') + \
                 self.material.get_value(coord, 'delta_so') / 3
  
  vpot = MaterialPotential1D(elementary, gradient)

  # Create and initiate class for spin-orbit splitting
  class MaterialDeltaSo1D(Expression):
    def __init__(self, material, gradient):
      self.material = MatParamBasics(material, gradient)

    def eval(self, value, coord):
      # Return band parameter
      value[0] = self.material.get_value(coord, 'delta_so')
  
  delta_so = MaterialDeltaSo1D(elementary, gradient)

  # Create and initiate class for band parameters
  class MaterialBand1D(Expression):
    def __init__(self, material, parameter, gradient):
      self.material = MatParamBasics(material, gradient)
      self.parameter = parameter

    def eval(self, value, coord):
      # Return band parameter
      value[0] = self.material.get_value(coord, self.parameter)
  
  L = MaterialBand1D(elementary, 'L', gradient)
  M = MaterialBand1D(elementary, 'M', gradient)
  N = MaterialBand1D(elementary, 'N', gradient)
  Np = MaterialBand1D(elementary, 'Np', gradient)
  Nm = MaterialBand1D(elementary, 'Nm', gradient)
  kane = [L, M, N, Np, Nm]

  # Create and initiate class for strain tensor
  class MaterialEpsilon1D(Expression):
    def __init__(self, material, asub, direction, gradient):
      self.material = MatParamBasics(material, gradient)
      self.asub = asub
      self.direction = direction

    def eval(self, value, coord):
      eps_para = self.asub / self.material.get_value(coord, 'a0') - 1
        
      if self.direction == 'para':
        value[0] = eps_para
      else:
        cratio = self.material.get_value(coord, 'C12') / \
                   self.material.get_value(coord, 'C11')
        value[0] = -2 * cratio * eps_para
  
  if asub:
    eps_para = MaterialEpsilon1D(elementary, asub, 'para', gradient)
    eps_perp = MaterialEpsilon1D(elementary, asub, 'perp', gradient)
    epsilon = [eps_para, eps_perp]
  else:
    epsilon = None

  # Create and initiate class for deformation potentials
  class MaterialDeformation1D(Expression):
    def __init__(self, material, parameter, gradient):
      self.material = MatParamBasics(material, gradient)
      self.parameter = parameter

    def eval(self, value, coord):
      av = self.material.get_value(coord, 'av')
      bv = self.material.get_value(coord, 'bv')

      if self.parameter == 'l':
        value[0] = av + 2 * bv
      elif self.parameter == 'm':
        value[0] = av - bv
      elif self.parameter == 'n':
        value[0] = 0
  
  if asub:
    def_l = MaterialDeformation1D(elementary, 'l', gradient)
    def_m = MaterialDeformation1D(elementary, 'm', gradient)
    def_n = MaterialDeformation1D(elementary, 'n', gradient)
    defpot = [def_l, def_m, def_n]
  else:
    defpot = None
  
  return [vpot, delta_so, kane, epsilon, defpot]

def materialMEMA_profile_1d(material, band, gradient = None, asub = None):
  # Check if right material definition is given
  # NEEDS TO BE DONE! 

  # Check if given band identifier is allowed
  allowed_bands = ['hh', 'lh', 'so']
  if band not in allowed_bands:
    raise Exception("Given band '" + band + "' unknown.")

  # Create and initiate class for effective mass
  class MaterialMeff1D(Expression):
    def __init__(self, material, band, asub, gradient):
      self.material = MatParamBasics(material, gradient)
      self.band = band
      self.asub = asub

    def strain_tensor(self, coord):
      # Retrieve material parameters
      a0 = self.material.get_value(coord, 'a0')
      C11 = self.material.get_value(coord, 'C11')
      C12 = self.material.get_value(coord, 'C12')

      # Calculate pseudomorphic strain (if required)
      if self.asub:
        eps_xx = self.asub / a0 - 1
	eps_zz = -2. * C12 / C11 * eps_xx
      else:
        eps_xx = 0
	eps_zz = 0
      
      return (eps_xx, eps_zz)

    def fplusminus(self, coord, sign):
      # Calculate strain dependent part 
      (eps_xx, eps_zz) = self.strain_tensor(coord)
      
      Qeps = self.material.get_value(coord, 'bv') * (eps_zz - eps_xx)
      dso = self.material.get_value(coord, 'delta_so')
      x = Qeps / dso
      
      # Put full equation together
      xroot = sqrt(1 + 2 * x + 9 * x**2)
      fnom = 2 * x * (1 + 3./2. * (x - 1 + sign * xroot)) + 6 * x**2
      fden = 3./4. * (x - 1 + sign * xroot)**2 + x - 1 - sign * xroot - 3 * x**2

      return fnom / fden

    def eval(self, value, coord):
      # Convert LK to gamma parameters
      L = self.material.get_value(coord, 'L')
      M = self.material.get_value(coord, 'M')
      N = self.material.get_value(coord, 'N')

      (g1, g2, g3) = LK2gamma((L, M, N))
      
      # Calculate effective mass
      if self.band == 'hh':
        value[0] = 1 / (g1 - 2 * g2)
      elif self.band == 'lh':
        fp = self.fplusminus(coord, +1)
	value[0] = 1 / (g1 + 2 * fp * g2)
      elif self.band == 'so':
        fp = self.fplusminus(coord, -1)
	value[0] = 1 / (g1 + 2 * fp * g2)
      else:
        raise Exception("Given band '" + band + "' unknown.")

  meff = MaterialMeff1D(material, band, asub, gradient)

  # Create class for potential function
  class MaterialVpot1D(Expression):
    def __init__(self, material, band, asub, gradient):
      self.material = MatParamBasics(material, gradient)
      self.band = band
      self.asub = asub

    def strain_tensor(self, coord):
      # Retrieve material parameters
      a0 = self.material.get_value(coord, 'a0')
      C11 = self.material.get_value(coord, 'C11')
      C12 = self.material.get_value(coord, 'C12')

      # Calculate pseudomorphic strain (if required)
      if self.asub:
        eps_xx = self.asub / a0 - 1
	eps_zz = -2. * C12 / C11 * eps_xx
      else:
        eps_xx = 0
	eps_zz = 0
      
      return (eps_xx, eps_zz)

    def bandh(self, (eps_xx, eps_zz), coord,):
      # Calculate energy shift due to hydrostatic and sheer strain
      E001hyv = self.material.get_value(coord, 'av') * (2. * eps_xx + eps_zz)
      E001shv = self.material.get_value(coord, 'bv') * 2. * (eps_zz - eps_xx)

      dso = self.material.get_value(coord, 'delta_so')

      # Set potential according to hole type
      if self.band == 'hh':
        vpot = self.material.get_value(coord, 'Eav') + \
               E001hyv - E001shv / 2. + dso / 3.
      elif self.band == 'lh':
        vpot = self.material.get_value(coord, 'Eav') + \
	       E001hyv - dso / 6. + E001shv / 4. + \
	       0.5 * sqrt(dso**2 + dso * E001shv + 9. / 4. * E001shv**2)
      elif self.band == 'so':
        vpot = self.material.get_value(coord, 'Eav') + \
               E001hyv - dso / 6. + E001shv / 4. - \
               0.5 * sqrt(dso**2 + dso * E001shv + 9. / 4. * E001shv**2)
 
      return vpot

    def eval(self, value, coord):
      # Calculate strain tensor
      epsilon = self.strain_tensor(coord)

      # Calculate band energy
      if self.band == 'hh' or self.band == 'lh' or self.band == 'so':
        value[0] = self.bandh(epsilon, coord)
      else:
        raise Exception("Given band '" + band + "' unknown.")

 
  vpot = MaterialVpot1D(material, band, asub, gradient)

  return [vpot, meff]

def materialEMA_profile_1d(material, band, gradient = None, asub = None):
  
  # Check if given band identifier is allowed
  allowed_bands = ['G', 'L', 'X100', 'X010', 'X001', 'hh', 'lh', 'so']
  if band not in allowed_bands:
    raise Exception("Given band '" + band + "' unknown.")

  # Create and initiate class for effective mass
  class MaterialMeff1D(Expression):
    def __init__(self, material, band, gradient, direction=None):
      self.material = MatParamBasics(material, gradient)
      self.band = band
      self.direction = direction

    def eval(self, value, coord):
      # Calculate effective mass
      direction_table = {'l':0, 't': 0}
      
      if self.band[0] == 'X' or self.band == 'L':
        index = direction_table[self.direction]
        value[0] = self.material.get_value(coord, 'meff' + self.band[0])[index]
      else:
        value[0] = self.material.get_value(coord, 'meff' + self.band)
  
  if band[0] == 'X' or band == 'L':
    meff_l = MaterialMeff1D(material, band, gradient, 'l')
    meff_t = MaterialMeff1D(material, band, gradient, 't')

    meff = [meff_l, meff_t]
  else:
    meff = MaterialMeff1D(material, band, gradient)

  # Create class for potential function
  class MaterialVpot1D(Expression):
    def __init__(self, material, band, asub, gradient):
      self.material = MatParamBasics(material, gradient)
      self.band = band
      self.asub = asub

    def strain_tensor(self, coord):
      # Retrieve material parameters
      a0 = self.material.get_value(coord, 'a0')
      C11 = self.material.get_value(coord, 'C11')
      C12 = self.material.get_value(coord, 'C12')

      # Calculate pseudomorphic strain (if required)
      if self.asub:
        eps_xx = self.asub / a0 - 1
	eps_zz = -2. * C12 / C11 * eps_xx
      else:
        eps_xx = 0
	eps_zz = 0
      
      return (eps_xx, eps_zz)

    def bandGL(self, (eps_xx, eps_zz), coord):
      vpot = self.material.get_value(coord, 'E' + self.band) + \
             self.material.get_value(coord, 'a' + self.band) * (2 * eps_xx + eps_zz)
      
      return vpot
 
    def bandX(self, (eps_xx, eps_zz), coord):
      vpot = self.material.get_value(coord, 'EX') + \
             self.material.get_value(coord, 'aX') * (2 * eps_xx + eps_zz)
 
      if self.band == 'X001':
        vpot = vpot + \
  	       2. / 3. * self.material.get_value(coord, 'xiX') * (eps_zz - eps_xx)
      else:
        vpot = vpot - \
	       1. / 3. * self.material.get_value(coord, 'xiX') * (eps_zz - eps_xx)

      return vpot

    def bandh(self, (eps_xx, eps_zz), coord,):
      # Calculate energy shift due to hydrostatic and sheer strain
      E001hyv = self.material.get_value(coord, 'av') * (2. * eps_xx + eps_zz)
      E001shv = self.material.get_value(coord, 'bv') * 2. * (eps_zz - eps_xx)

      dso = self.material.get_value(coord, 'delta_so')

      # Set potential according to hole type
      if self.band == 'hh':
        vpot = self.material.get_value(coord, 'Eav') + \
               E001hyv - E001shv / 2. + dso / 3.
      elif self.band == 'lh':
        vpot = self.material.get_value(coord, 'Eav') + \
	       E001hyv - dso / 6. + E001shv / 4. + \
	       0.5 * sqrt(dso**2 + dso * E001shv + 9. / 4. * E001shv**2)
      elif self.band == 'so':
        vpot = self.material.get_value(coord, 'Eav') + \
               E001hyv - dso / 6. + E001shv / 4. - \
               0.5 * sqrt(dso**2 + dso * E001shv + 9. / 4. * E001shv**2)
 
      return vpot

    def eval(self, value, coord):
      # Calculate strain tensor
      epsilon = self.strain_tensor(coord)

      # Calculate band energy
      if self.band == 'G' or self.band == 'L':
        value[0] = self.bandGL(epsilon, coord)
      elif self.band == 'X100' or self.band == 'X010' or self.band == 'X001':
        value[0] = self.bandX(epsilon, coord)
      elif self.band == 'hh' or self.band == 'lh' or self.band == 'so':
        value[0] = self.bandh(epsilon, coord)
      else:
        raise Exception("Given band '" + band + "' unknown.")

 
  vpot = MaterialVpot1D(material, band, asub, gradient)

  return [vpot, meff]
   

