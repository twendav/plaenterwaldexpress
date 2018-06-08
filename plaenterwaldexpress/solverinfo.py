class SolverInfo:
  
  def __init__(self, num_ev = None):
    # Set default solver properties
    self.num_ev = num_ev
    self.solver_type = 'arnoldi'
    self.spectrum = 'smallest magnitude'
    self.max_iterations = 1000
    self.dolfin_version = None
    self.plaenterwald_version = None
 
  def set_tolerance(self, tolerance):
    self.tolerance = tolerance

  def set_spectral_shift(self, spectral_shift):
    self.spectral_shift = spectral_shift

  def set_shift_and_invert(self, shift_and_invert):
    self.set_shift_and_invert = shift_and_invert

  def set_spectral_transform(self, spectral_transform):
    self.spectral_transform = spectral_transform

  def set_spectrum(self, spectrum):
    self.spectrum = spectrum

  def set_max_iterations(self, max_iterations):
    self.max_iterations = max_iterations

  def set_solver_type(self, solver_type):
    self.solver_type = solver_type

  def get_solver_type(self):
    return self.solver_type

  def set_runtime(self, runtime):
    self.runtime = runtime

  def get_runtime(self):
    return self.runtime

  def set_required_iterations(self, num_iterations):
    self.num_iterations = num_iterations

  def get_required_iterations(self):
    return self.num_iterations

  def set_version(self, dolfin_version, plaenterwald_version):
    self.dolfin_version = dolfin_version
    self.plaenterwald_version = plaenterwald_version
  
