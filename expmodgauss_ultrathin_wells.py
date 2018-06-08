import numpy as np

from plaenterwaldexpress.solverinfo import SolverInfo
from plaenterwaldexpress.emasolver import ema_solver_1d
from plaenterwaldexpress.kp6solver import kp6_solver_1d
from plaenterwaldexpress.materials import *
from plaenterwaldexpress.materialprofiles import *
from plaenterwaldexpress.alloygradients import *

from ge_sn_definition_Tonkikh2013 import *


########## Run Calculations ###################################################

### Set well parameter
layer_thickness = [0.1,30.0]
efield = 0.0
alloy_param = [14.0] + [0.28651, 21.74742, 0.76027, 0.31335, 0.53637]
num_ev = 10

### Set numerics parameter
elem_per_nm = 40
elem_per_layer = [int(elem_per_nm * k) for k in layer_thickness]
porder = 1
max_iter = 10000


### Calculate conduction bands

bands = ['G', 'L']

# Set up material models
Ge_meff = dict(Ge_basics.items() + Ge_ema.items())
Sn_meff = dict(Sn_basics.items() + Sn_ema.items())

Ge_elem = ElementaryEMA(Ge_meff)
Sn_elem = ElementaryEMA(Sn_meff)
  
GeSn_binary = Binary(Sn_elem, Ge_elem, GeSn_bowing)

substrate = Ge_elem
barrier = GeSn_binary.make_elementary(0.0029)
well = GeSn_binary
 
for i in range(len(bands)):
  # Barrier material profile 
  [vpot_b, meff_b] = materialEMA_profile_1d(barrier, bands[i], asub=substrate.get_value('a0'))
    
  # Well material profile
  xgrad_w = alloy_gradients_1d('expmodgauss', alloy_param)
  
  [vpot_w, meff_w] = materialEMA_profile_1d(well, bands[i], xgrad_w, asub=substrate.get_value('a0'))
  
  # Assemble well and barrier
  vpot = [vpot_b, vpot_w, vpot_b]
  meff = [meff_b, meff_w, meff_b]
 
  # Set up solver information object
  solver_info = SolverInfo(num_ev)
  solver_info.set_max_iterations(max_iter)
    
  # Run calculation
  solution = ema_solver_1d(layer_thickness, vpot, meff, bands[i], efield, \
                             elem_per_layer, porder, solver_info)
  
  # Print first eigenvalues
  print 'First eigenvalues from ema calculation of band ' + bands[i] + ':'
  for n in range(num_ev):
    try:
      eigval = solution.get_eigenstate(n).get_eigenvalue()
      print str(n) + ': ' + str(eigval)
    except:
      print str(n) + ': not converged'

  # Write solution to file
  fname = 'GeSn_' + bands[i] + '.h5'
  solution.export_solution(fname)

### Calculate valence bands
bands = ['hh', 'lh']
 
# Set up material models
Ge_meff = dict(Ge_basics.items() + Ge_kp6.items())
Sn_meff = dict(Sn_basics.items() + Sn_kp6.items())

Ge_elem = ElementaryKP6(Ge_meff)
Sn_elem = ElementaryKP6(Sn_meff)
  
GeSn_binary = Binary(Sn_elem, Ge_elem, GeSn_bowing)
  
barrier = Ge_elem
well = GeSn_binary
 
for i in range(len(bands)):
  # Set up material profile for barrier
  [vpot_b, meff_b] = materialMEMA_profile_1d(barrier, bands[i])
    
  # Create Fermi profile. Parameters: zstart, thickness
  xgrad_w = alloy_gradients_1d('expmodgauss', alloy_param)

  [vpot_w, meff_w] = materialMEMA_profile_1d(well, bands[i], xgrad_w, asub=barrier.get_value('a0'))
    
  vpot = [vpot_b, vpot_w, vpot_b]
  meff = [meff_b, meff_w, meff_b]
 
  # Set up solver information object
  solver_info = SolverInfo(num_ev)
  solver_info.set_max_iterations(max_iter)
    
  # Run calculation
  solution = ema_solver_1d(layer_thickness, vpot, meff, bands[i], efield, \
                           elem_per_layer, porder, solver_info)
  
  # Print first eigenvalues
  print 'First eigenvalues from ema calculation of band ' + bands[i] + ':'
  for n in range(num_ev):
    try:
      eigval = solution.get_eigenstate(n).get_eigenvalue()
      print str(n) + ': ' + str(eigval)
    except:
      print str(n) + ': not converged'
  
  # Write solution to file
  fname = 'GeSn_' + bands[i] + '.h5'
  solution.export_solution(fname)

