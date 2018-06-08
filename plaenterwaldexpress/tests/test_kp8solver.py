import numpy as np
import csv

from ..solverinfo import SolverInfo
from ..kp8solver import kp8_solver_1d
from ..alloyprofiles import *
from ..materials import *
from ..materialprofiles import *

def algaas_wells(porder=1, elem_per_nm=[1,2,4], num_iterations=1000):
  """Convergence test for a GaAs well on a Al(0.5)Ga(0.5)As barrier."""
  # Define test scenario ###############################################
  
  # Define parallel k-vector
  kvec = [0.00, 0.05]
  
  # Define material parameters (well = GaAs, barrier = Al(0.5)Ga(0.5)As)
  mat01 = {'EG':2.520, 'EL':1., 'EX':1., 'Eav':0.88667, 'delta_so':0.340, \
           'a0':0.5653, 'C11':153.1, 'C12':58.76, 'C44':0., \
           'av':1.972, 'bv':-2.404, 'aG':0.,'aL':0.,'aX':0., 'xiX':0.} 

  mat02 = {'EG':3.880, 'EL':1., 'EX':1., 'Eav':0.75667, 'delta_so':0.280, \
           'a0':0.5731, 'C11':167.5, 'C12':65.0, 'C44':0., \
	   'av':2.46, 'bv':-2.10, 'aG':0.,'aL':0.,'aX':0., 'xiX':0.}
  
  mat01['Ac'] = 2.5600;     mat02['Ac'] = 0.8900
  mat01['P'] = 0.828;       mat02['P'] = 0.828
  mat01['Lp'] = -2.0800;    mat02['Lp'] = -0.2300
  mat01['M'] = -2.6500;     mat02['M'] = -2.0900
  mat01['Np'] = -4.2300;    mat02['Np'] = -1.8000 
  mat01['Npp'] = -0.5800;   mat02['Npp'] = 1.2900 
  mat01['Nm'] = -3.6500;    mat02['Nm'] = -3.0900

  well = ElementaryKP8(mat01)
  barrier = ElementaryKP8(mat02)
  
  [cpot_w, vpot_w, dso_w, kane_w, eps_w, defpot_w] = \
    materialKP8_profile_1d(well)
  [cpot_b, vpot_b, dso_b, kane_b, eps_b, defpot_b] = \
    materialKP8_profile_1d(barrier)

  # Define quantum well properties
  layer_thickness = [8.0, 7.0, 8.0]
 
  cpot = [cpot_b, cpot_w, cpot_b]
  vpot = [vpot_b, vpot_w, vpot_b]
  delta_so = [dso_b, dso_w, dso_b]

  kane = [kane_b, kane_w, kane_b]
  defpot = [defpot_b, defpot_w, defpot_b]
  epsilon = None
  efield = None

  # Setting up solver 
  solver_info = SolverInfo()
  solver_info.set_max_iterations(num_iterations)
  solver_info.set_spectrum('target magnitude')
  solver_info.set_spectral_shift(1.9)
  solver_info.set_spectral_transform('shift-and-invert')
  
  # Run test calculation  #####################################################
  
  print 'TEST: k.p solver with GaAs as well and Al(0.5)Ga(0.5)As as barrier.' 
  print 'Convergence test started...'
  
  eigval = []
  runtime = []
  num_iter = []
  num_eigenstates = []
  
  for i in range(len(elem_per_nm)):
    elem_per_layer = [int(elem_per_nm[i] * x) for x in layer_thickness]
  
    eigenpair = kp8_solver_1d(layer_thickness, cpot, vpot, \
                              kane, delta_so, \
                              defpot, epsilon, \
                              kvec, efield, \
                              elem_per_layer, porder, \
                              solver_info)
  
    # Get first 6 eigenvalues
    ev_current_run = []
    for k in range(1,32,4):
      try:
        nth_eigenstate = eigenpair.get_eigenstate(k)
        ev_current_run.append(nth_eigenstate.get_eigenvalue())
      except:
        print 'Too few eigenvalues have converged.'
  
    eigval.append(ev_current_run)
  
    # Save runtime, number iterations, and number eigenvalues as well
    runtime.append(solver_info.get_runtime())
    num_iter.append(solver_info.get_required_iterations())
    num_eigenstates.append(eigenpair.get_num_eigenstates())
  
  
  # OUTPUT: Show results ###############################################
  print 'Showing results...'
  
  # Compare numerical values with nextnano3 values. 
  test_comparison = [2.5918, 0.9866, \
                     2.7792, 0.9718, \
                     0.9491, 0.8935, \
		     0.8897, -1.0]
  
  # Reminder for test size
  print 'Elements per nm: ' + str(elem_per_nm)
  
  for n in range(len(eigval[0])):
    for i in range(len(elem_per_nm)):
      try:
        if i == 0:
          print str(n + 1) + ': ' + str(eigval[0][n])
        else:
          print '   ' + str(eigval[i][n])
      except:
        print '   not converged'
    print '  --------------------'
    print '   ' + str(test_comparison[n])
    print ' '
  
  # Show general benchmarks
  print 'time per ev:'
  print [x[0] / x[1] for x in zip(runtime, num_eigenstates)]
  print 'iterations per ev:'
  print [1.0 * x[0] / x[1] for x in zip(num_iter, num_eigenstates)]
  
  # Show general benchmarks
  print 'time per ev:'
  print [x[0] / x[1] for x in zip(runtime, num_eigenstates)]
  print 'iterations per ev:' 
  print [1.0 * x[0] / x[1] for x in zip(num_iter, num_eigenstates)]
  
def algaas_strain_wells(porder=1, elem_per_nm=[1,2,4], num_iterations=1000):
  """Convergence test for a GaAs well on a Al(0.5)Ga(0.5)As barrier (strain included)."""
  
  # Define test scenario ###############################################
 
  # Define parallel k-vector
  kvec = [0.00, 0.00]
  
  # Define material parameters (well = GaAs, barrier = Al(0.5)Ga(0.5)As)
  mat01 = {'EG':2.520, 'EL':1., 'EX':1., 'Eav':0.88667, 'delta_so':0.340, \
           'a0':0.5653, 'C11':123.15, 'C12':43.2650, 'C44':0., \
           'av':-0.3785, 'bv':-2.8175, 'aG':-10.3145,'aL':0.,'aX':0., 'xiX':0.} 

  mat02 = {'EG':3.880, 'EL':1., 'EX':1., 'Eav':0.75667, 'delta_so':0.280, \
           'a0':0.5731, 'C11':123.15, 'C12':43.2650, 'C44':0., \
	   'av':-0.3785, 'bv':-2.8175, 'aG':-10.3145,'aL':0.,'aX':0., 'xiX':0.}
  
  mat01['Ac'] = 2.5600;     mat02['Ac'] = 0.8900
  mat01['P'] = 0.828;       mat02['P'] = 0.828
  mat01['Lp'] = -2.0800;    mat02['Lp'] = -0.2300
  mat01['M'] = -2.6500;     mat02['M'] = -2.0900
  mat01['Np'] = -4.2300;    mat02['Np'] = -1.8000 
  mat01['Npp'] = -0.5800;   mat02['Npp'] = 1.2900 
  mat01['Nm'] = -3.6500;    mat02['Nm'] = -3.0900

  a0_sub = 0.56579

  well = ElementaryKP8(mat01)
  barrier = ElementaryKP8(mat02)

  [cpot_w, vpot_w, dso_w, kane_w, eps_w, defpot_w] = \
    materialKP8_profile_1d(well, asub=a0_sub)
  [cpot_b, vpot_b, dso_b, kane_b, eps_b, defpot_b] = \
    materialKP8_profile_1d(barrier, asub=a0_sub)

  # Define quantum well properties
  layer_thickness = [8.0, 7.0, 8.0]
 
  cpot = [cpot_b, cpot_w, cpot_b]
  vpot = [vpot_b, vpot_w, vpot_b]
  delta_so = [dso_b, dso_w, dso_b]

  kane = [kane_b, kane_w, kane_b]
  defpot = [defpot_b, defpot_w, defpot_b]
  epsilon = [eps_b, eps_w, eps_b]
  efield = None

  # Setting up solver 
  solver_info = SolverInfo()
  solver_info.set_max_iterations(num_iterations)
  solver_info.set_spectrum('target magnitude')
  solver_info.set_spectral_shift(1.9)
  solver_info.set_spectral_transform('shift-and-invert')
  
  # Run test calculation  #####################################################
  
  print 'TEST: k.p solver with GaAs as well and Al(0.5)Ga(0.5)As as barrier.' 
  print 'Convergence test started...'
  
  eigval = []
  runtime = []
  num_iter = []
  num_eigenstates = []
  
  for i in range(len(elem_per_nm)):
    elem_per_layer = [int(elem_per_nm[i] * x) for x in layer_thickness]
  
    eigenpair = kp8_solver_1d(layer_thickness, cpot, vpot, \
                              kane, delta_so, \
                              defpot, epsilon, \
                              kvec, efield, \
                              elem_per_layer, porder, \
                              solver_info)
  
    # Get first 6 eigenvalues
    ev_current_run = []
    for k in range(1,32,4):
      try:
        nth_eigenstate = eigenpair.get_eigenstate(k)
        ev_current_run.append(nth_eigenstate.get_eigenvalue())
      except:
        print 'Too few eigenvalues have converged.'
  
    eigval.append(ev_current_run)
  
    # Save runtime, number iterations, and number eigenvalues as well
    runtime.append(solver_info.get_runtime())
    num_iter.append(solver_info.get_required_iterations())
    num_eigenstates.append(eigenpair.get_num_eigenstates())
  
  
  # OUTPUT: Show results ###############################################
  print 'Showing results...'
  
  # Compare numerical values with nextnano3 results 
  test_comparison = [2.5803, 0.9846, \
                     2.7712, 0.9721, \
                     0.9534, 0.9166, \
		     -1.0, -1.0]
  
  # Reminder for test size
  print 'Elements per nm: ' + str(elem_per_nm)
  
  for n in range(len(eigval[0])):
    for i in range(len(elem_per_nm)):
      try:
        if i == 0:
          print str(n + 1) + ': ' + str(eigval[0][n])
        else:
          print '   ' + str(eigval[i][n])
      except:
        print '   not converged'
    print '  --------------------'
    print '   ' + str(test_comparison[n])
    print ' '
  
  # Show general benchmarks
  print 'time per ev:'
  print [x[0] / x[1] for x in zip(runtime, num_eigenstates)]
  print 'iterations per ev:'
  print [1.0 * x[0] / x[1] for x in zip(num_iter, num_eigenstates)]
  
  # Show general benchmarks
  print 'time per ev:'
  print [x[0] / x[1] for x in zip(runtime, num_eigenstates)]
  print 'iterations per ev:' 
  print [1.0 * x[0] / x[1] for x in zip(num_iter, num_eigenstates)]

  eigenpair.export_solution('strain.h5')

def algaas_tilted_wells(porder=1, elem_per_nm=[1,2,4], num_iterations=1000):
  """Convergence test for a Al(x)Ga(1-x)As well (tilted) on Al(0.5)Ga(0.5)As barrier."""
  
  # Define test scenario ###############################################
  
  # Define parallel k-vector
  kvec = [0.00, 0.00]
  
  # Define material parameters (well = GaAs, barrier = Al(0.5)Ga(0.5)As)
  mat01 = {'EG':2.520, 'EL':1., 'EX':1., 'Eav':0.88667, 'delta_so':0.340, \
           'a0':0.5653, 'C11':153.1, 'C12':58.76, 'C44':0., \
           'av':1.972, 'bv':-2.404, 'aG':0.,'aL':0.,'aX':0., 'xiX':0.} 

  mat02 = {'EG':3.880, 'EL':1., 'EX':1., 'Eav':0.75667, 'delta_so':0.280, \
           'a0':0.5731, 'C11':167.5, 'C12':65.0, 'C44':0., \
	   'av':2.46, 'bv':-2.10, 'aG':0.,'aL':0.,'aX':0., 'xiX':0.}
  
  mat01['Ac'] = 2.5600;     mat02['Ac'] = 0.8900
  mat01['P'] = 0.828;       mat02['P'] = 0.828
  mat01['Lp'] = -2.0800;    mat02['Lp'] = -0.2300
  mat01['M'] = -2.6500;     mat02['M'] = -2.0900
  mat01['Np'] = -4.2300;    mat02['Np'] = -1.8000 
  mat01['Npp'] = -0.5800;   mat02['Npp'] = 1.2900 
  mat01['Nm'] = -3.6500;    mat02['Nm'] = -3.0900
  
  bowing = {'EG':0.305, 'EX': 0.055}
 
  elem01 = ElementaryKP8(mat01)
  elem02 = ElementaryKP8(mat02)
  binary0102 = Binary(elem01, elem02, bowing)

  lin_grad = alloy_gradients_1d('linear', (8.0, 7.0, 0., 0.5)) 

  [cpot_w, vpot_w, dso_w, kane_w, eps_w, defpot_w] = \
    materialKP8_profile_1d(binary0102, lin_grad)
  [cpot_b, vpot_b, dso_b, kane_b, eps_b, defpot_b] = \
    materialKP8_profile_1d(elem02)

  # Define quantum well properties
  layer_thickness = [8.0, 7.0, 8.0]
 
  cpot = [cpot_b, cpot_w, cpot_b]
  vpot = [vpot_b, vpot_w, vpot_b]
  delta_so = [dso_b, dso_w, dso_b]

  kane = [kane_b, kane_w, kane_b]
  defpot = [defpot_b, defpot_w, defpot_b]
  epsilon = None
  efield = None

  # Setting up solver 
  solver_info = SolverInfo()
  solver_info.set_max_iterations(num_iterations)
  solver_info.set_spectrum('target magnitude')
  solver_info.set_spectral_shift(2.5)
  solver_info.set_spectral_transform('shift-and-invert')
  
  # Run test calculation  #####################################################
  
  print 'TEST: k.p solver with GaAs as well and Al(0.5)Ga(0.5)As as barrier.' 
  print 'Convergence test started...'
  
  eigval = []
  runtime = []
  num_iter = []
  num_eigenstates = []
  
  for i in range(len(elem_per_nm)):
    elem_per_layer = [int(elem_per_nm[i] * x) for x in layer_thickness]
  
    eigenpair = kp8_solver_1d(layer_thickness, cpot, vpot, \
                              kane, delta_so, \
                              defpot, epsilon, \
                              kvec, efield, \
                              elem_per_layer, porder, \
                              solver_info)
  
    # Get first 6 eigenvalues
    ev_current_run = []
    for k in range(1,32,4):
      try:
        nth_eigenstate = eigenpair.get_eigenstate(k)
        ev_current_run.append(nth_eigenstate.get_eigenvalue())
      except:
        print 'Too few eigenvalues have converged.'
  
    eigval.append(ev_current_run)
  
    # Save runtime, number iterations, and number eigenvalues as well
    runtime.append(solver_info.get_runtime())
    num_iter.append(solver_info.get_required_iterations())
    num_eigenstates.append(eigenpair.get_num_eigenstates())
  
  
  # OUTPUT: Show results ###############################################
  print 'Showing results...'
  
  # Compare numerical values from nextnano3 calculation 
  test_comparison = [3.3973, 0.8877, \
                     3.6445, 0.8786, \
                     0.8543, 0.8481, \
		     0.8474, -1.0]
  
  # Reminder for test size
  print 'Elements per nm: ' + str(elem_per_nm)
  
  for n in range(len(eigval[0])):
    for i in range(len(elem_per_nm)):
      try:
        if i == 0:
          print str(n + 1) + ': ' + str(eigval[0][n])
        else:
          print '   ' + str(eigval[i][n])
      except:
        print '   not converged'
    print '  --------------------'
    print '   ' + str(test_comparison[n])
    print ' '
  
  # Show general benchmarks
  print 'time per ev:'
  print [x[0] / x[1] for x in zip(runtime, num_eigenstates)]
  print 'iterations per ev:'
  print [1.0 * x[0] / x[1] for x in zip(num_iter, num_eigenstates)]
  
  # Show general benchmarks
  print 'time per ev:'
  print [x[0] / x[1] for x in zip(runtime, num_eigenstates)]
  print 'iterations per ev:' 
  print [1.0 * x[0] / x[1] for x in zip(num_iter, num_eigenstates)]

