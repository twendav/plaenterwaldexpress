import numpy as np
import csv

from ..solverinfo import SolverInfo
from ..emasolver import ema_solver_1d
from ..materialprofiles import *
from ..alloygradients import *
from ..materials import *

def matlab_simple_qw(porder=1, elem_per_nm=[1,2,4,8,16], num_iterations=200):
  r"""Convergence test for a simple QW structure."""
  # Define test scenario ######################################################

  # Define material parameters
  mat01 = {'EG':0., 'EL':1., 'EX':1., 'Eav':-0.1, 'delta_so':0.3, \
           'a0':0.5, 'C11':0., 'C12':0., 'C44':0., \
           'av':0., 'bv':0., 'aG':0.,'aL':0.,'aX':0., 'xiX':0., \
           'meffG':0.7, 'meffL':0.7, 'meffX':0.7, \
           'meffhh':0.7, 'mefflh':0.7, 'meffso':0.7}

  mat02 = {'EG':1., 'EL':1., 'EX':1., 'Eav':-1.1, 'delta_so':0.3, \
           'a0':0.7, 'C11':0., 'C12':0., 'C44':0., \
	   'av':0., 'bv':0., 'aG':0.,'aL':0.,'aX':0., 'xiX':0., \
	   'meffG':0.5, 'meffL':0.5, 'meffX':0.5, \
	   'meffhh':0.5, 'mefflh':0.5, 'meffso':0.5}
  
  # Set up well and barrier potentials
  well = ElementaryEMA(mat01)
  barrier = ElementaryEMA(mat02)

  # Define quantum well properties
  layer_thickness = [8.0, 1.0, 8.0]                  # [nm]
  efield = 0.0                                       # [V/nm]
  
  # Set up solver
  solver_info = SolverInfo(3)
  solver_info.set_max_iterations(num_iterations)

  # Define precomputed reference values #######################################

  # From: semianalytic results for a finite potential well 
  comparison = [[0.264332962309381, 0.870265536700407, -1.0], \
                [-0.264332962309381, -0.870265536700407, 1.0]]
  
  # Run test calculation ######################################################
  valley_types = ['G', 'hh']
  
  print 'TEST: Calculate EV of simple of QW with different effective masses.'
  print 'Computation started ...'
  
  computed = []
  runtime = []
  num_iter = []
  num_eigenstates = []
  
  for i in range(len(elem_per_nm)):
    for vid in range(len(valley_types)):
      [vpot_w, meff_w] = materialEMA_profile_1d(well, valley_types[vid])
      [vpot_b, meff_b] = materialEMA_profile_1d(barrier, valley_types[vid])
      vpot = [vpot_b, vpot_w, vpot_b]
      meff = [meff_b, meff_w, meff_b]

      elem_per_layer = [int(elem_per_nm[i] * x) for x in layer_thickness]
  
      solution = ema_solver_1d(layer_thickness, vpot, meff, \
                               valley_types[vid], efield, elem_per_layer, \
			       porder, solver_info)
 
      # Get first 3 eigenvalues
      ev_current_run = []
      for k in range(len(comparison[vid])):
        try:
          nth_eigenstate = solution.get_eigenstate(k)
          ev_current_run.append(nth_eigenstate.get_eigenvalue())
        except:
          print 'Too few eigenvalues have converged.'
  
      computed.append(ev_current_run)
  
    # Save runtime, number iterations, and number eigenvalues as well
    runtime.append(solver_info.get_runtime())
    num_iter.append(solver_info.get_required_iterations())
    num_eigenstates.append(solution.get_num_eigenstates())
  
  
  # OUTPUT: Show results ###############################################
  print 'Done. Showing results...'
  
  # Reminder for test size
  print 'Elements per nm: ' + str(elem_per_nm)
 
  for vid in range(len(valley_types)):
    print 'Valley type: ' + valley_types[vid]

    for evid in range(len(comparison[vid])):
      for rid in range(len(elem_per_nm)):
        try:
          print '   ' + str(computed[(vid + 2 * rid)][evid])
        except:
          print '   not converged'
      print '  --------------------'
      print '   ' + str(comparison[vid][evid])
      print ' '
  
    # Show general benchmarks
    print 'time per ev:'
    print [x[0] / x[1] for x in zip(runtime[vid::2], num_eigenstates[vid::2])]
    print 'iterations per ev:'
    print [1.0 * x[0] / x[1] for x in zip(num_iter[vid::2], num_eigenstates[vid::2])]

def matlab_efield_qw(porder=1, elem_per_nm=[1,2,4,8,16], num_iterations=200):
  r"""Convergence test for a simple QW structure with applied efield."""
  # Define test scenario #############################################
  
  # Define material parameters
  mat01 = {'EG':0., 'EL':0., 'EX':0., 'Eav':0., 'delta_so':-0.1, \
           'a0':0.5, 'C11':0., 'C12':0., 'C44':0., \
           'av':0., 'bv':0., 'aG':0.,'aL':0.,'aX':0., 'xiX':0., \
           'meffG':1., 'meffL':1., 'meffX':1., \
           'meffhh':1., 'mefflh':1., 'meffso':1.}

  mat02 = {'EG':1., 'EL':0., 'EX':0., 'Eav':0., 'delta_so':-0.1, \
           'a0':0.7, 'C11':0., 'C12':0., 'C44':0., \
           'av':0., 'bv':0., 'aG':0.,'aL':0.,'aX':0., 'xiX':0., \
           'meffG':1., 'meffL':1., 'meffX':1., \
           'meffhh':1., 'mefflh':1., 'meffso':1.}

  # Set up potential functions
  barrier = ElementaryEMA(mat02)
  well = ElementaryEMA(mat01)

  [vpot_b, meff_b] = materialEMA_profile_1d(barrier, 'G')
  [vpot_w, meff_w] = materialEMA_profile_1d(well, 'G')

  # Define quantum well properties
  layer_thickness = [3.0, 1.0, 3.0]                  # [nm]
  vpot = [vpot_b, vpot_w, vpot_b]
  meff = [meff_b, meff_w, meff_b]
  efield = 0.1                                       # [V/nm]
 
  # Set up solver  
  solver_info = SolverInfo()
  solver_info.set_max_iterations(num_iterations)
  
  # Define precomputed reference values #######################################

  # From TMM calculation for an applied electric field.  
  test_comparison = [0.189320449, 0.70060, -1.0]

  # Run test calculation ######################################################
  
  print 'TEST: Effective mass solver with an applied electric field.' 
  print 'Computation started...'
  
  eigval = []
  runtime = []
  num_iter = []
  num_eigenstates = []
  
  for i in range(len(elem_per_nm)):
    elem_per_layer = [int(elem_per_nm[i] * x) for x in layer_thickness]
  
    eigenpair = ema_solver_1d(layer_thickness, vpot, meff, \
                              'G', efield, elem_per_layer, porder, \
  			      solver_info)
  
    # Get first 3 eigenvalues
    ev_current_run = []
    for k in range(3):
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
  print 'Done. Showing results...'
  
 
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
  
def nextnano_tilted_qw(porder=1, elem_per_nm=[1,2,4,8,16], num_iterations=200):
  r"""Convergence test for a QW structure with tilted well."""
  # Define test scenario #############################################
  
  # Define material parameters
  mat01 = {'EG':0.5, 'EL':0., 'EX':0., 'Eav':0., 'delta_so':-0.1, \
           'a0':0.5, 'C11':0., 'C12':0., 'C44':0., \
           'av':0., 'bv':0., 'aG':0.,'aL':0.,'aX':0., 'xiX':0., \
           'meffG':1., 'meffL':(1.,1.), 'meffX':(1.,1.), \
           'meffhh':1., 'mefflh':1., 'meffso':1.}

  mat02 = {'EG':1., 'EL':0., 'EX':0., 'Eav':0., 'delta_so':-0.1, \
           'a0':0.5, 'C11':0., 'C12':0., 'C44':0., \
           'av':0., 'bv':0., 'aG':0.,'aL':0.,'aX':0., 'xiX':0., \
           'meffG':1., 'meffL':(1.,1.), 'meffX':(1.,1.), \
           'meffhh':1., 'mefflh':1., 'meffso':1.}

  # Set up alloy gradient
  linear_gradient = alloy_gradients_1d('linear', (8.0, 1.0, 0., 1.))

  # Set up potential functions
  elem01 = ElementaryEMA(mat01)
  elem02 = ElementaryEMA(mat02)
  binary0102 = Binary(elem02, elem01)

  [vpot_b, meff_b] = materialEMA_profile_1d(elem02, 'G')
  [vpot_w, meff_w] = materialEMA_profile_1d(binary0102, 'G', linear_gradient)

  # Define quantum well properties
  layer_thickness = [8.0, 1.0, 8.0]                  # [nm]
  vpot = [vpot_b, vpot_w, vpot_b]
  meff = [meff_b, meff_w, meff_b]
  efield = 0.                                       # [V/nm]
  
  # Set up solver  
  solver_info = SolverInfo()
  solver_info.set_max_iterations(num_iterations)

  # Define precomputed reference values #######################################

  # From nextnano3 calculations 
  test_comparison = [0.8357, -1.0]
   
  # Run test calculation ######################################################
  
  print 'TEST: Effective mass solver with a potential.' 
  print 'Computation started...'
  
  eigval = []
  runtime = []
  num_iter = []
  num_eigenstates = []
  
  for i in range(len(elem_per_nm)):
    elem_per_layer = [int(elem_per_nm[i] * x) for x in layer_thickness]
  
    eigenpair = ema_solver_1d(layer_thickness, vpot, meff, \
                              'G', efield, elem_per_layer, porder, \
                              solver_info)
  
    # Get first 3 eigenvalues
    ev_current_run = []
    for k in range(len(test_comparison)):
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
  print 'Done. Showing results...'
  
  # Reminder for test size
  print 'Elements per nm: ' + str(elem_per_nm)
  
  for n in range(len(eigval[0])):
    for i in range(len(elem_per_nm)):
      try:
        if i == 0:
          print str(n + 1) + ': ' + str(eigval[0][n] - \
  	                              layer_thickness[0] * \
  				      efield)
        else:
          print '   ' + str(eigval[i][n] - \
  	                  layer_thickness[0] * \
  			  efield)
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

def modexpgauss_qw(porder=1, elem_per_nm=[1,2,4,8,16], num_iterations=200):
  r"""Convergence test for a QW structure with tilted well."""
  # Define test scenario #############################################
  
  # Define material parameters
  mat01 = {'EG':0.5, 'EL':0., 'EX':0., 'Eav':0., 'delta_so':-0.1, \
           'a0':0.5, 'C11':0., 'C12':0., 'C44':0., \
           'av':0., 'bv':0., 'aG':0.,'aL':0.,'aX':0., 'xiX':0., \
           'meffG':1., 'meffL':(1.,1.), 'meffX':(1.,1.), \
           'meffhh':1., 'mefflh':1., 'meffso':1.}

  mat02 = {'EG':1., 'EL':0., 'EX':0., 'Eav':0., 'delta_so':-0.1, \
           'a0':0.5, 'C11':0., 'C12':0., 'C44':0., \
           'av':0., 'bv':0., 'aG':0.,'aL':0.,'aX':0., 'xiX':0., \
           'meffG':1., 'meffL':(1.,1.), 'meffX':(1.,1.), \
           'meffhh':1., 'mefflh':1., 'meffso':1.}

  # Set up alloy gradient (y0, A, xc, w, t0)
  linear_gradient = alloy_gradients_1d('expmodgauss', (8.0, 0.28651, 21.74742, 0.76027, 0.31335, 0.53670))

  # Set up potential functions
  elem01 = ElementaryEMA(mat01)
  elem02 = ElementaryEMA(mat02)
  binary0102 = Binary(elem02, elem01)

  [vpot_b, meff_b] = materialEMA_profile_1d(elem02, 'G')
  [vpot_w, meff_w] = materialEMA_profile_1d(binary0102, 'G', linear_gradient)

  # Define quantum well properties
  layer_thickness = [8.0, 10.0, 8.0]                  # [nm]
  vpot = [vpot_b, vpot_w, vpot_b]
  meff = [meff_b, meff_w, meff_b]
  efield = 0.                                         # [V/nm]
  
  # Set up solver  
  solver_info = SolverInfo()
  solver_info.set_max_iterations(num_iterations)

  # Define precomputed reference values #######################################

  # TODO: Test calculations still have to be performed! 
  test_comparison = [-1.0, -1.0]
   
  # Run test calculation ######################################################
  
  print 'TEST: Effective mass solver with a exp. mod. gaussian potential.' 
  print 'Computation started...'
  
  eigval = []
  runtime = []
  num_iter = []
  num_eigenstates = []
  
  for i in range(len(elem_per_nm)):
    elem_per_layer = [int(elem_per_nm[i] * x) for x in layer_thickness]
  
    eigenpair = ema_solver_1d(layer_thickness, vpot, meff, \
                              'G', efield, elem_per_layer, porder, \
                              solver_info)
  
    # Get first 3 eigenvalues
    ev_current_run = []
    for k in range(len(test_comparison)):
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
  print 'Done. Showing results...'
  
  # Reminder for test size
  print 'Elements per nm: ' + str(elem_per_nm)
  
  for n in range(len(eigval[0])):
    for i in range(len(elem_per_nm)):
      try:
        if i == 0:
          print str(n + 1) + ': ' + str(eigval[0][n] - \
  	                              layer_thickness[0] * \
  				      efield)
        else:
          print '   ' + str(eigval[i][n] - \
  	                  layer_thickness[0] * \
  			  efield)
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
 

def nextnano_sige_qw():
  r"""Calculates and compares eigenenergies of a SiGe/SiGe QW structure."""

  # Define test scenario ######################################################

  # Define material parameters (in accordence with nextnano3 parameters)
  Si = {'EG':4.49, 'EL':3.11466, 'EX':2.270, 'Eav':1.090, 'delta_so':0.044, \
              'a0':0.54304, 'C11':165.77, 'C12':63.93, 'C44':0., \
	      'av':2.05, 'bv':-2.10, 'aG':-10.39,'aL':-2.02,'aX':3.40, 'xiX':9.16, \
	      'meffG':0.156, 'meffL':(1.420,0.130), 'meffX':(0.916,0.190), \
	      'meffhh':0.537, 'mefflh':0.153, 'meffso':0.234}

  Ge = {'EG':2.6663, 'EL':2.5063, 'EX':2.6973, 'Eav':1.67, 'delta_so':0.289, \
              'a0':0.5658, 'C11':128.53, 'C12':48.26, 'C44':0., \
	      'av':-0.35, 'bv':-2.86, 'aG':-10.41,'aL':-4.35,'aX':0.14, 'xiX':9.42, \
	      'meffG':0.038, 'meffL':(1.57,0.0807), 'meffX':(1.350,0.290), \
	      'meffhh':0.316, 'mefflh':0.0424, 'meffso':0.095}
  
  SiGe_bowing = {'EX': 0.206}

  # Define quantum well structure
  SiGe = Binary(ElementaryEMA(Ge), ElementaryEMA(Si), SiGe_bowing)

  substrate = SiGe.make_elementary(0.90)
  barrier = SiGe.make_elementary(0.85)
  well = SiGe.make_elementary(1.00)

  layer_thickness = [10.0, 10.0, 10.0]  # in nm
  efield = 0.0                          # in V/nm

  # Set numerical parameters ##################################################
  porder = 1
  elem_per_nm = 10   # Remark: Convergence roughly up to 10^-4eV
  elem_per_layer = [int(w * elem_per_nm) for w in layer_thickness]
  num_iterations = 100

  solver_info = SolverInfo()
  solver_info.set_max_iterations(num_iterations)
  solver_info.set_spectrum('target magnitude')
  solver_info.set_spectral_transform('shift-and-invert')
  shift_and_invert = [2.710, 2.520, 2.610, 2.610, 2.580, 1.790, 1.751, 1.480]

  # Define precomputed values (nextnano) as comparison for results ############
  comparison = [[2.7598, 2.8737], \
                [2.5422, 2.5780], \
		[0], \
		[0], \
		[0], \
		[1.7799, 1.7555, 1.7168, 1.6693], \
		[1.7248], \
		[1.4602, 1.4202]]

  # Run test calculation ######################################################
  
  print 'TEST: Calculate potentials of Ge/Si QW well and barrier material.'
  print 'Computation started ...'
  
  valley_types = ['G','L','X100','X010','X001','hh','lh','so']

  eigenvalues = []

  for i in range(len(valley_types)):
    # Set up potentials and effective masses for each layer
    [vpot_w, meff_w] = materialEMA_profile_1d(well, valley_types[i], \
                                             asub=substrate.get_value('a0'))
    [vpot_b, meff_b] = materialEMA_profile_1d(barrier, valley_types[i], \
                                             asub=substrate.get_value('a0'))
    vpot = [vpot_b, vpot_w, vpot_b]
    meff = [meff_b, meff_w, meff_b]
    
    # Run calculation
    solver_info.set_spectral_shift(shift_and_invert[i])
    solution = ema_solver_1d(layer_thickness, vpot, meff, \
                             valley_types[i], efield, elem_per_layer, porder, \
	                     solver_info)

    # Extract eigenvalues
    ev_to_be_saved = []
    for n in range(len(comparison[i])):
      ev_to_be_saved.append(solution.get_eigenstate(n).get_eigenvalue())

    eigenvalues.append(ev_to_be_saved)

  # Output results and compare with precomputed values ########################

  print 'Done. Showing results ...'

  for i in range(len(valley_types)):
    print 'Valley: ' + valley_types[i]
    print '            lowest eigenvalues:'
    print 'Computed    ' + str(eigenvalues[i])
    print '             -----------------------------'
    print 'Comparison: ' + str(comparison[i])


def matlab_sige_qw():
  r"""Calculates and compares potentials of a strained SiGe/SiGe QW."""

  # Define test scenario ######################################################

  # Define material parameters (in accordence with nextnano3 parameters)
  Si = {'EG':4.49, 'EL':3.11466, 'EX':2.270, 'Eav':1.090, 'delta_so':0.044, \
        'a0':0.54304, 'C11':165.77, 'C12':63.93, 'C44':0., \
        'av':2.05, 'bv':-2.10, 'aG':-10.39,'aL':-2.02,'aX':3.40, 'xiX':9.16, \
        'meffG':0.156, 'meffL':(1.420,0.130), 'meffX':(0.916,0.190), \
        'meffhh':0.537, 'mefflh':0.153, 'meffso':0.234}

  Ge = {'EG':2.6663, 'EL':2.5063, 'EX':2.6973, 'Eav':1.67, 'delta_so':0.289, \
        'a0':0.5658, 'C11':128.53, 'C12':48.26, 'C44':0., \
        'av':-0.35, 'bv':-2.86, 'aG':-10.41,'aL':-4.35,'aX':0.14, 'xiX':9.42, \
        'meffG':0.038, 'meffL':(1.57,0.0807), 'meffX':(1.350,0.290), \
        'meffhh':0.316, 'mefflh':0.0424, 'meffso':0.095}

  SiGe_bowing = {'EX': 0.206}

  # Define quantum well structure
  SiGe = Binary(ElementaryEMA(Ge), ElementaryEMA(Si), SiGe_bowing)
  substrate = SiGe.make_elementary(0.90)
  barrier = SiGe.make_elementary(0.85)
  well = SiGe.make_elementary(1.00)

  layer_thickness = [10.0, 10.0, 10.0]  # in nm
  efield = 0.0                          # in V/nm

  # Set numerical parameters ##################################################
  porder = 1
  elem_per_nm = 10   # Remark: Convergence roughly up to 10^-4eV
  elem_per_layer = [int(w * elem_per_nm) for w in layer_thickness]
  num_iterations = 100
  solver_info = SolverInfo()
  solver_info.set_max_iterations(num_iterations)
  solver_info.set_spectrum('target magnitude')
  solver_info.set_spectral_transform('shift-and-invert')
  shift_and_invert = [2.710, 2.520, 2.610, 2.610, 2.580, 1.790, 1.751, 1.480]

  # Define precomputed values (nextnano) as comparison for results ############

  # Precomputed values (MATLAB)
  comparison = [[2.91363, 2.71860], \
                [2.58747, 2.52816], \
		[2.61963, 2.67448], \
		[2.61963, 2.67448], \
		[2.58632, 2.74083], \
                [1.65736, 1.78824], \
                [1.67758, 1.75093], \
                [1.41413, 1.47611]]      
  
  # Run test calculation ######################################################
  valley_types = ['G','L','X100','X010','X001','hh','lh','so']

  print 'TEST: Calculate potentials of Ge/Si QW well and barrier material.'
  print 'Computation started ...'
  
  potentials = np.zeros([len(valley_types),2])

  for i in range(len(valley_types)):
    [vpot_w, meff_w] = materialEMA_profile_1d(well, valley_types[i], \
                                             asub=substrate.get_value('a0'))
    [vpot_b, meff_b] = materialEMA_profile_1d(barrier, valley_types[i], \
                                             asub=substrate.get_value('a0'))
    vpot = [vpot_b, vpot_w, vpot_b]
    meff = [meff_b, meff_w, meff_b]
  
    solution = ema_solver_1d(layer_thickness, vpot, meff, \
                             valley_types[i], efield, elem_per_layer, porder, \
  			     solver_info)
    
    # Extract potential information
    potential = solution.get_potential()

    barrier_vpot = potential.values[0]
    potentials[i][0] = barrier_vpot
   
    well_vpot = potential.values[int(sum(layer_thickness) / 2) * elem_per_nm]
    potentials[i][1] = well_vpot
    
  
  # Output results and compare with precomputed values ########################

  print 'Done. Showing results ...'
  
  for i in range(len(valley_types)):
    print 'Valley: ' + valley_types[i]
    print '               barrier           well     '
    print 'Computed    ' + str(potentials[i][0]) + ' ' + str(potentials[i][1])
    print '             -----------------------------'
    print 'Comparison: ' + str(comparison[i][0]) + ' ' + str(comparison[i][1])

def mod_sige_wells(porder=1, elem_per_nm=[1,2,4], num_iterations=1000):
  r"""Convergence test for Si/Ge well on Si barrier material including strain."""
  # Define test scenario ###############################################
  
  # Define material parameters (well = Si(0.6)Ge(0.4), barrier = Si)
  # grown on a Si(0.7)Ge(0.3) substrate
  # Parameter according to Harrison, Quantum Wells, Wires, and Dots
  # (pg. 350)
  mat01 = {'EG':0., 'EL':1., 'EX':1., 'Eav':-0.047467, 'delta_so':0.1424, \
           'a0':0.55214, 'C11':153.1, 'C12':58.76, 'C44':0., \
           'av':1.972, 'bv':-2.404, 'aG':0.,'aL':0.,'aX':0., 'xiX':0.} 

  mat02 = {'EG':1., 'EL':1., 'EX':1., 'Eav':-0.574667, 'delta_so':0.044, \
           'a0':0.5431, 'C11':167.5, 'C12':65.0, 'C44':0., \
	   'av':2.46, 'bv':-2.10, 'aG':0.,'aL':0.,'aX':0., 'xiX':0.}
  
  # Conversion to LK parameter set (for 6-band k.p)
  gamma_well = (7.892, 1.934, 3.14)
  [L_w, M_w, N_w, Np_w, Nm_w] = gamma2LK(gamma_well)
  gamma_barrier = (4.22, 0.39, 1.44)
  [L_b, M_b, N_b, Np_b, Nm_b] = gamma2LK(gamma_barrier) 

  mat01['L'] = L_w; mat02['L'] = L_b 
  mat01['M'] = M_w; mat02['M'] = M_b 
  mat01['N'] = N_w; mat02['N'] = N_b 
  mat01['Np'] = Np_w; mat02['Np'] = Np_b 
  mat01['Nm'] = Nm_w; mat02['Nm'] = Nm_b 

  well = ElementaryKP6(mat01)
  barrier = ElementaryKP6(mat02)
  a0_sub = 0.54988
  
   
  # Define quantum well properties
  layer_thickness = [16.6, 8.3, 16.6]
  efield = 0.0

  # Setting up solver 
  solver_info = SolverInfo()
  solver_info.set_max_iterations(num_iterations)
  
  # Define precomputed reference values #######################################

  # From: semianalytic results for a finite potential well 
  comparison = [[0.009609845262025, 0.060602841965926, 0.144924950957342, \
                 0.261127085199461, 0.405422389353317, 0.558944994397931], \
                [0.053367154764179, 0.143363162169284, 0.284220121750633, \
		 0.436881567311465], \
		[0.184958777237626, 0.270426023411146, 0.414903657722452, \
                 0.596153799148395]]
  
  # Compare numerical values with exact numerics TMM results (no re- 
  # normalization of energy level)
  #test_comparison = [0.009609176542656, 0.050297509665262, \
  #                   0.060601929415464, 0.123926866154943, \
  #                   0.144923595785290, 0.182115436816121]
  
  # If strain is turned off...
  #test_comparison = [0.016995562285266, 0.034073043564595, \
  #                   0.067804273096485, 0.122281860154654, \
  #                   0.151791414122053, 0.167047977422029]

  # Run test calculation ######################################################
  valley_types = ['hh', 'lh', 'so']
  
  print 'TEST: Calculate EV of a SiGe QW using the modified effective mass approximation.'
  print 'Computation started ...'
  
  computed = []
  runtime = []
  num_iter = []
  num_eigenstates = []

  for i in range(len(elem_per_nm)):
    for vid in range(len(valley_types)):
      [vpot_w, meff_w] = materialMEMA_profile_1d(well, valley_types[vid], \
                                                 asub=a0_sub)
      [vpot_b, meff_b] = materialMEMA_profile_1d(barrier, valley_types[vid], \
                                                 asub=a0_sub)
  
      vpot = [vpot_b, vpot_w, vpot_b]
      meff = [meff_b, meff_w, meff_b]

      elem_per_layer = [int(elem_per_nm[i] * x) for x in layer_thickness]
  
      solution = ema_solver_1d(layer_thickness, vpot, meff, \
                               valley_types[vid], efield, elem_per_layer, \
			       porder, solver_info)
 
      # Extract and save eigenvalues
      ev_current_run = []
      for k in range(len(comparison[vid])):
        try:
          nth_eigenstate = solution.get_eigenstate(k)
          ev_current_run.append(nth_eigenstate.get_eigenvalue())
        except:
          print 'Too few eigenvalues have converged.'
  
      computed.append(ev_current_run)

      # Save runtime, number iterations, and number eigenvalues as well
      runtime.append(solver_info.get_runtime())
      num_iter.append(solver_info.get_required_iterations())
      num_eigenstates.append(solution.get_num_eigenstates())
  
  
  # OUTPUT: Show results ###############################################
  print 'Done. Showing results...'
  
  # Reminder for test size
  print 'Elements per nm: ' + str(elem_per_nm)
 
  for vid in range(len(valley_types)):
    print 'Valley type: ' + valley_types[vid]

    for evid in range(len(comparison[vid])):
      for rid in range(len(elem_per_nm)):
        try:
          print '   ' + str(computed[(vid + len(valley_types) * rid)][evid])
        except:
          print '   not converged'
      print '  --------------------'
      print '   ' + str(comparison[vid][evid])
      print ' '
  
    # Show general benchmarks
    print 'time per ev:'
    print [x[0] / x[1] for x in zip(runtime[vid::3], num_eigenstates[vid::3])]
    print 'iterations per ev:'
    print [1.0 * x[0] / x[1] for x in zip(num_iter[vid::3], num_eigenstates[vid::3])]

