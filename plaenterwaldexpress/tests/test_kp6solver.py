import numpy as np
import csv

from ..solverinfo import SolverInfo
from ..kp6solver import kp6_solver_1d
from ..materialprofiles import *
from ..materials import *
from ..utilities import *

def sige_wells(porder=1, elem_per_nm=[1,2,4], num_iterations=1000):
  r"""Convergence test for Si/Ge well on Si barrier material including strain."""
  # Define test scenario ###############################################
  
  # Define parallel k-vector
  kvec = [0.0, 0.0]
  
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
  
  [vpot_well, dso_well, kane_well, eps_well, defpot_well] = \
    materialKP6_profile_1d(well, asub=a0_sub)
  [vpot_barrier, dso_barrier, kane_barrier, eps_barrier, defpot_barrier] = \
    materialKP6_profile_1d(barrier, asub=a0_sub)

  # Define quantum well properties
  layer_thickness = [16.6, 8.3, 16.6]
  
  vpot = [vpot_barrier, vpot_well, vpot_barrier]
  delta_so = [dso_barrier, dso_well, dso_barrier]

  kane = [kane_barrier, kane_well, kane_barrier]
  defpot = [defpot_barrier, defpot_well, defpot_barrier]
  epsilon = [eps_barrier, eps_well, eps_barrier]
  efield = None

  # Setting up solver 
  solver_info = SolverInfo()
  solver_info.set_max_iterations(num_iterations)
  
  # TEST: Same effective masses (simple example) #######################
  
  print 'TEST: k.p solver with Si as barriers and Si(0.6)Ge(0.4) as ' \
        'wells grown on Si(0.7)Ge(0.3) substrate for different number ' \
        'of elements per nm.'
  print 'Convergence test started...'
  
  eigval = []
  runtime = []
  num_iter = []
  num_eigenstates = []
  
  for i in range(len(elem_per_nm)):
    elem_per_layer = [int(elem_per_nm[i] * x) for x in layer_thickness]
  
    eigenpair = kp6_solver_1d(layer_thickness, vpot, \
                              kane, delta_so, \
                              defpot, epsilon, \
                              kvec, efield, \
                              elem_per_layer, porder, \
                              solver_info)
  
    # Get first 6 eigenvalues
    ev_current_run = []
    for k in range(1,24,4):
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
  
  # Compare numerical values with exact numerics TMM results (no re- 
  # normalization of energy level)
  test_comparison = [0.009609176542656, 0.050297509665262, \
                     0.060601929415464, 0.123926866154943, \
                     0.144923595785290, 0.182115436816121]
  
  # If strain is turned off...
  #test_comparison = [0.016995562285266, 0.034073043564595, \
  #                   0.067804273096485, 0.122281860154654, \
  #                   0.151791414122053, 0.167047977422029]
  
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
  
def algaas_wells(porder=1, elem_per_nm=[1,2,4], num_iterations=100):
  r"""Convergence test for a GaAs well on Al(0.5)Ga(0.5)As barrier (no strain)."""
  # Define test scenario ######################################################

  # Define material parameters (well = GaAs, barrier = Al(0.5)Ga(0.5)As)
  mat01 = {'EG':0., 'EL':1., 'EX':1., 'Eav':-0.11367, 'delta_so':0.341, \
           'a0':0., 'C11':0., 'C12':0., 'C44':0., \
           'av':0., 'bv':0., 'aG':0.,'aL':0.,'aX':0., 'xiX':0.} 

  mat02 = {'EG':1., 'EL':1., 'EX':1., 'Eav':-0.36850, 'delta_so':0.3105, \
           'a0':0., 'C11':0., 'C12':0., 'C44':0., \
	   'av':0., 'bv':0., 'aG':0.,'aL':0.,'aX':0., 'xiX':0.}
  
  gamma_well = (6.98, 2.06, 2.93)
  [L_w, M_w, N_w, Np_w, Nm_w] = gamma2LK(gamma_well)
  gamma_barrier = (5.37, 1.44, 2.18)
  [L_b, M_b, N_b, Np_b, Nm_b] = gamma2LK(gamma_barrier) 

  mat01['L'] = L_w; mat02['L'] = L_b 
  mat01['M'] = M_w; mat02['M'] = M_b 
  mat01['N'] = N_w; mat02['N'] = N_b 
  mat01['Np'] = Np_w; mat02['Np'] = Np_b 
  mat01['Nm'] = Nm_w; mat02['Nm'] = Nm_b 

  well = ElementaryKP6(mat01)
  barrier = ElementaryKP6(mat02)

  [vpot_well, dso_well, kane_well, eps_well, defpot_well] = \
    materialKP6_profile_1d(well)
  [vpot_barrier, dso_barrier, kane_barrier, eps_barrier, defpot_barrier] = \
    materialKP6_profile_1d(barrier)

  # Define parallel k-vector
  kvec = [0.05, 0.05]

  # Define quantum well properties ############################################
  layer_thickness = [8.5, 8.5, 8.5]

  vpot = [vpot_barrier, vpot_well, vpot_barrier]
  delta_so = [dso_barrier, dso_well, dso_barrier]

  kane = [kane_barrier, kane_well, kane_barrier]

  defpot = [defpot_barrier, defpot_well, defpot_barrier]
  epsilon = None

  efield = None

  # Setting up solver
  solver_info = SolverInfo()
  solver_info.set_max_iterations(num_iterations)
 
  # TEST: Same effective masses (simple example) #######################
    
  print 'TEST: k.p solver with Al(0.5)Ga(0.5)As as barriers and GaAs ' \
        'as wells with no strain for different number of elements per nm.'
  print 'Convergence test started...'
		
  eigval = []
  runtime = []
  num_iter = []
  num_eigenstates = []
			
  for i in range(len(elem_per_nm)):
    elem_per_layer = [int(elem_per_nm[i] * x) for x in layer_thickness]
  
    eigenpair = kp6_solver_1d(layer_thickness, vpot, \
                              kane, delta_so, \
                              defpot, epsilon, \
                              kvec, efield, \
                              elem_per_layer, porder, \
                              solver_info)
																			      
    # Get first 6 eigenvalues
    ev_current_run = []
    for k in range(1,24,4):
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

  # Compare numerical values with semianalytic results for a finite 
  # potential well 
  test_comparison = [-0.012094879698699, -0.029665641942436, -0.046219503090893, \
                     -0.095549250545684, -0.116954220930596, -0.169892811954494]

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

def algaas_wells_efield(porder=1, elem_per_nm=[1,2,4], num_iterations=100):
  r"""Convergence test for a GaAs well on Al(0.5)Ga(0.5)As barrier and E-field."""
  # Define test scenario ######################################################

  # Define material parameters (well = GaAs, barrier = Al(0.5)Ga(0.5)As)
  mat01 = {'EG':0., 'EL':1., 'EX':1., 'Eav':-0.11367, 'delta_so':0.341, \
           'a0':0., 'C11':0., 'C12':0., 'C44':0., \
           'av':0., 'bv':0., 'aG':0.,'aL':0.,'aX':0., 'xiX':0.} 

  mat02 = {'EG':1., 'EL':1., 'EX':1., 'Eav':-0.36850, 'delta_so':0.3105, \
           'a0':0., 'C11':0., 'C12':0., 'C44':0., \
	   'av':0., 'bv':0., 'aG':0.,'aL':0.,'aX':0., 'xiX':0.}
  
  gamma_well = (6.98, 2.06, 2.93)
  [L_w, M_w, N_w, Np_w, Nm_w] = gamma2LK(gamma_well)
  gamma_barrier = (5.37, 1.44, 2.18)
  [L_b, M_b, N_b, Np_b, Nm_b] = gamma2LK(gamma_barrier) 

  mat01['L'] = L_w; mat02['L'] = L_b 
  mat01['M'] = M_w; mat02['M'] = M_b 
  mat01['N'] = N_w; mat02['N'] = N_b 
  mat01['Np'] = Np_w; mat02['Np'] = Np_b 
  mat01['Nm'] = Nm_w; mat02['Nm'] = Nm_b 

  well = ElementaryKP6(mat01)
  barrier = ElementaryKP6(mat02)

  [vpot_well, dso_well, kane_well, eps_well, defpot_well] = \
    materialKP6_profile_1d(well)
  [vpot_barrier, dso_barrier, kane_barrier, eps_barrier, defpot_barrier] = \
    materialKP6_profile_1d(barrier)

  # Define parallel k-vector
  kvec = [0., 0.]

  # Define quantum well properties ############################################
  layer_thickness = [8.5, 8.5, 8.5]

  vpot = [vpot_barrier, vpot_well, vpot_barrier]
  delta_so = [dso_barrier, dso_well, dso_barrier]

  kane = [kane_barrier, kane_well, kane_barrier]

  defpot = [defpot_barrier, defpot_well, defpot_barrier]
  epsilon = None

  efield = 0.002   # in V/nm

  # Setting up solver
  solver_info = SolverInfo()
  solver_info.set_max_iterations(num_iterations)

  # Define precomputed reference values #######################################

  # Compare numerical values with TMM results 
  test_comparison = [0.010626621180084, 0.030379725706051, \
                     0.043960514829070, 0.097825700452811, \
                     0.113694709501860, 0.170640593382472]

  # Run test calculation ######################################################
    
  print 'TEST: k.p solver with Al(0.5)Ga(0.5)As as barriers and GaAs ' \
        'as wells (no strain) with an applied efield.'
  print 'Convergence test started...'
		
  eigval = []
  runtime = []
  num_iter = []
  num_eigenstates = []
			
  for i in range(len(elem_per_nm)):
    elem_per_layer = [int(elem_per_nm[i] * x) for x in layer_thickness]
  
    eigenpair = kp6_solver_1d(layer_thickness, vpot, \
                              kane, delta_so, \
                              defpot, epsilon, \
                              kvec, efield, \
                              elem_per_layer, porder, \
                              solver_info)
																			      
    # Get first 6 eigenvalues
    ev_current_run = []
    for k in range(1,4 * len(test_comparison),4):
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

  eigenpair.export_solution('solution_kp6.h5')
  
def algaas_wells_tilted(porder=1, elem_per_nm=[1,2,4], num_iterations=100):
  r"""Convergence test for a Al(x)Ga(1-x)As well on Al(0.5)Ga(0.5)As barrier."""
  # Define test scenario ######################################################

  # Define material parameters (well = GaAs, barrier = Al(0.5)Ga(0.5)As)
  mat01 = {'EG':0., 'EL':1., 'EX':1., 'Eav':-0.11367, 'delta_so':0.341, \
           'a0':0., 'C11':0., 'C12':0., 'C44':0., \
           'av':0., 'bv':0., 'aG':0.,'aL':0.,'aX':0., 'xiX':0.} 

  mat02 = {'EG':1., 'EL':1., 'EX':1., 'Eav':-0.36850, 'delta_so':0.3105, \
           'a0':0., 'C11':0., 'C12':0., 'C44':0., \
	   'av':0., 'bv':0., 'aG':0.,'aL':0.,'aX':0., 'xiX':0.}
  
  gamma_well = (6.98, 2.06, 2.93)
  [L_w, M_w, N_w, Np_w, Nm_w] = gamma2LK(gamma_well)
  gamma_barrier = (5.37, 1.44, 2.18)
  [L_b, M_b, N_b, Np_b, Nm_b] = gamma2LK(gamma_barrier) 

  mat01['L'] = L_w; mat02['L'] = L_b 
  mat01['M'] = M_w; mat02['M'] = M_b 
  mat01['N'] = N_w; mat02['N'] = N_b 
  mat01['Np'] = Np_w; mat02['Np'] = Np_b 
  mat01['Nm'] = Nm_w; mat02['Nm'] = Nm_b 

  elem01 = ElementaryKP6(mat01)
  elem02 = ElementaryKP6(mat02)
  
  lin_grad = alloy_gradients_1d('linear', (8.5, 8.5, 0., 0.5)) 
  binary0102 = Binary(elem01, elem02)

  [vpot_well, dso_well, kane_well, eps_well, defpot_well] = \
    materialKP6_profile_1d(binary0102, lin_grad)
  [vpot_barrier, dso_barrier, kane_barrier, eps_barrier, defpot_barrier] = \
    materialKP6_profile_1d(elem02)

  # Define parallel k-vector
  kvec = [0., 0.]

  # Define quantum well properties ############################################
  layer_thickness = [8.5, 8.5, 8.5]

  vpot = [vpot_barrier, vpot_well, vpot_barrier]
  delta_so = [dso_barrier, dso_well, dso_barrier]

  kane = [kane_barrier, kane_well, kane_barrier]

  defpot = [defpot_barrier, defpot_well, defpot_barrier]
  epsilon = None

  efield = None   # in V/nm

  # Setting up solver
  solver_info = SolverInfo()
  solver_info.set_max_iterations(num_iterations)

  # Define precalculated reference values #####################################

  # Compare numerical values with results from nextnano3 
  test_comparison = [-0.1854, -0.2050, \
                     -0.2341, -0.2649]
   # Run test calculation #####################################################
    
  print 'TEST: k.p solver with Al(0.5)Ga(0.5)As as barriers and Al(x)Ga(1-x)As ' \
        'as wells (tilted potential) with no strain.'
  print 'Convergence test started...'
		
  eigval = []
  runtime = []
  num_iter = []
  num_eigenstates = []
			
  for i in range(len(elem_per_nm)):
    elem_per_layer = [int(elem_per_nm[i] * x) for x in layer_thickness]
  
    eigenpair = kp6_solver_1d(layer_thickness, vpot, \
                              kane, delta_so, \
                              defpot, epsilon, \
                              kvec, efield, \
                              elem_per_layer, porder, \
                              solver_info)
																			      
    # Get first 6 eigenvalues
    ev_current_run = []
    for k in range(1,4 * len(test_comparison),4):
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

