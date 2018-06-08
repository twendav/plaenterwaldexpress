from dolfin import *
import numpy as np
import timeit as timeit 
import pkg_resources as pkg

from mesher import mesh_gen_1d
from materialinfo import *
from solutioncontainer import SolutionContainer

def ema_solver_1d(layer_thickness, vpot, meff, valley_type, efield, \
               elem_per_layer, porder, solver_info):
  
  # Testing FEniCS capabilities ######################################
 
  # Test for PETSc and SLEPc
  if not has_linear_algebra_backend("PETSc"):
    print "DOLFIN has not been configured with PETSc. Exiting."
    exit()
 
  if not has_slepc():
    print "DOLFIN has not been configured with SLEPc. Exiting."
    exit()
 
  # Define physical constants ########################################
  HB2M = Constant(0.03809982)

  # Calculate additional parameter from input ########################
  layer_startpos = np.cumsum([0] + layer_thickness)
  xmid = sum(layer_thickness) / 2.0

  # Generate mesh ####################################################
  mesh = mesh_gen_1d(layer_thickness, elem_per_layer)

  # Define function spaces ###########################################
  Vp = FunctionSpace(mesh, 'CG', porder)

  # Define boundary conditions #######################################

  def u0_boundary(x, on_boundary):
    return on_boundary

  bcs = DirichletBC(Vp, 0.0, u0_boundary)


  # Define subdomains ################################################

  class Layer(SubDomain):
    def __init__(self, interval):
      SubDomain.__init__(self)
      self.interval = interval
    def inside(self, x, on_boundary):
      return between(x[0], self.interval)

  # Create layers
  layers = []
  curpos = 0

  for k in range(len(layer_thickness)):
    layers.append(Layer((curpos, curpos + layer_thickness[k])))
    curpos = curpos + layer_thickness[k]

  subdomains = CellFunction('size_t', mesh)
  for k in range(len(layer_thickness)):
    layers[k].mark(subdomains, k)

  # Define test and trial function ###################################
  u = TrialFunction(Vp)
  v = TestFunction(Vp)
  
  # Set up problem ###################################################
  
  dx = Measure('dx')[subdomains]
 
  # Prepare stiffness and mass matrix
  stiffness = 0
  mass = 0
  
  for subdomain_id in range(len(layer_thickness)):
    if valley_type == 'G':
      mass_term = 1. / meff[subdomain_id]
      pot_term = vpot[subdomain_id]
    elif valley_type == 'L':
      mass_term = 1. / 3. / meff[subdomain_id][0] + 2. / 3. / meff[subdomain_id][1]
      pot_term = vpot[subdomain_id]
    elif valley_type == 'X100' or valley_type == 'X010':
      mass_term = 1. / meff[subdomain_id][1]
      pot_term = vpot[subdomain_id]
    elif valley_type == 'X001':
      mass_term = 1. / meff[subdomain_id][0]
      pot_term = vpot[subdomain_id]
    elif valley_type == 'hh' or valley_type == 'lh' or valley_type == 'so':
      mass_term = 1. / meff[subdomain_id]
      pot_term = -vpot[subdomain_id]
    else:
      raise Exception("Valley identifier '" + valley_type + "' unknown.")
    
    stiffness += HB2M * mass_term * inner(grad(u),grad(v)) * dx(subdomain_id) + \
                 pot_term * inner(u,v) * dx(subdomain_id)
    mass += inner(u,v) * dx(subdomain_id)

    if efield > 0 or efield < 0:
      stiffness += efield * Expression('x[0] - xm', xm=xmid) * \
                   inner(u,v) * dx(subdomain_id)
  
  # Assemble eigenvalue problem
  S = PETScMatrix()
  M = PETScMatrix()
  
  L = sum(map(lambda k: Constant(0.) * v * dx(k), \
              range(len(layer_thickness))))
  
  assemble_system(stiffness,L,bcs,A_tensor=S)
  assemble_system(mass,L,A_tensor=M)
  
  # Create eigensolver for generalized EVP from given settings
  eigensolver = SLEPcEigenSolver(S, M)

  eigensolver.parameters['solver'] = \
    solver_info.solver_type
  eigensolver.parameters['spectrum'] = \
    solver_info.spectrum
  eigensolver.parameters['maximum_iterations'] = \
    solver_info.max_iterations

  if hasattr(solver_info, 'tolerance'):
    eigensolver.parameters['tolerance'] = solver_info.tolerance

  if hasattr(solver_info, 'spectral_transform'):
    eigensolver.parameters['spectral_transform'] = \
        solver_info.spectral_transform

  if hasattr(solver_info, 'shift_and_invert'):
    eigensolver.parameters['shift-and-invert'] = \
        solver_info.shift_and_invert

  if hasattr(solver_info, 'spectral_shift'):
    eigensolver.parameters['spectral_shift'] = \
        solver_info.spectral_shift


  # Compute generalized eigenvalue decomposition and time process
  print "Computing eigenvalue decomposition. This may take a while."
  start_solver = timeit.default_timer()
  if solver_info.num_ev:
    eigensolver.solve(solver_info.num_ev)
  else:
    eigensolver.solve()
  stop_solver = timeit.default_timer()

  # Save solution to Solution object #################################

  # Check if enough eigenvalues have converged
  num_ev_conv = eigensolver.get_number_converged()
  
  print 'Eigenvalue solver: ' + str(num_ev_conv) + \
      ' eigenvalues have converged.'

  # Get mesh coordinates 
  dof_coordinates = Vp.dofmap().tabulate_all_coordinates(mesh)
  space_dim = mesh.geometry().dim()

  # Prepare Solution object for usage
  solution = SolutionContainer(dof_coordinates, space_dim)

  solver_info.set_runtime(stop_solver - start_solver)
  solver_info.set_required_iterations(eigensolver.get_iteration_number())

  vdolfin = dolfin_version()
  try:
    vplaenterwald = pkg.require("plaenterwaldexpress")[0].version
  except:
    vplaenterwald = "dev"
  solver_info.set_version(vdolfin, vplaenterwald)

  solution.add_solver_info(solver_info)

  # Set up function space for extracting solution
  eigenfct = Function(Vp)

  for i in range(num_ev_conv):
    ev_re, ev_im, ef_re, ef_im = eigensolver.get_eigenpair(i)
    
    # Save eigenvalue
    if valley_type == 'hh' or valley_type == 'lh' or valley_type == 'so':
      ev_complex = -(ev_re + 1j * ev_im)
    else:
      ev_complex = ev_re + 1j * ev_im

    # Extract and save eigenfunction
    u_re = Function(Vp)
    u_im = Function(Vp)
    u_re.vector()[:] = ef_re
    u_im.vector()[:] = ef_im
    ef_complex = np.array([a + 1j * b for a, b in \
                 zip(u_re.vector().array(), u_im.vector().array())])

    solution.add_eigenstate(ev_complex, ef_complex)

  # Add potential information to solution object
  solution.add_material_info(potential_to_materialinfo(subdomains, \
                             vpot, efield, xmid), valley_type)

  return solution

