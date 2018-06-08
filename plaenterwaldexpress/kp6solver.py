from dolfin import *
import numpy as np
import timeit as timeit

from mesher import mesh_gen_1d
from materialinfo import *
from solutioncontainer import SolutionContainer

def kp6_solver_1d(layer_thickness, vpot, kane, delta_so, \
                  def_pot, epsilon, kvec, efield, elem_per_layer, \
		  porder, solver_info):

  # Define numerics parameters #########################################
  num_bands = 6
  
  # Define physical constants ##########################################
  HB2M = 0.03809982
 
  # Renormalize Kane parameter
  tkane = [[row[i] for row in kane] for i in range(5)]
  Lp = [HB2M * elem for elem in tkane[0]]
  M = [HB2M * elem for elem in tkane[1]]
  Np = [HB2M * elem for elem in tkane[2]]
  Npp = [HB2M * elem for elem in tkane[3]]
  Nm = [HB2M * elem for elem in tkane[4]]

  # Redefine k-vector
  kx = kvec[0]
  ky = kvec[1]

  # Generate mesh ######################################################
  mesh = mesh_gen_1d(layer_thickness, elem_per_layer)
  xmid = sum(layer_thickness) / 2.0

  # Set function spaces ################################################
  VS = VectorFunctionSpace(mesh, 'CG', porder, num_bands)
  MS = MixedFunctionSpace([VS, VS])
  
  # Define domains and boundaries ######################################
  
  # Define boundary conditions
  def u0_boundary(x, on_boundary):
    return on_boundary
  
  bcs = DirichletBC(MS, (0.0,) * 2 * num_bands, u0_boundary)
  
  # Define subdomains
  class Layers(SubDomain):
    def __init__(self, interval):
      SubDomain.__init__(self)
      self.interval = interval
    def inside(self, x, on_boundary):
      return between(x[0], self.interval)
  
  # Create layers
  layers = []
  curpos = 0
  
  for k in range(len(layer_thickness)):
    layers.append(Layers((curpos, curpos + layer_thickness[k])))
    curpos = curpos + layer_thickness[k]
  
  subdomains = MeshFunction('size_t', mesh, 1)
  for k in range(len(layer_thickness)):
    layers[k].mark(subdomains, k)
  
  # Define test and trail function ############################################
  (ur, uc) = TrialFunctions(MS)
  (vr, vc) = TestFunctions(MS)

  # Set up problem ############################################################

  dx = Measure('dx')[subdomains]

  # Prepare stiffness and mass matrix
  stiffness = 0
  mass = 0

  for subdomain_id in range(len(layer_thickness)):
    # Set up Re(H0)
    h0_r = as_matrix( \
             [[vpot[subdomain_id] - delta_so[subdomain_id] / 3 + Lp[subdomain_id] * kx**2 + M[subdomain_id] * ky**2, \
	       Np[subdomain_id] * kx * ky, 0, 0, 0, delta_so[subdomain_id] / 3], \
	      [Np[subdomain_id] * kx * ky, \
	       vpot[subdomain_id] - delta_so[subdomain_id] / 3 + Lp[subdomain_id] * ky**2 +  M[subdomain_id] * kx**2, \
               0, 0, 0, 0], \
	      [0, 0, vpot[subdomain_id] - delta_so[subdomain_id] / 3 + M[subdomain_id] * (kx**2 + ky**2), \
	       -delta_so[subdomain_id] / 3, 0, 0], \
              [0, 0, -delta_so[subdomain_id] / 3, \
	       vpot[subdomain_id] - delta_so[subdomain_id] / 3 + Lp[subdomain_id] * kx**2 + M[subdomain_id] * ky**2, \
	       Np[subdomain_id] * kx * ky, 0], \
	      [0, 0, 0, Np[subdomain_id] * kx * ky, \
	       vpot[subdomain_id] - delta_so[subdomain_id] / 3 + Lp[subdomain_id] * ky**2 +  M[subdomain_id] * kx**2, 0], \
	      [delta_so[subdomain_id] / 3, 0, 0, 0, 0, \
	       vpot[subdomain_id] - delta_so[subdomain_id] / 3 + M[subdomain_id] * (kx**2 + ky**2)]])
     
    # Set up Im(H0)
    h0_c = as_matrix( \
             [[0, -delta_so[subdomain_id] / 3, 0, 0, 0, 0], \
              [delta_so[subdomain_id] / 3, 0, 0, 0, 0, -delta_so[subdomain_id] / 3], \
	      [0, 0, 0, 0, delta_so[subdomain_id] / 3, 0], \
              [0, 0, 0, 0, delta_so[subdomain_id] / 3, 0], \
	      [0, 0, -delta_so[subdomain_id] / 3, -delta_so[subdomain_id] / 3, 0, 0], \
	      [0, delta_so[subdomain_id] / 3, 0, 0, 0, 0]])

    # Set up Re(H1R)
    h1r_r = as_matrix( \
              [[0, 0, Nm[subdomain_id] * kx, 0, 0, 0], \
	       [0, 0, Nm[subdomain_id] * ky, 0, 0, 0], \
	       [Npp[subdomain_id] * kx, Npp[subdomain_id] * ky, 0, 0, 0, 0], \
	       [0, 0, 0, 0, 0, Nm[subdomain_id] * kx], \
	       [0, 0, 0, 0, 0, Nm[subdomain_id] * ky], \
	       [0, 0, 0, Npp[subdomain_id] * kx, Npp[subdomain_id] * ky, 0]])
    
    # Set up Re(H1L)
    h1l_r = as_matrix(
              [[0, 0, kx * Npp[subdomain_id], 0, 0, 0], \
               [0, 0, ky * Npp[subdomain_id], 0, 0, 0], \
               [kx * Nm[subdomain_id], ky * Nm[subdomain_id], 0, 0, 0, 0], \
               [0, 0, 0, 0, 0, kx * Npp[subdomain_id]], \
	       [0, 0, 0, 0, 0, ky * Npp[subdomain_id]], \
	       [0, 0, 0, kx * Nm[subdomain_id], ky * Nm[subdomain_id], 0]])
    
    # Set up H2
    h2 = as_matrix( \
           [[M[subdomain_id], 0, 0, 0, 0, 0], \
	    [0, M[subdomain_id], 0, 0, 0, 0], \
	    [0, 0, Lp[subdomain_id], 0, 0, 0], \
	    [0, 0, 0, M[subdomain_id], 0, 0], \
	    [0, 0, 0, 0, M[subdomain_id], 0], \
	    [0, 0, 0, 0, 0, Lp[subdomain_id]]])

    # Set up Hs (strain matrix) if strain information are given
    if epsilon:
      eps_para = epsilon[subdomain_id][0]
      eps_perp = epsilon[subdomain_id][1]
      def_l = def_pot[subdomain_id][0]
      def_m = def_pot[subdomain_id][1]

      hs = as_matrix( \
           [[def_l * eps_para + def_m * (eps_para + eps_perp), 0, 0, 0, 0, 0], \
	    [0, def_l * eps_para + def_m * (eps_para + eps_perp), 0, 0, 0, 0], \
            [0, 0, def_l * eps_perp + def_m * 2 * eps_para, 0, 0, 0], \
            [0, 0, 0, def_l * eps_para + def_m * (eps_para + eps_perp), 0, 0], \
	    [0, 0, 0, 0, def_l * eps_para + def_m * (eps_para + eps_perp), 0], \
            [0, 0, 0, 0, 0, def_l * eps_perp + def_m * 2 * eps_para]])
    else:
      hs = None


    i, j = indices(2)
    hamiltonian_rrr = vr[i]*h0_r[i,j]*ur[j]*dx(subdomain_id) + \
      Dx(vr[i],0)*h2[i,j]*Dx(ur[j],0)*dx(subdomain_id) 
    if hs:
      hamiltonian_rrr += vr[i]*hs[i,j]*ur[j]*dx(subdomain_id) 
 
    hamiltonian_rcc = vr[i]*h0_c[i,j]*uc[j]*dx(subdomain_id) 
    if kx != 0.0 or ky != 0.0:
      hamiltonian_rcc += Dx(vr[i],0)*h1r_r[i,j]*uc[j]*dx(subdomain_id) - \
                         vr[i]*h1l_r[i,j]*Dx(uc[j],0)*dx(subdomain_id)
  
    hamiltonian_ccr = vc[i]*h0_c[i,j]*ur[j]*dx(subdomain_id) 
    if kx != 0.0 or ky != 0.0:
      hamiltonian_ccr += Dx(vc[i],0)*h1r_r[i,j]*ur[j]*dx(subdomain_id) - \
                        vc[i]*h1l_r[i,j]*Dx(ur[j],0)*dx(subdomain_id)

    hamiltonian_crc = vc[i]*h0_r[i,j]*uc[j]*dx(subdomain_id) + \
      Dx(vc[i],0)*h2[i,j]*Dx(uc[j],0)*dx(subdomain_id)
    if hs:
      hamiltonian_crc += vc[i]*hs[i,j]*uc[j]*dx(subdomain_id)

    stiffness += (hamiltonian_rrr - hamiltonian_rcc + \
      hamiltonian_ccr + hamiltonian_crc)
    
    if efield and efield != 0.0:
      stiffness += vr[i]*Expression('x[0] - xm', xm=xmid)* \
                   efield*ur[i]*dx(subdomain_id) + \
                   vc[i]*Expression('x[0] - xm', xm=xmid)* \
                   efield*uc[i]*dx(subdomain_id) 

    mass += ur[i]*vr[i]*dx(subdomain_id) + uc[i]*vc[i]*dx(subdomain_id)
 
  # Assemle eigenvalue problem
  S = PETScMatrix()
  M = PETScMatrix()
  
  L = sum(map(lambda k: dot(Constant((0.0,) * num_bands), vr) * dx(k) + \
        dot(Constant((0.0,) * num_bands), vc) * dx(k), range(len(layer_thickness))))
  
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

  if hasattr(solver_info, 'spectral_transform'):
    eigensolver.parameters['spectral_transform'] = \
        solver_info.spectral_transform

  if hasattr(solver_info, 'spectral_shift'):
    eigensolver.parameters['spectral_shift'] = \
        solver_info.spectral_shift

  # Compute generalized eigenvalue decomposition
  print "Computing eigenvalue decomposition. This may take a while."
  start_solver = timeit.default_timer()
  if solver_info.num_ev:
    eigensolver.solve(solver_info.num_ev)
  else:
    eigensolver.solve()
  stop_solver =  timeit.default_timer()

  # Save solution to Solution object #################################

  # Check if enough eigenvalues have converged
  num_ev_conv = eigensolver.get_number_converged()
  
  print 'Eigenvalue solver: ' + str(num_ev_conv) + \
      ' eigenvalues have converged.'

  # Get mesh coordinates
  dof_coordinates = VS.dofmap().tabulate_all_coordinates(mesh)
  space_dim = mesh.geometry().dim()

  # Prepare Solution object for usage
  dofs_0 = VS.sub(0).dofmap().dofs()
  dofs_1 = VS.sub(1).dofmap().dofs()
  dofs_2 = VS.sub(2).dofmap().dofs()
  dofs_3 = VS.sub(3).dofmap().dofs()
  dofs_4 = VS.sub(4).dofmap().dofs()
  dofs_5 = VS.sub(5).dofmap().dofs()
  
  solution = SolutionContainer(dof_coordinates, \
                               space_dim, \
                               dofs_0, dofs_1, dofs_2, \
			       dofs_3, dofs_4, dofs_5)

  solver_info.set_runtime(stop_solver - start_solver)
  solver_info.set_required_iterations(eigensolver.get_iteration_number()) 

  vdolfin = dolfin_version()
  try:
    vplaenterwald = pkg.require("plaenterwaldexpress")[0].version
  except:
    vplaenterwald = "dev"
  solver_info.set_version(vdolfin, vplaenterwald)

  solution.add_solver_info(solver_info)

  # Extract eigenvalues and eigenfunctions
  for i in range(num_ev_conv):
    # Get information of current eigenstate from solver
    ev_re, ev_im, ef_re, ef_im = eigensolver.get_eigenpair(i)

    # Save eigenvalue
    ev_complex = ev_re + 1j * ev_im

    # Extract and save eigenfunction
    u = Function(MS)
    u.vector()[:] = ef_re

    u_re, u_im = u.split(True)
    
    # TESTING
    if i==0:
      mesh2 = UnitSquareMesh(1,1)
      #plot(mesh2, interactive=True)
      #plot(u_im, interactive=True)
    # TESTING

    ef_complex = np.array([a + 1j * b for a, b in zip(u_re.vector().array(), \
                                             u_im.vector().array())])

    solution.add_eigenstate(ev_complex, ef_complex)
  
  # Add potential information to solution object
  solution.add_material_info(potential_to_materialinfo(subdomains, \
                             vpot, efield, xmid), 'vband')

  return solution

