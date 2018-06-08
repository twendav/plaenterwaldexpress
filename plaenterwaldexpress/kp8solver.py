from dolfin import *
import numpy as np
import timeit as timeit

from mesher import mesh_gen_1d
from materialinfo import *
from solutioncontainer import SolutionContainer

def kp8_solver_1d(layer_thickness, cpot, vpot, kane, delta_so, \
                  defpot, epsilon, kvec, efield, elem_per_layer, \
                  porder, solver_info):

  # Define numerics parameters ################################################
  num_bands = 8

  # Define physical constants ##########################################
  HB2M = 0.03809982
  
  # Renormalize Kane parameter
  tkane = [[row[i] for row in kane] for i in range(7)]
  Ac = [HB2M * elem for elem in tkane[0]]
  P = tkane[1]
  Lp = [HB2M * elem for elem in tkane[2]]
  M = [HB2M * elem for elem in tkane[3]]
  Np = [HB2M * elem for elem in tkane[4]]
  Npp = [HB2M * elem for elem in tkane[5]]
  Nm = [HB2M * elem for elem in tkane[6]]

  # Renormalize vpot (spin-orbit correction)
  vpotp = [(vp - dso / 3.0) for (vp, dso) in zip(vpot, delta_so)]

  # Redefine k-vector
  kx = kvec[0]
  ky = kvec[1]

  # Generate mesh #############################################################
  mesh = mesh_gen_1d(layer_thickness, elem_per_layer)
  xmid = sum(layer_thickness) / 2.0

  # Set function spaces #######################################################
  VS = VectorFunctionSpace(mesh, 'CG', porder, num_bands)
  MS = MixedFunctionSpace([VS, VS])

  # Define domains and boundaries #############################################
  
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

  # Define discontinuous Hamiltonians #########################################
  dx = Measure('dx')[subdomains]
  
  # Defining stiffness and mass matrix
  stiffness = 0
  mass = 0

  for subdomain_id in range(len(layer_thickness)):
    # Set up Re(H0)
    h0_r = as_matrix( \
             [[cpot[subdomain_id] + Ac[subdomain_id] * (kx**2 + ky**2), 0, 0, 0, 0, 0, 0, 0], \
	      [0, vpotp[subdomain_id] + Lp[subdomain_id] * kx**2 + M[subdomain_id] * ky**2, \
	       Np[subdomain_id] * kx * ky, 0, 0, 0, 0, 0], \
	      [0, Np[subdomain_id] * kx * ky, \
	       vpotp[subdomain_id] + Lp[subdomain_id] * ky**2 +  M[subdomain_id] * kx**2, \
               0, 0, 0, 0, 0], \
	      [0, 0, 0, vpotp[subdomain_id] + M[subdomain_id] * (kx**2 + ky**2), \
	       0, 0, 0, 0], \
              [0, 0, 0, 0, \
	       cpot[subdomain_id] + Ac[subdomain_id] * (kx**2 + ky**2), 0, 0, 0], \
	      [0, 0, 0, 0, 0, \
	       vpotp[subdomain_id] + Lp[subdomain_id] * kx**2 + M[subdomain_id] * ky**2, \
	       Np[subdomain_id] * kx * ky, 0], \
	      [0, 0, 0, 0, 0, Np[subdomain_id] * kx * ky, \
	       vpotp[subdomain_id] + Lp[subdomain_id] * ky**2 + M[subdomain_id] * kx**2, 0], \
	      [0, 0, 0, 0, 0, 0, 0, \
	       vpotp[subdomain_id] + M[subdomain_id] * (kx**2 + ky**2)]]) + \
	    delta_so[subdomain_id] / 3.0 * as_matrix( \
	      [[0, 0, 0, 0, 0, 0, 0, 0], \
	       [0, 0, 0, 0, 0, 0, 0, 1], \
	       [0, 0, 0, 0, 0, 0, 0, 0], \
	       [0, 0, 0, 0, 0,-1, 0, 0], \
	       [0, 0, 0, 0, 0, 0, 0, 0], \
	       [0, 0, 0,-1, 0, 0, 0, 0], \
	       [0, 0, 0, 0, 0, 0, 0, 0], \
	       [0, 1, 0, 0, 0, 0, 0, 0]])

    #h0_r.vector()[index_list[0]] = cband[subdomain_id] + \
    #  Ac[subdomain_id] * kpara2
    #h0_r.vector()[index_list[9]] = vband[subdomain_id] - \
    #  delta_so[subdomain_id] / 3 + Lp[subdomain_id] * kx**2 + \
    #  M[subdomain_id] * ky**2
    #h0_r.vector()[index_list[10]] = Np[subdomain_id] * kx * ky
    #h0_r.vector()[index_list[17]] = Np[subdomain_id] * kx * ky
    #h0_r.vector()[index_list[18]] = vband[subdomain_id] - \
    #  delta_so[subdomain_id] / 3 + Lp[subdomain_id] * ky**2 + \
    #  M[subdomain_id] * kx**2
    #h0_r.vector()[index_list[27]] = vband[subdomain_id] - \
    #  delta_so[subdomain_id] / 3 + M[subdomain_id] * kpara2
    #h0_r.vector()[index_list[0 + 36]] = h0_r.vector()[index_list[0]]
    #h0_r.vector()[index_list[9 + 36]] = h0_r.vector()[index_list[9]]
    #h0_r.vector()[index_list[10 + 36]] = h0_r.vector()[index_list[10]]
    #h0_r.vector()[index_list[17 + 36]] = h0_r.vector()[index_list[17]]
    #h0_r.vector()[index_list[18 + 36]] = h0_r.vector()[index_list[18]]
    #h0_r.vector()[index_list[27 + 36]] = h0_r.vector()[index_list[27]]
    #h0_r.vector()[index_list[15]] = delta_so[subdomain_id] / 3
    #h0_r.vector()[index_list[29]] = -delta_so[subdomain_id] / 3
    #h0_r.vector()[index_list[43]] = -delta_so[subdomain_id] / 3
    #h0_r.vector()[index_list[57]] = delta_so[subdomain_id] / 3

    # Set up Im(H0)
    h0_c = as_matrix(\
             [[0, P[subdomain_id] * kx, P[subdomain_id] * ky, 0, 0, 0, 0, 0], \
	      [-kx * P[subdomain_id], 0, 0, 0, 0, 0, 0, 0], \
	      [-ky * P[subdomain_id], 0, 0, 0, 0, 0, 0, 0], \
	      [0, 0, 0, 0, 0, 0, 0, 0], \
	      [0, 0, 0, 0, 0, P[subdomain_id] * kx, P[subdomain_id] * ky, 0], \
	      [0, 0, 0, 0, -kx * P[subdomain_id], 0, 0, 0], \
	      [0, 0, 0, 0, -ky * P[subdomain_id], 0, 0, 0], \
	      [0, 0, 0, 0, 0, 0, 0, 0]]) + \
	   delta_so[subdomain_id] / 3.0 * as_matrix(\
	     [[0, 0, 0, 0, 0, 0, 0, 0], \
	      [0, 0,-1, 0, 0, 0, 0, 0], \
	      [0, 1, 0, 0, 0, 0, 0,-1], \
	      [0, 0, 0, 0, 0, 0, 1, 0], \
	      [0, 0, 0, 0, 0, 0, 0, 0], \
	      [0, 0, 0, 0, 0, 0, 1, 0], \
	      [0, 0, 0,-1, 0,-1, 0, 0], \
	      [0, 0, 1, 0, 0, 0, 0, 0]])
    #h0_c.vector()[index_list[1]] = P[subdomain_id] * kx
    #h0_c.vector()[index_list[2]] = P[subdomain_id] * ky
    #h0_c.vector()[index_list[8]] = -kx * P[subdomain_id]
    #h0_c.vector()[index_list[16]] = -ky * P[subdomain_id]
    #h0_c.vector()[index_list[1 + 36]] = h0_c.vector()[index_list[1]]
    #h0_c.vector()[index_list[2 + 36]] = h0_c.vector()[index_list[2]]
    #h0_c.vector()[index_list[8 + 36]] = h0_c.vector()[index_list[8]]
    #h0_c.vector()[index_list[16 + 36]] = h0_c.vector()[index_list[16]]
    #h0_c.vector()[index_list[10]] = delta_so[subdomain_id] / 3
    #h0_c.vector()[index_list[17]] = -delta_so[subdomain_id] / 3
    #h0_c.vector()[index_list[23]] = delta_so[subdomain_id] / 3
    #h0_c.vector()[index_list[30]] = -delta_so[subdomain_id] / 3
    #h0_c.vector()[index_list[46]] = -delta_so[subdomain_id] / 3
    #h0_c.vector()[index_list[51]] = delta_so[subdomain_id] / 3
    #h0_c.vector()[index_list[53]] = delta_so[subdomain_id] / 3
    #h0_c.vector()[index_list[58]] = -delta_so[subdomain_id] / 3

    # Set up Re(H1R)
    h1r_r = as_matrix( \
              [[0, 0, 0, 0, 0, 0, 0, 0], \
	       [0, 0, 0, Nm[subdomain_id] * kx, 0, 0, 0, 0], \
	       [0, 0, 0, Nm[subdomain_id] * ky, 0, 0, 0, 0], \
	       [0, Npp[subdomain_id] * kx, Npp[subdomain_id] * ky, 0, 0, 0, 0, 0], \
	       [0, 0, 0, 0, 0, 0, 0, 0], \
	       [0, 0, 0, 0, 0, 0, 0, Nm[subdomain_id] * kx], \
	       [0, 0, 0, 0, 0, 0, 0, Nm[subdomain_id] * ky], \
	       [0, 0, 0, 0, 0, Npp[subdomain_id] * kx, Npp[subdomain_id] * ky, 0]])
    
    #h1r_r.vector()[index_list[11]] = Nm[subdomain_id] * kx
    #h1r_r.vector()[index_list[19]] = Nm[subdomain_id] * ky
    #h1r_r.vector()[index_list[25]] = Npp[subdomain_id] * kx
    #h1r_r.vector()[index_list[26]] = Npp[subdomain_id] * ky
    #h1r_r.vector()[index_list[11 + 36]] = h1r_r.vector()[index_list[11]] 
    #h1r_r.vector()[index_list[19 + 36]] = h1r_r.vector()[index_list[19]] 
    #h1r_r.vector()[index_list[25 + 36]] = h1r_r.vector()[index_list[25]]
    #h1r_r.vector()[index_list[26 + 36]] = h1r_r.vector()[index_list[26]]
    
    # Set up Im(H1R)
    h1r_c = as_matrix(\
              [[0, 0, 0, 0, 0, 0, 0, 0], \
	       [0, 0, 0, 0, 0, 0, 0, 0], \
	       [0, 0, 0, 0, 0, 0, 0, 0], \
	       [-P[subdomain_id], 0, 0, 0, 0, 0, 0, 0], \
	       [0, 0, 0, 0, 0, 0, 0, 0], \
	       [0, 0, 0, 0, 0, 0, 0, 0], \
	       [0, 0, 0, 0, 0, 0, 0, 0], \
	       [0, 0, 0, 0, -P[subdomain_id], 0, 0, 0]])

    #h1r_c.vector()[index_list[3]] = -P[subdomain_id]
    #h1r_c.vector()[index_list[3 + 36]] =  h1r_c.vector()[index_list[3]]

    # Set up Re(H1L)
    h1l_r = as_matrix( \
              [[0, 0, 0, 0, 0, 0, 0, 0], \
	       [0, 0, 0, kx * Npp[subdomain_id], 0, 0, 0, 0], \
	       [0, 0, 0, ky * Npp[subdomain_id], 0, 0, 0, 0], \
	       [0, kx * Nm[subdomain_id], ky * Nm[subdomain_id], 0, 0, 0, 0, 0], \
	       [0, 0, 0, 0, 0, 0, 0, 0], \
	       [0, 0, 0, 0, 0, 0, 0, kx * Npp[subdomain_id]], \
	       [0, 0, 0, 0, 0, 0, 0, ky * Npp[subdomain_id]], \
	       [0, 0, 0, 0, 0, kx * Nm[subdomain_id], ky * Nm[subdomain_id], 0]])
     
    #h1l_r.vector()[index_list[11]] = kx * Npp[subdomain_id] 
    #h1l_r.vector()[index_list[19]] = ky * Npp[subdomain_id]
    #h1l_r.vector()[index_list[25]] = kx * Nm[subdomain_id] 
    #h1l_r.vector()[index_list[26]] = ky * Nm[subdomain_id] 
    #h1l_r.vector()[index_list[11 + 36]] = h1l_r.vector()[index_list[11]] 
    #h1l_r.vector()[index_list[19 + 36]] = h1l_r.vector()[index_list[19]] 
    #h1l_r.vector()[index_list[25 + 36]] = h1l_r.vector()[index_list[25]]
    #h1l_r.vector()[index_list[26 + 36]] = h1l_r.vector()[index_list[26]]
    
    # Set up Im(H1L)
    h1l_c = as_matrix(\
              [[0, 0, 0, P[subdomain_id], 0, 0, 0, 0], \
	       [0, 0, 0, 0, 0, 0, 0, 0], \
	       [0, 0, 0, 0, 0, 0, 0, 0], \
	       [0, 0, 0, 0, 0, 0, 0, 0], \
	       [0, 0, 0, 0, 0, 0, 0, P[subdomain_id]], \
	       [0, 0, 0, 0, 0, 0, 0, 0], \
	       [0, 0, 0, 0, 0, 0, 0, 0], \
	       [0, 0, 0, 0, 0, 0, 0, 0]])

    #h1l_c.vector()[index_list[24]] = P[subdomain_id]
    #h1l_c.vector()[index_list[24 + 36]] =  h1l_c.vector()[index_list[24]]

    # Set up H2
    h2 = as_matrix( \
           [[Ac[subdomain_id], 0, 0, 0, 0, 0, 0, 0], \
	    [0, M[subdomain_id], 0, 0, 0, 0, 0, 0], \
	    [0, 0, M[subdomain_id], 0, 0, 0, 0, 0], \
	    [0, 0, 0, Lp[subdomain_id], 0, 0, 0, 0], \
	    [0, 0, 0, 0, Ac[subdomain_id], 0, 0, 0], \
	    [0, 0, 0, 0, 0, M[subdomain_id], 0, 0], \
	    [0, 0, 0, 0, 0, 0, M[subdomain_id], 0], \
	    [0, 0, 0, 0, 0, 0, 0, Lp[subdomain_id]]])

    #h2.vector()[index_list[0]] = Ac[subdomain_id] 
    #h2.vector()[index_list[9]] = M[subdomain_id] 
    #h2.vector()[index_list[18]] = M[subdomain_id]
    #h2.vector()[index_list[27]] = Lp[subdomain_id]
    #h2.vector()[index_list[0 + 36]] = h2.vector()[index_list[0]] 
    #h2.vector()[index_list[9 + 36]] = h2.vector()[index_list[9]] 
    #h2.vector()[index_list[18 + 36]] = h2.vector()[index_list[18]]
    #h2.vector()[index_list[27 + 36]] = h2.vector()[index_list[27]]

    # Set up Hs (strain matrix) if strain information are given
    if epsilon:
      eps_para = epsilon[subdomain_id][0]
      eps_perp = epsilon[subdomain_id][1]
      def_aG = defpot[subdomain_id][0]
      def_l = defpot[subdomain_id][1]
      def_m = defpot[subdomain_id][2]
 
      hs = as_matrix( \
           [[def_aG * (2 * eps_para + eps_perp), 0, 0, 0, 0, 0, 0, 0], \
	    [0, def_l * eps_para + def_m * (eps_para + eps_perp), 0, 0, 0, 0, 0, 0], \
	    [0, 0, def_l * eps_para + def_m * (eps_para + eps_perp), 0, 0, 0, 0, 0], \
            [0, 0, 0, def_l * eps_perp + def_m * 2 * eps_para, 0, 0, 0, 0], \
	    [0, 0, 0, 0, def_aG * (2 * eps_para + eps_perp), 0, 0, 0], \
            [0, 0, 0, 0, 0, def_l * eps_para + def_m * (eps_para + eps_perp), 0, 0], \
	    [0, 0, 0, 0, 0, 0, def_l * eps_para + def_m * (eps_para + eps_perp), 0], \
            [0, 0, 0, 0, 0, 0, 0, def_l * eps_perp + def_m * 2 * eps_para]])
    else:
      hs = None

    # Weak formulation for current layer
    i, j = indices(2)
 
    hamiltonian_rrr = vr[i]*h0_r[i,j]*ur[j]*dx(subdomain_id) - \
                      Dx(vr[i],0)*h1r_c[i,j]*ur[j]*dx(subdomain_id) + \
                      vr[i]*h1l_c[i,j]*Dx(ur[j],0)*dx(subdomain_id) + \
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

    hamiltonian_crc = vc[i]*h0_r[i,j]*uc[j]*dx(subdomain_id) - \
                      Dx(vc[i],0)*h1r_c[i,j]*uc[j]*dx(subdomain_id) + \
                      vc[i]*h1l_c[i,j]*Dx(uc[j],0)*dx(subdomain_id) + \
                      Dx(vc[i],0)*h2[i,j]*Dx(uc[j],0)*dx(subdomain_id)
    if hs:
      hamiltonian_crc += vc[i]*hs[i,j]*uc[j]*dx(subdomain_id)

    stiffness += hamiltonian_rrr - hamiltonian_rcc + \
      hamiltonian_ccr + hamiltonian_crc 
    
    if efield and efield != 0.0:
      stiffness += vr[i]*Expression('x[0] - xm', xm=xmid)* \
                   efield*ur[i]*dx(subdomain_id) + \
                   vc[i]*Expression('x[0] - xm', xm=xmid)* \
                   efield*uc[i]*dx(subdomain_id) 

    mass += ur[i]*vr[i]*dx(subdomain_id) + uc[i]*vc[i]*dx(subdomain_id)
  
  # Assemle matrices
  S = PETScMatrix()
  M = PETScMatrix()
  
  L = dot(Constant((0.0,) * num_bands), vr) * dx(subdomain_id) + \
        dot(Constant((0.0,) * num_bands), vc) * dx(subdomain_id)
  
  assemble_system(stiffness,L,bcs,A_tensor=S)
  assemble_system(mass,L,A_tensor=M)
  
  # Create eigensolver for the generalized EVP
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
  dofs_6 = VS.sub(6).dofmap().dofs()
  dofs_7 = VS.sub(7).dofmap().dofs()

  solution = SolutionContainer(dof_coordinates, \
                               space_dim, \
                               dofs_0, dofs_1, dofs_2, dofs_3, \
			       dofs_4, dofs_5, dofs_6, dofs_7)

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
    ef_complex = np.array([a + 1j * b for a, b in zip(u_re.vector().array(), \
                                             u_im.vector().array())])

    solution.add_eigenstate(ev_complex, ef_complex)
  
  # Add potential information to solution object
  solution.add_material_info(potential_to_materialinfo(subdomains, \
                             vpot, efield, xmid), 'vband')
  solution.add_material_info(potential_to_materialinfo(subdomains, \
                             cpot, efield, xmid), 'cband')
  
  if epsilon:
    # If epsilon is 
    eps_para = [epsilon[i][0] for i in range(len(epsilon))]
    solution.add_material_info(property_to_materialinfo(subdomains, \
                               eps_para), 'eps_para')
  
    eps_perp = [epsilon[i][1] for i in range(len(epsilon))]
    solution.add_material_info(property_to_materialinfo(subdomains, \
                               eps_perp), 'eps_perp')
   
  return solution
