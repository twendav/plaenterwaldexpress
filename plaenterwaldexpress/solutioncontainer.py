import matplotlib.pyplot as pypl
import numpy as np
import h5py

class SolutionContainer:

  def __init__(self, dof_coords, space_dim, *args):
    self.material_info = {} 
    self.solver_info = []
    self.eigenstates = []
    self.space_dim = space_dim
    self.num_dofs = len(dof_coords) / self.space_dim

    # Prepare object for required solution dimension
    num_args = len(args)
    if len(args) > 0:
      self.solution_dim = num_args
      self.indices = args
    else:
      self.solution_dim = 1
      self.indices = [np.arange(self.num_dofs)]

    # Rearrange coordinates
    if self.space_dim > 1:
      dof_coords.reshape((self.num_dofs, self.space_dim))
    self.dof_coords = dof_coords
    
    # Create sorted master key for coordinates
    sorted_dof_coords = self.sort_dof_coords(dof_coords[self.indices[0]])
    self.dof_coords = sorted_dof_coords

    # Create lookup table for mapping of ... to ...
    self.dof_lookup = []
    for i in range(self.solution_dim):
      component_dofs = dof_coords[self.indices[i]]
      dof_dict = self.create_dof_dict(component_dofs)
      self.dof_lookup.append(self.create_dof_lookup(self.dof_coords, dof_dict))

  @staticmethod
  def create_dof_lookup(sorted_dofs, dof_dict):
    # Get space dimensions
    dof_shape = sorted_dofs.shape
    if len(dof_shape) == 1:
      space_dim = 1
      num_dofs = dof_shape[0]
    else:
      space_dim = dof_shape[0]
      num_dofs = dof_shape[1]

    # Create mapping from unsorted coords to sorted
    dof_lookup = []
    for i in range(num_dofs):
      coord = sorted_dofs[i]
      if space_dim > 1:
        coord = tuple(coord)
      
      position = dof_dict[coord]
      dof_lookup.append(position)

    return dof_lookup

  @staticmethod
  def sort_dof_coords(dof_coords):
    # Get space dimensions
    dof_shape = dof_coords.shape
    if len(dof_shape) == 1:
      space_dim = 1
    else:
      space_dim = dof_shape[0]

    # Sort list of tuples
    if space_dim == 1:
      sorted_dof_coords = np.sort(dof_coords)
    else:
      # Build list of keys for sorting
      lex_key = []
      for i in range(space_dim):
        lex_key.append(dof_coords[:,i])
      
      # Get indices of sorted array
      sorted_index = np.lexsort(lex_key)

      sorted_dof_coords = dof_coords[sorted_index]

    return sorted_dof_coords

  @staticmethod
  def create_dof_dict(dof_coords):
    # Get space dimension
    dof_shape = dof_coords.shape
    if len(dof_shape) == 1:
      space_dim = 1
      num_dofs = dof_shape[0]
    else:
      space_dim = dof_shape[0]
      num_dofs = dof_shape[1]

    # Create dict from list
    dof_dict = {}
    for i in range(num_dofs):
      if space_dim == 1:
        dof_dict[dof_coords[i]] = i
      else:
        dof_dict[tuple(dof_coords[i])] = i
    
    return dof_dict

  def add_eigenstate(self, eigen_nrg, eigen_fct):
    # Create empty eigenstate instance
    new_eigenstate = \
      Eigenstate(eigen_nrg, self.num_dofs, self.solution_dim)

    # Add eigenfunction componentwise
    for i in range(self.solution_dim):
      # Get eigenfunction for component i
      component = eigen_fct[self.indices[i]]
      
      # Order according to sorted (master) coordinates
      sorted_component = component[self.dof_lookup[i]]
      
      # Add component to new eigenstate
      new_eigenstate.add_component(sorted_component, i) 

    # Add eigenstate to list
    self.eigenstates.append(new_eigenstate)

  def __iter__(self):
    self.index = 0 
    return self

  def get_num_eigenstates(self):
    return len(self.eigenstates)

  def get_eigenstate(self, index):
    return self.eigenstates[index]

  def next(self):
    if self.index >= len(self.eigenstates):
      raise StopIteration

    self.index = self.index + 1  
    return self.eigenstates[self.index - 1]

  def add_material_info(self, material_property, name):
    self.material_info[name] = material_property

  def get_material_info(self, name):
    return self.material_info[name]

  def add_solver_info(self, solver):
    self.solver_info = solver

  def get_solver_info(self):
    return self.solver_info

#  def plot_eigenfunction(self, index, component=None):
#    eigenstate = self.get_eigenstate(index)
#
#    if index == None:
#      # Summation over all components at each node before plotting
#      sum_components = eigenstate.eigenfunction.sum(0)
#      plot_function = np.multiply(sum_components.conjugate(), \
#                             sum_components)
#    else:
#      plot_function = eigenstate.eigenfunction[index].conjugate() * \
#                             eigenstate.eigenfunction[index]
#
#    pypl.figure()
#    # Only plot real part to suppress numpy warning
#    pypl.plot(self.coordinates, plot_function.real)
#
#  def show_plots(self):
#    pypl.show()
  
  def export_solution(self, filename):
    # Open new hdf5 file
    fid = h5py.File(filename, "w")
    
    # Save dofs
    dset_data = self.dof_coords
    fid.create_dataset('dof_coords', data=dset_data)

    # Save eigenvalues and eigenfunctions
    for i in range(self.get_num_eigenstates()):
      # Create group for data
      grp_name = 'eigenstate' + str(i)
      gid = fid.create_group(grp_name)

      # Create dataset for eigenfunction
      dset_data = self.get_eigenstate(i).eigenfunction.transpose()
      dset = gid.create_dataset('eigfct_re', data=dset_data.real)
      dset = gid.create_dataset('eigfct_im', data=dset_data.imag)
       
      # Create dataset for eigenvalue
      dset_data = self.get_eigenstate(i).eigenvalue
      dset = gid.create_dataset('eigval_re', data=dset_data.real)
      dset = gid.create_dataset('eigval_im', data=dset_data.imag)

    # Attach material information 
    for property_name in self.material_info:
      # Create group for data
      grp_name = property_name
      gid = fid.create_group(grp_name)

      # Create dataset for coordinates
      dset_data = self.material_info[property_name].coordinates
      dset = gid.create_dataset('coords', data=dset_data)

      # Create dataset for function values 
      dset_data = self.material_info[property_name].values
      dset = gid.create_dataset('values', data=dset_data)
 
    # Attach meta-data to hdf5 file
    dset_data = self.solver_info.dolfin_version
    fid.create_dataset('dolfin_version', data=dset_data)

    dset_data = self.solver_info.plaenterwald_version
    fid.create_dataset('plaenterwald_version', data=dset_data)

    # Close file
    fid.close()

class Eigenstate:
  
  def __init__(self, eigenvalue, num_dofs, solution_dim, data_type=complex):
    self.eigenvalue = eigenvalue
    self.solution_dim = solution_dim
    self.num_coords = num_dofs / self.solution_dim
    self.eigenfunction = np.zeros((self.solution_dim, self.num_coords), \
                                  dtype=data_type)

  def add_eigenvalue(self, eigenvalue):
    self.eigenvalue = eigenvalue

  def add_component(self, values, index=0):
    self.eigenfunction[index] = values

  def get_eigenvalue(self):
    return self.eigenvalue
  
  def get_eigenfunction(self, index=0):
    return self.eigenfunction[index]



