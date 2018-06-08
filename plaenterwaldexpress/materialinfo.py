class MaterialInfo:
  """Container for material property as a function of space coordinate."""  
  def __init__(self, coords, values):
    self.coordinates = coords
    self.values = values

def potential_to_materialinfo(subdomains, potential, efield, middle_coord):
  # If efield is set to None, set to zero
  if efield == None:
    efield = 0.0
  
  # Get mesh coordinates
  subd_array = subdomains.array()
  mesh = subdomains.mesh()
  mesh_coord = mesh.coordinates().tolist()

  # Loop through mesh coordinates and calculate potential
  potential_coord = []
  potential_values = []
  band = [0, ]
  
  # For coordinate somewhere in the middle
  for i in range(len(subd_array)):
    # Get corresponding coordinates
    coord_before = mesh_coord[i][0]
    coord_after = mesh_coord[i + 1][0]

    # Get current and previous subdomain id
    subd_id_current = subd_array[i]
    try:
      subd_id_before = subd_array[i - 1]
    except:
      # Only for first coordinate
      subd_id_before = -1

    if subd_id_before != subd_id_current:
      # Subdomain interface: add two values for same coordinate
      psi = efield * (coord_before - middle_coord)
      potential[subd_id_current].eval(band, [coord_before, ])
      potential_values.append(band[0] + psi)
      potential_coord.append(coord_before)
    
    psi = efield * (coord_after - middle_coord)
    potential[subd_id_current].eval(band, [coord_after, ])
    potential_values.append(band[0] + psi)
    potential_coord.append(coord_after)
  
  return MaterialInfo(potential_coord, potential_values)

def property_to_materialinfo(subdomains, function):
  # Get mesh coordinates
  subd_array = subdomains.array()
  mesh = subdomains.mesh()
  mesh_coord = mesh.coordinates().tolist()

  # Loop through mesh coordinates and calculate potential
  function_coord = []
  function_values = []
  fvalue = [0, ]
  
  # For coordinate somewhere in the middle
  for i in range(len(subd_array)):
    # Get corresponding coordinates
    coord_before = mesh_coord[i][0]
    coord_after = mesh_coord[i + 1][0]

    # Get current and previous subdomain id
    subd_id_current = subd_array[i]
    try:
      subd_id_before = subd_array[i - 1]
    except:
      # Only for first coordinate
      subd_id_before = -1

    if subd_id_before != subd_id_current:
      # Subdomain interface: add two values for same coordinate
      function[subd_id_current].eval(fvalue, [coord_before, ])
      function_values.append(fvalue[0])
      function_coord.append(coord_before)

    function[subd_id_current].eval(fvalue, [coord_after, ])
    function_values.append(fvalue[0])
    function_coord.append(coord_after)

  return MaterialInfo(function_coord, function_values)

