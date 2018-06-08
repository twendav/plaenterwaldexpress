from dolfin import MeshEditor
from dolfin import Mesh

import numpy as np

def mesh_gen_1d(layer_thickness, elem_per_layer):
  editor = MeshEditor()
  mesh = Mesh()

  editor.open(mesh, 1, 1)  

  num_cells = sum(elem_per_layer)
  num_vertices = num_cells + 1

  editor.init_vertices(num_vertices)  # number of vertices
  editor.init_cells(num_cells)     # number of cells

  editor.add_vertex(0, np.array([0.0]))
  vertex_index = 0

  # Add vertices
  layer_startpos = np.cumsum([0] + layer_thickness)
  
  for i in range(len(layer_thickness)):
    dist = layer_thickness[i] / elem_per_layer[i]
    offset = layer_startpos[i]

    for k in range(elem_per_layer[i]):
      vertex_index += 1
      editor.add_vertex(vertex_index, \
                        np.array([offset + dist * (k + 1)]))

  # Add cells
  for i in range(num_cells):
    editor.add_cell(i, np.array([i, i + 1], dtype=np.uintp))
    
  editor.close()

  return mesh

