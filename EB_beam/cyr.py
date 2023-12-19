import h5py
import meshio
import numpy as np
import os

dirname = os.path.dirname(__file__)
points = [
    [0.0, 0.0],
    [1.0, 0.0],
    [0.0, 1.0],
    [1.0, 1.0],
    [2.0, 0.0],
    [2.0, 1.0],
]
cells = [
    ("triangle", [[0, 1, 2], [1, 3, 2]]),
    ("quad", [[1, 4, 5, 3]]),
]

mesh = meshio.Mesh(
    points,
    cells,
    # Optionally provide extra data on points, cells, etc.
    point_data={"T": [0.3, -1.2, 0.5, 0.7, 0.0, -3.0]},
    # Each item in cell data must match the cells array
    cell_data={"a": [[0.1, 0.2], [0.4]]},
)


print(mesh.cells)
filename = os.path.join(dirname, "results/test_motion.xdmf")
with meshio.xdmf.TimeSeriesWriter(filename) as writer:
    writer.write_points_cells(mesh.points,  mesh.cells_dict)
    
    for t in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]:
        writer.write_data(t, point_data={"disp": t*np.array([0.3, -1.2, 0.5, 0.7, 0.0, -3.0])})
