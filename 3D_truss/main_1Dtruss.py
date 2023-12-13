# from MATLAB codes for Finite Element Analysis
# problem2.m
import numpy as np
import scipy.linalg as LA
import meshio 
import matplotlib
from matplotlib import  pyplot as plt
import sys
import time
import os

import utilities1D as utilities
np.set_printoptions(precision=2, linewidth=300, threshold=sys.maxsize)
# GEOMETRY + MESH READING
# Читаем сетку в vtk формате, внимательно следим за номерами, которые присвоены физическим группам (по умолчанию -1)
dirname = os.path.dirname(__file__)
mesh = meshio.read(os.path.join(dirname, 'mesh/beam.vtk'))   # string, os.PathLike, or a buffer/open file)

xx = mesh.points[: , 0]
yy = mesh.points[: , 1]
nodeCoordinates = np.array(xx)
numberNodes = xx.size
print(numberNodes)

# 2 point elements are used!
LineElementNodesIds = mesh.cells_dict['line']
VertexElementNodesIds = mesh.cells_dict['vertex']
LineCellData = mesh.cell_data_dict['CellEntityIds']['line']
VertexCellData = mesh.cell_data_dict['CellEntityIds']['vertex']

GDof = numberNodes
numberElements = np.size(LineElementNodesIds,0)
print("Степеней свободы:", GDof)
print("Кол-во элементов:", numberElements)
print("Кол-во узлов:", numberNodes)

# Создаём лист индексов для задания краевых условий
# TO DO: автоматизировать создаваемые объекты для записи в списки loads, bcs
LeftCellsIds  = []
RightCellsIds = []
LoadCellsIds  = []

LineLoadElements = []
# Идём по элементам типа line3 и выбираем те, значение которых равно нужным нам иднексам физическиз групп
for i in range( len(mesh.cell_data_dict['CellEntityIds']['line']) ):
    if mesh.cell_data_dict['CellEntityIds']['line'][i] == 3: #!!! 
        LoadCellsIds.append(i)
        LineLoadElements.append(LineElementNodesIds[i])


for i in range( len(VertexCellData ) ):
    if VertexCellData [i] == 1: #!!! 
      LeftCellsIds.append(i)
    if VertexCellData [i] == 2: #!!! 
      RightCellsIds.append(i)

LeftNodeIds  = np.array(utilities.extractNodeIdsFromCells(VertexElementNodesIds[LeftCellsIds])) 
RightNodeIds = np.array(utilities.extractNodeIdsFromCells(VertexElementNodesIds[RightCellsIds]))
LoadNodeIds =  np.array(utilities.extractNodeIdsFromCells(LineElementNodesIds[LoadCellsIds])) 

bcs = [LeftNodeIds, RightNodeIds]
loads = [LoadNodeIds]

print(mesh.cell_data_dict['CellEntityIds']['vertex'])
print(mesh.cells_dict['vertex'])
print(bcs, loads)


# MATERIAL PARAMETERS
E = 1; A = 1.; EA = E*A

# AXIAL DISTRIBUTED LOAD
P = 0.

displacements = np.zeros(numberNodes)
force = np.zeros(numberNodes)
stiffness = np.zeros([numberNodes, numberNodes])
    
stiffness, force = utilities.formStiffnessBar2pt(GDof,numberElements, LineElementNodesIds,numberNodes, nodeCoordinates, xx, EA, P)

#print(stiffness)
#print(force)

#BCs
prescribedU1= np.array(LeftNodeIds)
prescribedU2= np.array(RightNodeIds)
print(prescribedU1, prescribedU2)
prescribedDof = [prescribedU1]
pointLoad = [ [prescribedU2, 1] ]

displacements = utilities.solution(GDof, prescribedDof, pointLoad, stiffness, force)
print("displacements:", displacements)

axial_stress = utilities.computeStress(displacements, GDof,numberElements, LineElementNodesIds,numberNodes, nodeCoordinates, xx, EA)
print("stresses:", axial_stress)

s_axial = np.reshape( axial_stress, (axial_stress.size, 1))
null_vertex_data = np.zeros((mesh.cell_data_dict['CellEntityIds']['vertex'].size, 1))

# SAVE RESULTS
res_mesh = meshio.Mesh(
    mesh.points,
    mesh.cells_dict,
    # Optionally provide extra data on points, cells, etc.
    point_data={ "disp": displacements.T},
    cell_data={ "s_axial": [null_vertex_data, s_axial]}
)

res_mesh.write(os.path.join(dirname,"results/results_1D_Bar.vtk"),  # str, os.PathLike, or buffer/open file
    # file_format="vtk",  # optional if first argument is a path; inferred from extension
)

'''
def exact_solution(x):
    u = p*L*x / (2*E*A)*(1 - x / L)
    s_x =  p*L/A*(0.5 - x/L)
    return u, s_x
'''
'''
x_ = np.linspace(0, L, 100)
ex_ans = exact_solution(x_)
plt.plot(x_, ex_ans[0], label = "Ex")
plt.plot(nodeCoordinates, displacements, label ="Num")
plt.scatter(nodeCoordinates, displacements)
plt.legend()
plt.show()
'''