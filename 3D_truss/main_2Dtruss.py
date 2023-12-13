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

import utilities2D as utilities
np.set_printoptions(precision=2, linewidth=300, threshold=sys.maxsize)

# GEOMETRY + MESH READING
# Читаем сетку в vtk формате, внимательно следим за номерами, которые присвоены физическим группам (по умолчанию -1)
dirname = os.path.dirname(__file__)
mesh = meshio.read(os.path.join(dirname, 'mesh/2DTruss.vtk'))   # string, os.PathLike, or a buffer/open file)

xx = mesh.points[: , 0]
yy = mesh.points[: , 1]
nodeCoordinates = np.column_stack((xx, yy))
numberNodes = xx.size

# 2 point elements are used!
LineElementNodesIds = mesh.cells_dict['line']
VertexElementNodesIds = mesh.cells_dict['vertex']
LineCellData = mesh.cell_data_dict['CellEntityIds']['line']
VertexCellData = mesh.cell_data_dict['CellEntityIds']['vertex']

GDof = numberNodes*2
numberElements = np.size(LineElementNodesIds,0)
print("Степеней свободы:", GDof)
print("Кол-во элементов:", numberElements)
print("Кол-во узлов:", numberNodes)


# Создаём лист индексов для задания краевых условий
# TO DO: автоматизировать создаваемые объекты для записи в списки loads, bcs
LeftCellsIds  = []
RightCellsIds = []
LoadCellsIds  = []
LoadPointsCellsIds = []

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
    if VertexCellData [i] == 3: #!!! 
      LoadPointsCellsIds.append(i)

#print(LoadPointsCellsIds)
LeftNodeIds  = np.array(utilities.extractNodeIdsFromCells(VertexElementNodesIds[LeftCellsIds])) 
RightNodeIds = np.array(utilities.extractNodeIdsFromCells(VertexElementNodesIds[RightCellsIds]))
LoadPointsIds = np.array(utilities.extractNodeIdsFromCells(VertexElementNodesIds[LoadPointsCellsIds]))
LoadNodeIds =  np.array(utilities.extractNodeIdsFromCells(LineElementNodesIds[LoadCellsIds])) 

bcs = [LeftNodeIds, RightNodeIds]
loads = [LoadNodeIds]

#print(mesh.cell_data_dict['CellEntityIds']['vertex'])
#print(mesh.cells_dict['vertex'])
#print(bcs, loads)


# MATERIAL PARAMETERS
E = 1; A = 1.; EA = E*A

# AXIAL DISTRIBUTED LOAD
P = 0.

stiffness, force = utilities.formStiffnessBar2pt(GDof,numberElements, LineElementNodesIds,numberNodes, nodeCoordinates, xx, EA, P)

#print(stiffness)
#print(force)

#BCs
prescribedUx1= np.array(LeftNodeIds)
prescribedUy1= np.array(LeftNodeIds) + numberNodes
prescribedUx2= np.array(RightNodeIds)
prescribedUy2= np.array(RightNodeIds) + numberNodes
prescribedDof = [prescribedUx1, prescribedUy1,   prescribedUy2]

prescribedPoint1Fx = np.array(LoadPointsIds)
prescribedPoint2Fy = np.array(LoadPointsIds)  + numberNodes
pointLoad = [ [prescribedPoint1Fx, 0], [prescribedPoint2Fy, 1] ]

#print(prescribedDof)
old_stiff = stiffness.copy()
displacements = utilities.solution(GDof, prescribedDof, pointLoad, stiffness, force)
print("displacements:", displacements)

#print(old_stiff)
print("Reaction:", old_stiff @ np.concatenate((displacements[0],  displacements[1]), axis=0))

axial_stress = utilities.computeStress(displacements, GDof,numberElements, LineElementNodesIds,numberNodes, nodeCoordinates, xx, EA)
print("stresses:", axial_stress)
s_axial = np.reshape( axial_stress, (axial_stress.size, 1))
null_vertex_data = np.zeros((mesh.cell_data_dict['CellEntityIds']['vertex'].size, 1))

# SAVE RESULTS
res_mesh = meshio.Mesh(
    mesh.points,
    mesh.cells_dict,
    # Optionally provide extra data on points, cells, etc.
    point_data={ "deflection": np.array(displacements).T},
    cell_data={ "s_axial": [null_vertex_data, s_axial]}
)

res_mesh.write(os.path.join(dirname,"results/results_2DTruss_ex1.vtk"),  # str, os.PathLike, or buffer/open file
    # file_format="vtk",  # optional if first argument is a path; inferred from extension
)