# from MATLAB codes for Finite Element Analysis
# 2 noded EB beam element

import numpy as np
import scipy.linalg as LA
import meshio 
import matplotlib
from matplotlib import  pyplot as plt
import sys
import time
import os
import math
import utilities as utilities
np.set_printoptions(precision=2, linewidth=300, threshold=sys.maxsize)

'''
nodes = np.array([[0,0,0], [0,0,1], [1,1,1], [1,1,0]])
elems = np.array([[0,1],[1,2],[2,3]])
def_beam_vec = np.array([0,0,1]) # задает ориентацию n1, сонаправленную с указанным вектором

for index,i in enumerate(elems):
    a = i[0]
    b = i[1]
    x1,y1,z1 = nodes[a]
    x2,y2,z2 = nodes[b]
    length = math.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
    tangent = np.array( [x2-x1, y2-y1, z2-z1])/length

    cross = np.cross(tangent, def_beam_vec) 
    n2 = cross/ np.linalg.norm(cross)
    n1 = np.cross( n2, tangent) 
    print(length)
    print(tangent, n1, n2)

# class material
# class beam_property
# class cross_section
'''


# GEOMETRY + MESH READING
# Читаем сетку в vtk формате, внимательно следим за номерами, которые присвоены физическим группам (по умолчанию -1)
dirname = os.path.dirname(__file__)
mesh = meshio.read(os.path.join(dirname, 'mesh/2el_test.vtk'))   # string, os.PathLike, or a buffer/open file)

xx = mesh.points[: , 0]
yy = mesh.points[: , 1]
zz = mesh.points[: , 2]

nodeCoordinates = np.column_stack((xx, yy, zz))
numberNodes = xx.size

# 2 point elements are used!
LineElementNodesIds = mesh.cells_dict['line']
VertexElementNodesIds = mesh.cells_dict['vertex']
LineCellData = mesh.cell_data_dict['CellEntityIds']['line']
VertexCellData = mesh.cell_data_dict['CellEntityIds']['vertex']

GDof = numberNodes*6
numberElements = np.size(LineElementNodesIds,0)
print("Степеней свободы:", GDof)
print("Кол-во элементов:", numberElements)
print("Кол-во узлов:", numberNodes)


LeftCellsIds  = []
RightCellsIds = []
LoadCellsIds  = []
LoadPointsCellsIds = []
LineLoadElements = []
# Идём по элементам типа line3 и выбираем те, значение которых равно нужным нам иднексам физическиз групп
for i in range( len(mesh.cell_data_dict['CellEntityIds']['line']) ):
    if mesh.cell_data_dict['CellEntityIds']['line'][i] == 1: #!!! 
        LoadCellsIds.append(i)
        LineLoadElements.append(LineElementNodesIds[i])


for i in range( len(VertexCellData ) ):
    if VertexCellData [i] == 2: #!!! 
      LeftCellsIds.append(i)
    if VertexCellData [i] == 3: #!!! 
      RightCellsIds.append(i)
    if VertexCellData [i] == 1: #!!! 
      LoadPointsCellsIds.append(i)

#print(LoadPointsCellsIds)
LeftNodeIds  = np.array(utilities.extractNodeIdsFromCells(VertexElementNodesIds[LeftCellsIds])) 
RightNodeIds = np.array(utilities.extractNodeIdsFromCells(VertexElementNodesIds[RightCellsIds]))
LoadPointsIds = np.array(utilities.extractNodeIdsFromCells(VertexElementNodesIds[LoadPointsCellsIds]))
LoadNodeIds =  np.array(utilities.extractNodeIdsFromCells(LineElementNodesIds[LoadCellsIds])) 

bcs = [LeftNodeIds]
loads = [LoadPointsIds]

print(bcs, loads)

# MATERIAL PARAMETERS
E = 1; A = 1.; EA = E*A;  Iy = 2; Iz = 2; G = E/(2*(1+0.3)); J = 5; 
mat_data = np.array([E, G, A, Iy, Iz, J])
stiffness, force = utilities.formStiffness3D_EBbeam_2pt(GDof,numberElements, LineElementNodesIds,numberNodes, nodeCoordinates, mat_data)

#BCs
prescribedUx1= np.array(LeftNodeIds)
prescribedUy1= np.array(LeftNodeIds) + numberNodes
prescribedUz1= np.array(LeftNodeIds) + 2*numberNodes
prescribedRx1= np.array(LeftNodeIds)+  3*numberNodes
prescribedRy1= np.array(LeftNodeIds) + 4*numberNodes
prescribedRz1= np.array(LeftNodeIds) + 5*numberNodes

prescribedDof = [prescribedUx1, prescribedUy1,   prescribedUz1, prescribedRx1, prescribedRy1, prescribedRz1  ]
print(prescribedDof)
prescribedPoint1Fz = np.array(LoadPointsIds)  + 2*numberNodes
pointLoad = [ [prescribedPoint1Fz, 1] ]

displacements = utilities.solution(GDof, prescribedDof, pointLoad, stiffness, force)
