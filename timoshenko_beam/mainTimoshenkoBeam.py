import numpy as np
import scipy.linalg as LA
from scipy import sparse
import meshio 
import matplotlib
from matplotlib import  pyplot as plt
import sys
import time
import os

import utilities

np.set_printoptions(precision=2, linewidth=300, threshold=sys.maxsize)


def shapeFunctionL2(xi):
  shape = np.array([ 0.5*(1-xi), 0.5*(1+xi)]).T
  naturalDerivatives = 0.5* np.array([-1,1]);
  
  return shape, naturalDerivatives

def formStiffnessMassTimoshenkoBeam(GDof,numberElements, elementNodes,numberNodes,xx,C,P,rho,I,thickness):
  # global
  stiffness = np.zeros( (GDof, GDof))
  force = np.zeros( (GDof,1) )

  # 2x2 Gauss quadrature
  gaussLocations = np.array([0.577350269189626 , -0.577350269189626])
  gaussWeights = np.ones((2))

  ########### Bending part ###################
  for e in range(numberElements):
    # индексы текущего элемента
    indice = elementNodes[e,:]
    # степени свободы текущего элемента
    elementDof = np.concatenate((indice,  indice+numberNodes), axis=0)

    #print(elementDof)

    ndof = indice.size # плохое название!!! Причем тут степени свободы?? Скорее количество узлов в элементе
    length_element = xx[ int(indice[1]) ]-xx[ int(indice[0]) ]
    detJacobian = length_element/2
    invJacobian=1/detJacobian

    for q in range(gaussWeights.size):
      # coordinate
      pt = gaussLocations[q]
      shape,naturalDerivatives = shapeFunctionL2( pt )
      Xderivatives = naturalDerivatives*invJacobian
      # deformation matrix
      B = np.zeros( (2,2*ndof))
      B[0, ndof:2*ndof] = Xderivatives.T

      add = B.T @ B *gaussWeights[q]*detJacobian*C[0,0]
      #print(add)

      for i in range(elementDof.size):
        for j in range(elementDof.size):
          stiffness[ int(elementDof[i]), int(elementDof[j]) ] += add[i,j]

      #print(stiffness)

      for i in range(indice.size):
        force[int(indice[i])] += shape[i]*P*detJacobian*gaussWeights[q]

  ############ Shear part ##########################

  ############ Reduced integration #################
  gaussLocations = 0.
  gaussWeights = 2.
  for e in range(numberElements):
    # индексы текущего элемента
    indice = elementNodes[e,:]
    # степени свободы текущего элемента
    elementDof = np.concatenate((indice,  indice+numberNodes), axis=0)
    ndof = indice.size # плохое название!!!
    length_element = xx[ int(indice[1]) ]-xx[ int(indice[0]) ]
    detJacobian = length_element/2
    invJacobian=1/detJacobian

    # coordinate
    pt = gaussLocations
    shape, naturalDerivatives = shapeFunctionL2( pt )
    Xderivatives = naturalDerivatives*invJacobian
    # deformation matrix
    B = np.zeros( (2,2*ndof))
    B[1, 0:ndof] = Xderivatives.T
    B[1, ndof:2*ndof] = shape

    add = B.T @ B *gaussWeights*detJacobian*C[1,1]
    #print(add)

    for i in range(elementDof.size):
        for j in range(elementDof.size):
          stiffness[ int(elementDof[i]), int(elementDof[j]) ] += add[i,j]

  return stiffness, force

def computeStress(solution, GDof,numberElements, elementNodes,numberNodes,xx,C,P,rho,I,thickness):
  # global
  solution = np.concatenate((solution[0],  solution[1]), axis=0)
  bend_stress = np.zeros( (numberElements,1) )  # вычисляем в узлах
  shear_stress = np.zeros( (numberElements,1) )  # вычисляем в узлах

  ########### Bending part ###################
  for e in range(numberElements):
    # индексы текущего элемента
    indice = elementNodes[e,:]
    # степени свободы текущего элемента
    elementDof = np.concatenate((indice,  indice+numberNodes), axis=0)
    ndof = indice.size # плохое название!!! Причем тут степени свободы?? Скорее количество узлов в элементе
    length_element = xx[ int(indice[1]) ]-xx[ int(indice[0]) ]
    detJacobian = length_element/2
    invJacobian=1/detJacobian

    shape,naturalDerivatives = shapeFunctionL2( 0.0 )
    Xderivatives = naturalDerivatives * invJacobian
    # deformation matrix
    B = np.zeros( (1,2*ndof))
    B[0, ndof:2*ndof] = Xderivatives.T
    #print(B)
    #print(solution[ elementDof.astype(int) ])
    bend_stress[int(e)] = E* thickness/2 * B @ solution[ elementDof.astype(int) ]

  for e in range(numberElements):
    # индексы текущего элемента
    indice = elementNodes[e,:]
    # степени свободы текущего элемента
    elementDof = np.concatenate((indice,  indice+numberNodes), axis=0)
    ndof = indice.size # плохое название!!!
    length_element = xx[ int(indice[1]) ]-xx[ int(indice[0]) ]
    detJacobian = length_element/2
    invJacobian=1/detJacobian

    shape, naturalDerivatives = shapeFunctionL2( 0.0 )
    Xderivatives = naturalDerivatives * invJacobian
    # deformation matrix
    B = np.zeros( (1,2*ndof))
    B[0, 0:ndof] = Xderivatives.T
    B[0, ndof:2*ndof] = shape
    shear_stress[int(e)] = C[1,1]/ thickness * B @ solution[ elementDof.astype(int) ]
  return bend_stress, shear_stress

def solution(GDof,prescribedDof,stiffness,force):
  # activeDof = np.setdiff1d(np.arange(GDof), prescribedDof)
  # Exclude prescribed DOFs, solve the system
  for i in prescribedDof:
    stiffness[i,:] = 0.0
    stiffness[i, i] = 1.0
    force[i] = 0.

  start = time.time()
  #disp = LA.solve(stiffness,force)
  disp = sparse.linalg.spsolve(stiffness,force)
  stop = time.time()
  print("Время расчёта:", stop-start)  
  disp_1 = disp[: int(GDof/2)]
  disp_2 = disp[int(GDof/2) : ]
  return [disp_1, disp_2]



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

GDof = 2*numberNodes
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

E = 1e8
poisson = 0.30
L = 1
thickness = 0.1
I = thickness**3/12;
EI = E*I
kapa = 5./6
P = -1
G = E/2/(1+poisson)
C = np.array([ [EI, 0], [0, kapa*thickness*G ]])

stiffness,force = formStiffnessMassTimoshenkoBeam(GDof,numberElements, LineElementNodesIds,numberNodes,xx, C, P, 1, I, thickness)

#print(stiffness)
#print(force)

#BCs
prescribedU1= np.array(LeftNodeIds)
prescribedRot1 = np.array(LeftNodeIds) + numberNodes
prescribedU2= np.array(RightNodeIds)
prescribedRot2 = np.array(RightNodeIds) + numberNodes
prescribedDof = [prescribedU1,  prescribedU2]

#solution
displacements = solution(GDof,prescribedDof,stiffness,force)
sigma_z, sigma_xz = computeStress(displacements, GDof,numberElements, LineElementNodesIds,numberNodes, xx, C, P, 1, I, thickness)

s_bend = np.reshape( sigma_z, (sigma_z.size, 1))
s_shear = np.reshape( sigma_xz, (sigma_xz.size, 1))
null_vertex_data = np.zeros((mesh.cell_data_dict['CellEntityIds']['vertex'].size, 1))

# SAVE RESULTS
res_mesh = meshio.Mesh(
    mesh.points,
    mesh.cells_dict,
    # Optionally provide extra data on points, cells, etc.
    point_data={ "deflection": np.array(displacements[0]).T, "rotation": np.array(displacements[1]).T },
    cell_data={ "s_bending": [null_vertex_data, s_bend], \
                "s_shear":   [null_vertex_data, s_shear]}
)


res_mesh.write(os.path.join(dirname,"results/resultsTB.vtk"),  # str, os.PathLike, or buffer/open file
    # file_format="vtk",  # optional if first argument is a path; inferred from extension
)

from matplotlib import  pyplot as plt
plt.scatter(nodeCoordinates, displacements[0], label ="Num")
plt.legend()
plt.show()