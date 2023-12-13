import numpy as np
import time
from scipy import sparse

# Идем по элементу и возвращаем вектор узлов входящих в него
# Нужно для передачи нагрузок и закреплений по узлам
def extractNodeIdsFromCells( CellNodes ):
    NodeIds = []
    for Cell in CellNodes:
        for Node in Cell:
            NodeIds.append(Node)
    NodesIdsUnique = list(set(NodeIds))
    return NodesIdsUnique


def shapeFunctionL2(xi):
    my_shape = np.array([[(1-xi)/2,(1+xi)/2 ]])
    my_Xderivatives = np.array([[ -1./2 , 1./2 ]])
    return my_shape.T, my_Xderivatives





def solution(GDof, prescribedDof, pointLoad,  stiffness, force):
  # activeDof = np.setdiff1d(np.arange(GDof), prescribedDof)
  # Exclude prescribed DOFs, solve the system
  for i in prescribedDof:
    stiffness[i,:] = 0.0
    stiffness[i, i] = 1.0
    force[i] = 0.
  for i in pointLoad:
    dofs = i[0]
    value = i[1]
    force[dofs] += value
  
  print(force)
  start = time.time()
  #disp = LA.solve(stiffness,force)
  disp = sparse.linalg.spsolve(stiffness,force)
  stop = time.time()
  print("Время расчёта:", stop-start)  
  return disp


def formStiffnessBar2pt(GDof, numberElements, elementNodes, numberNodes, nodeCoordinates, xx,EA,P):
  # global
  stiffness = np.zeros( (GDof, GDof))
  force = np.zeros( (GDof,1) )

  # 2x2 Gauss quadrature
  gaussLocations = np.array([0])
  gaussWeights = np.array([2])

  for e in range(numberElements):
    # индексы текущего элемента
    indice = elementNodes[e,:]
    # степени свободы текущего элемента
    elementDof = np.array(indice)

    ndof = indice.size #
    length_element = abs(nodeCoordinates[elementDof[1]] - nodeCoordinates[elementDof[0]])
    detJacobian = length_element/2
    invJacobian = 1./detJacobian

    for q in range(gaussWeights.size):
      # coordinate
      pt = gaussLocations[q]
      shape,naturalDerivatives = shapeFunctionL2( pt )
      Xderivatives = naturalDerivatives*invJacobian
      # deformation matrix
      B = Xderivatives.copy()
      add =  gaussWeights[q]*(B.T@ B) * detJacobian * EA

      for i in range(elementDof.size):
        for j in range(elementDof.size):
          stiffness[ int(elementDof[i]), int(elementDof[j]) ] += add[i,j]
      
      # axial distributed load
      for i in range(indice.size):
        force[int(indice[i])] += shape[i]*P*detJacobian*gaussWeights[q]

  return stiffness, force



def computeStress(solution, GDof,numberElements, elementNodes,numberNodes, nodeCoordinates, xx,EA):
  # global
  solution = np.array(solution)
  axial_stress = np.zeros( (numberElements,1) )  # вычисляем в узлах

  for e in range(numberElements):
    # индексы текущего элемента
    indice = elementNodes[e,:]
    # степени свободы текущего элемента
    elementDof = np.array(indice)
    print(elementDof)
    ndof = indice.size # 
    length_element = abs(nodeCoordinates[elementDof[1]] - nodeCoordinates[elementDof[0]])
    detJacobian = length_element/2
    invJacobian=1/detJacobian

    shape,naturalDerivatives = shapeFunctionL2( 0.0 )
    Xderivatives = naturalDerivatives * invJacobian
    # deformation matrix
    B = Xderivatives.copy()
    axial_stress[int(e)] = EA * B @ solution[ elementDof.astype(int) ]

  return axial_stress
