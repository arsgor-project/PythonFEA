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
    elementDof = np.concatenate((indice,  indice+numberNodes), axis=0)
    x = nodeCoordinates[indice][:,0]
    y = nodeCoordinates[indice][:,1]
    length_element = np.sqrt((x[0]-x[1])**2 + (y[0]-y[1])**2)
    
    l = (x[1]-x[0])/length_element # cos theta
    m = (y[1]-y[0])/length_element # sin theta

    L = np.array([[l,0,m,0], [0,l,0,m]])
    #print(L)
    ndof = indice.size #
    detJacobian = length_element/2
    invJacobian = 1./detJacobian

    for q in range(gaussWeights.size):
      # coordinate
      pt = gaussLocations[q]
      shape,naturalDerivatives = shapeFunctionL2( pt )
      Xderivatives = naturalDerivatives*invJacobian
      # deformation matrix
      B = Xderivatives.copy()
      add =  gaussWeights[q]* detJacobian * EA * L.T @ (B.T@ B) @ L
      #print(add)
      for i in range(elementDof.size):
        for j in range(elementDof.size):
          stiffness[ int(elementDof[i]), int(elementDof[j]) ] += add[i,j]
      
      fadd = shape.T @ L *P*detJacobian*gaussWeights[q]
      # print(fadd)
      # axial distributed load
      for i in range(indice.size):
        force[int(indice[i])] += fadd[0,i]

  return stiffness, force


def solution(GDof, prescribedDof, pointLoad,  stiffness, force):
    
  # activeDof = np.setdiff1d(np.arange(GDof), prescribedDof)
  # Exclude prescribed DOFs, solve the system
  for i in prescribedDof:
    stiffness[i,:] = 0.0
    stiffness[i, i] = 1.0
    force[i] = 0.
  
  for i in pointLoad:
    
    dofs = i[0]
    #print(dofs)
    value = i[1]
    force[dofs] += value
  
  
  # check for zero rows
  for i in range(GDof):
    is_all_zero = np.all((stiffness[i,:] == 0))
    if is_all_zero:
       stiffness[i,i] = 1
       force[i] = 0

  #print(force)
  #print(stiffness)
  start = time.time()
  #disp = LA.solve(stiffness,force)
  disp = sparse.linalg.spsolve(stiffness,force)
  stop = time.time()
  print("Время расчёта:", stop-start)  
  disp_1 = disp[: int(GDof/2)]
  disp_2 = disp[int(GDof/2) : ]

  
  return [disp_1, disp_2]



def computeStress(solution, GDof,numberElements, elementNodes,numberNodes, nodeCoordinates, xx,EA):
  # global
  solution = np.concatenate((solution[0],  solution[1]), axis=0)
  axial_stress = np.zeros( (numberElements,1) )  # вычисляем в узлах

  for e in range(numberElements):
    # индексы текущего элемента
    indice = elementNodes[e,:]
    elementDof = np.concatenate((indice,  indice+numberNodes), axis=0)
    x = nodeCoordinates[indice][:,0]
    y = nodeCoordinates[indice][:,1]
    length_element = np.sqrt((x[0]-x[1])**2 + (y[0]-y[1])**2)
    
    l = (x[1]-x[0])/length_element # cos theta
    m = (y[1]-y[0])/length_element # sin theta

    L = np.array([[l,0,m,0], [0,l,0,m]])
    ndof = indice.size # 
    
    detJacobian = length_element/2
    invJacobian = 1./detJacobian
    shape,naturalDerivatives = shapeFunctionL2( 0.0 )
    Xderivatives = naturalDerivatives * invJacobian
    # deformation matrix
    B = Xderivatives.copy()
    axial_stress[int(e)] = EA * B @ L @ solution[ elementDof.astype(int) ]

  return axial_stress
