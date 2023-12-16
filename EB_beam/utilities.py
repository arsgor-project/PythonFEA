import numpy as np
import time
from scipy import sparse
import math

def extractNodeIdsFromCells( CellNodes ):
    NodeIds = []
    for Cell in CellNodes:
        for Node in Cell:
            NodeIds.append(Node)
    NodesIdsUnique = list(set(NodeIds))
    return NodesIdsUnique

def formStiffness3D_EBbeam_2pt(GDof,numberElements, elementNodes ,numberNodes, nodeCoordinates, mat_data):
    # global
    stiffness = np.zeros( (GDof, GDof))
    force = np.zeros( (GDof,1) )
    
    E, G, A, Iy, Iz, J = mat_data
    
    
    for e in range(numberElements):
    #for e in range(1):
        # индексы текущего элемента
        indice = elementNodes[e,:]
        # степени свободы текущего элемента
        elementDof = np.concatenate((indice,  indice+numberNodes, indice + 2*numberNodes, \
            indice + 3*numberNodes, indice + 4*numberNodes, indice + 5*numberNodes), axis=0)
        
        #elementDof= np.array([0,1,2,3,4,5,6,7,8,9,10,11])
        #print(elementDof)
        x = nodeCoordinates[indice][:,0]
        y = nodeCoordinates[indice][:,1]
        z = nodeCoordinates[indice][:,2]
        
        length_element = np.sqrt((x[0]-x[1])**2 + (y[0]-y[1])**2 + (z[0]-z[1])**2)
        L = length_element

        if (x[0] == x[1]) and (y[0] == y[1]):
            if (z[1]>z[0]):
                r = np.array([[0,0,1], [0,1,0], [-1,0,0]])
            else:
                r = np.array([[0,0,-1], [0,1,0], [1,0,0]])
        else:
            CXx = (x[1] - x[0])/L
            CYx = (y[1] - y[0])/L
            CZx = (z[1] - z[0])/L
            D = math.sqrt(CXx*CXx + CYx*CYx)
            CXy = -CYx/D 
            CYy = CXx/D 
            CZy = 0
            CXz =  -CXx*CZx/D
            CYz = -CYx*CZx/D
            CZz = D
            r = np.array([[CXx, CYx, CZx] , [ CXy, CYy , CZy] , [CXz,  CYz, CZz]])

        # матрица перехода
        #r = np.array([[1,1,1],[1,1,1],[1,1,1]])
        zero = np.zeros((6,6))

        first_block = np.array([[r[0,0], 0, r[0,1], 0, r[0,2], 0], \
                                [0, r[0,0], 0, r[0,1], 0, r[0,2]], \
                                [r[1,0], 0, r[1,1], 0, r[1,2], 0], \
                                [0, r[1,0], 0, r[1,1], 0,  r[1,2]], \
                                [r[2,0], 0, r[2,1], 0, r[2,2], 0], \
                                [0, r[2,0], 0, r[2,1], 0, r[2,2]]    ])
        
        R = np.block([[first_block, zero], [zero, first_block ] ])
        #print(r)

                
        k1 = E*A/L
        k2 = 12*E*Iz/(L*L*L)
        k3 = 6*E*Iz/(L*L)
        k4 = 4*E*Iz/L
        k5 = 2*E*Iz/L
        k6 = 12*E*Iy/L**3
        k7 = 6*E*Iy/L**2
        k8 = 4*E*Iy/L 
        k9 = 2*E*Iy/L 
        k10 = G*J/L

        '''                    
        LocalK = np.array( [[k1, 0, 0, 0, 0, 0, -k1, 0, 0, 0, 0, 0], \
                            [0, k2, 0, 0, 0, k3, 0, -k2, 0, 0, 0, k3], \
                            [0, 0, k6, 0, -k7, 0, 0, 0, -k6, 0, -k7, 0], \
                            [0, 0, 0, k10, 0, 0, 0, 0, 0, -k10, 0, 0], \
                            [0, 0, -k7, 0, k8, 0, 0, 0, k7, 0, k9, 0], \
                            [0, k3, 0, 0, 0, k4, 0, -k3, 0, 0, 0, k5], \
                            [-k1, 0, 0, 0, 0, 0, k1, 0, 0, 0, 0, 0], \
                            [0, -k2, 0, 0, 0, -k3, 0, k2, 0, 0, 0, -k3], \
                            [0, 0, -k6, 0, k7, 0, 0, 0, k6, 0, k7, 0], \
                            [0, 0, 0, -k10, 0, 0, 0, 0, 0, k10, 0, 0], \
                            [0, 0, -k7, 0, k9, 0, 0, 0, k7, 0, k8, 0], \
                            [0, k3, 0, 0, 0, k5, 0, -k3, 0, 0, 0, k4]  ])
        '''
        '''
        LocalK = np.array( [[k1, 0, 0, 0, 0, 0, -k1, 0, 0, 0, 0, 0], \
                            [-k1, 0, 0, 0, 0, 0, k1, 0, 0, 0, 0, 0], \
                            [0, k2, 0, 0, 0, k3, 0, -k2, 0, 0, 0, k3], \
                            [0, -k2, 0, 0, 0, -k3, 0, k2, 0, 0, 0, -k3], \
                            [0, 0, k6, 0, -k7, 0, 0, 0, -k6, 0, -k7, 0], \
                            [0, 0, -k6, 0, k7, 0, 0, 0, k6, 0, k7, 0], \
                            [0, 0, 0, k10, 0, 0, 0, 0, 0, -k10, 0, 0], \
                            [0, 0, 0, -k10, 0, 0, 0, 0, 0, k10, 0, 0], \
                            [0, 0, -k7, 0, k8, 0, 0, 0, k7, 0, k9, 0], \
                             [0, 0, -k7, 0, k9, 0, 0, 0, k7, 0, k8, 0], \
                            [0, k3, 0, 0, 0, k4, 0, -k3, 0, 0, 0, k5], \
                            [0, k3, 0, 0, 0, k5, 0, -k3, 0, 0, 0, k4]  ])
        '''
        LocalK = np.array( [[k1, -k1, 0, 0,    0,    0,    0, 0,     0, 0,   0, 0], \
                    [-k1, k1, 0, 0,    0,    0,    0, 0,     0, 0,   0, 0], \
                    [0, 0,    k2, -k2, 0,    0,    0, 0,     0, 0,   k3, k3], \
                    [0, 0,   -k2,  k2, 0,    0,    0, 0,     0, 0,  -k3, -k3], \
                    [0, 0,     0, 0,   k6, -k6,    0, 0,   -k7, -k7, 0, 0], \
                    [0, 0,     0, 0,  -k6,  k6,    0, 0,    k7, k7,  0, 0], \
                    [0, 0,     0, 0,    0,   0,  k10, -k10, 0,   0, 0, 0], \
                    [0, 0,     0, 0,    0,   0, -k10,  k10, 0,   0, 0, 0], \
                    [0, 0,     0, 0,   - k7,  k7,    0,    0, k8, k9, 0, 0], \
                    [0, 0,     0, 0,   - k7,  k7,    0,    0, k9, k8, 0, 0], \
                    [0, 0,    k3, -k3,    0,   0,    0,    0,  0,  0, k4, k5], \
                    [0, 0,    k3, -k3,    0,   0,    0,    0,  0,  0, k5, k4]  ])
                            

        add =  R.T @ (LocalK) @ R
        #print(add)
        for i in range(elementDof.size):
            for j in range(elementDof.size):
                stiffness[ int(elementDof[i]), int(elementDof[j]) ] += add[i,j]        
        #print(stiffness)
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
    #print(np.linalg.cond(stiffness))
    start = time.time()
    #disp = LA.solve(stiffness,force)
    disp = sparse.linalg.spsolve(stiffness,force)
    stop = time.time()
    print("Время расчёта:", stop-start)  
    disp_1 = disp[: int(GDof/6)]
    disp_2 = disp[int(GDof/6) : int(2*GDof/6)]
    disp_3 = disp[int(2*GDof/6) : int(3*GDof/6)]
    disp_4 = disp[int(3*GDof/6) : int(4*GDof/6)]
    disp_5 = disp[int(4*GDof/6) : int(5*GDof/6)]
    disp_6 = disp[int(5*GDof/6) : int(6*GDof/6)]
    
    return [disp_1, disp_2, disp_3, disp_4, disp_5, disp_6]