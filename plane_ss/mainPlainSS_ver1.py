import numpy as np
import scipy.linalg as LA
from scipy import sparse
import meshio 
import matplotlib
from matplotlib import  pyplot as plt
import sys
import time
import os
np.set_printoptions(precision=2, linewidth=300, threshold=sys.maxsize)

#  Квадратура Гаусса
def gaussQuadrature(type):
    if type == "complete":
        locations = np.array([ [-0.577350269189626, -0.577350269189626], [0.577350269189626, -0.577350269189626], [ 0.577350269189626, 0.577350269189626], [-0.577350269189626, 0.577350269189626] ])
        weights = np.array([1,1,1,1])

    elif type == "third":
        locations = np.array([ [-0.774596669241483, -0.774596669241483], \
                               [0.,  -0.774596669241483],  \
                               [ 0.774596669241483, -0.774596669241483],\
                                [-0.774596669241483, 0.], \
                                [0., 0.],
                                [0.774596669241483, 0.], \
                                [-0.774596669241483, 0.774596669241483],
                                [0. ,  0.774596669241483],
                                [0.774596669241483, 0.774596669241483]    ])
        
        weights = np.array([0.555555555555556*0.555555555555556, 0.555555555555556*0.888888888888889, \
                            0.555555555555556*0.555555555555556, 0.888888888888889*0.555555555555556, \
                            0.888888888888889*0.888888888888889, 0.555555555555556*0.888888888888889, \
                            0.555555555555556*0.555555555555556, 0.555555555555556*0.888888888888889, \
                            0.555555555555556*0.555555555555556 ]) 
    
    elif type == "reduced":
        locations = np.array([ [0.,0.]])
        weights = np.array([4.])
    
    else:
        print("Non-defined quadrature type")

    return weights, locations

def shapeFunctionL2(xi, eta):
    shape = 1./4* np.array([ xi*eta*(xi-1)*(eta-1), xi*eta*(xi+1)*(eta-1), xi*eta*(xi+1)*(eta+1), xi*eta*(xi-1)*(eta+1),\
                           -2*eta*(xi*xi-1)*(eta-1),  -2*xi*(xi+1)*(eta*eta-1), -2*eta*(xi*xi-1)*(eta+1), -2*xi*(xi-1)*(eta*eta-1), 4*(xi*xi-1)*(eta*eta-1) ]).T
  
    naturalDerivatives = 1./4 * np.array( [ [eta*(2*xi-1)*(eta-1),  xi*(xi-1)*(2*eta-1)],  \
                                           [eta*(2*xi+1)*(eta-1), xi*(xi+1)*(2*eta-1)], \
                                            [eta*(2*xi+1)*(eta+1),  xi*(xi+1)*(2*eta+1)], \
                                                 [eta*(2*xi-1)*(eta+1), xi*(xi-1)*(2*eta+1)], \
                                                      [-4*xi*eta*(eta-1),  -2*(xi+1)*(xi-1)*(2*eta-1)], \
                                                          [-2*(2*xi+1)*(eta+1)*(eta-1),   -4*xi*eta*(xi+1)], \
                                                              [-4*xi*eta*(eta+1),  -2*(xi+1)*(xi-1)*(2*eta+1)], \
                                                                  [-2*(2*xi-1)*(eta+1)*(eta-1),   -4*xi*eta*(xi-1) ], \
                                                                      [8*xi*(eta**2-1), 8*eta*(xi**2-1)] ])
  
    return shape, naturalDerivatives

def shapeFuncLineLinear(xi):
    shape = np.array([ xi, 1- xi]).T
    naturalDerivatives = np.array( [ [1,  -1] ])
    return shape, naturalDerivatives

def Jacobian(nodeCoords, naturalDerivatives):

    JacobianMatrix = nodeCoords.T @ naturalDerivatives
    invJacobian = LA.inv(JacobianMatrix)
    XYDerivatives = naturalDerivatives @ invJacobian

    return JacobianMatrix, invJacobian, XYDerivatives

def formStiffnessAxisymmetric2D(GDof,numberElements, QuadElementNodesIds, numberNodes, nodeCoords, C, th, quadType):
    # global
    # TO DO: implement sparse matrix realization
    stiffness = np.zeros( (GDof, GDof))
    
    # 2x2 Gauss quadrature
    gaussWeights, gaussLocations = gaussQuadrature(quadType)

    for e in range(numberElements):
        # индексы текущего элемента
        indice = np.array(QuadElementNodesIds[e,:])
        # трансформируем индексы в матрицу
        # степени свободы текущего элемента
        elementDof = np.concatenate((indice,  indice+numberNodes), axis=0)
        #print(elementDof)
        ndof = indice.size # количество узлов в элементе

        for q in range(gaussWeights.size):
            # coordinate
            pt = gaussLocations[q]
            xi = pt[0]
            eta =  pt[1]

            shapeFunction, naturalDerivatives = shapeFunctionL2( xi, eta )
            Jacob, invJacobian, XYDerivatives = Jacobian( nodeCoords[indice], naturalDerivatives )

            #B deformation matrix
            B = np.zeros( (3,2*ndof))
            B[0, 0:ndof] = XYDerivatives[:,0].T
            B[1, ndof:2*ndof] = XYDerivatives[:,1].T
            B[2, 0:ndof] = XYDerivatives[:,1].T
            B[2, ndof:2*ndof] = XYDerivatives[:,0].T

            add = B.T @ C @ B *gaussWeights[q]*LA.det(Jacob)
            #print(add)
            # сформировали локальную МЖ

            #
            for i in range(elementDof.size):
                for j in range(elementDof.size):
                    stiffness[ int(elementDof[i]), int(elementDof[j]) ] += add[i,j]
   
    return stiffness

def solution(GDof,prescribedDof,stiffness,force):
  # Exclude prescribed DOFs, solve the system
  for i in prescribedDof:
    stiffness[i,:] = 0.0
    stiffness[i, i] = 1.0
    force[i] = 0.
  #spStiffness = sparse.csr_matrix(stiffness)
  #spForce = sparse.csr_matrix(force)
  #print("Число обусловленнности:", np.linalg.cond(stiffness) )
  start = time.time()
  #disp = LA.solve(stiffness,force)
  disp = sparse.linalg.spsolve(stiffness,force)
  stop = time.time()
  
  print("Время расчёта:", stop-start)
  disp_1 = disp[: int(GDof/2)]
  disp_2 = disp[int(GDof/2) : ]
  return [disp_1, disp_2]


# временная функция для отображения напряжений в середине элемента через CellData
def computeStressAtMiddle(solution, GDof, numberElements, QuadElementNodesId,numberNodes,nodeCoords , C):
  # global
  solution = np.concatenate((solution[0],  solution[1]), axis=0)
  stress_tensor = np.zeros( (numberElements, 3) )  # вычисляем в узлах
  quadType = "reduced"
  # 2x2 Gauss quadrature
  gaussWeights, gaussLocations = gaussQuadrature(quadType)
  
  for e in range(numberElements):
    # индексы текущего элемента
    indice = np.array(QuadElementNodesIds[e,:])
    # степени свободы текущего элемента
    elementDof = np.concatenate((indice,  indice+numberNodes), axis=0)
    ndof = indice.size

    #  заходим ровно один раз
    for q in range(gaussWeights.size):
        # coordinate
        pt = gaussLocations[q]
        xi = pt[0]
        eta =  pt[1]
        shapeFunction, naturalDerivatives = shapeFunctionL2( xi, eta )
        Jacob, invJacobian, XYDerivatives = Jacobian( nodeCoords[indice], naturalDerivatives )
        # deformation matrix
        #B deformation matrix
        B = np.zeros( (3,2*ndof))
        B[0, 0:ndof] = XYDerivatives[:,0].T
        B[1, ndof:2*ndof] = XYDerivatives[:,1].T
        B[2, 0:ndof] = XYDerivatives[:,1].T
        B[2, ndof:2*ndof] = XYDerivatives[:,0].T
        stress_tensor[int(e)] = C @ B @ solution[ elementDof.astype(int) ]
  
  return stress_tensor

# функция позволяет отрисовывать quad9 элементы и постоянные значения на элементе
# есть опция отрисовки физических групп
def showMeshPlot(nodes, elements, values, bcs, loads):

    y = nodes[:,0]
    z = nodes[:,1]

    elements = np.delete(elements, [4,5,6,7,8] , axis=1)
    # cut quad9 into quad4 

    def quatplot(y,z, quatrangles, values, ax=None, **kwargs):

        if not ax: ax=plt.gca()
        yz = np.c_[y,z]
        verts= yz[quatrangles]
        #print(verts[0])
        pc = matplotlib.collections.PolyCollection(verts, **kwargs)
        pc.set_array(values)
        ax.add_collection(pc)
        ax.autoscale()
        return pc

    fig, ax = plt.subplots()
    ax.set_aspect('equal')

    pc = quatplot(y,z, np.asarray(elements), values, ax=ax, 
             edgecolor="crimson", cmap="rainbow")
    fig.colorbar(pc, ax=ax)        
    ax.plot(y,z, marker="o", markersize= 2, ls ="", color="black")
    for k in bcs:
        ax.plot(y[k], z[k], marker='x', markersize= 10, ls ="", color='red')
    for k in loads:
        ax.plot(y[k], z[k], marker='d', markersize= 10, ls ="", color='blue')

    ax.set(title='Mesh', xlabel='R axis', ylabel='Z axis')

    plt.show()

# plot deform solution
def showSolution(nodes, elements, values, displacements):
    y = nodes[:,0]
    z = nodes[:,1]

    elements = np.delete(elements, [4,5,6,7,8] , axis=1)
    # cut quad9 into quad4 

    def quatplot(y,z, quatrangles, values, ax=None, **kwargs):
        if not ax: ax=plt.gca()
        yz = np.c_[y,z]
        verts= yz[quatrangles]
        pc = matplotlib.collections.PolyCollection(verts, **kwargs)
        pc.set_array(values)
        ax.add_collection(pc)
        ax.autoscale()
        return pc

    fig, ax = plt.subplots()
    ax.set_aspect('equal')

    pc = quatplot(y,z, np.asarray(elements), values, ax=ax, 
             edgecolor="crimson", cmap="rainbow")
    fig.colorbar(pc, ax=ax)        
    
    ax.plot(y,z, marker="o", markersize= 2, ls ="", color="black")

    scale = 0.1
    plot =ax.plot(y+ scale*displacements[0],z + scale*displacements[1], marker="o", markersize= 2, ls ="", color="purple")
    #plot = ax.contour(y,z, scale*displacements[0])
    pc1 = quatplot(y+ scale*displacements[0],z+ scale*displacements[1], np.asarray(elements), values, ax=ax, cmap="Greys", alpha= 0.2)
    fig.colorbar(pc1, ax=ax)

    ax.set(title='Mesh', xlabel='R axis', ylabel='Z axis')
    #fig.colorbar(plot, ax=ax)
    plt.show()

# Идем по элементу и возвращаем вектор узлов входящих в него
# Нужно для передачи нагрузок и закреплений по узлам
def extractNodeIdsFromCells( CellNodes ):
    NodeIds = []
    for Cell in CellNodes:
        for Node in Cell:
            NodeIds.append(Node)
    NodesIdsUnique = list(set(NodeIds))
    return NodesIdsUnique

def apply_traction_load_on_line(f_vector, LineElemNodeIds, nodeCoords, traction_vector):
    # В каждом элементе два отрезка, для каждого вычисляем криволинейный интеграл и добавляем в вектор правой части
    numberNodes = nodeCoords[0,:].size
    for line in LineElemNodeIds:

        left, right, mid = line
        x1, y1  = nodeCoords[left]
        x2, y2 = nodeCoords[right]
        xmid, ymid = nodeCoords[mid]

        ds1 = np.sqrt((xmid-x1)**2 + (ymid-y1)**2)
        cos1 = (x1-xmid)/ds1
        sin1 = (y1-ymid)/ds1

        ds2 = np.sqrt((x2-xmid)**2 + (y2-ymid)**2)
        cos2 = -(x2-xmid)/ds2
        sin2 = -(y2-ymid)/ds2

        f_vector[left]  += 2*ds1/3*traction_vector[0]
        f_vector[right] += 2*ds1/3*traction_vector[0]
        f_vector[mid]   += 2*ds1/3*traction_vector[0]

        f_vector[left + numberNodes] += 2*ds1/3*traction_vector[1]
        f_vector[right + numberNodes] += 2*ds1/3*traction_vector[1]
        f_vector[mid + numberNodes] += 2*ds1/3*traction_vector[1]
    
    return



# GEOMETRY + MESH GENERATION

# TO DO: gmsh+meshio support
# Quad9

# Читаем сетку в vtk формате, внимательно следим за номерами, которые присвоены физическим группам (по умолчанию -1)
dirname = os.path.dirname(__file__)
mesh = meshio.read(os.path.join(dirname, 'mesh/rectangle_w_c_struct_r0.1.vtk'))   # string, os.PathLike, or a buffer/open file)

xx = mesh.points[:, 0]
yy = mesh.points[:, 1]
nodeCoords = np.column_stack((xx, yy))
numberNodes = xx.size

#print(nodeCoords)
QuadElementNodesIds = mesh.cells_dict['quad9']
AllNodesQuadMesh = np.array( list(set(QuadElementNodesIds.reshape(-1))) )
QuadNodesRenumbered = np.arange( AllNodesQuadMesh.size )

#print(QuadElementNodesIds)
#print(max(AllNodesQuadMesh))

numberElements = np.size(mesh.cells_dict['quad9'], 0)
LineElementNodesIds = mesh.cells_dict['line3']
QuadCellData = mesh.cell_data_dict['CellEntityIds']['quad9']
QuadCellData = np.reshape(QuadCellData, len(QuadCellData)) # reshape to values for plotting

#print(mesh.cell_data_dict)

# Создаём лист индексов для задания краевых условий
# TO DO: автоматизировать создаваемые объекты для записи в списки loads, bcs
LeftCellsIds  = []
RightCellsIds = []
LoadCellsIds  = []


LineLoadElements = []
# Идём по элементам типа line3 и выбираем те, значение которых равно нужным нам иднексам физическиз групп
for i in range( len(mesh.cell_data_dict['CellEntityIds']['line3']) ):
    if mesh.cell_data_dict['CellEntityIds']['line3'][i] == 5: #!!!
        LeftCellsIds.append(i)
    if mesh.cell_data_dict['CellEntityIds']['line3'][i] == 1: #!!! 
        LoadCellsIds.append(i)
        LineLoadElements.append(LineElementNodesIds[i])
    if mesh.cell_data_dict['CellEntityIds']['line3'][i] == 7: #!!! 
        RightCellsIds.append(i)



LeftNodeIds  = np.array(extractNodeIdsFromCells(LineElementNodesIds[LeftCellsIds])) 
RightNodeIds = np.array(extractNodeIdsFromCells(LineElementNodesIds[RightCellsIds]))
LoadNodeIds =  np.array(extractNodeIdsFromCells(LineElementNodesIds[LoadCellsIds])) 

bcs = [LeftNodeIds] #, RightNodeIds]
loads = [LoadNodeIds]

# MATERIAL INPUT CARD

# 1 global material only supported
# TO DO: property class support

# Material elastic input
E = 1
rho = 1
th = 1
poisson = 0.30
C = th*E/(1-poisson**2)*np.array([ [1, poisson, 0], [poisson, 1, 0 ], [0 ,0 , (1-poisson)/2] ])
#C = th*E*np.array( [ [1, poisson, 0], [poisson, 1, 0 ], [0 ,0 , 1./(2*(1+poisson))] ] )

# LOAD CARD
pres = 1

# PLOT BLOCK
# drawing mesh with physical groups

showMeshPlot(nodeCoords, QuadElementNodesIds, QuadCellData, bcs, loads)

# ASSEMBLE BLOCK
GDof = 2*numberNodes
print(GDof)

stiffness = formStiffnessAxisymmetric2D(GDof, numberElements, QuadElementNodesIds, numberNodes, nodeCoords, C, th, "third")
#print(stiffness)

prescribedDofUr = np.array(LeftNodeIds)
prescribedDofSymmetry = np.array(RightNodeIds) + numberNodes 
prescribedDofUz = np.array(LeftNodeIds) + numberNodes 
prescribedDofUz2 = np.array(LoadNodeIds) + numberNodes
prescribedDof = [prescribedDofUr, prescribedDofSymmetry, 3, 3 + numberNodes ]
#print(prescribedDof)

force = np.zeros(GDof)
apply_traction_load_on_line(force, LineLoadElements, nodeCoords, np.array([0, pres]))

#force[ np.asarray(LoadNodeIds)] = -pres
#force[8] = 0
#force[8 + numberNodes ] = 0
#print(LeftNodeIds)
#print(LoadNodeIds)
print(sum(force))
#force[ np.asarray(LoadNodeIds) + numberNodes ] = pres
#print(prescribedDofUz)

#SOLVE
displacements = solution(GDof,prescribedDof,stiffness,force)
#print(displacements)

showSolution(nodeCoords, QuadElementNodesIds, QuadCellData, displacements)
stress_tensor = computeStressAtMiddle(displacements, GDof, numberElements, QuadElementNodesIds, numberNodes,nodeCoords , C)
#print("stress tensor",  mesh.cell_data_dict['CellEntityIds']['vertex'], stress_tensor)
sx = np.reshape( stress_tensor[:, 0], (stress_tensor[:, 0].size, 1))
sy = np.reshape( stress_tensor[:, 1], (stress_tensor[:, 1].size, 1))
sxy = np.reshape(stress_tensor[:, 2], (stress_tensor[:, 2].size, 1))
null_vertex_data = np.zeros((mesh.cell_data_dict['CellEntityIds']['vertex'].size, 1))
null_line_data = np.zeros((mesh.cell_data_dict['CellEntityIds']['line3'].size, 1))
print("maximum stresses:", max(sx))

# SAVE RESULTS
res_mesh = meshio.Mesh(
    mesh.points,
    mesh.cells_dict,
    # Optionally provide extra data on points, cells, etc.
    point_data={ "disp": np.array(displacements).T},
    cell_data={ "sx": [null_vertex_data, null_line_data, sx], \
               "sy": [null_vertex_data, null_line_data, sy], \
                      "sxy": [null_vertex_data, null_line_data, sxy]}
    # Each item in cell data must match the cells array
)
res_mesh.write(os.path.join(dirname,"results/resultsPlaneSS.vtk"),  # str, os.PathLike, or buffer/open file
    # file_format="vtk",  # optional if first argument is a path; inferred from extension
)