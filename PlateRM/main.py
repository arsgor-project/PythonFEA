import h5py
import numpy as np
import scipy.linalg as LA
import meshio 
import matplotlib
from matplotlib import  pyplot as plt
import sys
import time
import os
import math
np.set_printoptions(precision=2, linewidth=300, threshold=sys.maxsize)


class FEmodel:
    def __init__(self, mesh= None, loads = None, constraints = None, solver_params = None):
        self.mesh = mesh
        self.loads = loads
        self.constraints = constraints
        self.solver_params = solver_params


class Mesh:

    meshVTKisAssigned = False
    
    def __init__(self):
        pass
    
    def read_vtk(self, filename):
        self.vtkMesh = meshio.read(filename)
        self.xx = self.vtkMesh.points[:, 0]
        self.yy = self.vtkMesh.points[:, 1]
        self.zz = self.vtkMesh.points[:, 2]
        self.nodeCoords = np.column_stack((self.xx, self.yy, self.zz))
        self.numberNodes = self.xx.size
        self.cells = self.vtkMesh.cells_dict
        self.meshVTKisAssigned = True
        
        # типы поддерживаемых элементов
        self.vertexNodeIds = None
        self.line3NodeIds = None
        self.quad9NodeIds = None

        # Read cells that are supported
        for key in self.cells:
            if key == "vertex":
                self.vertexNodeIds = self.cells[key]
            elif key == "line3":
                self.line3NodeIds = self.cells[key]
            elif key == "quad9":
                self.quad9NodeIds =  self.cells[key]
            else:
                raise Exception("Неверно указан тип элемента или такие элементы не поддерживаются")  

    def instance_info(self):
        if self.meshVTKisAssigned:
            print("Количество узлов:", self.numberNodes)
            #print("Элементы:", self.cells)
            for key in self.cells:
                print(key, self.cells[key])
        else:
            raise Exception("Объекту не назначена сетка, используйте метод read_vtk")    

    @staticmethod
    def extractNodeIdsFromCells( CellNodes ):
        NodeIds = []
        for Cell in CellNodes:
            for Node in Cell:
                NodeIds.append(Node)
        NodesIdsUnique = list(set(NodeIds))
        return NodesIdsUnique


my_model = FEmodel()
my_mesh = Mesh()
dirname = os.path.dirname(__file__)

my_mesh.read_vtk(os.path.join(dirname, 'mesh/rectangle_quad_2x2.vtk'))
my_mesh.instance_info()


# Global Matrix Object
# Element Object
'''
ForEachElementInTheModelDo {
elementObj->ComputeStiffness( stiffmatrixObj );
globalmatrixObj->AssembleInMatrix( stiffmatrixObj );
}
'''

# How mesh is related to element and property?

'''

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

'''
