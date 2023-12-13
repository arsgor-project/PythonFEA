import numpy as np

# Идем по элементу и возвращаем вектор узлов входящих в него
# Нужно для передачи нагрузок и закреплений по узлам
def extractNodeIdsFromCells( CellNodes ):
    NodeIds = []
    for Cell in CellNodes:
        for Node in Cell:
            NodeIds.append(Node)
    NodesIdsUnique = list(set(NodeIds))
    return NodesIdsUnique
