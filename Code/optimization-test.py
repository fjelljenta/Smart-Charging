from pyscipopt import Model, quicksum, multidict

def gpp(V,E):
    """gpp -- model for the graph partitioning problem
    Parameters:
        - V: set/list of nodes in the graph
        - E: set/list of edges in the graph
    Returns a model, ready to be solved.
    """
    model = Model("gpp")
    x = {}
    y = {}
    for i in V:
        x[i] = model.addVar(vtype="B", name="x(%s)"%i)
    for (i,j) in E:
        y[i,j] = model.addVar(vtype="B", name="y(%s,%s)"%(i,j))

    model.addCons(quicksum(x[i] for i in V) == len(V)/2, "Partition")

    for (i,j) in E:
        model.addCons(x[i] - x[j] <= y[i,j], "Edge(%s,%s)"%(i,j))
        model.addCons(x[j] - x[i] <= y[i,j], "Edge(%s,%s)"%(j,i))

    model.setObjective(quicksum(y[i,j] for (i,j) in E), "minimize")
    model.data = x
    return model

model = gpp([0,1,2],[(0,1),(1,2)])
model.optimize()
print(model.getObjVal())