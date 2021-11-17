from data_processing import generate_data, generate_data_from_file, generate_data_mis, generate_data_mis_from_file
import networkx as nx
import matplotlib.pyplot as plt

"""
G = generate_data(3,2)
pos = nx.spring_layout(G)
nx.draw(G, pos, with_labels=True)
nx.draw_networkx_edge_labels(G, pos)
plt.show()


G = generate_data_from_file("2021-11-16_02-56.json")
pos = nx.spring_layout(G)
nx.draw(G, pos, with_labels=True)
nx.draw_networkx_edge_labels(G, pos)
plt.show()

G = generate_data_mis(4,10,500,100)
pos = nx.spring_layout(G)
nx.draw(G, pos, with_labels=True)
plt.show()

G = generate_data_mis_from_file("2021-11-16_03-15_mis.json")
pos = nx.spring_layout(G)
nx.draw(G, pos, with_labels=True)
plt.show()

G, passed_time = benchmarking(generate_data, 5, 5)
print(passed_time)
values,  passed_time = benchmarking(Monte_Carlo_solver, G, 4)
print(values, passed_time)
values, passed_time = benchmarking(brut_force, G, 4)
print(values, passed_time)
#benchmarking(generate_data_mis_from_file, "2021-11-16_03-15_mis.json")
"""

def graphSets(graph):
      
    # Base Case - Given Graph 
    # has no nodes
    if(len(graph) == 0):
        g = nx.Graph()
        return g
     
    # Base Case - Given Graph
    # has 1 node
    if(len(graph) == 1):
        g = nx.Graph()
        g.add_node(list(graph)[0])
        return g
      
    # Select a vertex from the graph
    vCurrent = list(graph)[0]
      
    # Case 1 - Proceed removing
    # the selected vertex
    # from the Maximal Set
    graph2 = graph.copy()
      
    # Delete current vertex 
    # from the Graph
    graph2.remove_node(vCurrent)
      
    # Recursive call - Gets 
    # Maximal Set,
    # assuming current Vertex 
    # not selected
    res1 = graphSets(graph2)
      
    # Case 2 - Proceed considering
    # the selected vertex as part
    # of the Maximal Set
    # Loop through its neighbours
    for v in list(graph.neighbors(vCurrent)):
          
        # Delete neighbor from 
        # the current subgraph
        if(v in graph2):
            graph2.remove_node(v)
      
    # This result set contains VFirst,
    # and the result of recursive
    # call assuming neighbors of vFirst
    # are not selected
    res2 = graphSets(graph2)
    res2.add_node(vCurrent)
    # Our final result is the one 
    # which is bigger, return it
    if(len(res1) > len(res2)):
        return res1
    return res2

G = generate_data_mis(10,50,500,100)
pos = nx.spring_layout(G)
nx.draw(G, pos, with_labels=True)
plt.show()

maximalIndependentSet = graphSets(G)
  
# Prints the Result 
for i in maximalIndependentSet:
    print(i, end =" ")