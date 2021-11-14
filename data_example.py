from data_processing import generate_data, generate_data_mis
from networkx.algorithms.approximation import one_exchange, randomized_partitioning
from Max_k_cut_functions import brut_force_k4_N5
import networkx as nx
import matplotlib.pyplot as plt

"""
G = generate_data(3,2)
for i in range(5):
    cut_value, partition = one_exchange(G)
    print(cut_value)
    print(partition)
for i in range(5):
    cut_value, partition = randomized_partitioning(G)
    print(cut_value)
    print(partition)
pos = nx.spring_layout(G)
c_opt, p = brut_force_k4_N5(G)
print(c_opt, p)
nx.draw(G, pos, with_labels=True)
nx.draw_networkx_edge_labels(G, pos)
#plt.show()
"""

G = generate_data_mis(4,10,500,100)
pos = nx.spring_layout(G)
nx.draw(G, pos, with_labels=True)
plt.show()