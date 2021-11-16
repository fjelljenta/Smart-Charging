from data_processing import generate_data, generate_data_from_file, generate_data_mis, generate_data_mis_from_file
from networkx.algorithms.approximation import one_exchange, randomized_partitioning
from Max_k_cut_functions import *
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
"""
G, passed_time = benchmarking(generate_data, 5, 5)
print(passed_time)
values,  passed_time = benchmarking(Monte_Carlo_solver, G, 4)
print(values, passed_time)
values, passed_time = benchmarking(brut_force, G, 4)
print(values, passed_time)
#benchmarking(generate_data_mis_from_file, "2021-11-16_03-15_mis.json")
