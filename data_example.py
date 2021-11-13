from data_processing import generate_data
import networkx as nx
import matplotlib.pyplot as plt

G = generate_data(3,2)
pos = nx.spring_layout(G)
nx.draw(G, pos, with_labels=True)
nx.draw_networkx_edge_labels(G, pos)
plt.show()