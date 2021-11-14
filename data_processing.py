import networkx as nx
import random
import matplotlib.pyplot as plt
import numpy as np

def generate_data(num_normal_vehicles, num_important_vehicles):
    """Generate vehicle data

    Args:
        num_normal_vehicles (int): Number of normal vehicles
        num_important_vehicles (int): Number of priorized vehicles

    Returns:
        nx.Graph: Returns a fully connected graph with the random assigned weights
    """
    total_number_of_vehicles = num_normal_vehicles + num_important_vehicles
    # 
    vehicles = {}
    # Time for charging in minutes
    min_charging = 15
    max_charging = 360
    # Create vehicles
    for i in range(total_number_of_vehicles):
        vehicles[str(i)] = {}
        if i < num_normal_vehicles:
            vehicles[str(i)]["weight"] = random.uniform(0.1,5) 
            vehicles[str(i)]["charging_time"] = random.randint(min_charging, max_charging)
        else:
            vehicles[str(i)]["weight"] = random.uniform(5,10)
            vehicles[str(i)]["charging_time"] = random.randint(min_charging, max_charging)
    #print(vehicles)

    G = nx.complete_graph(total_number_of_vehicles)
    for i in range(total_number_of_vehicles):
        for j in range(i):
            weight_1 = round(vehicles[str(i)]["weight"]*vehicles[str(j)]["charging_time"],2)
            weight_2 = round(vehicles[str(i)]["charging_time"]*vehicles[str(j)]["weight"],2)
            weight = int(min(weight_1, weight_2))
            #print(j,i,weight, weight_1, weight_2)
            G.add_edge(i,j,weight=weight)

    return G

def get_weight_matrix(G):
    """
    :param G: nx.Graph: a fully connected graph with the random assigned weights
    :return: numpy array: a matrix containing weights of edges
    """
    matrix = []
    for i in range(len(G)):
        row = []
        for j in range(len(G)):
            if j != i:
                row.append(G._adj[i][j]['weight'])
            else:
                row.append(0)
        matrix.append(row)
    return np.array(matrix)

if __name__ == "__main__":
    G = generate_data(3,2)
    pos = nx.spring_layout(G)
    nx.draw(G, pos, with_labels=True)
    nx.draw_networkx_edge_labels(G, pos)
    plt.show()