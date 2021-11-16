import networkx as nx
import random
import matplotlib.pyplot as plt
import numpy as np


def generate_data(num_normal_vehicles, num_important_vehicles):
    """Generate vehicle data. Use a number of normal and important vehicles to create a fully connected graph

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
        # TODO int weights
        if i < num_normal_vehicles:
            vehicles[str(i)]["weight"] = random.uniform(0.1, 5)
            vehicles[str(i)]["charging_time"] = random.randint(
                min_charging, max_charging)
        else:
            vehicles[str(i)]["weight"] = random.uniform(5, 10)
            vehicles[str(i)]["charging_time"] = random.randint(
                min_charging, max_charging)
    # print(vehicles)

    G = nx.complete_graph(total_number_of_vehicles)
    for i in range(total_number_of_vehicles):
        for j in range(i):
            weight_1 = round(vehicles[str(i)]["weight"]
                             * vehicles[str(j)]["charging_time"], 2)
            weight_2 = round(
                vehicles[str(i)]["charging_time"]*vehicles[str(j)]["weight"], 2)
            weight = int(min(weight_1, weight_2))
            #print(j,i,weight, weight_1, weight_2)
            G.add_edge(i, j, weight=weight)
    # TODO store graph data (json?)
    return G


def generate_data_mis(number_of_groups, number_of_vehicles, max_time, max_charging_time=360):
    """Generate data for the MIS. Specify the numbers of groups, the numbers of vehicles, the maximum time scale and the maximum charing time.

    Args:
        number_of_groups (int): Number of groups
        number_of_vehicles (int): Number of vehicles
        max_time (int): Start time charging limit (say, after 100 minutes, every car started at least charging)
        max_charging_time (int): maximum time a car can charge (like 6h)
    """
    G = nx.Graph()
    vehicles = {}
    edge_list = []
    for i in range(number_of_vehicles):
        vehicles[str(i)] = {}
        vehicles[str(i)]["group"] = random.randint(1, number_of_groups)
        vehicles[str(i)]["start_charging"] = random.randint(0, max_time)
        vehicles[str(i)]["end_charging"] = random.randint(vehicles[str(
            i)]["start_charging"], max_charging_time+vehicles[str(i)]["start_charging"])

    G.add_nodes_from(range(number_of_vehicles))
    for i in range(number_of_vehicles):
        for j in range(i):
            added_flag = False
            if vehicles[str(i)]["group"] == vehicles[str(j)]["group"]:
                edge_list.append((i, j))
            if not added_flag:
                if vehicles[str(i)]["start_charging"] <= vehicles[str(j)]["start_charging"] <= vehicles[str(i)]["end_charging"]:
                    edge_list.append((i, j))
                elif vehicles[str(i)]["start_charging"] <= vehicles[str(j)]["end_charging"] <= vehicles[str(i)]["end_charging"]:
                    edge_list.append((i, j))
                elif vehicles[str(j)]["start_charging"] <= vehicles[str(i)]["start_charging"] <= vehicles[str(j)]["end_charging"]:
                    edge_list.append((i, j))
                else:
                    # No connection
                    pass

    G.add_edges_from(edge_list)
    print(vehicles)
    print(edge_list)
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


def compute_normalization_scale(w):
    """
    a function to compute the normalization factor for the weight matrix given by the paper, R = max(w)*N^2/4
    :param w:  the original weight matrix
    :return:  the rescaled weight matrix
    """
    return np.max(w) * len(w)**2 / 4


if __name__ == "__main__":
    G = generate_data(3, 2)
    pos = nx.spring_layout(G)
    nx.draw(G, pos, with_labels=True)
    nx.draw_networkx_edge_labels(G, pos)
    plt.show()
