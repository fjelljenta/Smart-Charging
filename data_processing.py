import networkx as nx
import random
import matplotlib.pyplot as plt

def generate_data(num_normal_vehicles, num_important_vehicles):
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

if __name__ == "__main__":
    G = generate_data(3,2)
    pos = nx.spring_layout(G)
    nx.draw(G, pos, with_labels=True)
    nx.draw_networkx_edge_labels(G, pos)
    plt.show()