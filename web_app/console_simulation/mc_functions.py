import random
from console_simulation.data_processing import *
# from console_simulation.Max_k_cut_classical_functions import *
from console_simulation.Max_k_cut_quantum_functions import *
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

# Set parameters here

p = 3  # the depth of QAOA

l = 2  # the number of qubits for each node, k = 2^l where k is the number of subgraphs we'd like to cut,
       # or to say, the k in 'max-k-cut'. l = 2 for a max-4-cut example

init_params = [np.pi/8, np.pi]*p

local_optimization_method = 'Powell'

k = 2**l

def init_vehicles(n):
    d={"vehicles": [{"importance": round(0.1+0.9*random.random(), 3), "charging_time": round(20*(0.1+0.9*random.random()), 2)} for _ in range(n)]}
    return d

def get_mc(jobs):
    weights=[]
    for i in range(len(jobs)):
        row=[]
        for j in range(len(jobs)):
            row+=[min(float(jobs[i]["importance"])*float(jobs[j]["charging_time"]), float(jobs[j]["importance"])*float(jobs[i]["charging_time"]))]
        weights+=[row]
    R = compute_normalization_scale(weights)
    rescaled_weights = weights / R
    print("get_mc called")
    circ = make_full_circuit(5,l,rescaled_weights,p)
    simulator='qasm_simulator'
    backend = Aer.get_backend(simulator)
    backend.shots = 128
    circ_w_param = circ.bind_parameters([ 2.95085482e+00, -1.25155070e-01, -4.59385172e-03,  7.15978099e+00, 1.03754388e+01,  2.32943974e+00])
    transpiled_circ = transpile(circ_w_param, backend)
    counts = backend.run(transpiled_circ, shots=128).result().get_counts()
    measurement=max(counts, key=counts.get)
    return [int(measurement[i:i + l], 2) for i in range(0, len(measurement), l)]
