from Code.Max_k_cut_classical_functions import *
from Code.Max_k_cut_quantum_functions import *
from itertools import combinations


def make_cost_block_MIS(n, layer):
    """
    :param n: The number of nodes
    :param layer: The subscript to be given to the variable gamma, the variable would be named to 'gamma_i'
    :return: The QuantumCircuit object corresponding to the mixing block
    """
    globals()['gamma%s' % layer] = Parameter('\u03B3' + str(layer))
    cost_circ = QuantumCircuit(n)
    for q in range(0, n):
        cost_circ.rz((-1) * globals()['gamma%s' % layer], q)
    return cost_circ


def make_mixing_block_MIS(n, w, layer):
    """
    :param n: The number of nodes
    :param w: The adjacency matrix of the graph
    :param layer: The subscript to be given to the variable beta, the variable would be named to 'beta_i'
    :return: The QuantumCircuit object corresponding to the cost block
    """
    globals()['beta%s' % layer] = Parameter('\u03B2' + str(layer))
    mix_circ = QuantumCircuit(n)
    for node in range(n):  # node corresponding to term in Hc expression
        adj_nodes = [i for i in range(n) if w[node][i] == 1]
        mix_circ.h(node)  # apply the starting h gate
        for s in powerset(adj_nodes):  # generate a power set so that we can perform cx-rz-cx, cx-cx-rz-
            # cx-cx, etc. sub-blocks one by one
            seq = list(s)  # get list of all the qubits involved for a sub-block (s is a tuple)
            seq.sort()  # we want to apply cx gates in order
            seq += [node]
            for i in range(len(seq) - 1):
                mix_circ.cx(seq[i], seq[i + 1])  # apply cx gates in ascending order of index
            mix_circ.rz(globals()['beta%s' % layer] / 2**(len(s)-1), seq[-1])  # apply the rz gate
            for i in range(len(seq) - 1):
                mix_circ.cx(seq[-(i + 2)], seq[-(i + 1)])  # apply cx gates in descending order of index
            mix_circ.barrier()
        mix_circ.h(node)  # apply the ending h gate
        #mix_circ.barrier()
    return mix_circ


def make_full_circuit_MIS(n, w, p):
    """
    :param n: The number of nodes
    :param w: The adjacency matrix of the graph
    :param p: The number of layers
    :return: The QuantumCircuit object corresponding to the full circuit
    """
    circ = QuantumCircuit(n)
    for layer in range(p):
        circ.append(make_cost_block_MIS(n, layer), [i for i in range(n)])
        circ.append(make_mixing_block_MIS(n, w, layer), [i for i in range(n)])
    circ.measure_all()
    return circ


def compute_cost_MIS(counts, w, U, n_counts=512):
    """
    :param counts:  dict{measurement result in string:count}
    :param l: The number of qubits representing a node
    :param w: The adjacency matrix for edges
    :param U: The size of punishment added to the cost due to including nodes of the same group in the MIS
    :return:  The averaged cost
    """
    total_cost = 0
    for measurement, count in counts.items():
        preprocessed_chosen_set = np.argwhere(np.array(list(measurement))=='1')
        if len(preprocessed_chosen_set) == 0:
            continue
        chosen_set = np.concatenate(preprocessed_chosen_set)
        total_cost += -1 * len(chosen_set) * count
        if len(chosen_set) < 2:
            continue
        else:
            for edge in combinations(chosen_set, 2):
                if w[edge[0]][edge[1]] == 1:
                    total_cost += U * count
    average_cost = total_cost / n_counts
    return average_cost


def func_to_optimize_wrapper_MIS(circ, w, U, nshots=512, simulator='qasm_simulator'):
    """
    :param circ:  The QuantumCircuit object corresponding to the full circuit, without feeding parameters or transpiling
    :param w: The adjacency matrix for edges
    :param nshots:  The number of shots
    :param simulator:  The simulator
    :return:  The function to optimize, corresponding to the circuit, to be fed to scipy.optimization.minimize
    """
    backend = Aer.get_backend(simulator)
    backend.shots = nshots

    def func_to_optimize(param_list):
         circ_w_param = circ.bind_parameters(param_list)
         transpiled_circ = transpile(circ_w_param, backend)
         counts = backend.run(transpiled_circ, shots=nshots).result().get_counts()
         return compute_cost(counts, w, U, nshots)

    return func_to_optimize