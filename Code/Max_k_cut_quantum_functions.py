from Code.Max_k_cut_classical_functions import *

from qiskit import QuantumCircuit
from qiskit import Aer
from qiskit.compiler import transpile
from qiskit.circuit import Parameter

from scipy.optimize import differential_evolution, minimize
from Code.data_processing import *
import numpy as np


def make_mixing_block(n, l, layer):
    """
    :param n: The number of nodes
    :param l: The number of qubits representing a node
    :param layer: The subscript to be given to the variable beta, the variable would be named to 'beta_layer'
    :return: The QuantumCircuit object corresponding to the mixing block
    """
    nq = n * l
    globals()['beta%s' % layer] = Parameter('\u03B2' + str(layer))
    mix_circ = QuantumCircuit(nq)
    for q in range(0, nq):
        mix_circ.rx(2 * globals()['beta%s' % layer], q)
    return mix_circ


def make_cost_block(n, l, w, layer):
    """
    :param n: The number of nodes
    :param l: The number of qubits representing a node
    :param w: The weight matrix for edges
    :param layer: The subscript to be given to the variable gamma, the variable would be named to 'gamma_layer'
    :return: The QuantumCircuit object corresponding to the cost block
    """
    nq = n * l
    globals()['gamma%s' % layer] = Parameter('\u03B3' + str(layer))
    cost_circ = QuantumCircuit(nq)
    for n1 in range(n):  # first node
        for n2 in range(n1 + 1, n):  # second node, which is always larger than first node by index to avoid repetition
            for s in powerset([i for i in range(l)]):  # generate a power set so that we can perform cx-rz-cx, cx-cx-rz-
                # cx-cx, etc. sub-blocks one by one
                if len(s) == 0:  # empty set, we could do a global phase rotation but unnecessary
                    pass
                else:
                    seq = [l * n1 + i for i in s] + [l * n2 + i for i in s]  # list out all the qubits involved for a
                    # sub-block
                    seq.sort()  # we want to apply cx gates in order
                    for i in range(len(seq) - 1):
                        cost_circ.cx(seq[i], seq[i + 1])  # apply cx gates in ascending order of index
                    cost_circ.rz(globals()['gamma%s' % layer] * w[n1][n2] / 2, seq[-1])  # apply the rz gate
                    for i in range(len(seq) - 1):
                        cost_circ.cx(seq[-(i + 2)], seq[-(i + 1)])  # apply cx gates in descending order of index
                #cost_circ.barrier()
    return cost_circ


def make_initial_block(n, l):
    """
    :param n: The number of nodes
    :param l: The number of qubits representing a node
    :return: The QuantumCircuit object corresponding to the initial state
    """
    nq = n * l
    init_circ = QuantumCircuit(nq)
    for i in range(0, nq):
        init_circ.h(i)
    return init_circ


def make_full_circuit(n, l, w, p):
    """
    :param n: The number of nodes
    :param l: The number of qubits representing a node
    :param w: The weight matrix for edges
    :param p: The number of layers
    :return: The QuantumCircuit object corresponding to the full circuit
    """
    nq = n * l
    circ = QuantumCircuit(nq)
    circ.append(make_initial_block(n, l), [i for i in range(nq)])
    for layer in range(p):
        circ.append(make_cost_block(n, l, w, layer), [i for i in range(nq)])
        circ.append(make_mixing_block(n, l, layer), [i for i in range(nq)])
    circ.measure_all()
    return circ


def run_circuit(circ, param_list, nshots=512, simulator='qasm_simulator'):
    """
    :param circ:  The QuantumCircuit object corresponding to the full circuit
    :param param_list:  A list containing the parameters for the circuit
    :param nshots:  The number of shots
    :param simulator:  The simulator
    :return:  counts = dict{measurement result in string:count},
              transpiled_circ = The QuantumCircuit object corresponding to the transpiled circuit
    """
    circ_w_param = circ.bind_parameters(param_list)
    backend = Aer.get_backend(simulator)
    transpiled_circ = transpile(circ_w_param, backend)
    counts = backend.run(transpiled_circ, shots=nshots).result().get_counts()
    return counts, transpiled_circ



def func_to_optimize_wrapper(circ, l, w, nshots=512, simulator='qasm_simulator'):
    """
    :param circ:  The QuantumCircuit object corresponding to the full circuit, without feeding parameters or transpiling
    :param l: The number of qubits representing a node
    :param w: The weight matrix for edges
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
         return -1 * compute_cost(counts, l, w, nshots)

    return func_to_optimize


def full_optimization_loop(n, l, w, p, bounds=[(-np.pi, np.pi), (0, 4*np.pi)], nshots=512, simulator='qasm_simulator',
                           local_optimization_method='BFGS', optimal_cost=None):
    """
    :param n: The number of nodes
    :param l: The number of qubits representing a node
    :param w: The weight matrix for edges
    :param p: The number of layers
    :param bounds: The bounds for parameters, in the form of [(param_1_low, param_1_high),
                                                              (param_2_low, param_2_high), ...]
    :param nshots: The number of shots
    :param simulator: The simulator
    :param local_optimization_method: The local optimizer to use. By default, BFGS.
                                      Other recommended choices include Powell, Nelder-Mead, COBYLA, ...
                                      When set to None, no local optimization will be performed
    :param optimal_cost: The optimal cost found by brute force. By default, None.
    :return: three lists showing optimization history results, containing respectively
             1. parameters in list format,
             2. cost in float
             3. ciruits run as Qiskit QuantumCircuit objects
    """
    ####################################################################################################################
    # Initialize
    backend = Aer.get_backend(simulator)
    backend.shots = nshots
    param_history = []
    cost_history = []
    circ_history = []
    ####################################################################################################################
    # Run the educated global guess (EGG) optimization for the first time
    circ = make_full_circuit(n, l, w, 1)
    circ_history.append(circ)
    func_to_optimize = func_to_optimize_wrapper(circ, l, w)
    result = differential_evolution(func_to_optimize, bounds)
    param, cost = result.x, -1 * result.fun
    print('1st params', param)
    print('1st cost', cost)
    param_history.append(param)
    cost_history.append(cost)
    ####################################################################################################################
    # If depth = 1, no need to continue
    if p < 2:
        return param_history, cost_history, circ_history
    ####################################################################################################################
    # Else continue
    for i in range(2, p+1):
        if i == 2:
            abbrev = 'nd'
        else:
            abbrev = 'th'
    ####################################################################################################################
    # Run the educated global guess (EGG) optimization for ith iteration
        circ = make_full_circuit(n, l, w, i)
        param_names = circ.parameters
        param_bind_dict = {}
        for j in range(i-1):
            param_prev = param_history[-1]
            param_bind_dict[param_names[j]] = param_prev[j]
            param_bind_dict[param_names[j + i]] = param_prev[j + i - 1]
        circ_w_param = circ.bind_parameters(param_bind_dict)
        circ_history.append(circ_w_param)
        func_to_optimize = func_to_optimize_wrapper(circ_w_param, l, w)
        result = differential_evolution(func_to_optimize, bounds)
        param, cost = result.x, -1 * result.fun
        complete_param = np.concatenate((param_prev[:i-1], np.array([param[0]]),
                                         param_prev[i-1:], np.array([param[1]])))
        print(str(i) + abbrev + ' iteration (EGG), params', complete_param)
        print(str(i) + abbrev + ' iteration (EGG),  cost', cost)
        param_history.append(complete_param)
        cost_history.append(cost)
    ####################################################################################################################
    # Run the local optimization of choice for ith iteration if neededs
        if local_optimization_method is not None:
            func_to_optimize = func_to_optimize_wrapper(circ, l, w)
            result = minimize(func_to_optimize, complete_param, method=local_optimization_method)
            param, cost = result.x, -1 * result.fun
            print(str(i) + abbrev + ' iteration (' + local_optimization_method + '), params', param)
            print(str(i) + abbrev + ' iteration (' + local_optimization_method + '), cost', cost)
            param_history.append(param)
            cost_history.append(cost)
            circ_history.append(circ)
    ####################################################################################################################
    # If optimal cost found by brute force is provided, compute the approximation ratio evolution
    if optimal_cost is not None:
        print('Approximation Ratio Evolution ', cost_history / optimal_cost)
    ####################################################################################################################
    # lists of parameters, costs and quantum circuits are returned
    return param_history, cost_history, circ_history





def qaoa_run(G, l, p, local_optimization_method='Powell', nshots=512):
    n = G.number_of_nodes()
    k = l ** 2

    weights = get_weight_matrix(G)
    R = compute_normalization_scale(weights)
    rescaled_weights = weights / R

    param_history, cost_history, circ_history = full_optimization_loop(n, l, rescaled_weights, p,
                                                                       local_optimization_method=local_optimization_method)

    circ = make_full_circuit(n, l, rescaled_weights, p)
    counts, transpiled_circ = run_circuit(circ, param_history[-1], nshots=nshots)

    return show_distribution(counts, l)



def qaoa_solver(G, k, p):
    l = np.log2(k)
    if l != int(l):
        print('Error: Qaoa only works if k=l**2 with integer-valued l')
        return 0
    else:
        l = int(l)

    distribution_qaoa = qaoa_run(G, l, p, local_optimization_method='Nelder-Mead')

    plot_distribution_diagramm(distribution_qaoa)  # Plot to check

    highest_counts = qaoa_highest_counts(distribution_qaoa, 10)

    C_opt_qaoa = 0
    for key, value in highest_counts.items():
        P = str_array_into_dict(key)
        C = cost(G, P)
        if C > C_opt_qaoa:
            C_opt_qaoa = C
            P_opt_qaoa = P.copy()

    return C_opt_qaoa, P_opt_qaoa