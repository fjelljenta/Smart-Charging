from Code.Max_k_cut_quantum_functions import *


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
        measurement=measurement[::-1]
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
         return compute_cost_MIS(counts, w, U, nshots)

    return func_to_optimize


def full_optimization_loop_MIS(n, w, U, p, bounds=[(-np.pi, np.pi), (0, 4*np.pi)], nshots=512,
                               simulator='qasm_simulator', local_optimization_method='BFGS', optimal_cost=None):
    """
    :param n: The number of nodes
    :param w: The adjacency matrix for edges
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
    circ = make_full_circuit_MIS(n, w, 1)
    circ_history.append(circ)
    func_to_optimize = func_to_optimize_wrapper_MIS(circ, w, U, nshots=nshots, simulator=simulator)
    result = differential_evolution(func_to_optimize, bounds)
    param, cost = result.x, result.fun
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
        circ = make_full_circuit_MIS(n, w, i)
        param_names = circ.parameters
        param_bind_dict = {}
        for j in range(i-1):
            param_prev = param_history[-1]
            param_bind_dict[param_names[j]] = param_prev[j]
            param_bind_dict[param_names[j + i]] = param_prev[j + i - 1]
        circ_w_param = circ.bind_parameters(param_bind_dict)
        circ_history.append(circ_w_param)
        func_to_optimize = func_to_optimize_wrapper_MIS(circ_w_param, w, U, nshots=nshots, simulator=simulator)
        result = differential_evolution(func_to_optimize, bounds)
        param, cost = result.x, result.fun
        complete_param = np.concatenate((param_prev[:i-1], np.array([param[0]]),
                                         param_prev[i-1:], np.array([param[1]])))
        print(str(i) + abbrev + ' iteration (EGG), params', complete_param)
        print(str(i) + abbrev + ' iteration (EGG),  cost', cost)
        param_history.append(complete_param)
        cost_history.append(cost)
    ####################################################################################################################
    # Run the local optimization of choice for ith iteration if needed
        if local_optimization_method is not None:
            func_to_optimize = func_to_optimize_wrapper_MIS(circ_w_param, w, U, nshots=nshots, simulator=simulator)
            result = minimize(func_to_optimize, complete_param, method=local_optimization_method)
            param, cost = result.x, result.fun
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
