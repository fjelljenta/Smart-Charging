from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister
from qiskit import Aer, execute
from qiskit.circuit import Parameter
from itertools import combinations

from data_processing import *
import numpy as np


def make_mixing_block(n, l, layer):
    """
    :param n: The number of nodes
    :param l: The number of qubits representing a node
    :param i: The subscript to be given to the variable beta, the varibale would be named to 'beta_i'
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
    :param i: The subscript to be given to the variable gamma, the varibale would be named to 'gamma_i'
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
                cost_circ.barrier()
    return cost_circ


def powerset(l):
    """
    :param l: a list containing numbers usually
    :return: the powerset, or to say a set containing all the subsets, of the set of the list, in list format
    """
    res = []
    n = len(l)
    for i in range(0, n + 1):
        for element in combinations(l, i):
            res.append(element)
    return res


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
    nq = n * l
    circ = QuantumCircuit(nq)
    circ.append(make_initial_block(n, l), [i for i in range(nq)])
    for layer in range(p):
        circ.append(make_mixing_block(n, l, layer), [i for i in range(nq)])
        circ.append(make_cost_block(n, l, w, layer), [i for i in range(nq)])
    return circ


def brut_force_k4_N5(G):
    """Classical brut-force solution of the max-k-cut problem of Graph G
    Args:
        G (graph): Graph of the max-k-cut problem
        #TODO[k (int): Number of cuts (available loading stations)]
    Returns:
        C_opt (int): cost function of optimal max-k-cut solution
        P_opt (dict): dictionary with the number of the load station as keys 
            and array of jobs at a given loading station (still unordered!) as values
    Method: 
        compute every posible k-cut, save the corresponding cost function 
        and find the maximal index of the maximal cost
    """

    # TODO: implement for general number od Nodes N, at the moment only version for N=5!!!!
    # N=5

    # TODO: implememnt for general k, at the moment only version for k=4!!!!
    # k=4

    # TODO: at the moment only works if node names correspond to the number (int) of the node
    # jobs=list(G.nodes())

    P_opt = {'P1': [1], 'P2': [], 'P3': [], 'P4': []}
    C_opt = 0

    # PN0=0  # node o is always in Partition 0
    for PN1 in range(4):
        for PN2 in range(4):
            for PN3 in range(4):
                for PN4 in range(4):
                    P_try = {'P0': [0], 'P1': [], 'P2': [], 'P3': []}
                    P_try['P' + str(PN1)].append(1)
                    P_try['P' + str(PN2)].append(2)
                    P_try['P' + str(PN3)].append(3)
                    P_try['P' + str(PN4)].append(4)
                    C_try = cost(G, P_try)
                    if C_try > C_opt:
                        P = P_try.copy()
                        C_opt = C_try

    return C_opt, P


def cost(G, P):
    """Cost of a given partition P of Graph G
    Args:
        G (graph): graph the partition is defined on
        P (dict): dictionary with changing points as keys and array of jobs (nodes) for the charging point as values
    Returns:
        C (int): cost function of the given Partition
    """
    partitions = list(P.keys())
    C = 0
    for Pr in partitions:
        for Ps in partitions[partitions.index(Pr) + 1:]:
            C = C + cut_weight(G, P[Pr], P[Ps])
            # print('Pr=',Pr , ' Ps=', Ps, '   cut_weight=', cut_weight(G, P[Pr], P[Ps]), '\n')
    return C


def cut_weight(G, P1, P2):
    """weight of all edges connecting the subnodesets P1 and P2 in Graph G
    Args:
        G (graph): graph the partition is defined on
        P1 (dict): nodes of the first subset
        P2 (dict): nodes of the second subset
    Returns:
        w (int): sum of weights with one node in P1 and the other in P2
    """
    weights = get_weight_matrix(G)
    w = np.sum([weights[i, j] for i in P1 for j in P2])
    return w
