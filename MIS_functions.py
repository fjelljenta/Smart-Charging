from qiskit import QuantumCircuit
from qiskit import Aer
from qiskit.compiler import transpile
from qiskit.circuit import Parameter
from itertools import combinations
import datetime as dt
from Max_k_cut_functions import *

from data_processing import *
import numpy as np
import copy


def make_cost_block(n, layer):
    """
    :param n: The number of nodes
    :param i: The subscript to be given to the variable gamma, the varibale would be named to 'gamma_i'
    :return: The QuantumCircuit object corresponding to the mixing block
    """
    globals()['gamma%s' % layer] = Parameter('\u03B3' + str(layer))
    mix_circ = QuantumCircuit(n)
    for q in range(0, n):
        mix_circ.rz((-1) * globals()['gamma%s' % layer], q)
    return mix_circ


def make_mixing_block(n, w, layer):
    """
    :param n: The number of nodes
    :param w: The adjacency matrix of the graph
    :param i: The subscript to be given to the variable beta, the varibale would be named to 'beta_i'
    :return: The QuantumCircuit object corresponding to the cost block
    """
    globals()['beta%s' % layer] = Parameter('\u03B2' + str(layer))
    cost_circ = QuantumCircuit(n)
    for node in range(n):  # node corresponding to term in Hc expression
        adj_nodes=[i for i in range(n) if w[node][i]==1]
        cost_circ.h(node)  # apply the starting h gate
        for s in powerset(adj_nodes):  # generate a power set so that we can perform cx-rz-cx, cx-cx-rz-
            # cx-cx, etc. sub-blocks one by one
            seq = list(s)  # get list of all the qubits involved for a sub-block (s is a tuple)
            seq.sort()  # we want to apply cx gates in order
            seq+=[node]
            for i in range(len(seq) - 1):
                cost_circ.cx(seq[i], seq[i + 1])  # apply cx gates in ascending order of index
            cost_circ.rz(globals()['beta%s' % layer] / 2**(len(s)-1), seq[-1])  # apply the rz gate
            for i in range(len(seq) - 1):
                cost_circ.cx(seq[-(i + 2)], seq[-(i + 1)])  # apply cx gates in descending order of index
            cost_circ.barrier()
        cost_circ.h(node)  # apply the ending h gate
        #cost_circ.barrier()
    return cost_circ


def make_full_circuit(n, w, p):
    """
    :param n: The number of nodes
    :param w: The adjacency matrix of the graph
    :param p: The number of layers
    :return: The QuantumCircuit object corresponding to the full circuit
    """
    circ = QuantumCircuit(n)
    for layer in range(p):
        circ.append(make_mixing_block(n, w, layer), [i for i in range(n)])
        circ.append(make_cost_block(n, layer), [i for i in range(n)])
    circ.measure_all()
    return circ
