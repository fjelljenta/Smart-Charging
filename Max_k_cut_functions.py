from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister
from qiskit import Aer, execute
from qiskit.circuit import Parameter

from data_processing import *
import numpy as np

def make_mixing_block(nq):
    beta = Parameter('\u03B2')
    mix_circ = QuantumCircuit(nq)
    for q in range(0, nq):
        mix_circ.rx(2 * beta, q)
    return mix_circ


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
    
    #TODO: implement for general number od Nodes N, at the moment only version for N=5!!!!
    #N=5
    
    #TODO: implememnt for general k, at the moment only version for k=4!!!!
    #k=4
    
    #TODO: at the moment only works if node names correspond to the number (int) of the node
    #jobs=list(G.nodes())
    
    P_opt={'P1':[1],'P2':[],'P3':[],'P4':[]}
    C_opt = 0
    
    #PN0=0  # node o is always in Partition 0
    for PN1 in range(4):
        for PN2 in range(4):
            for PN3 in range(4):
                for PN4 in range(4):
                    P_try={'P0':[0], 'P1':[],'P2':[],'P3':[]}
                    P_try['P'+str(PN1)].append(1)
                    P_try['P'+str(PN2)].append(2)
                    P_try['P'+str(PN3)].append(3)
                    P_try['P'+str(PN4)].append(4)
                    C_try=cost(G,P_try)
                    if C_try>C_opt:
                        P=P_try.copy()
                        C_opt=C_try
    
    return C_opt, P 


def cost(G,P):
    """Cost of a given partition P of Graph G
    Args:
        G (graph): graph the partition is defined on
        P (dict): dictionary with changing points as keys and array of jobs (nodes) for the charging point as values
    Returns:
        C (int): cost function of the given Partition
    """
    partitions=list(P.keys())
    C=0
    for Pr in partitions:
        for Ps in partitions[partitions.index(Pr)+1:]:
            C= C+cut_weight(G, P[Pr], P[Ps])
            #print('Pr=',Pr , ' Ps=', Ps, '   cut_weight=', cut_weight(G, P[Pr], P[Ps]), '\n')
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
    w= np.sum([weights[i,j] for i in P1 for j in P2])
    return w
