from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister
from qiskit import Aer, execute
from qiskit.circuit import Parameter


def make_mixing_block(nq):
    beta = Parameter('\u03B2')
    mix_circ = QuantumCircuit(nq)
    for q in range(0, nq):
        mix_circ.rx(2 * beta, q)
    return mix_circ