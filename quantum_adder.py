from numpy import pi
from qiskit import QuantumCircuit, assemble, Aer
from qiskit.visualization import plot_histogram
from IPython.display import display

def binary_to_integer(num):
    encoding = 0
    for i in range(len(num)):
        encoding += (2**(len(num)-1-i))*int(num[i])
    return int(encoding)

def integer_to_binary(num):
    binary = ""
    while num > 1:
        remainder = num % 2
        num = num//2
        binary = str(remainder) + binary
    remainder = num % 2
    binary = str(remainder) + binary
    return binary

def transform(a, b):
    
    a = integer_to_binary(a)
    b = integer_to_binary(b)
    
    if len(a) > len(b):
        b = "0"*(len(a) - len(b)) + b
    elif len(b) > len(a):
        a = "0"*(len(b) - len(a)) + a
        
    a = "0" + a
    b = "0" + b
    
    return a, b
    
def initialize(qc, a, b):
    
    for i in range(len(a)):
        if a[i] == "1":
            qc.x(len(a) - i - 1)
            
    for i in range(len(b)):
        if b[i] == "1":
            qc.x(len(a) + len(b) - i - 1)
    
    return qc

def qft(qc, n, last_bit):
    if last_bit == n:
        return qc
    last_bit -= 1
    qc.h(last_bit)
    for qubit in range(n, last_bit):
        qc.cp(pi/2**(last_bit-qubit), qubit, last_bit)
    qft(qc, n, last_bit)
    
def inverse_qft(qc, n, last_bit):
    
    for tbit in range(n, last_bit):
        k = tbit - n
        for cbit in range(n, tbit):
            qc.cp(-pi/(2**(k)), cbit, tbit)
            k -= 1
        qc.h(tbit)
        
    return qc

def quantum_adder(qc):
            
    n = int(len(qc.qubits)/2)
    qc.barrier()
    qft(qc, n, 2*n)
    qc.barrier()
    
    for i in range(n):
        tbit = n + i
        k = i
        for cbit in range(i+1):
            qc.cp(pi/(2**k), cbit, tbit)
            k -= 1

    qc.barrier()
    inverse_qft(qc, n, 2*n)
    qc.barrier()
    qc.measure(range(n, 2*n), range(n))
    
    return qc

def qft_adder(a, b, display_circuit = False, display_counts = False):
            
    bin_a, bin_b = transform(a, b)
    n = len(bin_a)
    qc = QuantumCircuit(2*n, n)
    qc = initialize(qc, bin_a, bin_b)
    
    n = int(len(qc.qubits)/2)
    qc.barrier()
    qft(qc, n, 2*n)
    qc.barrier()
    
    for i in range(n):
        tbit = n + i
        k = i
        for cbit in range(i+1):
            qc.cp(pi/(2**k), cbit, tbit)
            k -= 1

    qc.barrier()
    inverse_qft(qc, n, 2*n)
    qc.barrier()
    qc.measure(range(n, 2*n), range(n))
    
    if display_circuit:
        display(qc.draw("mpl"))

    qobj = assemble(qc)
    sim = Aer.get_backend("aer_simulator")
    result = sim.run(qobj).result()
    counts = result.get_counts()
    if display_counts:
        display(plot_histogram(counts))
    max_count = max(zip(counts.values(), counts.keys()))[1]
    soln = binary_to_integer(max_count)

    return qc, soln
    