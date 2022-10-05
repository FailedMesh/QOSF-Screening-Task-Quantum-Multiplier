from numpy import pi
from qiskit import QuantumCircuit, assemble, Aer
from qiskit.visualization import plot_histogram
from IPython.display import display

#*****************************************************************************************************#
################################### BINARY TO INTEGER CONVERSION ######################################
#*****************************************************************************************************#

def binary_to_integer(binary):
    """
    Converts a binary string to an integer (int)
    Input : binary = String, binary number
    Output : integer = Int, binary number converted to decimal representation
    """
    integer = 0
    for i in range(len(binary)):
        integer += (2**(len(binary)-1-i))*int(binary[i])
    return integer

#*****************************************************************************************************#
################################### INTEGER TO BINARY CONVERSION ######################################
#*****************************************************************************************************#

def integer_to_binary(num):
    """
    Converts an integer (int) to binary string
    Input : num = Int, in decimal representation
    Output : binary = String, converted binary number
    """
    binary = ""
    while num > 1:
        remainder = num % 2
        num = num//2
        binary = str(remainder) + binary
    remainder = num % 2
    binary = str(remainder) + binary
    return binary

#*****************************************************************************************************#
############################### TRANSFORM NUMBERS TO BE MULTIPLIED ####################################
#*****************************************************************************************************#

def transform(a, b):
    """
    Converts the smaller number into a bit-string, then extends the bit-string with zeroes to match the size of
    the result of the multiplication.
    Input : a = Int, Integer
             b = Int, Integer
    Output : binary = String, extended bit-string of the same size as a*b
              integer = Int, the smaller of the 2 integers
    """
    
    bin_a = integer_to_binary(a)  #Convert a to bit-string
    bin_b = integer_to_binary(b)  #Convert b to bit-string
    
    #The number of bits required to represent the product of a and b will be the sum of the number of bits in
    #each of a and b:
    total_bits = len(bin_a) + len(bin_b)
    
    #Find the smaller of the numbers and extend the bitstring to match the size of the product bitstring.
    #We want to represent the smaller number in the circuit, so that lesser number of qubits are used, thus
    #reducing space complexity.
    if b < a:
        bin_b = "0"*(total_bits - len(bin_b)) + bin_b
        binary, integer =  bin_b, a
    else:
        bin_a = "0"*(total_bits - len(bin_a)) + bin_a
        binary, integer =  bin_a, b

    return binary, integer
    
#*****************************************************************************************************#
################################### INITIALIZE QUANTUM CIRCUIT ########################################
#*****************************************************************************************************#

def initialize(qc, a):
    """
    Initializes a given quantum circuit so that both the first n-qubits and the second n-qubits match the
    bit-string 'a'
    Input : qc = QuantumCircuit
            a = String, bit-string
    Output : qc = QuantumCircuit, modified quantum circuit
    """
    
    #Loop through qubits and apply an X-gate to the qubits that correspond to a "1" in the bit-string
    #The first n-qubits and second n-qubits are duplicated, where n = len(a)
    for i in range(len(a)):
        if a[i] == "1":
            qc.x(len(a) - i - 1)
            qc.x(2*len(a) - i - 1)
    
    return qc

#*****************************************************************************************************#
################################### QUANTUM FOURIER TRANSFORM #########################################
#*****************************************************************************************************#

def qft(qc, n, last_bit):
    """
    Applies a quantum fourier transform from the nth bit to the last_bit in a quantum circuit qc
    Input : qc = QuantumCircuit
            n = Int, starting bit for fourier transform
            last_bit = Int, last bit for fourier transform
    Output : qc = QuantumCircuit, modified quantum circuit
    """

    if last_bit == n:
        return qc
    last_bit -= 1
    qc.h(last_bit) #Change basis to |+> and |->
    for qubit in range(n, last_bit):
        qc.cp(pi/2**(last_bit-qubit), qubit, last_bit) #Controlled rotation around Z to match the phase required
    qft(qc, n, last_bit)

#*****************************************************************************************************#
############################### INVERSE QUANTUM FOURIER TRANSFORM #####################################
#*****************************************************************************************************#
    
def inverse_qft(qc, n, last_bit):
    """
    Applies an inverse quantum fourier transform from the nth bit to the last_bit in a quantum circuit qc
    Input : qc = QuantumCircuit
            n = Int, starting bit for inverse fourier transform
            last_bit = Int, last bit for inverse fourier transform
    Output : qc = QuantumCircuit, modified quantum circuit
    """

    #Loop from starting to last bit as control bit
    for tbit in range(n, last_bit):
        k = tbit - n
        for cbit in range(n, tbit):  #Loop through target bits
            qc.cp(-pi/(2**(k)), cbit, tbit) #Controlled reverse rotation around Z back to |+> or |-> 
            k -= 1
        qc.h(tbit) #Apply Hadamard to change basis to 0 and 1
        
    return qc

#*****************************************************************************************************#
##################################### QUANTUM MULTIPLIER ##############################################
#*****************************************************************************************************#

def quantum_multiplier(a, b, display_circuit = False, display_counts = False):
    """
    Multiplies two integers a and b using a quantum circuit, quantum fourier transform and quantum gates
    Input : a = Int, integer number
            b = Int, integer number
            display_circuit = Bool, whether to display the circuit or not
            display_counts = Bool, whether to display the measurement distribution or not
    Output : qc = QuantumCircuit, modified quantum circuit
             soln = Int, the product of a and b
    """

    multiplicant, multiplier = transform(a, b) #Transform a and b
    n = len(multiplicant)
    qc = QuantumCircuit(2*n, n)         #Create the quantum circuit
    qc = initialize(qc, multiplicant)   #Initialize the input qubits to the multiplicant
    
    n = int(len(qc.qubits)/2)
    qc.barrier()
    qft(qc, n, 2*n)         #Apply the fourier transform to the last 'n' qubits of the circuit
    qc.barrier()
    
    #REPEATED ADDITION:
    #Do a controlled rotation of each the last n-qubits which are in fourier basis, controlled by the first 
    #n-qubits which represent the same number but not converted to fourier basis.
    #Apply each controlled rotation "(multiplier - 1)" times, this performs repeated addition of the qubit-set
    #to itself "(multiplier - 1)" times, which gives us the product (multiplicant * multiplier) .
    for i in range(n):
        tbit = n + i
        k = i
        for cbit in range(i+1):
            qc.cp((multiplier-1)*pi/(2**k), cbit, tbit)
            k -= 1

    inverse_qft(qc, n, 2*n)  #Apply inverse fourier transform to convert the result back to 0-1 basis
    qc.barrier()
    qc.measure(range(n, 2*n), range(n))  #Measure the result and store in the classical register

    if display_circuit:
        display(qc.draw("mpl"))
    
    sim = Aer.get_backend("aer_simulator")
    qobj = assemble(qc)
    result = sim.run(qobj).result()
    counts = result.get_counts()   #Measure the counts stored in the classical register
    if display_counts:
        display(plot_histogram(counts))
    max_count = max(zip(counts.values(), counts.keys()))[1]  #Get the bitstring with highest counts
    soln = binary_to_integer(max_count) #Convert the highest probability bitstring to an integer to get the answer
    
    return qc, soln