import numpy as np
import math
import itertools

# P - parity check matrix in the form
# P = [A] / [I]
# Alphabet = {0,1}

P = np.array([[0, 0, 1, 1],
              [0, 1, 0, 1],
              [0, 1, 1, 0],
              [0, 1, 1, 1],
              [1, 0, 0, 1],
              [1, 0, 1, 0],
              [1, 0, 1, 1],
              [1, 1, 0, 0],
              [1, 0, 0, 0],
              [0, 1, 0, 0],
              [0, 0, 1, 0],
              [0, 0, 0, 1]])

# Message example

message = np.array([0, 0, 1, 1, 1, 0, 1, 0])
# encoded_m = np.array([0,0,1,0,0,1])

'''
def get_idenity_matrix(P):
    # I is the identity matrix
    I = P[-P.shape[1]:]
    print('The idenity matrix is ' + str(I))
    return I
'''
# Gets the A matrix from the parity check matrix


def get_A_matrix(P):
    n, k = P.shape
    print('P has ' + str(n) + ' rows and ' + str(k) + ' columns.', end='\n\n')
    # A matrix
    A = P[:(n-k), :]
    print('The A matrix is \n' + str(A), end='\n\n')
    return A

# Constructs the generator matrix G given the A matrix


def construct_generator_matrix(A):
    I = np.identity(A.shape[0], dtype=int)
    G = np.concatenate((I, A), axis=1)
    print('The generator matrix is \n' + str(G), end='\n\n')
    return G

# Constructs the syndrom table


def construct_syndrome_table(P):
    n = P.shape[0]
    syndrome_table = {}

    error = np.zeros(n, dtype=int)
    error[0] = 1

    for s in P:
        syndrome_table[str(s)] = str(error)
        error = np.roll(error, 1)

    print('The syndrome table is \n')
    for s, err in syndrome_table.items():
        print(" Syndrome: {} -> Error:{}".format(s, err))
    print()

    return syndrome_table

# Computes the minimum hamming distance given all the codewords
# Ref: https://stackoverflow.com/questions/42752610/python-how-to-generate-the-pairwise-hamming-distance-matrix


def compute_min_HammingDistance(X):
    hamming_d = (X[:, None, :] != X).sum(2)
    min_hamming_d = np.min(hamming_d[np.nonzero(hamming_d)])
    print('The minimum hamming distance is ' + str(min_hamming_d), end='\n\n')
    return min_hamming_d

# Generates the all the codewords for all messages


def generate_messages_codewords(G, nk):
    messages = np.array(list(itertools.product([0, 1], repeat=nk)))
    codewords = np.empty([messages.shape[0], G.shape[1]], dtype=int)
    for i in range(0, messages.shape[0]-1):
        m = encode(messages[i], G)
        codewords[i] = m

    # print('enodings' + str(encodings))

    messages_codewords = {}

    for m, e in zip(messages, codewords):
        messages_codewords[str(m)] = str(e)

    return (codewords, messages_codewords)

# Returns the maximum number of error that can be detected and covered


def can_detect(d):
    detected = d-1
    covered = math.floor((d-1)/2)
    return (detected, covered)

# Decodes a message given the parity check matrix


def decode(encoded_m, P):
    syndrome_table = construct_syndrome_table(P)

    n, k = P.shape
    m_P = np.dot(encoded_m, P)

    for i in range(0, m_P.size):
        if m_P[i] < 2:
            continue
        elif m_P[i] % 2 == 0:
            m_P[i] = 0
        else:
            m_P[i] = 1

    print('The result of multiplying the encoded message with the parity check matrix is ' + str(m_P), end='\n\n')

    encoded_without_error = np.array([])

    if np.sum(m_P) == 0:
        print('The code has no errors.\n')
    elif (np.sum(m_P) != 0 and (m_P in P)):
        print('The result is not 0 and it is present in the parity check matrix.\n')
        error = np.fromstring(
            syndrome_table[str(m_P)][1:-1], dtype=int, sep=' ')
        print(error)
        encoded_without_error = np.add(encoded_m, error)
        print('The encoded message without errors is ' +
              str(encoded_without_error), end='\n\n')
    else:
        print('Something went wrong. :)')
    return encoded_m[:n-k]

# Encodes a message given the Generator matrix


def encode(message, G):
    encoded_m = np.dot(message, G)

    for i in range(0, encoded_m.size):
        if encoded_m[i] < 2:
            continue
        elif encoded_m[i] % 2 == 0:
            encoded_m[i] = 0
        else:
            encoded_m[i] = 1

    return encoded_m


# Encoding
A = get_A_matrix(P)
G = construct_generator_matrix(A)


codewords, message_codewords = generate_messages_codewords(G, A.shape[0])
hamming_d = compute_min_HammingDistance(codewords)

detected, covered = can_detect(hamming_d)
print('The maximum number of errors that can be detected is ' +
      str(detected), end='\n\n')
print('The maximum number of errors that can be covered is ' +
      str(covered), end='\n\n')


encoded_m = encode(message, G)
print('The encoded message for ' + str(message) +
      ' is ' + str(encoded_m), end='\n\n')


# Decoding
decoded_m = decode(encoded_m, P)
print('The decoded message is ' + str(decoded_m), end='\n\n')
