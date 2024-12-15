import numpy as np
import time
import ngesh
import PyQt5
import random
from memory_profiler import memory_usage
import tracemalloc

def dict_to_matrix(D):
    # size = len(D)
    # matrix = np.zeros((size, size))
    # for i in D:
    #     for j in D[i]:
    #         matrix[i][j] = D[i][j]
    # return matrix
    max_key = max(max(D.keys()), max(max(d.keys()) for d in D.values()))
    size = max_key + 1  # Add 1 because indices start at 0
    
    matrix = np.zeros((size, size))
    for i in D:
        for j in D[i]:
            matrix[i][j] = D[i][j]
    return matrix

def matrix_to_dict(M):
    D = {i: {j: int(M[i][j]) for j in range(len(M[i]))} for i in range(len(M))}
    return D

def non_degenerate(mat):
    """
    Checks if the matrix contains any degenerate triple.
    """
    for i in range(len(mat)):
        for j in range(len(mat)):
            for k in range(len(mat)):
                if i != j and i != k and j != k:
                    if mat[i][j] + mat[j][k] == mat[i][k]:
                        return False
    return True

def triples(mat):
    """
    Finds the first degenerate triple in the matrix.
    Returns indices i, j, k.
    """
    n = len(mat)
    for i in range(n):
        for j in range(i + 1, n): 
            for k in range(j + 1, n):
                if abs(mat[i][j] + mat[j][k] - mat[i][k]) < 1e-10:  
                    return i, j, k
    return None

def compute_delta(mat):
    """
    Computes delta for a matrix. Delta is the minimal correction factor.
    """

    min_delta = float("inf")
    n = len(mat)
    
    for x in range(n):
        for a in range(x + 1, n):
            for b in range(a + 1, n):
                delta = 0.5 * (mat[a][b] - abs(mat[a][x] - mat[b][x]))
                if delta > 0:  # Only consider positive deltas
                    min_delta = min(delta, min_delta)
                    
    return min_delta

def update_tree(T, i, j, k, x):
    """
    Updates the tree T with the resolved triple (i, j, k) and edge weight x.
    """
    for vertex in (i, j, k):
        if vertex not in T:
            T[vertex] = {}
            
    # Add edges between i and k
    T[i][k] = x
    T[k][i] = x

def degenerative_triples(D):
    T = {}  # Tree to hold the result
    M = dict_to_matrix(D)  # Convert dictionary to matrix

    # Base case: If only two nodes remain
    if len(D.keys()) == 2:
        for node in D.keys():
            for other_node in D[node]:
                if node not in T:
                    T[node] = {}
                T[node][other_node] = D[node][other_node]
                if other_node not in T:
                    T[other_node] = {}
                T[other_node][node] = D[other_node][node]
        return T

    # If the matrix is non-degenerate
    if non_degenerate(M):
        delta = compute_delta(M)
        for i in D.keys():
            for j in D[i]:
                D[i][j] -= 2 * delta
                D[j][i] -= 2 * delta
    # Resolve a triple
    if triples(M) == None:
        return "No triple found"
    i, j, k = triples(M)  # Find a degenerate triple
    x = M[i][j]  # Distance between i and j

    # Update tree with resolved triple
    T = update_tree(T, i, j, k, x)

    # Remove j from the matrix and dictionary
    for node in list(D.keys()):
        if j in D[node]:
            del D[node][j]
    del D[j]

    degenerative_triples(T)
    # return T

D0 = [
    [0, 2, 8, 7],
    [2, 0, 6, 5],
    [8, 6, 0, 7],
    [7, 5, 7, 0]
]

def testing():
    D = matrix_to_dict(D0)
    T = degenerative_triples(D)
    print(T)

testing()