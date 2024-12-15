import numpy as np
import time
import ngesh
import PyQt5
import random
from memory_profiler import memory_usage
import tracemalloc


def is_additive(D):
    """
    Returns true if the square matrix D is additive

    :param: D is an nxn list of lists of ints
    :return: true if D is an additive distance matrix
    """
    n = len(D)
    for i in range(n - 3):
        for j in range(i + 1, n - 2):
            for k in range(j + 1, n - 1):
                for l in range(k + 1, n):
                    if D[i][j] + D[k][l] > D[i][k] + D[j][l]:
                        return False
    return True

def generate_matrix(length):
    # additive = False
    M = np.random.randint(0, 100, (length, length)); M = (M + M.T) // 2
    for i in range(length):
        M[i][i] = 0
    D = matrix_to_dict(M)
    return M

def check_tree(D,T):
    """
    Implement floyd-warshall algorithm on the graph defined in T

    """

    V = len(T)
    dT = [[float("inf") for i in range(V)] for j in range(V)]
    n = (V + 2)//2

    # Length check
    if n!=len(D):
        return False

    # fill in existing edges
    for k in T:
        for l in T[k]:
            dT[k][l] = T[k][l]
            dT[l][k] = T[l][k]

    # fill in the diagonal elements
    for i in range(len(dT)):
        dT[i][i] = 0.0

    # relax edges
    for k in range(V):
        for i in range(V):
            for j in range(V):
                if dT[i][j] > dT[i][k] + dT[k][j]:
                    dT[i][j] = dT[i][k] + dT[k][j]
                    dT[j][i] = dT[i][j]

    for i in range(n) :
        print(dT[i][:n])
    # Check each value in dT
    for i in range(n):
        for j in range(n):
            if D[i][j] != dT[i][j]:
                return False
    return True

def min_S_value(D, u):
    """
    returns the value (i,j) for which
    (m-2)*D[i][j] - u_i - u_j is minimum
    """
    m = len(D)
    min_S, min_i, min_j = float("inf"),-1,-1
    for k in D:
        for l in D[k]:
            if l!=k:
                crit = (m-2)*D[k][l] - u[k] - u[l]
                if crit < min_S:
                    min_S = crit
                    min_i = k
                    min_j = l
    return (min_i, min_j)


def neighbor_join(D):
    """
    Takes a distance matrix D, and returns the tree T
    consistent with the closest additive matrix D' to D.

    :param: D is a dict of dicts representing pairwise distances between leaves
    :return: a dict of dicts that contains all the edges with their weights in the tree defined by D'.
    """
    T = {}
    nodeCount = len(D)  # Counter for nodes

    # Initialize the tree with the original leaves
    for node in D:
        T[node] = {}

    while len(D) > 2:
        # Compute u_k

        u = {k: sum(D[k][j] for j in D[k] if j != k) for k in D}

        # Find i, j that minimizes S[i][j]
        i, j = min_S_value(D, u)

        # Create a new internal node
        r = nodeCount
        nodeCount += 1

        # Get edge weights for internal node
        T[r] = {}
        T[i][r] = 0.5 * (D[i][j] + (u[i] - u[j]) / (len(D) - 2))
        T[r][i] = T[i][r]
        D[i][r] = T[i][r]
        T[j][r] = 0.5 * (D[i][j] + (u[j] - u[i]) / (len(D) - 2))
        T[r][j] = T[j][r]
        D[j][r] = T[j][r]

        # Precompute and add new row and column in D for r with distances
        D[r] = {}
        for m in list(D.keys()):
            if m != i and m != j:
                D[r][m] = 0.5 * (D[i][m] + D[j][m] - D[i][j])
                D[m][r] = D[r][m]

        # remove the rows and columns for i and j from D
        del D[i]
        del D[j]
        for k in list(D.keys()):
          if i in D[k]:
            del D[k][i]
          if j in D[k]:
            del D[k][j]
    # Handle the last two remaining nodes
    if len(D) == 2:
        remaining_nodes = list(D.keys())
        dist = D[remaining_nodes[0]][remaining_nodes[1]]
        T[remaining_nodes[0]][remaining_nodes[1]] = dist
        T[remaining_nodes[1]][remaining_nodes[0]] = dist
    return T


def dict_to_matrix(D):
    size = len(D)
    matrix = np.zeros((size, size))
    for i in D:
        for j in D[i]:
            matrix[i][j] = D[i][j]
    return matrix

def matrix_to_dict(M):
    D = {i: {j: int(M[i][j]) for j in range(len(M[i]))} for i in range(len(M))}
    return D


def analyze_trees_nj():
    # generate tree
    tree = ngesh.gen_tree(1.0, 0.5, max_time=10, labels="enum", num_leaves=5, seed=20)
    print(tree)
    # create distance matrix dictonary 
    distance_matrix = {}
    # need to convert named "string" nodes to integer, this is to keep track of which node is which
    node_tracker = {}
    i = 0
    # go through each leaf and find its distance from every other leaf
    for node in tree.iter_leaves():
        distance_matrix[i] = {}
        node_tracker[i] = node.name
        j = 0
        for other_nodes in tree.iter_leaves():
            distance_matrix[i][j] = node.get_distance(other_nodes)
            j = j + 1
        i = i + 1
    print(distance_matrix)
    matrix = dict_to_matrix(distance_matrix)
    print("This distance matrix is additve: ", is_additive(matrix))
    start = time.time()
    tracemalloc.start()
    remade_tree = neighbor_join(distance_matrix)
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    end = time.time()
    print(f"Memory: {current / 1024} KBs")
    print(f"Time elapased: {end - start} seconds")
    print(remade_tree)

analyze_trees_nj()


            

