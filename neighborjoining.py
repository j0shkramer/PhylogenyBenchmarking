import numpy as np
import time
import ngesh
from memory_profiler import memory_usage
import tracemalloc

class NeighborJoining:
    def __init__(self):
        pass

    def is_additive(self, D):
        """
        Returns true if the square matrix D is additive.

        :param D: nxn list of lists of ints representing a distance matrix
        :return: True if D is an additive distance matrix
        """
        n = len(D)
        for i in range(n - 3):  # Only need to check up to n-3 elements
            for j in range(i + 1, n - 2):
                for k in range(j + 1, n - 1):
                    for l in range(k + 1, n):
                        if D[i][j] + D[k][l] > D[i][k] + D[j][l]:
                            return False
        return True

    def check_tree(self, D, T):
        """
        Implement the Floyd-Warshall algorithm on the graph defined in T to verify its consistency with D.

        :param D: Distance dictionary representing pairwise distances between leaves
        :param T: Dictionary of dictionaries representing the tree with edges and weights
        :return: True if the tree T is consistent with the distance matrix D, otherwise False
        """
        V = len(T)
        dT = [[float("inf")] * V for _ in range(V)]
        n = (V + 2) // 2  # Compute n as half of total number of nodes

        # Length check
        if n != len(D):
            return False

        # Initialize the distance matrix from existing edges
        for k in T:
            for l in T[k]:
                dT[k][l] = T[k][l]
                dT[l][k] = T[l][k]

        # Fill in the diagonal elements (distance to self)
        for i in range(len(dT)):
            dT[i][i] = 0.0

        # Apply Floyd-Warshall algorithm to find shortest paths
        for k in range(V):
            for i in range(V):
                for j in range(V):
                    if dT[i][j] > dT[i][k] + dT[k][j]:
                        dT[i][j] = dT[i][k] + dT[k][j]
                        dT[j][i] = dT[i][j]

        # Check if the reconstructed distances match the given distance matrix D
        for i in range(n):
            for j in range(n):
                if D[i][j] != dT[i][j]:
                    return False
        return True

    def min_S_value(self, D, u):
        m = len(D)
        min_S, min_i, min_j = float("inf"), -1, -1
        for k in D:
            for l in D[k]:
                if l != k:
                    crit = (m - 2) * D[k][l] - u[k] - u[l]
                    if crit < min_S:
                        min_S = crit
                        min_i = k
                        min_j = l
        return min_i, min_j

    def neighbor_join(self, D):
        T = {}
        nodeCount = len(D)

        # Initialize the tree with the original leaves
        for node in D:
            T[node] = {}

        while len(D) > 2:
            # Compute u_k
            u = {k: sum(D[k][j] for j in D[k] if j != k) for k in D}

            # Find i, j that minimizes S[i][j]
            i, j = self.min_S_value(D, u)

            # Create a new internal node
            r = nodeCount
            nodeCount += 1

            # Get edge weights for the internal node
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

            # Remove the rows and columns for i and j from D
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

    def dict_to_matrix(self, D):
        size = len(D)
        matrix = np.zeros((size, size))
        for i in D:
            for j in D[i]:
                matrix[i][j] = D[i][j]
        return matrix

    def matrix_to_dict(self, M):
        D = {i: {j: int(M[i][j]) for j in range(len(M[i]))} for i in range(len(M))}
        return D

def analyze_trees_nj():
    # Generate a tree
    tree = ngesh.gen_tree(1.0, 0.5, max_time=10, labels="enum", num_leaves=5, seed=20)
    print(tree)

    # Create distance matrix dictionary
    distance_matrix = {}
    node_tracker = {}
    # Go through each leaf and find its distance from every other leaf
    for i, node in enumerate(tree.iter_leaves()):
        node_tracker[i] = node.name
        distance_matrix[i] = {}
        for j, other_node in enumerate(tree.iter_leaves()):
            distance_matrix[i][j] = node.get_distance(other_node)

    # print(distance_matrix)
    matrix = NeighborJoining().dict_to_matrix(distance_matrix)
    print("This distance matrix is additive: ", NeighborJoining().is_additive(matrix))
    # compute running time and space complexity
    start = time.time()
    tracemalloc.start()
    remade_tree = NeighborJoining().neighbor_join(distance_matrix)
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    end = time.time()

    print(f"Memory: {current / 1024} KBs")
    print(f"Time elapsed: {end - start} seconds")
    print(remade_tree)

analyze_trees_nj()


            

