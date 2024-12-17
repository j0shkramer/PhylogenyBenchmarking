import numpy as np
import time
import ngesh
from memory_profiler import memory_usage
import tracemalloc
import queue

class AdditivePhylogeny:
    def __init__(self):
        pass

    def calculate_limb_length(self, distance_matrix, n, j):
        # Calculate the limb length for the leaf "j" in the distance matrix.
        limb_length = float('inf')
        if j > 0:
            i = j - 1
        else:
            j + 1

        for k in range(n):
            if k != i and k != j:
                curr_length = (distance_matrix[i][j] + distance_matrix[j][k] - distance_matrix[i][k]) // 2
                if curr_length < limb_length:
                    limb_length = curr_length
                    limb_indices = (i, k)
        
        return limb_length, limb_indices[0], limb_indices[1]

    def add_node_to_tree(self, adj, j, limb_length, i, k, x):
        # Adds a new node when needed during the phylogenetic tree construction.
        l = len(adj)
        dist = [float('inf')] * l
        parent = [-1] * l
        q = queue.Queue()
        dist[i] = 0
        q.put(i)

        while not q.empty():
            curr_node = q.get()
            # check all the neighbors of the current node
            for node, weight in adj[curr_node].items():
                if dist[node] == float('inf'):
                    dist[node] = dist[curr_node] + weight
                    parent[node] = curr_node
                    q.put(node)
                    # see if we have reached node "k"
                    if node == k:
                        prev_node = node
                        while x < dist[prev_node]:
                            curr_node = prev_node
                            prev_node = parent[curr_node]
                        # determine where to attach the new node "j"
                        if x == dist[prev_node]:
                            adj[prev_node][j] = limb_length
                            adj[j][prev_node] = limb_length
                        else:
                            # create a new node and update all the distances
                            adj.append({})
                            new_node = len(adj) - 1
                            adj[j][new_node] = limb_length
                            adj[new_node][j] = limb_length
                            del adj[prev_node][curr_node]
                            del adj[curr_node][prev_node]
                            adj[prev_node][new_node] = x - dist[prev_node]
                            adj[new_node][prev_node] = x - dist[prev_node]
                            adj[curr_node][new_node] = dist[curr_node] - x
                            adj[new_node][curr_node] = dist[curr_node] - x
                        return

    def reconstruct_phylogeny(self, D, n):
        # Reconstruct the phylogenetic tree from a distance matrix.
        adj = [{} for _ in range(n)]
        adj[0][1] = D[0][1]
        adj[1][0] = D[1][0]
        # iteratively add each leaf to the tree
        for j in range(2, n):
            limb_length, i, k = self.calculate_limb_length(D, j + 1, j)
            x = D[i][j] - limb_length
            self.add_node_to_tree(adj, j, limb_length, i, k, x)

        return adj

def dict_to_matrix(distance_dict):
    n = len(distance_dict)
    matrix = [[0] * n for _ in range(n)]
    for i, row in distance_dict.items():
        for j, value in row.items():
            matrix[i][j] = value
    return matrix

def analyze_trees_trips():
    # generate tree
    tree = ngesh.gen_tree(1.0, 0.5, max_time=10, labels="enum", num_leaves=100, seed=20)
    print(tree)

    # create distance matrix from generated tree
    distance_matrix = {}
    node_tracker = {}
    # Go through each leaf and find its distance from every other leaf
    for i, node in enumerate(tree.iter_leaves()):
        node_tracker[i] = node.name
        distance_matrix[i] = {}
        for j, other_node in enumerate(tree.iter_leaves()):
            distance_matrix[i][j] = node.get_distance(other_node)

    matrix = dict_to_matrix(distance_matrix)
    # compute running time and space complexity
    phylogeny = AdditivePhylogeny()
    start = time.time()
    tracemalloc.start()
    remade_tree = phylogeny.reconstruct_phylogeny(matrix, len(matrix))
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    end = time.time()
    print(f"Memory: {current / 1024} KBs")
    print(f"Time elapsed: {end - start} seconds")

analyze_trees_trips()
