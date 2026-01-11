#!/usr/bin/env python

import time
import tracemalloc
import ngesh
from degeneratetrips import AdditivePhylogeny
from Benchmarking.neighborjoining import NeighborJoining


def build_distance_dict(tree):
    """Build a dict-of-dicts distance matrix (float) from a tree's leaves."""
    distance_matrix = {}
    leaves = list(tree.iter_leaves())
    for i, node in enumerate(leaves):
        distance_matrix[i] = {}
        for j, other_node in enumerate(leaves):
            distance_matrix[i][j] = float(node.get_distance(other_node))
    return distance_matrix


def dict_to_matrix(D):
    """Convert dict-of-dicts to a list-of-lists matrix (floats)."""
    n = len(D)
    M = [[0.0] * n for _ in range(n)]
    for i, row in D.items():
        for j, val in row.items():
            M[i][j] = float(val)
    return M


def run_additive_phylogeny(D_matrix):
    print("\n=== Additive Phylogeny (degeneratetrips.py) ===")
    ap = AdditivePhylogeny()
    start = time.time()
    tracemalloc.start()

    # ORIGINAL signature: reconstruct_phylogeny(D, n)
    adj_ap = ap.reconstruct_phylogeny(D_matrix, len(D_matrix))

    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    elapsed = time.time() - start
    print(f"Current memory: {current / (1024 * 1024):.3f} MB")
    print(f"Peak memory:    {peak / (1024 * 1024):.3f} MB")
    print(f"Time elapsed:   {elapsed:.4f} s")
    return adj_ap


def run_neighbor_joining(D_dict):
    print("\n=== Neighbor Joining (neighborjoining.py) ===")
    nj = NeighborJoining()
    start = time.time()
    tracemalloc.start()
    T_nj = nj.neighbor_join(D_dict)
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    elapsed = time.time() - start
    print(f"Current memory: {current / (1024 * 1024):.3f} MB")
    print(f"Peak memory:    {peak / (1024 * 1024):.3f} MB")
    print(f"Time elapsed:   {elapsed:.4f} s")

    try:
        ok = nj.check_tree(D_dict, T_nj)
        print(f"Tree matches original distances: {ok}")
    except Exception as e:
        print(f"Verification skipped (error): {e}")
    return T_nj


def main():
    # Generate a single input tree to be shared by both algorithms
    tree = ngesh.gen_tree(1.0, 0.5, max_time=10, labels="enum", num_leaves=8, seed=20)
    print(tree)

    D_dict = build_distance_dict(tree)   # for Neighbor Joining
    D_matrix = dict_to_matrix(D_dict)    # for Additive Phylogeny

    # Run both algorithms on the SAME input data
    _ = run_additive_phylogeny(D_matrix)
    _ = run_neighbor_joining(D_dict)

if __name__ == "__main__":
    main()