import collections
import networkx as nx
import numpy as np
from strawberryfields.apps import sample, subgraph

# ---------- Graph and Utility Classes ----------

class Graph:
    def __init__(self, vertices, edges):
        self.vertices = set(vertices)
        self.edges = set()
        self.adj = collections.defaultdict(set)
        for u, v in edges:
            self.edges.add(tuple(sorted((u, v))))
            self.adj[u].add(v)
            self.adj[v].add(u)

    def __repr__(self):
        return f"Graph(V={len(self.vertices)}, E={len(self.edges)})"

    def get_induced_edges_count(self, S):
        S = set(S)
        count = 0
        for u in S:
            for v in self.adj[u]:
                if v in S and tuple(sorted((u, v))) in self.edges:
                    count += 1
        return count // 2

    def get_degree_in_subset(self, v, S):
        return sum(1 for neighbor in self.adj[v] if neighbor in S) if v in S else 0

# ---------- Metric Functions ----------

def calculate_density(graph, S):
    return graph.get_induced_edges_count(S) / len(S) if S else 0.0

def calculate_jaccard_index(S1, S2):
    return len(S1 & S2) / len(S1 | S2) if S1 | S2 else 1.0

def calculate_q(S_collection, graph_sequence, λ):
    k = len(S_collection)
    density = sum(calculate_density(graph_sequence[i], S_collection[i]) for i in range(k))
    jaccard = sum(calculate_jaccard_index(S_collection[i], S_collection[j])
                  for i in range(k) for j in range(i + 1, k))
    return density + λ * jaccard

# ---------- GBS Integration ----------

def build_graph_with_full_nodes(edges, nodes):
    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)
    return G

def gbs_find_densest_subgraph_single(adj_matrix, n_mean=4, samples=100, k_min=3, k_max=5):
    G = nx.Graph(adj_matrix)
    s = sample.sample(adj_matrix, n_mean, samples)
    s = sample.postselect(s, k_min, k_max)
    results = subgraph.search(s, G, k_min, k_max, max_count=1)

    best_S, max_d = set(), -1.0
    for k in results:
        for d, nodes in results[k]:
            if d > max_d:
                best_S, max_d = set(nodes), d
    return best_S

# ---------- ITR Algorithm ----------

def itr_algorithm(graph_sequence, λ):
    k = len(graph_sequence)
    if k == 0:
        return []

    V = graph_sequence[0].vertices
    print("\n--- Generating Candidate using GBS Samples ---")

    adj_matrices = [
        nx.to_numpy_array(build_graph_with_full_nodes(g.edges, g.vertices), nodelist=list(g.vertices))
        for g in graph_sequence
    ]
    S = [gbs_find_densest_subgraph_single(adj) for adj in adj_matrices]

    current_q = calculate_q(S, graph_sequence, λ)
    print("--- GBS Candidate Generation Complete ---")
    print(f"\nInitial Q = {current_q:.4f}")

    changed, iteration = True, 0
    max_iter = 50

    while changed and iteration < max_iter:
        changed = False
        iteration += 1
        for i in range(k):
            C = set(V)
            candidates = []

            for _ in range(len(V) - 1):
                best_remove, max_q = None, -1.0
                for v in C:
                    C_prime = C - {v}
                    S_temp = list(S)
                    S_temp[i] = C_prime
                    q_val = calculate_q(S_temp, graph_sequence, λ)
                    if q_val > max_q:
                        best_remove, max_q = v, q_val
                if best_remove:
                    C.remove(best_remove)
                    candidates.append((set(C), max_q))

            if candidates:
                best_C, best_q = max(candidates, key=lambda x: x[1])
                if best_q > current_q:
                    S[i] = best_C
                    current_q = best_q
                    changed = True

    return S

# ---------- Example Usage ----------

'''
nodes = [chr(i) for i in range(97, 109)]

G1 = Graph(nodes, [('a', 'b'), ('b', 'c'), ('c', 'd'), ('d', 'e'), ('e', 'f'), ('f', 'g'), ('g', 'h'), ('h', 'i'), ('i', 'j'), ('j', 'k'), ('k', 'l'), ('l', 'a')])
G2 = Graph(nodes, [('a', 'c'), ('b', 'd'), ('c', 'e'), ('d', 'f'), ('e', 'g'), ('f', 'h'), ('g', 'i'), ('h', 'j'), ('i', 'k'), ('j', 'l')])
G3 = Graph(nodes, [('a', 'b'), ('b', 'd'), ('c', 'e'), ('e', 'f'), ('f', 'g'), ('g', 'i'), ('h', 'j'), ('j', 'k'), ('k', 'l'), ('l', 'a')])
G4 = Graph(nodes, [('a', 'e'), ('b', 'f'), ('c', 'g'), ('d', 'h'), ('e', 'i'), ('f', 'j'), ('g', 'k'), ('h', 'l'), ('i', 'a'), ('j', 'b')])
G5 = Graph(nodes, [('a', 'd'), ('b', 'e'), ('c', 'f'), ('d', 'g'), ('e', 'h'), ('f', 'i'), ('g', 'j'), ('h', 'k'), ('i', 'l'), ('j', 'a')])

graph_sequence = [G1, G2, G3, G4, G5]
lambda_param = 0.3
'''

nodes = ['a','b','c','d','e','f','g','h','i','j']

G1 = Graph(nodes, [('a','b'), ('a','c'), ('a','d'), ('b','c'), ('b','d'), ('c','d'), ('a','e')])
G2 = Graph(nodes, [('e','f'), ('f','g'), ('g','e'), ('e','h')])
G3 = Graph(nodes, [('h','i'), ('i','j'), ('j','h'), ('i','a')])
G4 = Graph(nodes, [('c','d'), ('c','f'), ('f','g'), ('d','g')])
G5 = Graph(nodes, [('a','j'), ('j','e'), ('e','a')])

graph_sequence = [G1, G2, G3, G4, G5]
λ = 0.3

print("Running JWDS ITR algorithm...")
S_result = itr_algorithm(graph_sequence, λ)

print("\n--- ITR Algorithm Results ---")
for i, S_i in enumerate(S_result):
    d = calculate_density(graph_sequence[i], S_i)
    print(f"S_{i+1}: {sorted(S_i)}, Density: {d:.2f}")

final_q = calculate_q(S_result, graph_sequence, λ)
print(f"Final Q Score (lambda={λ}): {final_q:.4f}")

print("\nPairwise Jaccard Indices (Final S):")
for i in range(len(S_result)):
    for j in range(i + 1, len(S_result)):
        jaccard = calculate_jaccard_index(S_result[i], S_result[j])
        print(f"J(S_{i+1}, S_{j+1}): {jaccard:.2f}")
