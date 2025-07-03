import collections
import time

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

    def get_induced_edges_count(self, S_subset):
        S = set(S_subset)
        count = 0
        for u in S:
            for v in self.adj[u]:
                if v in S and tuple(sorted((u, v))) in self.edges:
                    count += 1
        return count // 2

    def get_degree_in_subset(self, vertex, S_subset):
        if vertex not in S_subset:
            return 0
        return sum(1 for neighbor in self.adj[vertex] if neighbor in S_subset)

# ---------- Metric Functions ----------

def calculate_density(graph, S_i):
    return graph.get_induced_edges_count(S_i) / len(S_i) if S_i else 0.0

def calculate_jaccard_index(set1, set2):
    return len(set1 & set2) / len(set1 | set2) if set1 | set2 else 1.0

def calculate_q(S_collection, graph_sequence, lambda_param):
    k = len(S_collection)
    density_sum = sum(calculate_density(graph_sequence[i], S_collection[i]) for i in range(k))
    jaccard_sum = sum(calculate_jaccard_index(S_collection[i], S_collection[j])
                      for i in range(k) for j in range(i + 1, k))
    return density_sum + lambda_param * jaccard_sum

# ---------- Densest Subgraph Initialization ----------

def find_densest_subgraph_single(graph):
    current = set(graph.vertices)
    candidates = [(set(current), calculate_density(graph, current))]

    while current:
        degrees = {v: graph.get_degree_in_subset(v, current) for v in current}
        v_remove = min(degrees, key=degrees.get)
        current.remove(v_remove)
        if current:
            d = calculate_density(graph, current)
            candidates.append((set(current), d))

    return max(candidates, key=lambda item: item[1])[0]

def find_densest_common_subgraph(graph_sequence):
    common = set.intersection(*(G.vertices for G in graph_sequence))
    current = set(common)
    candidates = [(set(current), sum(calculate_density(G, current) for G in graph_sequence))]

    while current:
        degrees = {v: sum(G.get_degree_in_subset(v, current) for G in graph_sequence)
                   for v in current}
        v_remove = min(degrees, key=degrees.get)
        current.remove(v_remove)
        if current:
            d_sum = sum(calculate_density(G, current) for G in graph_sequence)
            candidates.append((set(current), d_sum))

    return max(candidates, key=lambda item: item[1])[0]

# ---------- ITR Algorithm ----------

def itr_algorithm(graph_sequence, lambda_param):
    k = len(graph_sequence)
    if k == 0:
        return [], 0.0, 0

    V = graph_sequence[0].vertices

    print("\n--- Initializing with Classical Densest Subgraphs ---")
    start_time = time.time()

    common_init = find_densest_common_subgraph(graph_sequence)
    S_common = [set(common_init) for _ in range(k)]
    q_common = calculate_q(S_common, graph_sequence, lambda_param)

    separate_init = [find_densest_subgraph_single(G) for G in graph_sequence]
    S_separate = [set(s) for s in separate_init]
    q_separate = calculate_q(S_separate, graph_sequence, lambda_param)

    if q_common >= q_separate:
        print("Using common densest subgraph initialization (S)")
        S = S_common
        q_score = q_common
    else:
        print("Using separate densest subgraph initialization (S')")
        S = S_separate
        q_score = q_separate

    init_time = time.time() - start_time
    print(f"--- Classical Initialization Complete (Time: {init_time:.2f} sec) ---")
    print(f"\nInitial Q Score from Classical Initialization: Q = {q_score:.4f}")

    changed = True
    iteration = 0
    max_iterations = 50

    while changed and iteration < max_iterations:
        changed = False
        iteration += 1
        for i in range(k):
            current_C = set(V)
            best_C = set(S[i])
            best_q = calculate_q(S, graph_sequence, lambda_param)
            candidates = []

            for _ in range(len(V) - 1):
                best_remove = None
                max_q = -1.0
                for v in sorted(current_C):
                    C_prime = current_C - {v}
                    S_temp = list(S)
                    S_temp[i] = C_prime
                    q_val = calculate_q(S_temp, graph_sequence, lambda_param)
                    if q_val > max_q:
                        max_q = q_val
                        best_remove = v

                if best_remove:
                    current_C.remove(best_remove)
                    candidates.append((set(current_C), max_q))

            if candidates:
                best_C_peel, best_q_peel = max(candidates, key=lambda x: x[1])
                if best_q_peel > q_score:
                    S[i] = best_C_peel
                    q_score = best_q_peel
                    changed = True

    return S, init_time, iteration

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
lambda_param = 0.3

print("Running JWDS ITR algorithm with Classical initialization...")
S_itr, init_time, iterations = itr_algorithm(graph_sequence, lambda_param)

print("\n--- ITR Algorithm Results ---")
for i, S_i in enumerate(S_itr):
    d = calculate_density(graph_sequence[i], S_i)
    print(f"S_{i+1}: {sorted(S_i)}, Density: {d:.2f}")

final_q = calculate_q(S_itr, graph_sequence, lambda_param)
print(f"Final Q Score (lambda={lambda_param}): {final_q:.4f}")
print(f"Total Classical Initialization Time: {init_time:.2f} seconds")
print(f"Total Iterations Until Convergence: {iterations}")

print("\nPairwise Jaccard Indices (Final S):")
for i in range(len(S_itr)):
    for j in range(i + 1, len(S_itr)):
        jaccard = calculate_jaccard_index(S_itr[i], S_itr[j])
        print(f"J(S_{i+1}, S_{j+1}): {jaccard:.2f}")
