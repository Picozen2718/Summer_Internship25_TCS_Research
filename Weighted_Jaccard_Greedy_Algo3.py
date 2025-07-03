import collections

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
        """Calculates the number of induced edges within a given subset of vertices."""
        if not S_subset:
            return 0
        
        S_set = set(S_subset)
        induced_edges = 0
        for u in S_set:
            for v_neighbor in self.adj[u]:
                if v_neighbor in S_set:
                    # Count each edge once (e.g., if u-v is counted, don't count v-u)
                    if tuple(sorted((u, v_neighbor))) in self.edges:
                        induced_edges += 1
        return induced_edges // 2 # Each edge (u,v) is visited twice (u's neighbors, v's neighbors)

    def get_degree_in_subset(self, vertex, S_subset):
        """Calculates the degree of a vertex within an induced subgraph S_subset."""
        if vertex not in S_subset:
            return 0
        
        S_set = set(S_subset)
        degree = 0
        for neighbor in self.adj[vertex]:
            if neighbor in S_set:
                degree += 1
        return degree

def calculate_density(graph, S_i):
    """
    Calculates the density of a subgraph induced by vertex set S_i.
    Density d(S_i) = |E(S_i)| / |S_i|
    where E(S_i) is the set of edges with both endpoints in S_i.
    """
    if not S_i:
        return 0.0
    
    induced_edges = graph.get_induced_edges_count(S_i)
    return induced_edges / len(S_i)

def calculate_jaccard_index(set1, set2):
    """
    Calculates the Jaccard index between two sets.
    J(S, T) = |S intersect T| / |S union T|
    """
    if not set1 and not set2:
        return 1.0  # Both empty, considered perfectly similar
    
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    
    if union == 0:
        return 0.0  # Avoid division by zero if both sets are empty
    
    return intersection / union

def calculate_q(S_collection, graph_sequence, lambda_param):
    """
    Calculates the objective function q(S; lambda) for JWDS.
    q(S; lambda) = sum_i d(S_i) + lambda * sum_{i<j} J(S_i, S_j)
    """
    k = len(S_collection)
    sum_of_densities = 0.0
    for i in range(k):
        sum_of_densities += calculate_density(graph_sequence[i], S_collection[i])
    
    sum_of_jaccard_indices = 0.0
    for i in range(k):
        for j in range(i + 1, k):
            sum_of_jaccard_indices += calculate_jaccard_index(S_collection[i], S_collection[j])
            
    return sum_of_densities + lambda_param * sum_of_jaccard_indices

def grd_algorithm(graph_sequence, lambda_param):
    """
    Implements the GRD algorithm (Algorithm 3) for Jaccard-Weighted Densest Subgraphs (JWDS).
    This is a greedy peeling algorithm that iteratively removes the vertex (from any snapshot)
    that globally maximizes the q score upon its removal.

    Args:
        graph_sequence (list): A list of graph snapshots (Graph objects).
        lambda_param (float): The weighting parameter for Jaccard indices in the objective function.

    Returns:
        list: A list of vertex sets [S1, S2, ..., Sk] representing the collection
              of dense subgraphs that maximize the q score.
    """
    k = len(graph_sequence)
    if k == 0:
        return []

    V = graph_sequence[0].vertices # Assuming all graphs have the same vertex set

    # Step 1: Initialize S (Line 1 in Algo 3)
    # S_i = V for all i
    S = [set(V) for _ in range(k)] 
    
    # Initialize best_tested_S and its score
    best_overall_S = [set(s) for s in S] # Deep copy
    max_overall_q = calculate_q(best_overall_S, graph_sequence, lambda_param)
    
    # print(f"Starting GRD with initial Q = {max_overall_q:.4f}")

    # Step 2: While there are vertices (Line 2 in Algo 3)
    # This loop runs N * K times in the worst case (removing one vertex at a time)
    # The actual number of removals is sum(|S_i| starting with V), so N*K removals total
    
    # We maintain a list of all vertices that are still "active" in any S_i.
    # A more precise way to check "while there are vertices" would be if any S_i is non-empty.
    # The loop will implicitly terminate when all S_i become empty.
    
    all_vertices_present = True
    iteration = 0
    while all_vertices_present:
        all_vertices_present = False # Will be set to True if any S_i is not empty
        iteration += 1
        # print(f"\nIteration {iteration}. Current Q = {max_overall_q:.4f}")

        best_vertex_to_remove = None
        best_snapshot_index = -1
        max_q_if_removed = -1.0 # Initialize with a very low value

        # Line 3 in Algo 3: (u,j) <- arg max q(...)
        # Iterate through all snapshots and all vertices within them to find the best removal
        for j in range(k): # iterate over snapshots
            if not S[j]: # If S_j is already empty, no vertex to remove
                continue
            
            all_vertices_present = True # At least one snapshot has vertices

            for v in S[j]: # iterate over vertices currently in S_j
                # Create a hypothetical S' where S_j is replaced by S_j \ {v}
                hypothetical_S_prime_collection = list(S) # Create a copy of the current collection
                hypothetical_S_prime_collection[j] = S[j] - {v} # Update the specific S_j
                
                # Calculate the global q score for this hypothetical S'
                q_val_if_removed = calculate_q(hypothetical_S_prime_collection, graph_sequence, lambda_param)

                # Check if this removal yields the maximum q score so far
                if q_val_if_removed > max_q_if_removed:
                    max_q_if_removed = q_val_if_removed
                    best_vertex_to_remove = v
                    best_snapshot_index = j
        
        # If no vertex could be removed to improve the score (or if all sets are empty), break
        if best_vertex_to_remove is None or max_q_if_removed <= calculate_q(S, graph_sequence, lambda_param):
            break # No improvement possible or no more vertices to remove
            
        # Line 4 in Algo 3: S_j <- S_j \ {u}
        S[best_snapshot_index].remove(best_vertex_to_remove)
        # print(f"  Removed {best_vertex_to_remove} from S_{best_snapshot_index+1}. Hypothetical Q = {max_q_if_removed:.4f}")

        # Update best_overall_S if current_q is better
        # Note: The paper says "return best tested S". This implies we keep track
        # of the best S found at any point *during* the peeling process.
        current_q = calculate_q(S, graph_sequence, lambda_param)
        if current_q > max_overall_q:
            max_overall_q = current_q
            best_overall_S = [set(s) for s in S] # Deep copy

    # Line 5 in Algo 3: return best tested S
    return best_overall_S

# --- Example Usage ---
G1 = Graph(
    vertices=['a', 'b', 'c', 'd', 'e', 'f'],
    edges=[('a', 'b'), ('a', 'd'), ('b', 'd')]
)
G2 = Graph(
    vertices=['a', 'b', 'c', 'd', 'e', 'f'],
    edges=[('a', 'b'), ('a', 'c'), ('b', 'c'), ('b', 'f'), ('c', 'f'), ('a', 'f')]
)
G3 = Graph(
    vertices=['a', 'b', 'c', 'd', 'e', 'f'],
    edges=[('a', 'b'), ('a', 'd'), ('b', 'd'), ('b', 'e'), ('d', 'e'), ('a', 'e'), ('a', 'f'), ('e', 'f')]
)

graph_sequence_example = [G1, G2, G3]
lambda_example = 0.3

# Run the GRD algorithm
print("Running JWDS GRD algorithm...")
result_subgraphs_grd = grd_algorithm(graph_sequence_example, lambda_example)

print("\n--- GRD Algorithm Results ---")
for i, S_i in enumerate(result_subgraphs_grd):
    density = calculate_density(graph_sequence_example[i], S_i)
    print(f"S_{i+1}: {sorted(list(S_i))}, Density: {density:.2f}")

final_q_score_grd = calculate_q(result_subgraphs_grd, graph_sequence_example, lambda_example)
print(f"Final Q Score (lambda={lambda_example}): {final_q_score_grd:.4f}")

print("\nPairwise Jaccard Indices (Final S):")
k = len(result_subgraphs_grd)
for i in range(k):
    for j in range(i + 1, k):
        jaccard = calculate_jaccard_index(result_subgraphs_grd[i], result_subgraphs_grd[j])
        print(f"J(S_{i+1}, S_{j+1}): {jaccard:.2f}")

