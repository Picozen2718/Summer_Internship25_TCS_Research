import collections

def calculate_density(graph, S_i):
    """
    Calculates the density of a subgraph induced by vertex set S_i.
    Density d(S_i) = |E(S_i)| / |S_i|
    where E(S_i) is the set of edges with both endpoints in S_i.
    """
    if not S_i:
        return 0.0
    
    # Create a set for efficient lookup
    S_i_set = set(S_i)
    
    induced_edges = 0
    # Iterate through all edges in the graph
    for u, v in graph.edges:
        if u in S_i_set and v in S_i_set:
            induced_edges += 1
            
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
        return 0.0  # Avoid division by zero if both sets are empty but not a special case
    
    return intersection / union

def satisfies_jaccard_constraint(S, alpha):
    """
    Checks if all pairwise Jaccard indices in the collection of subgraphs S
    are at least alpha.
    S is a list of vertex sets [S1, S2, ..., Sk]
    """
    k = len(S)
    for i in range(k):
        for j in range(i + 1, k):
            if calculate_jaccard_index(S[i], S[j]) < alpha:
                return False
    return True

def hard_jcds_algorithm(graph_sequence, alpha):
    """
    Implements the HARD algorithm for Jaccard-Constrained Densest Subgraphs (JCDS).

    Args:
        graph_sequence (list): A list of graph snapshots, where each graph is
                               represented as an object with a 'vertices' attribute (list)
                               and an 'edges' attribute (list of (u, v) tuples).
                               For simplicity, we assume 'vertices' is the same for all graphs.
        alpha (float): The minimum required Jaccard index between any two subgraphs.

    Returns:
        list: A list of vertex sets [S1, S2, ..., Sk] representing the collection
              of dense subgraphs that satisfy the Jaccard constraint.
    """
    k = len(graph_sequence)
    if k == 0:
        return []

    V = graph_sequence[0].vertices # Assuming all graphs have the same vertex set

    # Step 1: Initialize S as the densest common subgraph (alpha = 1 equivalent)
    # This is a simplified initialization. In a full implementation, you'd find
    # the actual densest common subgraph. For this example, we'll start with all
    # vertices, and the iterative process will refine it.
    # A more robust DCS initialization would involve building a weighted graph
    # where edge weights are counts of occurrences across snapshots, then finding
    # its densest subgraph.
    # For now, let's assume a dummy DCS function or a simple starting point.
    
    # A basic starting point: all vertices in each snapshot.
    # This might not be optimal but allows the algorithm to start iterating.
    S = [set(V) for _ in range(k)] 
    
    # Calculate initial total density
    current_total_density = sum(calculate_density(graph_sequence[i], S[i]) for i in range(k))
    
    # Store initial Jaccard indices to avoid recalculating frequently
    # jaccard_matrix[i][j] stores J(S[i], S[j])
    jaccard_matrix = collections.defaultdict(lambda: collections.defaultdict(float))
    for i in range(k):
        for j in range(i + 1, k):
            jaccard_matrix[i][j] = calculate_jaccard_index(S[i], S[j])
            jaccard_matrix[j][i] = jaccard_matrix[i][j] # Symmetric

    # Step 2: Iterate while changes in the score occur
    changed = True
    iteration = 0
    max_iterations = 100 # Add a safeguard to prevent infinite loops in cases where convergence is slow

    while changed and iteration < max_iterations:
        changed = False
        iteration += 1
        # print(f"Iteration {iteration}: Current Total Density = {current_total_density:.4f}")

        # Iterate through each snapshot and each vertex
        for i in range(k): # foreach i = 1, ..., k do
            for v in V:    # foreach v in V do
                S_i_copy = set(S[i]) # Create a copy to modify
                original_S_i = set(S[i]) # Keep original for comparison

                # Option 1: Remove vertex v from S_i
                if v in original_S_i: # if v in S_i then
                    S_i_candidate = S_i_copy - {v} # S' with S_i replaced with S_i \ {v}
                    
                    # Create a temporary S' for evaluation
                    temp_S = list(S)
                    temp_S[i] = S_i_candidate

                    # Check Jaccard constraint and density gain
                    is_jaccard_satisfied = True
                    temp_jaccard_matrix = jaccard_matrix.copy() # Start with current matrix

                    # Only check Jaccard constraints involving the modified S_i_candidate
                    for j in range(k):
                        if i == j:
                            continue
                        new_jaccard = calculate_jaccard_index(S_i_candidate, temp_S[j])
                        if new_jaccard < alpha:
                            is_jaccard_satisfied = False
                            break
                        # Update temporary jaccard matrix for potential acceptance
                        temp_jaccard_matrix[i][j] = new_jaccard
                        temp_jaccard_matrix[j][i] = new_jaccard

                    if is_jaccard_satisfied:
                        new_total_density = sum(calculate_density(graph_sequence[idx], temp_S[idx]) for idx in range(k))
                        if new_total_density > current_total_density: # and d(S') > d(S)
                            S[i] = S_i_candidate # S <- S'
                            current_total_density = new_total_density
                            jaccard_matrix = temp_jaccard_matrix # Update Jaccard matrix permanently
                            changed = True
                            # print(f"  Removed {v} from G{i+1}. New density: {current_total_density:.4f}")
                            # Continue to next vertex after change
                            continue # Optimization: if a change occurred, restart inner loop to re-evaluate with updated S[i]

                # Option 2: Add vertex v to S_i
                if v not in original_S_i: # else (if v not in S_i)
                    S_i_candidate = S_i_copy | {v} # S' with S_i replaced with S_i U {v}
                    
                    # Create a temporary S' for evaluation
                    temp_S = list(S)
                    temp_S[i] = S_i_candidate

                    # Check Jaccard constraint and density gain
                    is_jaccard_satisfied = True
                    temp_jaccard_matrix = jaccard_matrix.copy() # Start with current matrix

                    # Only check Jaccard constraints involving the modified S_i_candidate
                    for j in range(k):
                        if i == j:
                            continue
                        new_jaccard = calculate_jaccard_index(S_i_candidate, temp_S[j])
                        if new_jaccard < alpha:
                            is_jaccard_satisfied = False
                            break
                        # Update temporary jaccard matrix for potential acceptance
                        temp_jaccard_matrix[i][j] = new_jaccard
                        temp_jaccard_matrix[j][i] = new_jaccard

                    if is_jaccard_satisfied:
                        new_total_density = sum(calculate_density(graph_sequence[idx], temp_S[idx]) for idx in range(k))
                        if new_total_density > current_total_density: # and d(S') > d(S)
                            S[i] = S_i_candidate # S <- S'
                            current_total_density = new_total_density
                            jaccard_matrix = temp_jaccard_matrix # Update Jaccard matrix permanently
                            changed = True
                            # print(f"  Added {v} to G{i+1}. New density: {current_total_density:.4f}")
                            # Continue to next vertex after change
                            continue # Optimization: if a change occurred, restart inner loop to re-evaluate with updated S[i]

    return S

# --- Example Usage ---

# Define a simple Graph class for demonstration
class Graph:
    def __init__(self, vertices, edges):
        self.vertices = set(vertices)
        self.edges = set()
        for u, v in edges:
            # Ensure edges are stored consistently (e.g., (min(u,v), max(u,v)))
            self.edges.add(tuple(sorted((u, v))))

    def __repr__(self):
        return f"Graph(V={len(self.vertices)}, E={len(self.edges)})"

# --- Example Usage ---
# G1 = ({a,b,c,d,e,f}, {(a,b), (a,d), (b,d), (b,e)})
G1 = Graph(
    vertices=['a', 'b', 'c', 'd', 'e', 'f'],
    edges=[('a', 'b'), ('a', 'd'), ('b', 'd')] # Only these for density 3/3 for (a,b,d)
)
# G2 = ({a,b,c,d,e,f}, {(a,b), (a,c), (b,c), (b,f), (c,f), (a,f)})
G2 = Graph(
    vertices=['a', 'b', 'c', 'd', 'e', 'f'],
    edges=[('a', 'b'), ('a', 'c'), ('b', 'c'), ('b', 'f'), ('c', 'f'), ('a', 'f')]
)
# G3 = ({a,b,c,d,e,f}, {(a,b), (a,d), (b,d), (b,e), (d,e), (a,e), (a,f), (e,f)})
G3 = Graph(
    vertices=['a', 'b', 'c', 'd', 'e', 'f'],
    edges=[('a', 'b'), ('a', 'd'), ('b', 'd'), ('b', 'e'), ('d', 'e'), ('a', 'e'), ('a', 'f'), ('e', 'f')]
)

graph_sequence_example = [G1, G2, G3]
alpha_example = 0.6

# Run the algorithm
print("Running JCDS HARD algorithm...")
result_subgraphs = hard_jcds_algorithm(graph_sequence_example, alpha_example)

print("\n--- Algorithm Results ---")
for i, S_i in enumerate(result_subgraphs):
    density = calculate_density(graph_sequence_example[i], S_i)
    print(f"S_{i+1}: {sorted(list(S_i))}, Density: {density:.2f}")

total_density_result = sum(calculate_density(graph_sequence_example[i], result_subgraphs[i]) for i in range(len(result_subgraphs)))
print(f"Total Density: {total_density_result:.2f}")

print("\nPairwise Jaccard Indices:")
k = len(result_subgraphs)
for i in range(k):
    for j in range(i + 1, k):
        jaccard = calculate_jaccard_index(result_subgraphs[i], result_subgraphs[j])
        print(f"J(S_{i+1}, S_{j+1}): {jaccard:.2f}")

is_satisfied = satisfies_jaccard_constraint(result_subgraphs, alpha_example)
print(f"Jaccard constraint (alpha={alpha_example}) satisfied: {is_satisfied}")

