import numpy as np

# --- LOAD COMPUTATION HELPERS ---
# Note: Using standard fixed-end moment formulas (Clockwise = Positive)

def get_point_load_fef(P, a, L):
    """Point load P at distance 'a' from left end."""
    b = L - a
    # Pre-calculate powers to keep the formulas clean
    L2 = L**2
    
    # Standard FEM formulas for fixed-end forces
    f1 = (P * b**2 * (3*a + b)) / L**3
    m1 = (P * a * b**2) / L2
    f2 = (P * a**2 * (3*b + a)) / L**3
    m2 = -(P * a**2 * b) / L2
    
    return np.array([f1, m1, f2, m2])

def get_udl_fef(w, L):
    """Standard UDL 'w' over full span L."""
    f = w * L / 2.0
    m = w * L**2 / 12.0
    # Returns [V1, M1, V2, M2]
    return np.array([f, m, f, -m])

# --- MAIN ASSEMBLY ---

def build_load_vector(spans, joint_loads):
    """
    spans: list of spans, each span is a dict with 'L' and 'loads'
    joint_loads: list of direct nodal forces [V1, M1, V2, M2...]
    """
    num_nodes = len(spans) + 1
    total_dof = 2 * num_nodes
    
    # Equivalent nodal forces (from member loads)
    ae = np.zeros(total_dof)
    
    for i, span in enumerate(spans):
        L = span['L']
        local_ae = np.zeros(4)
        
        # Sum up all loads on this specific member
        for load in span['loads']:
            l_type = load.get('type')
            if l_type == 'point':
                local_ae += get_point_load_fef(load['P'], load['a'], L)
            elif l_type == 'udl':
                local_ae += get_udl_fef(load['w'], L)
            else:
                print(f"Warning: Unknown load type '{l_type}' on span {i+1}")

        # Map to global vector (2 DOFs per node)
        idx = 2 * i
        # IMPORTANT: We subtract 'ae' from 'aj' later, or just -= here
        ae[idx : idx + 4] -= local_ae

    # Convert joint loads to numpy array for math
    aj = np.array(joint_loads, dtype=float)
    
    # Combined load vector: ac = aj + ae (nodal + equiv member loads)
    ac = aj + ae
    
    return ae, ac

# --- TEST CASE / RUNNER ---
if __name__ == "__main__":
    # Defining the structure data
    my_spans = [
        {
            'L': 5.0,
            'loads': [
                {'type': 'point', 'P': 10, 'a': 2.5}, # Mid-point load
                {'type': 'udl', 'w': 2}               # Background UDL
            ]
        },
        {
            'L': 4.0,
            'loads': [{'type': 'udl', 'w': 5}]
        }
    ]

    # No direct joint loads for now
    nodes = len(my_spans) + 1
    p_nodal = [0] * (2 * nodes)

    ae_vec, ac_vec = build_load_vector(my_spans, p_nodal)

    print("-" * 40)
    print("GLOBAL LOAD VECTOR (ac):")
    print(np.round(ac_vec, 3))
    print("-" * 40)
    # Quick sanity check: Total vertical load vs sum of Reactions
    print(f"Total applied vertical load: {10 + (2*5) + (5*4)} kN")
