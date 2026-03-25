import numpy as np

# Quick function for the 4x4 beam stiffness
# Based on standard Bernoulli beam theory
def beam_k(L, EI):
    # Added a quick check so we don't div by zero
    if L <= 0:
        raise ValueError("Length must be positive, man.")
        
    # Shorthand for the coefficients to save space
    c1 = 12 * EI / L**3
    c2 = 6 * EI / L**2
    c3 = 4 * EI / L
    c4 = 2 * EI / L
    
    return np.array([
        [c1,  c2, -c1,  c2],
        [c2,  c3, -c2,  c4],
        [-c1, -c2,  c1, -c2],
        [c2,  c4, -c2,  c3]
    ])

def run_analysis():
    # --- INPUTS ---
    # 2 spans (5m each), EI is constant at 1000
    L_list = [5.0, 5.0]
    EI_list = [1000.0, 1000.0]
    
    nodes = len(L_list) + 1
    total_dof = 2 * nodes
    
    # Initialize Global K
    K = np.zeros((total_dof, total_dof))

    # 1. Assembly loop
    for i in range(len(L_list)):
        k_local = beam_k(L_list[i], EI_list[i])
        idx = 2 * i
        # Splicing local into global
        K[idx:idx+4, idx:idx+4] += k_local

    # Store K_old for reaction calc later
    K_raw = K.copy()
    
    # 2. Boundary Conditions & Loading
    # Load: 10kN down at the middle node (Node 2, vertical)
    F = np.zeros(total_dof)
    F[2] = -10.0 
    
    # BCs: 0 = Fixed, 1 = Free
    # Current setup: Fixed-Pinned-Fixed
    # [V1, M1, V2, M2, V3, M3]
    supports = [0, 0, 0, 1, 0, 0] 

    # 3. Apply BCs using the big-number/penalty-ish approach 
    # Or just zeroing out rows/cols (Cleanest for small matrices)
    for i, state in enumerate(supports):
        if state == 0:
            K[i, :] = 0
            K[:, i] = 0
            K[i, i] = 1.0
            F[i] = 0

    # 4. Solve for 'd' (displacements)
    try:
        disps = np.linalg.solve(K, F)
    except np.linalg.LinAlgError:
        print("Matrix is singular! Check your supports.")
        return

    # 5. Get Reactions (R = K * d)
    reactions = K_raw @ disps

    # --- OUTPUT ---
    print("-" * 30)
    print("ANALYSIS DONE")
    print(f"Displacements: {np.round(disps, 5)}")
    print(f"Reactions: {np.round(reactions, 2)}")
    print("-" * 30)

if __name__ == "__main__":
    run_analysis()
