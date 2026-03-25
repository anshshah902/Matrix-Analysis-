import numpy as np

def get_beam_k(L, EI):
    """Calculates the 4x4 local stiffness matrix."""
    return (EI / L**3) * np.array([
        [12,    6*L,    -12,   6*L],
        [6*L,   4*L**2, -6*L,  2*L**2],
        [-12,  -6*L,    12,   -6*L],
        [6*L,   2*L**2, -6*L,  4*L**2]
    ])

def get_enl_udl(w, L):
    """Fixed-end forces for a UDL (w)."""
    return np.array([w*L/2, w*L**2/12, w*L/2, -w*L**2/12])

def run_interactive_analysis():
    print("--- BEAMSOLVER-FEA INTERACTIVE ---")
    
    try:
        nm = int(input("Enter number of beam members: "))
        L_list, EI_list = [], []
        
        for i in range(nm):
            L_list.append(float(input(f"Length of member {i+1}: ")))
            EI_list.append(float(input(f"EI of member {i+1}: ")))

        dof = 2 * (nm + 1)
        K_global = np.zeros((dof, dof))
        Pe_global = np.zeros(dof) # Equivalent nodal loads from span loads
        Pj_global = np.zeros(dof) # Direct joint loads

        # 1. Assembly & Member Loads
        for i in range(nm):
            L, EI = L_list[i], EI_list[i]
            K_global[2*i : 2*i+4, 2*i : 2*i+4] += get_beam_k(L, EI)
            
            has_udl = input(f"Does member {i+1} have a UDL? (y/n): ").lower()
            if has_udl == 'y':
                w = float(input(f"  Enter UDL value (w): "))
                Pe_global[2*i : 2*i+4] += get_enl_udl(w, L)

        # 2. Joint Loads
        print("\n--- Joint Loads (External Forces at Nodes) ---")
        for i in range(nm + 1):
            Pj_global[2*i] = float(input(f"Vertical Load at Node {i+1} (downward is negative): "))
            Pj_global[2*i+1] = float(input(f"Moment at Node {i+1} (Clockwise is positive): "))

        # Total Load Vector P = Pj - Pe
        P_total = Pj_global - Pe_global

        # 3. Boundary Conditions
        print("\n--- Boundary Conditions (0=Fixed, 1=Free) ---")
        supports = []
        for i in range(nm + 1):
            supports.append(int(input(f"Node {i+1} Vertical: ")))
            supports.append(int(input(f"Node {i+1} Rotation: ")))

        # Apply BCs (Elimination Method)
        K_mod = K_global.copy()
        P_mod = P_total.copy()
        for i, state in enumerate(supports):
            if state == 0:
                K_mod[i, :], K_mod[:, i] = 0, 0
                K_mod[i, i], P_mod[i] = 1.0, 0

        # 4. Solve
        disps = np.linalg.solve(K_mod, P_mod)
        
        # 5. Reactions: R = K*d + Pe - Pj (at constrained DOFs)
        reactions = K_global @ disps + Pe_global - Pj_global

        # --- FINAL REPORT ---
        print("\n" + "="*40)
        print("STRUCTURAL ANALYSIS RESULTS")
        print("="*40)
        print(f"Displacements (v, theta):\n{np.round(disps, 6)}")
        print(f"\nReactions (Force, Moment):\n{np.round(reactions, 3)}")
        print("="*40)

    except Exception as e:
        print(f"Error in input/calculation: {e}")

if __name__ == "__main__":
    run_interactive_analysis()300

    
