# Matrix-Analysis-
​BeamSolver-FEA: Python Structural Engine
​BeamSolver-FEA is a computational tool designed to solve multi-span beam problems using the Direct Stiffness Method (DSM). It automates the transition from complex physical loads to exact nodal displacements and support reactions, providing a transparent "whitebox" alternative to commercial FEA software.
​Core Capabilities
​Automated Load Processing: Converts Point, UDL, and UVL loads into Equivalent Nodal Force vectors using Fixed-End Moment (FEM) theory.
​Matrix Assembly: Generates 4 \times 4 local stiffness matrices based on Euler-Bernoulli beam theory and maps them into a 2n \times 2n Global Stiffness Matrix (K).
​Boundary Condition Solver: Employs the Elimination Method to handle various support configurations (Fixed, Pinned, Roller) and solves the system equilibrium equation Kd = F.
​Numerical Precision: Utilizes NumPy for optimized, high-performance linear algebra and LU-decomposition.
​Technical Workflow
​The engine discretizes a continuous beam into n elements. For each element, it calculates the flexural stiffness relative to its length (L) and rigidity (EI). After assembling the global matrix, the solver partitions the system based on user-defined boundary conditions to find unknown displacements and back-calculate reactions.
