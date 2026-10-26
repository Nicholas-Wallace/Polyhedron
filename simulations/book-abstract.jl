using Polyhedron 

A = [1 1; 0 1]
B = [2; 1;;]
C = [1 0]

X = [0.8 0;
     0 1;
     -1 0;
     0 -1;]
U = [1.2; -1.5;;]

sys = LTISystem(A, B, C)
result = is_pinvariant(sys, U, X) 
