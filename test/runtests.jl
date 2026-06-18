using Polyhedron 
using Test

A = [1 1; 0 1]
B = [2; 1]
X = [0.8 0; 0 1; -1 0; 0 -1]
U = [1.2; -1.5]
C = [1 0]

@test is_pinvariant(A, B, C, U, X) < 1
