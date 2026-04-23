using Pkg
Pkg.activate(".")
Pkg.instantiate()
include("Polyhedron.jl")
using .Poly

A = [0 1;
    -0.1 0.95] 
B = [0;
     1;;]

C = [0.025 0.025]

T = 0.4

delta = 0.5

A_exp = hcat(A, [0 0; 0 0])
A_exp = vcat(A_exp, [0 -T 1 T; 0 0 0 delta])

B_exp = vcat(B, [0;0])

E_exp = vcat(zeros(3), [1-delta])

Sx = [1/2 0
      0 1/2;
     -1/2 0;
      0 -1/2;]
Sv = [1/10;
      -1/10;]

Sw = [1/10;
      -1/10;]

Sx = vcat(Sx, zeros(4, 2))
Sv = vcat(vcat(zeros(4, 1), Sv) , zeros(2, 1))
Sw = vcat(zeros(6, 1), Sw)

S = hcat(hcat(Sx, Sv), Sw)

# Fr para a referência
R = [1/2;
      -1/2;;]

#testando para d = 1
d = 1

#primeiro teste
result_dis = Poly.finding_L_pinvariant_segref_delay(A_exp, B_exp, E_exp, S, R, d, time=30,lf=18,  lambda=0.99)

F = result_dis["F"]
G = result_dis["G"]
K = result_dis["K"]
H = result_dis["H"]
L = result_dis["L"]
M = result_dis["M"]
N = result_dis["N"]
P = result_dis["P"]
T = result_dis["T"]
J = result_dis["J"]

using DelimitedFiles

open("matrices_teste1.txt", "w") do f
	println(f, "teste 1")
    println(f, "F")
    writedlm(f, F)
    println(f, "\nG")
    writedlm(f, G)
    println(f, "\nK")
    writedlm(f, K)
    println(f, "\nH")
    writedlm(f, H)
    println(f, "\nL")
    writedlm(f, L)
    println(f, "\nM")
    writedlm(f, M)
    println(f, "\nN")
    writedlm(f, N)
    println(f, "\nP")
    writedlm(f, P)
    println(f, "\nT")
    writedlm(f, T)
    println(f, "\nJ")
    writedlm(f, J)
    println(f, "====================================")
end

F = result_dis["F"]
G = result_dis["G"]
K = result_dis["K"]
H = result_dis["H"]
L = result_dis["L"]
M = result_dis["M"]
N = result_dis["N"]
P = result_dis["P"]
T = result_dis["T"]
J = result_dis["J"]

result_dis = Poly.finding_L_pinvariant_segref_delay(A_exp, B_exp, E_exp, S, R, d, time=30,lf=18,  lambda=0.9)

open("matrices_teste1.txt", "w") do f
	println(f, "teste 2")
    println(f, "F")
    writedlm(f, F)
    println(f, "\nG")
    writedlm(f, G)
    println(f, "\nK")
    writedlm(f, K)
    println(f, "\nH")
    writedlm(f, H)
    println(f, "\nL")
    writedlm(f, L)
    println(f, "\nM")
    writedlm(f, M)
    println(f, "\nN")
    writedlm(f, N)
    println(f, "\nP")
    writedlm(f, P)
    println(f, "\nT")
    writedlm(f, T)
    println(f, "\nJ")
    writedlm(f, J)
    println(f, "====================================")
end