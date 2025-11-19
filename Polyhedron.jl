module Polyhedron

using LinearAlgebra, DelimitedFiles, JuMP, NEOSServer

struct Poly
    lin_ineq::Matrix
    vertices::Vector
end

# struct Poly 
#     lin_ineq::Matrix
#     vertices = calver(lin_ineq, ones(size(lin_ineq, 1)))
# end

function calver(G, w)
    vertex = Set{Tuple{Float64, Float64}}() 
    
    n_l = size(G, 1)

    # Loop through every unique pair of lines
    for i in 1 : n_l - 1
        for j in i + 1 : n_l # Start from i+1 to avoid self-intersection
            row1 = G[i, :]
            row2 = G[j, :]
            A = vcat(row1', row2')
            B = [w[i], w[j]]

            try
                solution = A \ B

                if satisfact_all(G, solution)
                    # Round the solution to a few decimal places to handle minor variations
                    # of the same point, then push to the set.
                    rounded_solution = Tuple(round.(solution, digits=8))
                    push!(vertex, rounded_solution)
                end

            catch e
                # This is expected for parallel lines, so no need to print unless debugging
                # println("Rows $i and $j: singular matrix, skipping")
                print(e)

            end
        end
    end

    # Convert the Set of vertices back to an Array
    P = Poly(G, reorderpoints(vertex))
    return P
end


# otimizar essa funcao (pense em pdroduto matricial)

function satisfact_all(G, x)

    tolerance = 1e-4 

    for v in eachrow(G)
        if v'*x > 1 + tolerance 
            return false
        end
    end
    return true
end

# function satisfact_all(G, x)
#     tolerance = 1e-4 

#     return G*x > ones(size(G)[1]) + tolerance 
# end



function reorderpoints(points)
    vet_in_order = Vector{Tuple{Float64, Float64}}() 
    # push!(vet_in_order, pop!(points))

    # print(vet_in_order)

    # finding the center of the Polyhedron

    pts = collect(points)

    cx = sum(p[1] for p in pts)/length(pts)
    cy = sum(p[2] for p in pts)/length(pts)

    vet_in_order = sort(pts, by = p -> atan(p[2] - cy, p[1] - cx)) # função anonima
    
    # sum 

    # soma_x = 0
    # soma_y = 0
    # for p in points
    #   soma_x += p[1]
    #   soma_y += p[2]
    # end
  
    # soma_x /= length(points)
    # soma_y /= length(points)

    # c = (soma_x , soma_y)


    # for i in range(1, length(points) - 1)

    #     menor = 2*pi # upper bound
    #     next_pt = Tuple{Float64, Float64}
    #     for p in points
    #       a = [p[1] - c[1]; p[2] - c[2]]
    #       b = vet_in_order[i] #[vet_in_order[i][1], vet_in_order[i][2]]
    #       d = acos(dot(a, b)/(norm(a)*norm(b)))
          
    #       if d < menor
    #           menor = d
    #           next_pt = p
    #       end
    #     end

    #     pop!(points, next_pt)

    #     push!(vet_in_order, next_pt)

    # end

    # push!(vet_in_order, pop!(points))



    return vet_in_order
end 

#conjunto de vetores igualmente espaçados 
function vet_eq_spc(n)
    angle = 2*pi/n

    #matriz de rotação
    m_rot = [cos(angle) -sin(angle);
             sin(angle) cos(angle)]
    
    vet_array = Array{Vector{Float64}}(undef, n)
    
    vet_array[1] = [0, 1]

    for i in range(2, n)
        vet_array[i] = m_rot*vet_array[i-1] 
    end

    return vet_array

end

function trajectory(x0, G, passos)
    for i in range(2, passos)
        try
            x0 = hcat(x0, G*x0[:, i-1])    
        catch e
            print(e)
        end
        
    end

    return x0   
end

# The Linear Constrained Regulation Problem
# verifica se um certo poliedro X é invariante
function is_pinvariant(A, B, C, U, X; SOF = false, d = 0) 
    
    model = Model() do
        return NEOSServer.Optimizer(; email = "wallace.lopes.162@ufrn.edu.br", solver = "Knitro")
    end

    #Parâmetros
    n = size(A, 1); #Ordem do sistema (Linhas de A)
    m = size(B, 2); #Número de Entradas (Colunas de B)
    p = size(C, 1); #Número de Saídas (Linhas de C)

    lx = size(X, 1); #Número de linhas de X
    lu = size(U, 1); #Número de linhas de U

    xi = ones(lx)
    phi = ones(lu);

    #Variáveis de decisão
    @variable(model, 0 <= H[1:lx,1:lx] <= 100)
    SOF ? @variable(model, K[1:m,1:p]) : @variable(model, F[1:m,1:n])
    @variable(model, 0 <= M[1:lu,1:lx] <= 100)
    @variable(model, 0 <= lambda)

    #Objetivo

    func_obj = lambda

    @objective(model, Min, func_obj)

    # tentando encontrar um L, caso bilinear


    if SOF # realimentação de saida 
        @constraint(model, contractiave_pi0, H*X == X*(A + B*K*C) ) 
        @constraint(model, u_included0, M*X == U*K*C )


    else # realimentação de estado
        @constraint(model, contractiave_pi0, H*X == X*(A + B*F) ) 
        @constraint(model, u_included0, M*X == U*F )

    end
    
    @constraint(model, contractiave_pi1, H*xi <= xi*lambda )
    @constraint(model, u_included1, M*xi <= phi)

    optimize!(model)

    lambda = value.(lambda)

    print(termination_status(model))

    return lambda

end

function finding_L_pinvariant(A, B, C, U, X; SOF = false, ll = 6, t = 8, pond = 0.02, d = 0)
    model = Model() do
        return NEOSServer.Optimizer(; email = "wallace.lopes.162@ufrn.edu.br", solver = "Knitro")
    end

    #Parâmetros
    #Parâmetros
    n = size(A, 1); #Ordem do sistema (Linhas de A)
    m = size(B, 2); #Número de Entradas (Colunas de B)
    p = size(C, 1); #Número de Saídas (Linhas de C)

    lx = size(X, 1); #Número de linhas de X
    lu = size(U, 1); #Número de linhas de U

    xi = ones(lx)
    phi = ones(lu);
    xl = ones(ll)

    vet = vet_array = Polyhedron.vet_eq_spc(t)
    xt = ones(t)

    #Variáveis de decisão
    @variable(model, 0 <= H[1:ll,1:ll] <= 100)
    @variable(model, 0 <= T[1:lx,1:ll] <= 100)
    @variable(model, L[1:ll,1:n])
    SOF ? @variable(model, K[1:m,1:p]) : @variable(model, F[1:m,1:n])
    @variable(model, 0 <= M[1:lu,1:ll] <= 100)
    @variable(model, 0 <= lambda <= 0.9999)
    @variable(model, J[1:n,1:ll])
    @variable(model, gamma[1:t] >= 0)


    #Objetivo
    func_obj = pond*lambda - (1-pond)*sum(gamma)/t

    @objective(model, Min, func_obj)

    # tentando encontrar um L, caso bilinear
    if SOF # realimentação de saida 
        @constraint(model, contractiave_pi0, H*L == L*(A + B*K*C) ) 
        @constraint(model, u_included0, M*L == U*K*C )


    else # realimentação de estado
        @constraint(model, contractiave_pi0, H*L == L*(A + B*F) ) 
        @constraint(model, u_included0, M*L == U*F )

    end

    @constraint(model, contractiave_pi1, H*xl <= xl*lambda )

    @constraint(model, x_included0, T*L == X)
    @constraint(model, x_included1, T*xl <= xi)

    @constraint(model, u_included1, M*xl <= phi)

    @constraint(model, rank_l, J*L == I(n))

    for i in range(1, t)
        @constraint(model, L*gamma[i]*vet[i] .<= xl)
    end

    optimize!(model)


    #Variáveis encontradas pelo modelo
    #V = value.(V)
    L = value.(L)
    SOF ? K = value.(K) :  F = value.(F)
    print("valor de lambda: ")
    print(value.(lambda))

    print(termination_status(model))

    return calver(L, xl)
end

function is_pinvariant_delay_simetric(A, Ad, F; d=0)
    model = Model() do
        return NEOSServer.Optimizer(; email = "wallace.lopes.162@ufrn.edu.br", solver = "Knitro")
    end

    # parâmetros

    n = size(A, 1) # ordem do sistema
    f = size(F, 1) # numero de linhas de f

    w = ones(f)

    @variable(model, K[1:n, 1:n])
    @variable(model, 0 <= H1[1:f, 1:f])
    @variable(model, 0 <= H2[1:f, 1:f])
    @variable(model, 0 <= L1[1:f, 1:f])
    @variable(model, 0 <= L2[1:f, 1:f])
    @variable(model, 0 <= M1[1:f, 1:f])
    @variable(model, 0 <= M2[1:f, 1:f])
    @variable(model, 0 <= N1[1:f, 1:f])
    @variable(model, 0 <= N2[1:f, 1:f])

    @constraint(model, (H1-H2)*F == F*(A + K))
    @constraint(model, (L1-L2)*F == F*(Ad - K))
    @constraint(model, (M1-M2)*F == -F*K*(A - I(n)))
    @constraint(model, (N1-N2)*F == -F*K*Ad)
    @constraint(model, ((H1 + H2) + (L1 + L2) + d*((M1 + M2) + (N1 + N2)))*w <= w)

    optimize!(model)

    K = value.(K)
    H = value.(H1) + value.(H2)
    L = value.(L1) + value.(L2)

    result = Dict("K" => K, "H" => H, "L" => L)

    print(termination_status(model))

    return result
end

function extended_space_matrix_A(A, Ad, d)
    n = size(A, 1)

    line1 = Matrix{Matrix{Float64}}(undef, 1, d)
    line1[1] = A

    for i in range(2, d)
        line1[i] = zeros(n,n) 
    end

    identity_dxd = Matrix{Matrix{Float64}}(undef, d, d)

    for i in range(1, d)
        for j in range(1, d)
            if i == j
                identity_dxd[i,j] = I(n)

            else
                identity_dxd[i,j] = zeros(n, n)
            end
        end    
    end

    es_matrix_A = vcat(line1, identity_dxd)

    columnd = Matrix{Matrix{Float64}}(undef, d+1, 1)
    columnd[1] = Ad
    for i in range(2, d+1)
        columnd[i] = zeros(n,n) 
    end

    es_matrix_A = hcat(es_matrix_A, columnd)

    return es_matrix_A

end

function extended_space_matrix_F(F, d)
    es_matrix_F = Matrix{Matrix{Float64}}(undef, d+1, d+1)

    n = size(F, 1)

    for i in range(1, d+1)
        for j in range(1, d+1)
            if i == j
                es_matrix_F[i, j] = F
            else
                es_matrix_F[i, j] = zeros(n, n)
            end
        end
    end
    return es_matrix_F
end

function cond_iniciais_adm(ext_F, ext_A, d)
    cond_iniciais_matrix = Matrix{Matrix{Matrix{Float64}}}(undef, d+1, 1)
    cond_iniciais_matrix[1] = ext_F
    for i in range(2, d+1)
        cond_iniciais_matrix[i] = ext_F*ext_A^(i-1)
    end
    return cond_iniciais_matrix
end

end

