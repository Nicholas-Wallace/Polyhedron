module Poly

using LinearAlgebra, DelimitedFiles, JuMP, HiGHS, NEOSServer, Polyhedra, CDDLib, Base.Iterators


# a cddlib retorna os vertices em um vetor, mas para as operações de plotar e calcular trajetórias 
# é melhor lidar com os vertices sendo um vetor de tuplas, em que cada tupla é uma cordenada no espaço de estados
struct Vertices{T <: Tuple} 
    points::Vector{T}
end

Base.getindex(v::Vertices, i::Int) = v.points[i]
Base.length(v::Vertices) = length(v.points)
Base.iterate(v::Vertices, state=1) = iterate(v.points, state)

function get_extreme_vertices(A, b, num_pontos=5)
    n_vars = size(A, 2)
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    @variable(model, x[1:n_vars])
    @constraint(model, A * x .<= b)
    
    lista_vertices = Vector{Vector{Float64}}()
    
    for i in 1:num_pontos
        c = randn(n_vars)
        @objective(model, Max, sum(c[j] * x[j] for j in 1:n_vars))
        
        optimize!(model)
        
        if termination_status(model) == MOI.OPTIMAL
            push!(lista_vertices, value.(x))
        end
    end
    return lista_vertices
end

# Não lembro se já ta tudo certinho, tem que revisar.
function get_vertices(F, w, init_cond, init_cond_w, index)
    hrep_ic = hrep(init_cond, init_cond_w)
    P = vrep(polyhedron(hrep_ic, CDDLib.Library()))

    ext_vertices = collect(points(P))
    ext_vertices_tuple = Tuple{Float64, Float64}[]

    for i in 1:2:(length(vertices[index]))
        push!(ext_vertices_tuple, (vertices[index][i], vertices[index][i+1]))
    end

    h_rep_F = hrep(vcat(F, -F), vcat(w, w))
    iter = points(vrep(polyhedron(h_rep_F, CDDLib.Library())))

    vertices_F = []
    for vet in iter
        push!(vertices_F, Tuple(vet))
    end

    return ext_vertices_tuple, vertices_F
end

# recebe como parâmetro um iterador de vector 
# iterador de vetor de vector é o que a vrep retorna ex: points(vrep(polyhedron(h_rep_F, CDDLib.Library()))) 
function vertices_tuple(vet, n)
    chunks = Iterators.partition(vet, n)
    
    return [Tuple(chunk) for chunk in chunks if length(chunk) == n]
end

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

function reorderpoints(points)
    vet_in_order = Vector{Tuple{Float64, Float64}}() 
    # push!(vet_in_order, pop!(points))

    # print(vet_in_order)

    # finding the center of the Polyhedron

    pts = collect(points)

    cx = sum(p[1] for p in pts)/length(pts)
    cy = sum(p[2] for p in pts)/length(pts)

    vet_in_order = sort(pts, by = p -> atan(p[2] - cy, p[1] - cx)) # função anonima
    return vet_in_order
end 

# cria um conjunto de vetores igualmente espaçados
# função auxiliar para o enl(.) dos problemas de otimização bilineares
function vet_eq_spc(n)
    angle = 2*pi/n

    #matriz de rotação
    m_rot = [cos(angle) -sin(angle);
             sin(angle) cos(angle)]
    
    vet_array = Array{Vector{Float64}}(undef, n)
    
    # por padrão iniciasse com um vetor unitário na direção do eixo y, e vai rotacionando ele  
    vet_array[1] = [0, 1]

    for i in range(2, n)
        vet_array[i] = m_rot*vet_array[i-1] 
    end
    return vet_array

end

# plota a tragetória de sistemas descritos por x(k+1) = A*x(k)
function trajectory(x0, A, passos)
    for i in range(2, passos)
        try
            x0 = hcat(x0, A*x0[:, i-1])    
        catch e
            print(e)
        end
        
    end
    return x0   
end

# plota a tragetória para sistemas com atraso descritos por x(k+1) = A*x(k) + A*x(k-d) (delay fixo)
#                                                    ou por x(k+1) = A*x(k) + A*x(k-d(k)) (delay variante)
# x0 é um vetor com todos as condições iniciais de [x(0), x(-1), ..., x(-d)]
function trajectory_delay(x0, A, Ad, passos, d; varying=false, reverse=true)
    # Como o vetor vem na forma [x[k]...x[k-d]] para plotar a trajetória
    # é melhor que esteja na ordem cronológica [x[k-d]...x[k]]
    x0_traj = copy(x0)

    if reverse
        reverse!(x0_traj)
    end
    
    if varying == false
        for i in 1:passos
            x_atual = collect(x0_traj[end])
            x_atrasado = collect(x0_traj[end - d])
            x_novo = A * x_atual + Ad * x_atrasado

            push!(x0_traj, Tuple(x_novo))   
        end
    else
        for i in 1:passos
            d_k = rand(1:d)

            x_atual = collect(x0_traj[end])
            x_atrasado = collect(x0_traj[end - d_k])
            x_novo = A * x_atual + Ad * x_atrasado

            push!(x0_traj, Tuple(x_novo))
        end
    end

    return x0_traj[d+1:end]
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

    # retorna o resultado da otimização, de forma que se lambda > 1 -> o poliedro não é PI
    return lambda

end

# Caso bilinear
# procuramos o maior poliedro dentro de X, que seja PI
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

    hrep_L = hrep(L, xl)

    # retorna um poliedro do tipo polyhedron da biblioteca Polyhedra (útil para transitar entre as hrep e vrep)
    return polyhedron(hrep_L, CDDLib.Library())
end

# Verifica se o poliedro F é PI, levando em consideração um atraso de tamanho d
# utilizasse o modelo transformado de x(k+1) = A*x(k) + Ad*x(k-d)
# Caso não simétrico e simétrico implementado com flag
function is_pinvariant_delay(A, Ad, F; d=0, symetric=true)
    model = Model() do
        return NEOSServer.Optimizer(; email = "wallace.lopes.162@ufrn.edu.br", solver = "Knitro")
    end

    # parâmetros

    n = size(A, 1) # ordem do sistema
    f = size(F, 1) # numero de linhas de f
    w = ones(f)

    @variable(model, 0 <= lambda <= 1)

    if symetric
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
        @constraint(model, ((H1 + H2) + (L1 + L2) + d*((M1 + M2) + (N1 + N2)))*w .<= lambda * w)

        @objective(model, Min, lambda)

        optimize!(model)
        
        lambda = value(lambda)
        K = value.(K)
        H = value.(H1) + value.(H2)
        L = value.(L1) + value.(L2)

        print(termination_status(model))

        # retorna um dicionário com lambda e as matrizes K, H e L
        result = Dict("lambda" => lambda, "K" => K, "H" => H, "L" => L)

        return result
    end
    
    @variable(model, K[1:n, 1:n])
    @variable(model, 0 <= H[1:f, 1:f])
    @variable(model, 0 <= L[1:f, 1:f])
    @variable(model, 0 <= M[1:f, 1:f])
    @variable(model, 0 <= N[1:f, 1:f])

    @constraint(model, H*F == F*(A + K))
    @constraint(model, L*F == F*(Ad - K))
    @constraint(model, M*F == -F*K*(A - I(n)))
    @constraint(model, N*F == -F*K*Ad)
    @constraint(model, (H + L + d*(M + N))*w .<= lambda * w)

    @objective(model, Min, lambda)

    optimize!(model)

    lambda = value(lambda)
    K = value.(K)
    H = value.(H)
    L = value.(L)

    print(termination_status(model))

    result = Dict("lambda" => lambda, "K" => K, "H" => H, "L" => L)

    return result
end

# Passo 1 da saturação
function step1_saturation(A, B, F, X, lambda=0.99; symetric=true)
    model = Model() do
        return NEOSServer.Optimizer(; email = "wallace.lopes.162@ufrn.edu.br", solver = "Knitro")
    end

    # parâmetros

    n = size(A, 1) # ordem do sistema
    f = size(F, 1) # numero de linhas de f
    x = size(X, 1)
    w = ones(f)
    q = ones(x)

    @variable(model, 0 <= d)

    if symetric
        @variable(model, G[1, 1:n])
        @variable(model, K[1:n, 1:n])

        @variable(model, 0 <= H1[1:f, 1:f])
        @variable(model, 0 <= H2[1:f, 1:f])
        @variable(model, 0 <= L1[1:f, 1:f])
        @variable(model, 0 <= L2[1:f, 1:f])
        @variable(model, 0 <= M1[1:f, 1:f])
        @variable(model, 0 <= M2[1:f, 1:f])
        @variable(model, 0 <= N1[1:f, 1:f])
        @variable(model, 0 <= N2[1:f, 1:f])
        @variable(model, 0 <= R1[1:x, 1:f])
        @variable(model, 0 <= R2[1:x, 1:f])

        @constraint(model, (H1-H2)*F == F*(A + K))
        @constraint(model, (L1-L2)*F == F*(B*G - K))
        @constraint(model, (M1-M2)*F == -F*K*(A - I(n)))
        @constraint(model, (N1-N2)*F == -F*K*B*G)
        @constraint(model, ((H1 + H2) + (L1 + L2) + d*((M1 + M2) + (N1 + N2)))*w .<= lambda * w)
        @constraint(model, (R1-R2)*F == X)
        @constraint(model, (R1+R2)*w .<= q)

        @objective(model, Max, d)

        optimize!(model)
        
        d = value(d)
        K = value.(K)
        H = value.(H1) + value.(H2)
        L = value.(L1) + value.(L2)

        print(termination_status(model))

        # retorna um dicionário com lambda e as matrizes K, H e L
        result = Dict("d" => d, "K" => K, "H" => H, "L" => L)

        return result
    end
    
    @variable(model, G[1:1, 1:n])
    @variable(model, K[1:n, 1:n])
    @variable(model, 0 <= H[1:f, 1:f])
    @variable(model, 0 <= L[1:f, 1:f])
    @variable(model, 0 <= M[1:f, 1:f])
    @variable(model, 0 <= N[1:f, 1:f])
    @variable(model, 0 <= R[1:x, 1:f])

    @constraint(model, H*F == F*(A + K))
    @constraint(model, L*F == F*(B*G - K))
    @constraint(model, M*F == -F*K*(A - I(n)))
    @constraint(model, N*F == -F*K*B*G)
    @constraint(model, (H + L + d*(M + N))*w .<= lambda * w)
    @constraint(model, R*F = X)
    @constraint(model, R*w .<= q)

    @objective(model, Max, d)

    optimize!(model)

    d = value(d)
    K = value.(K)
    H = value.(H)
    L = value.(L)

    print(termination_status(model))

    result = Dict("d" => d, "K" => K, "H" => H, "L" => L)

    return result
end

# FALTA IMPLEMENTAR
function step2_saturation(A, B, F, X, umax, umin, d, ll=6, lambda=0.99, v=8; symetric=true)
    model = Model() do
        return NEOSServer.Optimizer(; email = "wallace.lopes.162@ufrn.edu.br", solver = "Knitro")
    end

    # parâmetros

    n = size(A, 1) # ordem do sistema
    f = size(F, 1) # numero de linhas de f
    x = size(X, 1)
    w = ones(f)
    q = ones(x)

    vet_array = Polyhedron.vet_eq_spc(v)

    if symetric
        @variable(model, P[1:1, 1:n])
        @variable(model, G[1:1, 1:n])
        @variable(model, K[1:n, 1:n])

        @variable(model, L[1:ll, 1:n])
        @variable(model, 0 <= Hp[1:f, 1:f])
        @variable(model, 0 <= Hm[1:f, 1:f])
        @variable(model, 0 <= L1p[1:f, 1:f])
        @variable(model, 0 <= L1m[1:f, 1:f])
        @variable(model, 0 <= L2p[1:f, 1:f])
        @variable(model, 0 <= L2m[1:f, 1:f])
        @variable(model, 0 <= Mp[1:f, 1:f])
        @variable(model, 0 <= Mm[1:f, 1:f])
        @variable(model, 0 <= N1p[1:f, 1:f])
        @variable(model, 0 <= N1m[1:f, 1:f])
        @variable(model, 0 <= N2p[1:f, 1:f])
        @variable(model, 0 <= N2m[1:f, 1:f])
        @variable(model, 0 <= Rp[1:x, 1:f])
        @variable(model, 0 <= Rm[1:x, 1:f])
        @variable(model, 0 <= Qp[1:1, 1:f])
        @variable(model, 0 <= Qm[1:1, 1:f])
    
        @variable(model, gamma[1:t] >= 0)
        @objective(model, Max, gamma)
    
        @constraint(model, (Hp-Hm)*F == F*(A + K))
        @constraint(model, (L1p-L1m)*F == F*(B*G - K))
        @constraint(model, (L2p-L2m)*F == F*(B*P - K))
        @constraint(model, (Mp-Mm)*F == -F*K*(A - I(n)))
        @constraint(model, (N1p-N1m)*F == -F*K*B*G)
        @constraint(model, (N2p-N2m)*F == -F*K*B*P)
        @constraint(model, ((Hp + Hm) + (L1p + L1m) + d*((Mp + Mm) + (N1p + N1m)))*w .<= lambda * w)
        @constraint(model, ((Hp + Hm) + (L1p + L1m) + d*((Mp + Mm) + (N2p + N2m)))*w .<= lambda * w)
        @constraint(model, ((Hp + Hm) + (L2p + L2m) + d*((Mp + Mm) + (N1p + N1m)))*w .<= lambda * w)
        @constraint(model, ((Hp + Hm) + (L2p + L2m) + d*((Mp + Mm) + (N2p + N2m)))*w .<= lambda * w)

        @constraint(model, (Rp-Rm)*F == X)
        @constraint(model, (Rp+Rm)*w .<= q)

        @constraint(model, (Qp - Qm)*F == P)
        @constraint(model, (Qp + Qm)*w .<= umax)

        for i in range(1, t)
            @constraint(model, L*gamma[i]*vet[i] .<= xl)
        end

        @objective(model, Max, gamma)

        optimize!(model)
        
        d = value(d)
        K = value.(K)
        H = value.(H1) + value.(H2)
        L = value.(L1) + value.(L2)

        print(termination_status(model))

        # retorna um dicionário com lambda e as matrizes K, H e L
        result = Dict("d" => d, "K" => K, "H" => H, "L" => L)

        return result
    end
    
    @variable(model, P[1:1, 1:n])
    @variable(model, G[1:1, 1:n])
    @variable(model, K[1:n, 1:n])

    @variable(model, 0 <= H[1:f, 1:f])
    @variable(model, 0 <= L1[1:f, 1:f])
    @variable(model, 0 <= L2[1:f, 1:f])
    @variable(model, 0 <= M[1:f, 1:f])
    @variable(model, 0 <= N1[1:f, 1:f])
    @variable(model, 0 <= N2[1:f, 1:f])
    @variable(model, 0 <= R[1:x, 1:f])
    @variable(model, 0 <= Q1[1:1, 1:f])
    @variable(model, 0 <= Q2[1:1, 1:f])

    @constraint(model, H*F == F*(A + K))
    @constraint(model, L1*F == F*(B*G - K))
    @constraint(model, L2*F == F*(B*P - K))
    @constraint(model, M*F == -F*K*(A - I(n)))
    @constraint(model, N1*F == -F*K*B*G)
    @constraint(model, N2*F == -F*K*B*P)

    @constraint(model, (H + L1 + d*(M + N1))*w .<= lambda * w)
    @constraint(model, (H + L1 + d*(M + N2))*w .<= lambda * w)
    @constraint(model, (H + L2 + d*(M + N1))*w .<= lambda * w)
    @constraint(model, (H + L2 + d*(M + N2))*w .<= lambda * w)

    @constraint(model, R*F = X)
    @constraint(model, R*w .<= q)
    
    @constraint(model, Q1*F == P)
    @constraint(model, Q2*F == -P)
    @constraint(model, Q1*w .<= umax)
    @constraint(model, Q2*w .<= umin)

    @objective(model, Max, enl)

    optimize!(model)

    K = value.(K)
    H = value.(H)
    L = value.(L)

    print(termination_status(model))

    result = Dict("K" => K, "H" => H, "L" => L)

    return result
end


# funções auxiliares para gerar a matriz de condições iniciais admissiveis
# futuramente mudar o aproach de matriz de matrizes

function cond_iniciais_adm(ext_F, ext_A, d)
    cond_iniciais_matrix = Matrix{Matrix{Matrix{Float64}}}(undef, d+1, 1)
    cond_iniciais_matrix[1] = ext_F
    for i in range(2, d+1)
        cond_iniciais_matrix[i] = ext_F*ext_A^(i-1)
    end

    

    return cond_iniciais_matrix
end

function mat_cond_iniciais_adm(A, Ad, F, d; symetric=true)
  ext_A = Polyhedron.extended_A(A, Ad, d)
  ext_F = Polyhedron.extended_F(F, d)
  cond_iniciais_adm = Polyhedron.cond_iniciais_adm(ext_F, ext_A, d)
  m, n = size(cond_iniciais_adm)#mxn
  p, q = size(cond_iniciais_adm[1])#pxq
  r, s = size(cond_iniciais_adm[1][1])#rxs

  cia_final = cond_iniciais_adm[1, 1][1, 1]
  temp = cond_iniciais_adm[1, 1][1, 1]

  for i in 1:m
    for j in 1:n
      for k in 1:p
        for l in 1:q
          if l == 1
            temp = cond_iniciais_adm[i, j][k, l]
          else
            temp = hcat(temp, cond_iniciais_adm[i, j][k, l])
          end
        end
        if i == j == k == 1
          cia_final = temp
        else
          cia_final = vcat(cia_final, temp)
        end
      end
    end
  end
  if symetric
    return vcat(cia_final, -cia_final)
  else
    return cia_final
  end
  
end

function elimred(G, ro)
    h = hrep(G, ro)
    p = polyhedron(h, CDDLib.Library())
    removehredundancy!(p)
    h_clean = hrep(p)

    n_h = nhalfspaces(h_clean)
    dim = size(G, 2)

    Gn = Matrix{Float64}(undef, n_h, dim)
    ron = Vector{Float64}(undef, n_h)

    for (i, half_space) in enumerate(halfspaces(h_clean))
        Gn[i, :] = half_space.a
        ron[i] = half_space.β
    end

    return Gn, ron
end

# Cria lista de matrizes F e forma a matriz diagonal extended_F
function extended_F(F, d)
    blocks = [F for _ in 1:(d+1)]
    # Concatena na diagonal (concatena tanto verticalmente quanto horizontalmente)
    # Isso gera a matriz diagonal, os zeros são colocados automaticamente
    return cat(blocks..., dims=(1, 2))
end

function extended_A(A, Ad, d; dm=d)
    n = size(A, 1)
    N = n * (dm + 1)
    id = (d-1) * n + 1 : d * n

    row1 = zeros(n, N-n)
    row1[:, id] = Ad
    row1 = hcat(A, row1)

    identity = Float64.(I(n*dm))
    bottom_zero = zeros(n*dm, n)

    return vcat(row1, hcat(identity, bottom_zero))
end

function extended_A_Vector(A, Ad, dm)
    n = size(A, 1)
    N = n * (dm + 1)
    x = [zeros(n, N) for _ in 1:dm]
    
    for i in 1:dm
        x[i] = extended_A(A, Ad, i; dm=dm)
    end

    return x
end

function allPossibleComb(j, dm)
    # Cria lista de ranges 1:dm
    ranges = [1:dm for _ in 1:j]

    # Retorna o produto cartesiano desses ranges numa matriz
    return collect(Iterators.product(ranges...))
end

function admissable_initCond(A, Ad, F, dm, w; symetric=false, fixed_d=false)
    ext_F = extended_F(F, dm)
    ext_w = repeat(w, dm + 1)

    init_cond_F = ext_F
    init_cond_w = ext_w 

    n = size(A, 1)
    N = n * (dm + 1)

    if fixed_d == false
        A_array = extended_A_Vector(A, Ad, dm)
        
        for i in 1:dm
            indexes = allPossibleComb(i, dm)
            for index in indexes
                product = Matrix{Float64}(I, N, N)
                for elem in index
                    product = A_array[elem] * product
                end
                init_cond_F = vcat(init_cond_F, ext_F * product)
                init_cond_w = vcat(init_cond_w, ext_w)
            end

            init_cond_F, init_cond_w = elimred(init_cond_F, init_cond_w)
        end
    else 
        ext_A = extended_A(A, Ad, dm)
        current_A_pow = ext_F

        for i in 1:dm
            current_A_pow *= ext_A
            init_cond_F = vcat(init_cond_F, current_A_pow)
            init_cond_w = vcat(init_cond_w, ext_w)
            
            init_cond_F, init_cond_w = elimred(init_cond_F, init_cond_w)
        end
    end
        
    if symetric
        return vcat(init_cond_F, -init_cond_F), vcat(init_cond_w, init_cond_w)
    end

    return init_cond_F, init_cond_w
end

end