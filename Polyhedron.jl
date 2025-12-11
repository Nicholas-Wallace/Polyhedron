module Polyhedron

using LinearAlgebra, DelimitedFiles, JuMP, NEOSServer, Polyhedra, CDDLib, Base.Iterators


# a cddlib retorna os vertices em um vetor, mas para as operações de plotar e calcular trajetórias 
# é melhor lidar com os vertices sendo um vetor de tuplas, em que cada tupla é uma cordenada no espaço de estados
struct Vertices{T <: Tuple} 
    points::Vector{T}
end

Base.getindex(v::Vertices, i::Int) = v.points[i]
Base.length(v::Vertices) = length(v.points)
Base.iterate(v::Vertices, state=1) = iterate(v.points, state)

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

# plota a tragetória para sistemas com atraso descritos por x(k+1) = A*x(k) + A*X(k-d) 
# x0 é um vetor com todos as condições iniciais de [x(0), x(-1), ..., x(-d)]
function trajectory_delay(x0, A, Ad, passos, d)
    for i in 1:passos
        try 

            push!(x0, Tuple(A*collect(x0[d+i]) + Ad*collect(x0[i])))
        catch e
            print(e)
        end

    end

    return x0[d+1:end]
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
function is_pinvariant_delay_simetric(A, Ad, F; d=0) # falta o caso não simétrico (utilizar uma flag)
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

    # retorna um dicionário com as matrizes K, H, L
    return result
end

# funções auxiliares para gerar a matriz de condições iniciais admissiveis
# futuramente mudar o aproach de matriz de matrizes
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

function mat_cond_iniciais_adm(A, Ad, F, d; symetric=true)
  ext_A = Polyhedron.extended_space_matrix_A(A, Ad, d)
  ext_F = Polyhedron.extended_space_matrix_F(F, d)
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

#TALVEZ DE PRA USAR NO LUGAR DA EXTENDEDA_MATRIX ANTIGA
function ExtendedA_Matrix(A, Ad, d; dm=d)
    n = size(A, 1)
    N = n * (dm + 1)

    line1 = zeros(n, N-n)
    line1[:, (d-1) * n + 1 : d * n] = Ad
    line1 = hcat(A, line1)

    identity = Float64.(I(n*dm))
    bottom_zero = zeros(n*dm, n)

    return vcat(line1, hcat(identity, bottom_zero))
end

function ExtendedA_Vector(A, Ad, dm)
    n = size(A, 1)
    N = n * (dm + 1)
    x = [zeros(n, N) for _ in 1:dm]
    
    for i in 1:dm
        x[i] = ExtendedA_Matrix(A, Ad, i, dm)
    end

    return x
end

function AllPossibleComb(j, dm)
    ranges = [1:dm for _ in 1:j]
    return collect(Iterators.product(ranges...))
end

# FALTA TERMINAR AINDA
function VariantDelayAdmisInit(A, Ad, F, dm)

    A_array = ExtendedA_Vector(A, Ad, dm)
    extendedMatrix = []
    push!(extendedMatrix, F)

    for i in 1:dm
        indexes = AllPossibleComb(i, dm)
        for index in indexes
            for elem in index
                product *= A_array[elem]
            end
            !push(extendedMatrix, F * product)
        end
    end

    extendedMatrix = reduce(vcat, extendedMatrix)
    return extendedMatrix
end

end