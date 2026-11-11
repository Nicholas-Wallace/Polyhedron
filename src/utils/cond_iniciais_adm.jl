# utilizar o "!" e alterar a matrix diretamente, sem copia
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

"""
    admissable_initCond(A, Ad, F, dm, w; symetric=false, fixed_d=false) -> Matrix{Float64}, Vector{Float64}
    
    Gera a matriz de condições iniciais admissíveis para o sistema:

        x[k+1] = A * x[k] + Ad * x[k-d(k)]

    Onde d(k) é o delay, que pode ser fixo ou variável.

    E para o poliedro: 

        F * x <= w

    `dm` é o maior delay possível tal que o poliedro ainda seja invariante.
    A flag `symetric` controla se a matriz de condições iniciais deve ser simétrica, concatenando verticalmente ela
    com sua versão negativa.
    A flag `fixed_d` é para indicar se o sistema é com delay fixo ou variável

    Retorna a matriz de condições iniciais admissíveis e o vetor de limites correspondentes.
"""

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

function calcular_v(F, E, R)
    # Pré-computa a matriz M = F * E
    M = F * E
    
    n_linhas = size(M, 1)
    n_cols_R = size(R, 2)  # dimensão do vetor r
    
    ones_r = ones(size(R, 1))
    
    v = zeros(n_linhas)
    
    # Cria UM único modelo e resolve repetidamente (mais eficiente que criar N modelos)
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    
    @variable(model, r[1:n_cols_R])
    @constraint(model, restricoes, R * r .<= ones_r)
    
    for i in 1:n_linhas
        # Define o objetivo como a i-ésima linha de M
        @objective(model, Max, sum(M[i, j] * r[j] for j in 1:n_cols_R))
        
        optimize!(model)
        
        if termination_status(model) == OPTIMAL
            v[i] = objective_value(model)
        else
            error("LP infeasible or unbounded at row $i")
        end
    end
    
    return v
end
