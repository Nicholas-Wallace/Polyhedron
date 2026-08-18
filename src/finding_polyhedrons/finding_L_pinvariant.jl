"""
    finding_L_pinvariant(A, B, C, U, X; SOF = false, ll = 6, t = 8, pond = 0.02) -> Dict{String, Matrix}

Procuramos um poliedro de f linhas que seja invariante w.r.t

    x[k+1] = A*x(k) + B*u[k]

X é o poliedro de restrições 

Essa função suporta realimentação de estado e de saída

Retorna um dicionario com:
o poliedro encontrado, a matriz de ganho
lambda e todas as matrizes utilizadas no problema de otimização 
"""

function finding_L_pinvariant(A::Matrix, B::Matrix, C::Matrix, U::Matrix, X::Matrix; SOF = false, ll = 6, t = 8, pond = 0.02)
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

        @constraint(model, contractiave_pi1, H*xl <= xl*lambda )
        @constraint(model, x_included0, T*L == X)
        @constraint(model, x_included1, T*xl <= xi)
        @constraint(model, u_included1, M*xl <= phi)
        @constraint(model, rank_l, J*L == I(n))
        
        for i in range(1, t)
            @constraint(model, L*gamma[i]*vet[i] .<= xl)
        end

        optimize!(model)

        if termination_status(model) != MathOptInterface.LOCALLY_SOLVED
            error("Não foi possivel encontrar um poliedro com essas restrições")
        end

        L = value.(L)
        lambda = value(lambda)
        K = value.(K)    
        L = value.(L)
        result = Dict("lambda" => lambda, "L" => L, "K" => K)

        return result

    else # realimentação de estado
        @constraint(model, contractiave_pi0, H*L == L*(A + B*F) ) 
        @constraint(model, u_included0, M*L == U*F )

        @constraint(model, contractiave_pi1, H*xl <= xl*lambda )
        @constraint(model, x_included0, T*L == X)
        @constraint(model, x_included1, T*xl <= xi)
        @constraint(model, u_included1, M*xl <= phi)
        @constraint(model, rank_l, J*L == I(n))
        
        for i in range(1, t)
            @constraint(model, L*gamma[i]*vet[i] .<= xl)
        end

        optimize!(model)

        if termination_status(model) != MathOptInterface.LOCALLY_SOLVED
            error("Não foi possivel encontrar um poliedro com essas restrições")
        end
        
        L = value.(L)
        lambda = value(lambda)
        F = value.(F)    
        L = value.(L)
        result = Dict("lambda" => lambda, "L" => L, "F" => F)

        return result

    end
end

""" 
    finding_L_pinvariant_delay(A, B, X, f; d = 0) -> Dict{String, Matrix{Float}}

Procuramos um poliedro de f linhas que seja invariante w.r.t

    x[k+1] = A*x[k] + BG*x[k - d]

G é o ganho do controlador
X é o poliedro de restrições 

Retorna um dicionario com:
o poliedro encontrado, a matriz de ganho
lambda e todas as matrizes utilizadas no problema de otimização 
"""

function finding_L_pinvariant(A::Matrix, B::Matrix, X::Matrix; f::Int=8, d::Int = 0)
    model = Model() do
        return NEOSServer.Optimizer(; email = "wallace.lopes.162@ufrn.edu.br", solver = "Knitro")
    end

    # parâmetros

    n = size(A, 1) # ordem do sistema
    x = size(X, 1)
    b = size(B, 1)
    w = ones(f)
    q = ones(x)

    @variable(model, 0 <= lambda)

    @variable(model, F[1:f, 1:n])
    @variable(model, K[1:n, 1:n])
    @variable(model, G[1:1, 1:n])
    
    @variable(model, 0 <= H[1:f, 1:f] <= 200)
    @variable(model, 0 <= L[1:f, 1:f] <= 200)
    @variable(model, 0 <= M[1:f, 1:f] <= 200)
    @variable(model, 0 <= N[1:f, 1:f] <= 200)
    @variable(model, 0 <= R[1:x, 1:f] <= 200)
    @variable(model, BG[1:b, 1:n])

    @constraint(model, BG == B*G)
    @constraint(model, H*F == F*(A+K))
    @constraint(model, L*F == F*(BG - K))
    @constraint(model, M*F == -F*K*(A - I(n)))
    @constraint(model, N*F == -F*K*(BG))
    @constraint(model, (H + L + d*(M + N))*w .<= lambda * w)
    @constraint(model, R*F == X)
    @constraint(model, R*w .<= q)


    @objective(model, Min, lambda)
    optimize!(model)

    lambda = value(lambda)
    G = value.(G)
    F = value.(F)
    K = value.(K)
    H = value.(H)
    L = value.(L)
    M = value.(M)
    N = value.(N)
    R = value.(R)

    println(termination_status(model))

    result = Dict("lambda" => lambda, "F" => F, "G" => G, "K" => K, "H" => H, "L" => L, "M" => M, "N" => N, "R" => R)

    return result
end

