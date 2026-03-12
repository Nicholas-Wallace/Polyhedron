function is_pinvariant_seg_ref(A, B, E, S, R; lambda=0.99, ll=6) 
    
    model = Model() do
        return NEOSServer.Optimizer(; email = "wallace.lopes.162@ufrn.edu.br", solver = "Knitro")
    end

    #Parâmetros
    n = size(A, 1); #Ordem do sistema (Linhas de A)
    m = size(B, 2); #Número de Entradas (Colunas de B)
    r = size(R, 1); #Número de linhas das restrições em r[k]
    s = size(S, 1); #Número de linhas da matriz das restrições de X[k]

    l = ones(ll)
    phi = ones(r)
    w = ones(s)

    upperb = 300
    LKb = 500

    @variable(model, -LKb <= L[1:ll, 1:n] <= LKb) #Poliedro que queremos achar
    @variable(model, 0 <= H[1:ll,1:ll] <= upper) #Farkas para A + BK
    @variable(model, 0 <= P[1:ll,1:r] <= upper) #Farkas para E
    @variable(model, 0 <= M[1:s, 1:ll] <= upper) #Pra fazer a inclusão de L em S
    @variable(model, J[1:n, 1:ll]) #Pseudo inversa
    @variable(model, -LKb <= K[1:1, 1:n] <= LKb) #Ganhos do controlador

    @variable(model, igamma >= 1)
    @variable(model, 0 <= Z[1:ll, 1:s] <= upper)

    @objective(model, Min, igamma)

    @constraint(model, H*L == L*(A + B*K)) 
    @constraint(model, P*R == L*E)
    @constraint(model, M*L == S)
    @constraint(model, J*L == I(n))

    # Restrições do igamma
    @constraint(model, Z*S == L)
    @constraint(model, Z*w .<= l*igamma)
    
    @constraint(model, H*l + P*phi .<= l*lambda)
    @constraint(model, M*l .<= w)

    #0  Auto (Knitro escolhe sozinho).
    #1: Active Set 
    #2: Interior Point / CG
    #3: Interior Point / Direct
    #4: SQP (Programação Quadrática Sequencial)

    set_optimizer_attribute(model, "outlev", 2)
    set_optimizer_attribute(model, "algorithm", 1)
    set_optimizer_attribute(model, "maxtime", 3600.0)
    set_optimizer_attribute(model, "maxit", 10000)

    optimize!(model)

    L = value.(L)
    K = value.(K)

    result = Dict("L" => L, "K" => K)

    print(termination_status(model))

    return result
end
