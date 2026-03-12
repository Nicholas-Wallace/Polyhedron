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

    @variable(model, 0 <= L[1:ll, 1:n] <= 100) #Poliedro que queremos achar
    @variable(model, 0 <= H[1:ll,1:ll] <= 100) #Farkas para A + BK
    @variable(model, 0 <= P[1:ll,1:r] <= 100) #Farkas para E
    @variable(model, 0 <= M[1:s, 1:ll] <= 100) #Pra fazer a inclusão de L em S
    @variable(model, 0 <= J[1:n, 1:ll]) #Pseudo inversa
    @variable(model, K[1:1, 1:n]) #Ganhos do controlador

    @variable(model, igamma >= 1)
    @variable(model, 0 <= r[1:ll, 1:s] <= 100)

    @objective(model, Min, igamma)

    @constraint(model, H*L == L*(A + B*K)) 
    @constraint(model, P*R == L*E)
    @constraint(model, M*L == S)
    @constraint(model, J*L == I(n))

    # Restrições do igamma
    @constraint(model, r*S == L)
    @constraint(model, r*w .<= l*igamma)
    
    @constraint(model, H*l + P*phi .<= l*lambda)
    @constraint(model, M*l .<= w)

    optimize!(model)

    L = value.(L)
    K = value.(K)

    result = Dict("L" => L, "K" => K)

    print(termination_status(model))

    return result
end