function step1_ref_tracking(A, B, E, Rref, X, f; lambda=0.99, symetric=false)
        model = Model() do
        return NEOSServer.Optimizer(; email = "wallace.lopes.162@ufrn.edu.br", solver = "Knitro")
    end

    # parâmetros

    n = size(A, 1) # ordem do sistema
    x = size(X, 1)
    r_ref = size(Rref, 1)
    w = ones(f)
    q = ones(x)

    @variable(model, 0 <= d <= 100)

    @variable(model, G[1:1, 1:n])
    @variable(model, K[1:n, 1:n])
    @variable(model, F[1:f, 1:n])
    @variable(model, Y[1:f, 1:n])
    @variable(model, J[1:n, 1:f])

    @constraint(model, Y == F*K)
    @constraint(model, J*F == I(n))

    @objective(model, Max, d)

    if symetric
        @variable(model, 0 <= H1[1:f, 1:f] <= 100)
        @variable(model, 0 <= H2[1:f, 1:f] <= 100)
        @variable(model, 0 <= L1[1:f, 1:f] <= 100)
        @variable(model, 0 <= L2[1:f, 1:f] <= 100)
        @variable(model, 0 <= M1[1:f, 1:f] <= 100)
        @variable(model, 0 <= M2[1:f, 1:f] <= 100)
        @variable(model, 0 <= N1[1:f, 1:f] <= 100)
        @variable(model, 0 <= N2[1:f, 1:f] <= 100)
        @variable(model, 0 <= R1[1:x, 1:f] <= 100)
        @variable(model, 0 <= R2[1:x, 1:f] <= 100)

        @constraint(model, (H1-H2)*F == F*(A + K))
        @constraint(model, (L1-L2)*F == F*(B*G - K))
        @constraint(model, (M1-M2)*F == -Y*(A - I(n)))
        @constraint(model, (N1-N2)*F == -Y*(B*G))
        @constraint(model, ((H1 + H2) + (L1 + L2) + d*((M1 + M2) + (N1 + N2)))*w .<= lambda * w)
        @constraint(model, (R1-R2)*F == X)
        @constraint(model, (R1+R2)*w .<= q)

        optimize!(model)
        
        d = value(d)
        F = value.(F)
        K = value.(K)
        H = value.(H1) - value.(H2)
        L = value.(L1) - value.(L2)
        M = value.(M1) - value.(M2)
        N = value.(N1) - value.(N2)
        R = value.(R1) - value.(R2)
        G = value.(G)

        print(termination_status(model))

        # retorna um dicionário com d e as matrizes
        result = Dict("d" => d, "F" => F, "G" => G, "K" => K, "H" => H, "L" => L, "M" => M, "N" => N, "R" => R)

        return result
    end
    
    @variable(model, 0 <= H[1:f, 1:f] <= 200)
    @variable(model, 0 <= L[1:f, 1:f] <= 200)
    @variable(model, 0 <= M[1:f, 1:f] <= 200)
    @variable(model, 0 <= N[1:f, 1:f] <= 200)
    @variable(model, 0 <= R[1:x, 1:f] <= 200)
    @variable(model, 0 <= Jref[1:f, 1:r_ref] <= 200)

    @constraint(model, H*F == F*(A + K))
    @constraint(model, L*F == F*(B*G - K))
    @constraint(model, M*F == -Y*(A - I(n)))
    @constraint(model, N*F == -Y*(B*G))
    @constraint(model, (H + L + d*(M + N))*w .<= lambda * w)
    @constraint(model, R*F == X)
    @constraint(model, R*w .<= q)

    optimize!(model)

    d = value(d)
    G = value.(G)
    F = value.(F)
    K = value.(K)
    H = value.(H)
    L = value.(L)
    M = value.(M)
    N = value.(N)
    R = value.(R)

    print(termination_status(model))

    result = Dict("d" => d, "F" => F, "G" => G, "K" => K, "H" => H, "L" => L, "M" => M, "N" => N, "R" => R)

    return result
end