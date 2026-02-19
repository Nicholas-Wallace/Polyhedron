function finding_L_pinvariant_delay(A, B, X, f; d = 0)
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