
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