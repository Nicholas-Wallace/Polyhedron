function enl_simples(S, X)

    # S contido em X

    s = size(S, 1); #Número de linhas da matriz das restrições de X[k]
    x = size(X, 1)

    model = Model() do
        return NEOSServer.Optimizer(; email = "wallace.lopes.162@ufrn.edu.br", solver = "Knitro")
    end

    @variable(model, 0.001 <= igamma <= 1000)
    @variable(model, 0 <= Z[1:x, 1:s] <= 30000)

    @constraint(model, Z*S == X)
    @constraint(model, Z*ones(s) .<= ones(x)*igamma)

    @objective(model, Min, igamma)

    optimize!(model)


    A = [1 1; 0 1]
    B = [2; 1;;]
    C = [1 0]
    system = SLIT(A, B, C)
    ip_int = InvariantPolyhedron(X*(value.(igamma)^(-1)), ones(x), 1, system)
    ip_ext = InvariantPolyhedron(S, ones(s), 1, system)
    
    @show value.(igamma)

    return ip_int, ip_ext
end
