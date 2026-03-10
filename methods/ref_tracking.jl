function is_pinvariant_seg_ref(A, B, E, Fx, Fr) 
    
    model = Model() do
        return NEOSServer.Optimizer(; email = "wallace.lopes.162@ufrn.edu.br", solver = "Knitro")
    end

    #Parâmetros
    n = size(A, 1); #Ordem do sistema (Linhas de A)
    m = size(B, 2); #Número de Entradas (Colunas de B)

    fx = size(Fx, 1)
    fr = size(Fr, 1)

    @variable(model, 0 <= H[1:fx,1:fx] <= 100)
    @variable(model, 0 <= J[1:fx,1:fr] <= 100)
    @variable(model, 0 <= lambda <= 1)
    @variable(model, K[1:1, 1:4])

    #Objetivo

    func_obj = lambda

    @objective(model, Min, func_obj)

    # Não lembro se é Fx * E ou Fr * E, ta ruim de ler na foto que eu tirei da folha
    @constraint(model, H*Fx == Fx*(A + B*K)) 
    @constraint(model, J*Fr == Fx*E )
    
    @constraint(model, H*ones(fx) + J*ones(fr)  .<= ones(fx)*lambda)

    optimize!(model)

    lambda = value.(lambda)
    K = value.(K)

    result = Dict("lambda" => lambda, "K" => K)

    print(termination_status(model))

    # retorna o resultado da otimização, de forma que se lambda > 1 -> o poliedro não é PI
    return result

end