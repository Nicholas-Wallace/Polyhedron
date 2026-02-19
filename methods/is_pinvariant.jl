
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