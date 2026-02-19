
function finding_L_pinvariant(A, B, C, U, X; SOF = false, ll = 6, t = 8, pond = 0.02, d = 0)
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


    else # realimentação de estado
        @constraint(model, contractiave_pi0, H*L == L*(A + B*F) ) 
        @constraint(model, u_included0, M*L == U*F )

    end

    @constraint(model, contractiave_pi1, H*xl <= xl*lambda )

    @constraint(model, x_included0, T*L == X)
    @constraint(model, x_included1, T*xl <= xi)

    @constraint(model, u_included1, M*xl <= phi)

    @constraint(model, rank_l, J*L == I(n))

    for i in range(1, t)
        @constraint(model, L*gamma[i]*vet[i] .<= xl)
    end

    optimize!(model)


    #Variáveis encontradas pelo modelo
    #V = value.(V)
    L = value.(L)
    SOF ? K = value.(K) :  F = value.(F)
    print("valor de lambda: ")
    print(value.(lambda))

    print(termination_status(model))

    hrep_L = hrep(L, xl)

    # retorna um poliedro do tipo polyhedron da biblioteca Polyhedra (útil para transitar entre as hrep e vrep)
    return polyhedron(hrep_L, CDDLib.Library())
end

