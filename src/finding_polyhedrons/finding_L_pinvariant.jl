"""
    finding_L_pinvariant(system::SLIT, U, X; ll = 6, t = 8, pond = 0.02, verbose = false) -> InvariantPolyhedron

Procuramos um poliedro de f linhas que seja invariante w.r.t

    x[k+1] = system.A*x[k] + system.B*u[k]
    y[k+1] = system.C*x[k]

X é o poliedro de restrições 

Essa função suporta realimentação de estado e de saída

Retorna um `InvariantPolyhedron` com o poliedro encontrado (`L x <= ones(ll)`),
o ganho encontrado (escrito em `system.G`) e o valor de lambda.

Com `verbose=true` retorna uma `NamedTuple` `(poly=..., matrices=...)`, onde
`matrices` é um `Dict{String,Any}` com as matrizes do programa linear
(problema: `A`, `B`, `C`, `U`, `X`; variáveis: `H`, `T`, `L`, `M`, `J`,
`gamma`, o ganho `K`/`F`, `lambda`; e `status`).
"""

function finding_L_pinvariant(system::SLIT, U::Matrix, X::Matrix; ll=6, t=8, pond=0.02, verbose::Bool=false)
  model = Model() do
    return NEOSServer.Optimizer(; email="wallace.lopes.162@ufrn.edu.br", solver="Knitro")
  end

  #Parâmetros
  n, m, p = checkdims(system.A, system.B, system.C)

  lx = size(X, 1); #Número de linhas de X
  lu = size(U, 1); #Número de linhas de U

  xi = ones(lx)
  phi = ones(lu);
  xl = ones(ll)

  vet = vet_array = Polyhedron.vet_eq_spc(t)
  xt = ones(t)

  #Variáveis de decisão
  @variable(model, 0 <= H[1:ll, 1:ll] <= 100)
  @variable(model, 0 <= T[1:lx, 1:ll] <= 100)
  @variable(model, L[1:ll, 1:n])
  system.feedback == SOF ? @variable(model, K[1:m, 1:p]) : @variable(model, F[1:m, 1:n])
  @variable(model, 0 <= M[1:lu, 1:ll] <= 100)
  @variable(model, 0 <= lambda <= 0.9999)
  @variable(model, J[1:n, 1:ll])
  @variable(model, gamma[1:t] >= 0)


  #Objetivo
  func_obj = pond*lambda - (1-pond)*sum(gamma)/t

  @objective(model, Min, func_obj)

  # tentando encontrar um L, caso bilinear
  if system.feedback == SOF # realimentação de saida 
    @constraint(model, contractiave_pi0, H*L == L*(system.A + system.B*K*system.C))
    @constraint(model, u_included0, M*L == U*K*system.C)

    @constraint(model, contractiave_pi1, H*xl <= xl*lambda)
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

    lambda = value(lambda)
    K = value.(K)
    L = value.(L)
    H = value.(H)
    T = value.(T)
    M = value.(M)
    J = value.(J)
    gamma = value.(gamma)

    system.G = K

    ip = InvariantPolyhedron(L, xl, lambda, system)

    matrices = Dict{String,Any}(
      "A" => system.A, "B" => system.B, "C" => system.C,
      "U" => U, "X" => X,
      "H" => H, "T" => T, "L" => L, "M" => M, "J" => J, "gamma" => gamma,
      "K" => K, "lambda" => lambda, "status" => string(termination_status(model)),
    )

    return verbose ? (poly=ip, matrices=matrices) : ip

  else # realimentação de estado
    @constraint(model, contractiave_pi0, H*L == L*(system.A + system.B*F))
    @constraint(model, u_included0, M*L == U*F)

    @constraint(model, contractiave_pi1, H*xl <= xl*lambda)
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
    H = value.(H)
    T = value.(T)
    M = value.(M)
    J = value.(J)
    gamma = value.(gamma)

    system.G = F

    ip = InvariantPolyhedron(L, xl, lambda, system)

    matrices = Dict{String,Any}(
      "A" => system.A, "B" => system.B, "C" => system.C,
      "U" => U, "X" => X,
      "H" => H, "T" => T, "L" => L, "M" => M, "J" => J, "gamma" => gamma,
      "F" => F, "lambda" => lambda, "status" => string(termination_status(model)),
    )

    return verbose ? (poly=ip, matrices=matrices) : ip

  end
end

""" 
  finding_L_pinvariant(system::SystemDelayed, U::Matrix, X::Matrix; lambda=0.99, time=10, lf=10)  -> Dict{String, Matrix{Float}}

Procuramos um poliedro de f linhas que seja invariante w.r.t

    x[k+1] = Ax[k] + Adu[k-d]
    y[k+1] = Cx[k]

Respeitando as restrições de estado e saída respectivamente:
    
    Xx <= ones(size(X, 1))
    Uu <= ones(size(X, 1))

Retorna um dicionario com:
o poliedro encontrado, a matriz de ganho
lambda e todas as matrizes utilizadas no problema de otimização 
"""
function finding_L_pinvariant(system::SystemDelayed, U::Matrix, X::Matrix; time=10, lf=10) 
   
    model = Model() do
        return NEOSServer.Optimizer(; email = "wallace.lopes.162@ufrn.edu.br", solver = "Knitro")
    end

    #Parâmetros
    n = size(A, 1); #Ordem do sistema (Linhas de A)
    m = size(B, 2); #Número de Entradas (Colunas de B)
    r = size(R, 1); #Número de linhas das restrições em r[k]
    s = size(S, 1); #Número de linhas da matriz das restrições de X[k]

    ones_f = ones(lf)
    ones_r = ones(r)
    ones_s = ones(s)


    @variable(model, K[1:n, 1:n]) #Variavel auxiliar trans model
    @variable(model, G[1:m, 1:n]) #Ganho do controlador 
    @variable(model, -500 <= F[1:lf,1:n] <= 500) #Poliedro que queremos achar
    @variable(model, 0 <= H[1:lf,1:lf] <= 300)
    @variable(model, 0 <= L[1:lf,1:lf] <= 300) #Farkas para A
    @variable(model, 0 <= M[1:lf,1:lf] <= 300) #Farkas para A
    @variable(model, 0 <= N[1:lf,1:lf] <= 300) #Farkas para A
    @variable(model, 0 <= P[1:lf,1:r] <= 300) #Farkas para E
    @variable(model, 0 <= T[1:s, 1:lf] <= 300) #Pra fazer a inclusão de L em S
    @variable(model, J[1:n, 1:lf]) #Pseudo inversa

    @variable(model, 0.001 <= igamma <= 1000)
    @variable(model, 0 <= Z[1:lf, 1:s] <= 300)

    @objective(model, Min, igamma)

    @constraint(model, H*F == F*(A+K)) 
    @constraint(model, L*F == F*(B*G - K))
    @constraint(model, M*F == -F*K*(A-I(n)))
    @constraint(model, N*F == -F*K*B*G)
    @constraint(model, P*R == F*E)
    @constraint(model, T*F == S)
    @constraint(model, J*F == I(n))

    # Restrições do igamma
    @constraint(model, Z*S == F)
    @constraint(model, Z*ones_s .<= ones_f*igamma)
    
    @constraint(model, (H + L + d*(M + N))*ones_f + P*ones_r .<= ones_f*lambda)
    @constraint(model, T*ones_f .<= ones_s)

    #0  Auto (Knitro escolhe sozinho).
    #1: Active Set 
    #2: Interior Point / CG
    #3: Interior Point / Direct
    #4: SQP (Programação Quadrática Sequencial)

    set_optimizer_attribute(model, "outlev", 2)
    set_optimizer_attribute(model, "algorithm", 2)
    set_optimizer_attribute(model, "maxtime", time*60.0)
    set_optimizer_attribute(model, "maxit", 10000)
    set_optimizer_attribute(model, "ms_enable", 1)
    set_optimizer_attribute(model, "ms_maxsolves", 100)

    optimize!(model)

    F = value.(F)
    G = value.(G)
    K = value.(K)
    H = value.(H)
    L = value.(L)
    M = value.(M)
    N = value.(N)
    P = value.(P)
    T = value.(T)
    Z = value.(Z)
    J = value.(J)

    result = Dict(
    "F" => F,
    "G" => G,
    "K" => K,
    "H" => H,
    "L" => L, 
    "M" => M,
    "N" => N,
    "P" => P, 
    "T" => T,
    "Z" => Z
    )

    # Salvar matrizes em arquivo
    open("results_pinvariant_segref_delay.txt", "w") do io
        println(io, "=== Resultados da otimização ===")
        println(io, "Status: $(termination_status(model))\n")
        
        println(io, "=== Matriz F (Poliedro) ===")
        writedlm(io, F)
        println(io, "")
        
        println(io, "=== Matriz G (Ganho do controlador) ===")
        writedlm(io, G)
        println(io, "")
        
        println(io, "=== Matriz K (Variável auxiliar) ===")
        writedlm(io, K)
        println(io, "")
        
        println(io, "=== Matriz H (Farkas) ===")
        writedlm(io, H)
        println(io, "")
        
        println(io, "=== Matriz L (Farkas) ===")
        writedlm(io, L)
        println(io, "")
        
        println(io, "=== Matriz M (Farkas) ===")
        writedlm(io, M)
        println(io, "")
        
        println(io, "=== Matriz N (Farkas) ===")
        writedlm(io, N)
        println(io, "")
        
        println(io, "=== Matriz P (Farkas para E) ===")
        writedlm(io, P)
        println(io, "")
        
        println(io, "=== Matriz T (Inclusão em S) ===")
        writedlm(io, T)
        println(io, "")
        
        println(io, "=== Matriz Z ===")
        writedlm(io, Z)
        println(io, "")
    end

    print(termination_status(model))

    return result
end
