"""
    finding_L_pinvariant(system::SLIT, U, X; ll = 6, t = 8, pond = 0.02, verbose = false) -> InvariantPolyhedron

`U` pode ser uma `Matrix` ou `nothing` (nesse último caso, também se pode chamar
`finding_L_pinvariant(system, X; ...)`).

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

function finding_L_pinvariant(system::SLIT, U::Union{Matrix,Nothing}, X::Matrix; ll=6, t=8, pond=0.02, verbose::Bool=false)
  model = Model() do
    return NEOSServer.Optimizer(; email="wallace.lopes.162@ufrn.edu.br", solver="Knitro")
  end

  # TODO: tratar U === nothing (pular M, phi e as restrições u_included0/u_included1)
  if U === nothing
    error("finding_L_pinvariant: caso sem U ainda não implementado (TODO)")
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
    finding_L_pinvariant(system::SLIT, X::Matrix; ll = 6, t = 8, pond = 0.02, verbose = false) -> InvariantPolyhedron

Versão sem restrição de entrada: equivale a chamar com `U = nothing`.
"""
function finding_L_pinvariant(system::SLIT, X::Matrix; ll=6, t=8, pond=0.02, verbose::Bool=false)
  return finding_L_pinvariant(system, nothing, X; ll=ll, t=t, pond=pond, verbose=verbose)
end

""" 
  finding_L_pinvariant(system::SystemDelayed, U::Matrix, X::Matrix; lambda=0.99, time=10, lf=10)  -> Dict{String, Matrix{Float}}

`U` pode ser uma `Matrix` ou `nothing` (nesse último caso, também se pode chamar
`finding_L_pinvariant(system, X; ...)`).

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
function finding_L_pinvariant(system::SystemDelayed, U::Matrix, X::Matrix; time=10, lf=14, pond=0.2, t=8) 

    model = Model() do
        return NEOSServer.Optimizer(; email = "wallace.lopes.162@ufrn.edu.br", solver = "Knitro")
    end

    # TODO: tratar U === nothing (usar U nas restrições de entrada u_included0/u_included1)
    if U === nothing
        error("finding_L_pinvariant: caso sem U ainda não implementado (TODO)")
    end

    #Parâmetros
    n = size(system.A, 1); #Ordem do sistema (Linhas de A)
    m = size(system.Ad, 2); #Número de Entradas (Colunas de Ad)
    u = size(U, 1) # Número de linhas da Restrição de Entrada
    x = size(X, 1) # Número de linhas da Restrição de Estado 

   vet = vet_array = Polyhedron.vet_eq_spc(t)
   xt = ones(t)


    ones_f = ones(lf)
    ones_x = ones(x)
    ones_u = ones(u)


    @variable(model, K[1:n, 1:n]) #Variavel auxiliar trans model
    @variable(model, G[1:m, 1:n]) #Ganho do controlador 
    @variable(model, -500 <= F[1:lf,1:n] <= 500) #Poliedro que queremos achar
    @variable(model, 0 <= H[1:lf,1:lf] <= 300)
    @variable(model, 0 <= L[1:lf,1:lf] <= 300) #Farkas para A
    @variable(model, 0 <= M[1:lf,1:lf] <= 300) #Farkas para A
    @variable(model, 0 <= N[1:lf,1:lf] <= 300) #Farkas para A
#    @variable(model, 0 <= P[1:lf,1:r] <= 300) #Farkas para E
    @variable(model, 0 <= R[1:u,1:lf] <=300) # Farkas V
    @variable(model, 0 <= T[1:x, 1:lf] <= 300) #Pra fazer a inclusão de L em X
    @variable(model, J[1:n, 1:lf]) #Pseudo inversa

    @variable(model, gamma[1:t] >= 0)
    @variable(model, 0.001 <= lambda <= 1)
    @variable(model, 0 <= Z[1:lf, 1:x] <= 300)

    func_obj = pond*lambda - (1-pond)*sum(gamma)/t

    @objective(model, Min, func_obj)

    @constraint(model, H*F == F*(system.A+K)) 
    @constraint(model, L*F == F*(system.Ad*G - K))
    @constraint(model, M*F == -F*K*(system.A-I(n)))
    @constraint(model, N*F == -F*K*system.Ad*G)
    @constraint(model, T*F == X)
    @constraint(model, R*F == U*G)
    @constraint(model, J*F == I(n))

    # Restrições do igamma
    for i in range(1, t)
      @constraint(model, F*gamma[i]*vet[i] .<= ones_f)
    end
#    @constraint(model, Z*X == F)
#    @constraint(model, Z*ones_x .<= ones_f*igamma)
    
    @constraint(model, (H + L + system.d*(M + N))*ones_f  .<= ones_f*lambda)
    @constraint(model, T*ones_f .<= ones_x)
    @constraint(model, R*ones_f .<= ones_u)

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

    system.G = value.(G) 

    ip = InvariantPolyhedron(value.(F), ones_f, value.(lambda), system)
    
    return ip 
end

"""
    finding_L_pinvariant(system::SystemDelayed, X::Matrix; time=10, lf=10) -> Dict{String, Matrix{Float64}}

Versão sem restrição de entrada: equivale a chamar com `U = nothing`.
"""
function finding_L_pinvariant(system::SystemDelayed, X::Matrix; time=10, lf=14)
    return finding_L_pinvariant(system, nothing, X; time=time, lf=lf)
end

"""
    finding_L_pinvariant(system::LTIDisturbanceDelayedInput, S, R, V; lambda=0.99, time=10, lf=10) -> InvariantPolyhedron

Procura um poliedro `F x <= ones(lf)` invariante w.r.t.

    x[k+1] = system.A*x[k] + system.B*u[k-d] + system.E*w[k]
    y[k] = system.C*x[k]

com realimentação de estado (`SSF`), isto é, `u[k] = G*x[k]`.

`S`, `R` e `V` são as matrizes das restrições de estado, distúrbio e entrada:

    S*x <= ones(size(S, 1))
    R*w <= ones(size(R, 1))
    V*u <= ones(size(V, 1))

A formulação SOF ainda não está implementada. Retorna um `InvariantPolyhedron`
e escreve o ganho encontrado em `system.G`.
"""
function finding_L_pinvariant(system::LTIDisturbanceDelayedInput, S::AbstractMatrix, R::AbstractMatrix, V::AbstractMatrix; lambda=0.99, time=10, lf=10)
    if system.feedback == SOF
        error("finding_L_pinvariant(::LTIDisturbanceDelayedInput): formulação SOF ainda não implementada (TODO)")
    elseif system.feedback == OL
        error("finding_L_pinvariant(::LTIDisturbanceDelayedInput): realimentação OL não possui ganho a otimizar")
    end

    n = size(system.A, 1)
    m = size(system.B, 2)
    q = size(system.E, 2)

    size(S, 2) == n ||
        throw(DimensionMismatch("S deve ter $n colunas (s×n), recebido $(size(S))"))
    size(R, 2) == q ||
        throw(DimensionMismatch("R deve ter $q colunas (r×q), recebido $(size(R))"))
    size(V, 2) == m ||
        throw(DimensionMismatch("V deve ter $m colunas (v×m), recebido $(size(V))"))

    result = _solve_pinvariant_segref_delay_sim2(
        system.A, system.B, system.E, S, R, V, system.d;
        lambda=lambda, time=time, lf=lf,
    )

    status = result["status"]
    status == MathOptInterface.LOCALLY_SOLVED ||
        error("Não foi possível encontrar um poliedro com essas restrições (status: $status)")

    F = result["F"]
    G = result["G"]
    system.G = G

    return InvariantPolyhedron(F, ones(lf), lambda, system)
end
