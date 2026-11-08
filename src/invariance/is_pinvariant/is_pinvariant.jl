
"""
    is_pinvariant(system::SLIT, X; U=nothing, return_model=false) -> InvariantPolyhedron

Verifica se o poliedro X é p-invariante w.r.t

    x[k+1] = system.A*x[k] + system.B*u[k]
    y[k+1] = system.C*x[k]

O tipo de realimentação é lido de `system.feedback`:
  SOF -> u[k] = G*y[k] = G*system.C*x[k]
  SSF -> u[k] = G*x[k]
  OL  -> não suportado (não há ganho a otimizar)

As restrições são adicionadas pelas funções auxiliares
`invariance_constraints!` e `input_constraints!` (ambas em modo `find=false`).
Se `U` for `nothing`, as restrições de entrada são omitidas.

Após a solução do Knitro o ganho encontrado é escrito em `system.G`.
Retorna um `InvariantPolyhedron` com o valor de lambda e as representações do
poliedro.

Com `return_model=true` o modelo JuMP é retornado antes da otimização.

X é o poliedro a ser verificado
U é a restrição na entrada (opcional)
"""
function is_pinvariant(system::SLIT, X::Matrix;
                       U::Union{Matrix,Nothing}=nothing,
                       return_model::Bool=false)

  model = Model() do
    return NEOSServer.Optimizer(; email="wallace.lopes.162@ufrn.edu.br", solver="Knitro")
  end

  system.feedback == OL && error("realimentação OL não possui ganho a otimizar")

  @variable(model, 0 <= lambda <= 1)

  invariance_constraints!(model, system; find=false, X=X)

  if U !== nothing
    input_constraints!(model, system, U, size(system.A, 1); find=false, X=X)
  end

  @objective(model, Min, lambda)

  return_model && return model

  optimize!(model)
  @show termination_status(model)

  termination_status(model) == MathOptInterface.LOCALLY_SOLVED ||
    error("Poliedro não é invariante")

  system.G = value.(model[:G])

  return InvariantPolyhedron(X, ones(size(X, 1)), value(lambda), system)
end


"""
    is_pinvariant(system::SystemDelayed, X; symetric=true, return_model=false) -> InvariantPolyhedron

Verifica se o poliedro X é p-invariante w.r.t

    x[k+1] = system.A*x[k] + system.Ad*x[k - system.d]

O tipo de realimentação é lido de `system.feedback`:

  SSF -> u[k] = K*x[k]
  SOF -> u[k] = G*system.C*x[k] (formulação a ser implementada)

Após a solução do Knitro o ganho encontrado é escrito em `system.G`.
Retorna um `InvariantPolyhedron` com todas as representações relevantes e o
lambda.

Com `return_model=true` retorna um JuMP.Model resultante do problema de otimização

Onde system.d é o tamanho do atraso
"""
function is_pinvariant(system::SystemDelayed, X; symetric=true, return_model::Bool=false)
  model = Model() do
    return NEOSServer.Optimizer(; email="wallace.lopes.162@ufrn.edu.br", solver="Knitro")
  end

  # parâmetros

  n = size(system.A, 1) # ordem do sistema
  x = size(X, 1) # numero de linhas de x
  w = ones(x)

  if system.feedback == SOF # realimentação de saida
    # TODO: formular o problema de otimização para SOF
    error("Formulação SOF para SystemDelayed ainda não implementada")
  end

  # realimentação de estado (SSF)
  @variable(model, 0 <= lambda <= 1)

  if symetric
    @variable(model, K[1:n, 1:n])
    @variable(model, 0 <= H1[1:x, 1:x])
    @variable(model, 0 <= H2[1:x, 1:x])
    @variable(model, 0 <= L1[1:x, 1:x])
    @variable(model, 0 <= L2[1:x, 1:x])
    @variable(model, 0 <= M1[1:x, 1:x])
    @variable(model, 0 <= M2[1:x, 1:x])
    @variable(model, 0 <= N1[1:x, 1:x])
    @variable(model, 0 <= N2[1:x, 1:x])

    @constraint(model, (H1-H2)*X == X*(system.A + K))
    @constraint(model, (L1-L2)*X == X*(system.Ad - K))
    @constraint(model, (M1-M2)*X == -X*K*(system.A - I(n)))
    @constraint(model, (N1-N2)*X == -X*K*system.Ad)
    @constraint(model, ((H1 + H2) + (L1 + L2) + system.d*((M1 + M2) + (N1 + N2)))*w .<= lambda * w)

    @objective(model, Min, lambda)
    optimize!(model)


    @show termination_status(model)

    if termination_status(model) != MathOptInterface.LOCALLY_SOLVED
      error("polyhedro não é invariante")
    end

    lambda = value(lambda)
    K = value.(K)
    H1 = value.(H1); H2 = value.(H2)
    L1 = value.(L1); L2 = value.(L2)
    M1 = value.(M1); M2 = value.(M2)
    N1 = value.(N1); N2 = value.(N2)

    system.G = K

    ip = InvariantPolyhedron(X, w, lambda, system)

    return return_model ? (poly=ip, model=model) : ip
  end

  @variable(model, K[1:n, 1:n])
  @variable(model, 0 <= H[1:x, 1:x])
  @variable(model, 0 <= L[1:x, 1:x])
  @variable(model, 0 <= M[1:x, 1:x])
  @variable(model, 0 <= N[1:x, 1:x])

  @constraint(model, H*X == X*(system.A + K))
  @constraint(model, L*X == X*(system.Ad - K))
  @constraint(model, M*X == -X*K*(system.A - I(n)))
  @constraint(model, N*X == -X*K*system.Ad)
  @constraint(model, (H + L + system.d*(M + N))*w .<= lambda * w)

  @objective(model, Min, lambda)
  optimize!(model)


  @show termination_status(model)

  if termination_status(model) != MathOptInterface.LOCALLY_SOLVED
    error("polyhedro não é invariante")
  end

  system.G = value.(K)

  ip = InvariantPolyhedron(X, w, lambda, system)

  return model ? (poly=ip, model=model) : ip
end
