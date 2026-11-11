
"""
    is_pinvariant(system::SLIT, X; U=nothing, time=10, ms_seed=0, return_model=false) -> InvariantPolyhedron

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
                       time::Real=10,
                       ms_seed::Int=0,
                       return_model::Bool=false)

  model = Model() do
    return NEOSServer.Optimizer(; email="wallace.lopes.162@ufrn.edu.br", solver="Knitro")
  end

  _set_knitro_options!(model; time=time, ms_seed=ms_seed)

  system.feedback == OL && error("realimentação OL não possui ganho a otimizar")

  @variable(model, 0 <= lambda <= 1)

  invariance_constraints!(model, system, X; find=false)

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
    is_pinvariant(system::SystemDelayed, X; U=nothing, time=10, ms_seed=0, return_model=false) -> InvariantPolyhedron

Verifica se o poliedro X é p-invariante w.r.t

    x[k+1] = system.A*x[k] + system.Ad*x[k - system.d]

O tipo de realimentação é lido de `system.feedback`:

  SOF -> u[k] = G*system.C*x[k]
  SSF -> u[k] = G*x[k]
  OL  -> não suportado (não há ganho a otimizar)

As restrições são adicionadas pelas funções auxiliares
`invariance_constraints!` e `input_constraints!` (ambas em modo `find=false`).
Se `U` for `nothing`, as restrições de entrada são omitidas.

Após a solução do Knitro o ganho encontrado é escrito em `system.G`.
Retorna um `InvariantPolyhedron` com todas as representações relevantes e o
lambda.

Com `return_model=true` retorna um JuMP.Model resultante do problema de otimização

Onde system.d é o tamanho do atraso
"""
function is_pinvariant(system::SystemDelayed, X::Matrix;
                       U::Union{Matrix,Nothing}=nothing,
                       time::Real=10,
                       ms_seed::Int=0,
                       return_model::Bool=false)

  model = Model() do
    return NEOSServer.Optimizer(; email="wallace.lopes.162@ufrn.edu.br", solver="Knitro")
  end

  _set_knitro_options!(model; time=time, ms_seed=ms_seed)

  system.feedback == OL && error("realimentação OL não possui ganho a otimizar")

  @variable(model, 0 <= lambda <= 1)

  invariance_constraints!(model, system, X; find=false)

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
    is_pinvariant(system::LTIDisturbanceDelayedInput, X, D; U=nothing, time=10, ms_seed=0, return_model=false) -> InvariantPolyhedron

Verifica se o poliedro X é p-invariante w.r.t

    x[k+1] = system.A*x[k] + system.B*u[k] + system.E*w[k]
    y[k+1] = system.C*x[k]

com entrada atrasada de `system.d` amostras.

O tipo de realimentação é lido de `system.feedback`:

  SOF -> u[k] = G*system.C*x[k]
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
D é o poliedro da perturbação
U é a restrição na entrada (opcional)
"""
function is_pinvariant(system::LTIDisturbanceDelayedInput, X::Matrix, D::Matrix;
                       U::Union{Matrix,Nothing}=nothing,
                       time::Real=10,
                       ms_seed::Int=0,
                       return_model::Bool=false)

  model = Model() do
    return NEOSServer.Optimizer(; email="wallace.lopes.162@ufrn.edu.br", solver="Knitro")
  end

  _set_knitro_options!(model; time=time, ms_seed=ms_seed)

  system.feedback == OL && error("realimentação OL não possui ganho a otimizar")

  @variable(model, 0 <= lambda <= 1)

  invariance_constraints!(model, system, X, D; find=false)

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
