"""
    finding_L_pinvariant(system::SLIT, U, X; ll = 6, t = 8, pond = 0.02, return_model = false) -> InvariantPolyhedron

`U` pode ser uma `Matrix` ou `nothing` (nesse último caso, também se pode chamar
`finding_L_pinvariant(system, X; ...)`).

Procuramos um poliedro de f linhas que seja invariante w.r.t

    x[k+1] = system.A*x[k] + system.B*u[k]
    y[k+1] = system.C*x[k]

X é o poliedro de restrições 

Essa função suporta realimentação de estado e de saída

Retorna um `InvariantPolyhedron` com o poliedro encontrado (`L x <= ones(ll)`),
o ganho encontrado (escrito em `system.G`) e o valor de lambda.

Com o return_model = true, a função retorna o JuMP.Model resultante, que pode ser utilizado para testar as restrições com as variáveis auxiliares criadas 
"""

function pinvariant_poly(system::SLIT, X::Matrix; kwargs...)

  # coletando kwargs
  U            = get(kwargs, :U, nothing)
  lf           = get(kwargs, :lf, 6)
  t            = get(kwargs, :t, 8)
  pond         = get(kwargs, :pond, 0.02)
  time         = get(kwargs, :time, 10)
  ms_seed      = get(kwargs, :ms_seed, 0)
  F_init       = get(kwargs, :F_init, nothing)
  return_model = get(kwargs, :return_model, false)

  # creating model
  model = Model() do
    return NEOSServer.Optimizer(; email="wallace.lopes.162@ufrn.edu.br", solver="Knitro")
  end

  _set_knitro_options!(model; time=time, ms_seed=ms_seed)

  # Create F Polyhedron and Lambda
  n = size(system.A, 1)
  ones_f = ones(lf)
  @variable(model, -500 <= F[1:lf, 1:n] <= 500)
  if F_init !== nothing
    size(F_init) == (lf, n) ||
      throw(DimensionMismatch("F_init must be $(lf)×$(n), got $(size(F_init))"))
    set_start_value.(F, F_init) # solução inicial
    set_optimizer_attribute(model, "ms_enable", 0) # desativando multistart

    @assert all(has_start_value, F)

  end
  @variable(model, 0.001 <= lambda <= 1)

  # Problem constraints
  invariance_constraints!(model, system, X; find=true)
  state_constraints!(model, system, X)

  if U !== nothing
    input_constraints!(model, system, U, size(system.A, 1); find=true)
  end

  # por padrao enl é eq_spaced_vet_enl!
  enl = get(kwargs, :enl, eq_spaced_vet_enl!)
  # TODO: test if enl is a valid function
  # its important to let the user create you own function
  enl(model; X=X, kwargs...)

  if return_model
    return model
  end

  optimize!(model)
  @show termination_status(model)

  if termination_status(model) != MathOptInterface.LOCALLY_SOLVED
    if return_model
      # return the model even if the problem was not solved
      return model
    else
      error("Não foi possivel encontrar um poliedro com essas restrições")
    end
  end

  system.G = value.(model[:G]) # Update gains
  ip = InvariantPolyhedron(value.(F), ones_f, value.(lambda), system)

  return return_model ? model : ip
end


function pinvariant_poly(system::SystemDelayed, X::Matrix; kwargs...)

  # coletando kwargs
  U            = get(kwargs, :U, nothing)
  lf           = get(kwargs, :lf, 6)
  t            = get(kwargs, :t, 8)
  pond         = get(kwargs, :pond, 0.02)
  time         = get(kwargs, :time, 10)
  ms_seed      = get(kwargs, :ms_seed, 0)
  F_init       = get(kwargs, :F_init, nothing)
  return_model = get(kwargs, :return_model, false)

  # creating model
  model = Model() do
    return NEOSServer.Optimizer(; email="wallace.lopes.162@ufrn.edu.br", solver="Knitro")
  end

  _set_knitro_options!(model; time=time, ms_seed=ms_seed)

  # Create F Polyhedron and Lambda
  n = size(system.A, 1)
  ones_f = ones(lf)
  @variable(model, -500 <= F[1:lf, 1:n] <= 500)
  if F_init !== nothing
    size(F_init) == (lf, n) ||
      throw(DimensionMismatch("F_init must be $(lf)×$(n), got $(size(F_init))"))
    set_start_value.(F, F_init) # solução inicial
    set_optimizer_attribute(model, "ms_enable", 0) # desativando multistart

    @assert all(has_start_value, F)

  end
  @variable(model, 0.001 <= lambda <= 1)

  # Problem constraints
  invariance_constraints!(model, system, X)
  state_constraints!(model, system, X)

  if U !== nothing
    input_constraints!(model, system, U, size(system.A, 1); find=true)
  end

  # por padrao enl é eq_spaced_vet_enl!
  enl = get(kwargs, :enl, eq_spaced_vet_enl!)
  # TODO: test if enl is a valid function
  # its important to let the user create you own function
  enl(model; X=X, kwargs...)

  if return_model
    return model
  end

  optimize!(model)
  @show termination_status(model)

  if termination_status(model) != MathOptInterface.LOCALLY_SOLVED
    error("Não foi possivel encontrar um poliedro com essas restrições")
  end

  system.G = value.(model[:G]) # Update gains
  ip = InvariantPolyhedron(value.(F), ones_f, value.(lambda), system)

  return ip
end


function pinvariant_poly(system::LTIDisturbanceDelayedInput, X::Matrix, D::Matrix; kwargs...)

  # coletando kwargs
  U            = get(kwargs, :U, nothing)
  lf           = get(kwargs, :lf, 6)
  t            = get(kwargs, :t, 8)
  pond         = get(kwargs, :pond, 0.02)
  time         = get(kwargs, :time, 10)
  ms_seed      = get(kwargs, :ms_seed, 0)
  F_init       = get(kwargs, :F_init, nothing)
  return_model = get(kwargs, :return_model, false)

  # creating model
  model = Model() do
    return NEOSServer.Optimizer(; email="wallace.lopes.162@ufrn.edu.br", solver="Knitro")
  end

  _set_knitro_options!(model; time=time, ms_seed=ms_seed)

  # Create F Polyhedron and Lambda
  n = size(system.A, 1)
  ones_f = ones(lf)
  @variable(model, -500 <= F[1:lf, 1:n] <= 500)
  if F_init !== nothing
    size(F_init) == (lf, n) ||
      throw(DimensionMismatch("F_init must be $(lf)×$(n), got $(size(F_init))"))
    set_start_value.(F, F_init) # solução inicial
    set_optimizer_attribute(model, "ms_enable", 0) # desativando multistart

    @assert all(has_start_value, F)
    
  end
  @variable(model, 0.001 <= lambda <= 1)

  # Problem constraints
  invariance_constraints!(model, system, X, D)
  state_constraints!(model, system, X)

  if U !== nothing
    input_constraints!(model, system, U, size(system.A, 1); find=true)
  end

  # por padrao enl é eq_spaced_vet_enl!
  enl = get(kwargs, :enl, eq_spaced_vet_enl!)
  # TODO: test if enl is a valid function
  # its important to let the user create you own function
  enl(model; X=X, kwargs...)

  if return_model
    return model
  end

  optimize!(model)
  @show termination_status(model)

  if termination_status(model) != MathOptInterface.LOCALLY_SOLVED
    error("Não foi possivel encontrar um poliedro com essas restrições")
  end

  system.G = value.(model[:G]) # Update gains
  ip = InvariantPolyhedron(value.(F), ones_f, value.(lambda), system)

  return ip
end

function _set_knitro_options!(model; time::Real=10, ms_seed::Int=0)
  set_optimizer_attribute(model, "outlev", 2)
  set_optimizer_attribute(model, "algorithm", 2)
  set_optimizer_attribute(model, "maxtime", time * 60.0)
  set_optimizer_attribute(model, "maxit", 10000)
  set_optimizer_attribute(model, "ms_enable", 1)
  set_optimizer_attribute(model, "ms_maxsolves", 100)
  set_optimizer_attribute(model, "ms_seed", ms_seed)
  set_optimizer_attribute(model, "numthreads", 4)
  return model
end


