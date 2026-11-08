"""
    state_constraint_enl!(model; X::Matrix, pond=0.02)

Cria variáveis e restrições para garantir que `X ⊆ igamma * model[:F]`
(isto é, `model[:F]` contém as restrições de estado a menos de escala) e
equilibra a contração `model[:lambda]` com a escala `igamma`.
"""
function state_constraint_enl!(model; X::Matrix, pond::Real=0.02, kwargs...)

  # Parametros
  lf = size(model[:F], 1); #Número de linhas da matriz das restrições de X[k]
  lx = size(X, 1)

  ones_f = ones(lf)
  ones_x = ones(lx)

  @variable(model, 0.001 <= igamma <= 1000)
  @variable(model, 0 <= Z[1:lx, 1:lf] <= 30000)

  @constraint(model, state_enl0, Z * model[:F] == X)
  @constraint(model, state_enl1, Z * ones_f <= ones_x * igamma)

  @objective(model, Min, pond*model[:lambda] + (1-pond)*igamma)

  return nothing
end
