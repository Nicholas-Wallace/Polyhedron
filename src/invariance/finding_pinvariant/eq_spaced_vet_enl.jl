"""
    eq_spaced_vet_enl!(model; t=8, pond=0.02)

Cria variáveis e restrições para forçar o alargamento do poliedro `model[:F]`
ao longo de `t` direções unitárias igualmente espaçadas.

Requer `n == 2` (dimensão do estado), pois `vet_eq_spc` só gera vetores 2D.
"""
function eq_spaced_vet_enl!(model; t::Int=8, pond::Real=0.02, kwargs...)

  # Parametros
  lf = size(model[:F], 1)
  n = size(model[:F], 2)

  n == 2 || error("eq_spaced_vet_enl!: vet_eq_spc só suporta n == 2 (recebido n=$n)")

  ones_f = ones(lf)

  vet = vet_eq_spc(t)

  @variable(model, gamma[1:t] >= 0)

  @constraint(model, [i=1:t], model[:F] * gamma[i]*vet[i] <= ones_f, base_name="enl")

  #Objetivo
  func_obj = pond*model[:lambda] - (1-pond)*sum(gamma)/t
  @objective(model, Min, func_obj)

  return nothing
end
