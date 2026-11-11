"""
  input_constraints!(model, system::System, U::Matrix; find::Bool=true, X=nothing)
Adicionando restrições ao problema para garantir que o poliedro nao viole as restrições de entrada.

Com `find=false`, o poliedro candidato `X` substitui `F`.

Essa restrição depende do system.feedback, forma que:

  se for realmentação de estado (SSF)

    M*F == U*G
    M*ones_f == ones_u

  se for realimentação de saída (SOF)

    M*F == U*G*C
    M*ones_f == ones_u
"""

function input_constraints!(model, system::System, U::Matrix, lx::Int; find::Bool=true, X::Union{Matrix,Nothing}=nothing)

  # parametros
  lu = size(U, 1)
  m = size(U, 2)

  if !find && isnothing(X)
    error("input_constraints!: X é obrigatório quando find=false")
  end

  # diferencia verificar invariancia de achar poliedro
  poly = find ? model[:F] : X

  lf = size(poly, 1)
  n = size(poly, 2)

  ones_u = ones(lu)
  ones_f = ones(lf)

  @variable(model, 0 <= M[1:lu, 1:lf] <= 100)

  # para realimentacao de saida
  if system.feedback == SOF

    p = size(system.C, 1)

    @constraint(model, u_included0, M * poly == U * model[:G] * system.C)

  # para realimentacao de estado
  elseif system.feedback == SSF

    @constraint(model, u_included0, M * poly == U * model[:G])
  else
    error("realimentação OL não possui ganho a otimizar")
  end

  @constraint(model, u_included1, M * ones_f <= ones_u)

  return nothing
end

