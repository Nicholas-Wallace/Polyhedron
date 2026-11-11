"""
  state_constraints!(model, system::System, X::Matrix)
Adicionando restrições ao problema para garantir que o poliedro nao viole as restrições de estado.

Com `find=false`, o poliedro candidato `X` substitui `F`.

  T*F == X
  T*ones_f == ones_x

Formulação do lema de farkas
"""

function state_constraints!(model, system::System, X::Matrix)
 
  # parametros
  lx = size(X, 1)
  lf = size(model[:F], 1)

  ones_x = ones(lx)
  ones_f = ones(lf)

  @variable(model, 0 <= T[1:lx, 1:lf])

  @constraint(model, x_included0, T * model[:F] == X)
  @constraint(model, x_included1, T * ones_f <= ones_x)

  return nothing
end


