
"""
    is_pinvariant(system::LTISystem, U, X; SOF = false) -> Float

Verifica se o poliedro X é p-invariante w.r.t

    x[k+1] = system.A*x[k] + system.B*u[k]
    y[k+1] = system.C*x[k]

E retorna o valor de lambda

Essa função suporta realimentação de estado e de saída

X é o poliedro a ser verificado
U é a restrição na entrada
"""

function is_pinvariant(system::LTISystem, U::Matrix{Float64}, X::Matrix{Float64}; SOF::Bool=false)

  model = Model() do
    return NEOSServer.Optimizer(; email="wallace.lopes.162@ufrn.edu.br", solver="Knitro")
  end

  #Parâmetros
  n = size(system.A, 1); #Ordem do sistema (Linhas de A)
  m = size(system.B, 2); #Número de Entradas (Colunas de B)
  p = size(system.C, 1); #Número de Saídas (Linhas de C)

  lx = size(X, 1); #Número de linhas de X
  lu = size(U, 1); #Número de linhas de U

  xi = ones(lx)
  phi = ones(lu);


  if SOF # realimentação de saida 
    @variable(model, 0 <= H[1:lx, 1:lx] <= 100)
    @variable(model, K[1:m, 1:p])
    @variable(model, 0 <= M[1:lu, 1:lx] <= 100)
    @variable(model, 0 <= lambda <= 1)

    @constraint(model, contractiave_pi0, H*X == X*(system.A + system.B*K*system.C))
    @constraint(model, u_included0, M*X == U*K*system.C)
    @constraint(model, contractiave_pi1, H*xi <= xi*lambda)
    @constraint(model, u_included1, M*xi <= phi)

    @objective(model, Min, lambda)

    optimize!(model)

    @show termination_status(model)
    if termination_status(model) != MathOptInterface.LOCALLY_SOLVED
      error("Poliedro não é invariante")
    end


    H = value.(H)
    K = value.(K)
    M = value.(M)
    lambda = value.(lambda)

    result = Dict("H" => H, "K" => K, "M" => M, "lambda" => lambda)

    return result

  else # realimentação de estado

    @variable(model, 0 <= H[1:lx, 1:lx] <= 100)
    @variable(model, F[1:m, 1:n])
    @variable(model, 0 <= M[1:lu, 1:lx] <= 100)
    @variable(model, 0 <= lambda <= 1)


    @constraint(model, contractiave_pi0, H*X == X*(system.A + (system.B)*F))
    @constraint(model, u_included0, M*X == U*F)
    @constraint(model, contractiave_pi1, H*xi <= xi*lambda)
    @constraint(model, u_included1, M*xi <= phi)

    @objective(model, Min, lambda)
    optimize!(model)


    @show termination_status(model)

    if termination_status(model) != MathOptInterface.LOCALLY_SOLVED
      error("polyhedro não é invariante")
    end


    H = value.(H)
    F = value.(F)
    M = value.(M)
    lambda = value.(lambda)
    result = Dict("H" => H, "F" => F, "M" => M, "lambda" => lambda)

    return result

  end
end


"""
    is_pinvariant(system::SystemDelayed, X; symetric=true) -> Dict{String, Float}

Verifica se o poliedro X é p-invariante w.r.t

    x[k+1] = system.A*x[k] + system.Ad*x[k - system.d]

E retorna um dicionario com todas as matrizes relevantes e o lambda

Onde system.d é o tamanho do atraso
"""
function is_pinvariant(system::SystemDelayed, X; symetric=true)
  model = Model() do
    return NEOSServer.Optimizer(; email="wallace.lopes.162@ufrn.edu.br", solver="Knitro")
  end

  # parâmetros

  n = size(system.A, 1) # ordem do sistema
  x = size(X, 1) # numero de linhas de x
  w = ones(x)

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
    H = value.(H1) + value.(H2)
    L = value.(L1) + value.(L2)

    # retorna um dicionário com lambda e as matrizes K, H e L
    result = Dict("lambda" => lambda, "K" => K, "H" => H, "L" => L)

    return result
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

  lambda = value(lambda)
  K = value.(K)
  H = value.(H)
  L = value.(L)

  result = Dict("lambda" => lambda, "K" => K, "H" => H, "L" => L)

  return result
end
