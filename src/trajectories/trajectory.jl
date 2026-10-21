function _check_state_dim(M, n)
    size(M, 1) == n ||
        throw(DimensionMismatch(
            "cada estado de x0 deve ter $n componentes (dimensão de A), recebido $(size(M, 1))"))
    return M
end

"""
    trajectory(system::SLIT, x0, passos) -> Matrix{Float64}

Gera a trajetória dos estados do sistema `SLIT` em malha fechada. A dinâmica
depende da realimentação escolhida em `system.feedback`:

    SOF -> x[k+1] = (A + B*G*C)*x[k]   (u[k] = G*C*x[k])
    SSF -> x[k+1] = (A + B*G)  *x[k]   (u[k] = G*x[k])
    OL  -> x[k+1] =  A         *x[k]   (sem realimentação)

onde `A = system.A`, `B = system.B`, `C = system.C` e `G = system.G`.

`passos` é o número de iterações. `x0` pode ser:
- uma matriz `n × L`, com cada coluna sendo um estado; ou
- um vetor de estados (ex.: `[x_init]` ou uma lista de estados), ou ainda um
  único estado (vetor de números).

A simulação parte do último estado de `x0`. Cada estado deve ter `n`
componentes, onde `n = size(A, 1)`.

Retorna uma matriz `n × K`, onde cada coluna é um estado (componentes nas
linhas) e `K = length(x0) + passos`; as condições iniciais são mantidas.
"""
function trajectory(system::SLIT, x0::Matrix, passos)
    _check_state_dim(x0, size(system.A, 1))
    traj = Matrix{Float64}(x0)

    if system.feedback == SOF
        for _ in 1:passos
            traj = hcat(traj, system.A * traj[:, end] + system.B * system.G * system.C * traj[:, end])
        end
    elseif system.feedback == SSF
        for _ in 1:passos
            traj = hcat(traj, system.A * traj[:, end] + system.B * system.G * traj[:, end])
        end
    else
        for _ in 1:passos
            traj = hcat(traj, system.A * traj[:, end])
        end
    end

    return traj
end

function trajectory(system::SLIT, x0::AbstractVector, passos)
    traj = eltype(x0) <: Number ? reshape(Float64.(collect(x0)), :, 1) :
                                  hcat([collect(x) for x in x0]...)
    return trajectory(system, traj, passos)
end


"""
    trajectory(system::SystemDelayed, x0, passos) -> Matrix{Float64}

Gera a trajetória dos estados do sistema linear com atraso. A dinâmica
depende da realimentação escolhida em `system.feedback`:

    SOF -> x[k+1] = A*x[k] + Ad*G*C*x[k-d]  (u[k] = G*C*x[k])
    SSF -> x[k+1] = A*x[k] + Ad*G*x[k-d]    (u[k] = G*x[k])
    OL  -> x[k+1] = A*x[k]                  (sem realimentação)

onde `A = system.A`, `Ad = system.Ad`, `C = system.C` e `G = system.G`.

`passos` é o número de iterações. `x0` pode ser:
- uma matriz `n × L`, com cada coluna sendo um estado; ou
- um vetor de estados (ex.: `[x_init]` ou uma lista de estados), ou ainda um
  único estado (vetor de números).

A simulação parte do último estado de `x0`. Cada estado deve ter `n`
componentes, onde `n = size(A, 1)`, e `x0` deve conter pelo menos `d + 1`
estados para que o atraso `d = system.d` seja acessível.

Retorna uma matriz `n × K`, onde cada coluna é um estado (componentes nas
linhas) e `K = length(x0) + passos`; as condições iniciais são mantidas.
"""
function trajectory(system::SystemDelayed, x0::Matrix, passos)
    n = size(system.A, 1)
    d = system.d
    _check_state_dim(x0, n)

    size(x0, 2) >= d + 1 ||
        throw(ArgumentError(
            "x0 deve ter pelo menos $(d + 1) estados (colunas) para atraso d = $d, recebido $(size(x0, 2))"))

    traj = Matrix{Float64}(x0)

    if system.feedback == SOF
        for _ in 1:passos
            traj = hcat(traj, system.A * traj[:, end] + system.Ad * system.G * system.C * traj[:, end - d])
        end
    elseif system.feedback == SSF
        for _ in 1:passos
            traj = hcat(traj, system.A * traj[:, end] + system.Ad * system.G * traj[:, end - d])
        end
    else
        for _ in 1:passos
            traj = hcat(traj, system.A * traj[:, end])
        end
    end

    return traj
end

function trajectory(system::SystemDelayed, x0::AbstractVector, passos)
    traj = eltype(x0) <: Number ? reshape(Float64.(collect(x0)), :, 1) :
                                  hcat([collect(x) for x in x0]...)
    return trajectory(system, traj, passos)
end

function _check_disturbance_dim(W, q, passos)
    size(W, 1) == q ||
        throw(DimensionMismatch(
            "cada distúrbio de W deve ter $q componentes (dimensão de E), recebido $(size(W, 1))"))
    size(W, 2) >= passos ||
        throw(DimensionMismatch(
            "W deve ter pelo menos $passos colunas (uma por passo), recebido $(size(W, 2))"))
    return W
end

"""
    trajectory(system::LTIDisturbance, x0, passos) -> Matrix{Float64}
    trajectory(system::LTIDisturbance, x0, W, passos) -> Matrix{Float64}

Gera a trajetória dos estados do sistema linear com distúrbio `w(k)`. A dinâmica
depende da realimentação escolhida em `system.feedback`:

    SOF -> x[k+1] = (A + B*G*C)*x[k] + E*w[k]   (u[k] = G*C*x[k])
    SSF -> x[k+1] = (A + B*G)  *x[k] + E*w[k]   (u[k] = G*x[k])
    OL  -> x[k+1] =  A         *x[k] + E*w[k]   (sem realimentação)

onde `A = system.A`, `B = system.B`, `E = system.E`, `C = system.C` e
`G = system.G`.

`W` é uma matriz `q × passos`, onde `q = size(E, 2)`; a coluna `W[:, i]` é o
distúrbio usado na i-ésima transição simulada. Se `W` for omitido, considera-se
distúrbio nulo.

`passos` é o número de iterações. `x0` pode ser:
- uma matriz `n × L`, com cada coluna sendo um estado; ou
- um vetor de estados (ex.: `[x_init]` ou uma lista de estados), ou ainda um
  único estado (vetor de números).

A simulação parte do último estado de `x0`. Cada estado deve ter `n`
componentes, onde `n = size(A, 1)`.

Retorna uma matriz `n × K`, onde cada coluna é um estado (componentes nas
linhas) e `K = size(x0, 2) + passos`; as condições iniciais são mantidas.
"""
function trajectory(system::LTIDisturbance, x0::Matrix, W::AbstractMatrix, passos)
    n = size(system.A, 1)
    _check_state_dim(x0, n)
    _check_disturbance_dim(W, size(system.E, 2), passos)

    traj = Matrix{Float64}(x0)

    if system.feedback == SOF
        for i in 1:passos
            x = traj[:, end]
            traj = hcat(traj, system.A * x + system.B * system.G * system.C * x + system.E * W[:, i])
        end
    elseif system.feedback == SSF
        for i in 1:passos
            x = traj[:, end]
            traj = hcat(traj, system.A * x + system.B * system.G * x + system.E * W[:, i])
        end
    else
        for i in 1:passos
            x = traj[:, end]
            traj = hcat(traj, system.A * x + system.E * W[:, i])
        end
    end

    return traj
end

function trajectory(system::LTIDisturbance, x0::AbstractVector, W::AbstractMatrix, passos)
    traj = eltype(x0) <: Number ? reshape(Float64.(collect(x0)), :, 1) :
                                  hcat([collect(x) for x in x0]...)
    return trajectory(system, traj, W, passos)
end

function trajectory(system::LTIDisturbance, x0, passos)
    W = zeros(size(system.E, 2), passos)
    return trajectory(system, x0, W, passos)
end

"""
    trajectory(system::LTIDisturbanceDelayedInput, x0, passos) -> Matrix{Float64}
    trajectory(system::LTIDisturbanceDelayedInput, x0, W, passos) -> Matrix{Float64}

Gera a trajetória dos estados do sistema linear com distúrbio `w(k)` e entrada
atrasada. A dinâmica depende da realimentação escolhida em `system.feedback`:

    SOF -> x[k+1] = A*x[k] + B*G*C*x[k-d] + E*w[k]  (u[k] = G*C*x[k])
    SSF -> x[k+1] = A*x[k] + B*G*x[k-d]   + E*w[k]  (u[k] = G*x[k])
    OL  -> x[k+1] = A*x[k]                + E*w[k]  (sem realimentação)

onde `A = system.A`, `B = system.B`, `E = system.E`, `C = system.C`,
`G = system.G` e `d = system.d`.

`W` é uma matriz `q × passos`, onde `q = size(E, 2)`; a coluna `W[:, i]` é o
distúrbio usado na i-ésima transição simulada. Se `W` for omitido, considera-se
distúrbio nulo.

`passos` é o número de iterações. `x0` pode ser:
- uma matriz `n × L`, com cada coluna sendo um estado; ou
- um vetor de estados (ex.: `[x_init]` ou uma lista de estados), ou ainda um
  único estado (vetor de números).

A simulação parte do último estado de `x0`. Cada estado deve ter `n`
componentes, onde `n = size(A, 1)`, e `x0` deve conter pelo menos `d + 1`
estados para que o atraso `d = system.d` seja acessível.

Retorna uma matriz `n × K`, onde cada coluna é um estado (componentes nas
linhas) e `K = size(x0, 2) + passos`; as condições iniciais são mantidas.
"""
function trajectory(system::LTIDisturbanceDelayedInput, x0::Matrix, W::AbstractMatrix, passos)
    n = size(system.A, 1)
    d = system.d
    _check_state_dim(x0, n)

    size(x0, 2) >= d + 1 ||
        throw(ArgumentError(
            "x0 deve ter pelo menos $(d + 1) estados (colunas) para atraso d = $d, recebido $(size(x0, 2))"))

    _check_disturbance_dim(W, size(system.E, 2), passos)

    traj = Matrix{Float64}(x0)

    if system.feedback == SOF
        for i in 1:passos
            x = traj[:, end]
            xd = traj[:, end - d]
            traj = hcat(traj, system.A * x + system.B * system.G * system.C * xd + system.E * W[:, i])
        end
    elseif system.feedback == SSF
        for i in 1:passos
            x = traj[:, end]
            xd = traj[:, end - d]
            traj = hcat(traj, system.A * x + system.B * system.G * xd + system.E * W[:, i])
        end
    else
        for i in 1:passos
            x = traj[:, end]
            traj = hcat(traj, system.A * x + system.E * W[:, i])
        end
    end

    return traj
end

function trajectory(system::LTIDisturbanceDelayedInput, x0::AbstractVector, W::AbstractMatrix, passos)
    traj = eltype(x0) <: Number ? reshape(Float64.(collect(x0)), :, 1) :
                                  hcat([collect(x) for x in x0]...)
    return trajectory(system, traj, W, passos)
end

function trajectory(system::LTIDisturbanceDelayedInput, x0, passos)
    W = zeros(size(system.E, 2), passos)
    return trajectory(system, x0, W, passos)
end
