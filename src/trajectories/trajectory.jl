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
