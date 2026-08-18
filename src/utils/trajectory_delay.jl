"""
    trajectory_delay(x0, A, B, G, passos, d; varying=false, reverse=true) -> Vector{Tuple{Float64}}

    Gera a trajetoria dos estados de um sistema com delay do tipo:

        x[k+1] = A*x[k] + B*G*x[k-d]

    `x0` é o vetor de condições iniciais e deve estar como Vector{Tuple{Float64}}, e na ordem
    [x[k-d]...x[k]]
    `G` é a matriz de ganho do controlador
    `varying` é uma flag para indicar se o delay é variável ou fixo
    `reverse` é uma flag para indicar se o vetor de condições iniciais deve ser invertido
    para ficar na ordem correta.

    Retorna a trajetória sem as condições iniciais.
"""

function trajectory_delay(x0, A, Ad, passos, d; varying=false, reverse=true)
    # Como o vetor vem na forma [x[k]...x[k-d]] para plotar a trajetória
    # é melhor que esteja na ordem cronológica [x[k-d]...x[k]]
    x0_traj = copy(x0)

    if reverse
        reverse!(x0_traj)
    end
    
    if varying == false
        for i in 1:passos
            x_atual = collect(x0_traj[end])
            x_atrasado = collect(x0_traj[end - d])
            x_novo = A * x_atual + Ad * x_atrasado

            push!(x0_traj, Tuple(x_novo))   
        end
    else
        for i in 1:passos
            d_k = rand(1:d)

            # porque tem esse "collect()"?
            x_atual = collect(x0_traj[end])
            x_atrasado = collect(x0_traj[end - d_k])
            x_novo = A * x_atual + Ad * x_atrasado

            push!(x0_traj, Tuple(x_novo))
        end
    end

    return x0_traj[d+1:end]
end

# leva em conta a saturação e a trajetoria + um vetor das entradas [sat(x[d]), sat(x[d-1]), ..., sat(x[passos])]
# w.r.t. x(k+1) = Ax(k) + B*sat(Gx(k-d))
function  trajectory_delay_sat(x0, A, Ad, passos, d, umax, umin)
    x0_traj = copy(x0)

    for i in 1:passos
        d_k = rand(1:d)

    # porque tem esse "collect()"?
        x_atual = collect(x0_traj[end])
        x_atrasado = collect(x0_traj[end - d_k])
        x_novo = A * x_atual + Ad * x_atrasado

        push!(x0_traj, Tuple(x_novo))
    end

    return x0_traj[d+1:end], sat(x_traj[1:end-d])
end

function sat(x, umax, umin)
    if x > umax return umax
    elseif x < umin return umin
    else return x
    end
end

