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

            x_atual = collect(x0_traj[end])
            x_atrasado = collect(x0_traj[end - d_k])
            x_novo = A * x_atual + Ad * x_atrasado

            push!(x0_traj, Tuple(x_novo))
        end
    end

    return x0_traj[d+1:end]
end