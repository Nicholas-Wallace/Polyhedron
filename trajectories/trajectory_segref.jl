#############################################
### essa função é referente a um sistema: ###
###     x[k+1] = A*x(k) + Er(k)     ### 
#############################################

function trajectory_segref(x0, A, E, r, passos)
    for i in range(2, passos)
        try
            x0 = hcat(x0, A*x0[:, i-1] + E*r[i-1])    
        catch e
            print(e)
        end
        
    end
    return x0   
end

function trajectory_segref_delay(x0, A, BG, E, r, passos, d; varying=false, reverse=true)
    # Como o vetor vem na forma [x[k]...x[k-d]] para plotar a trajetória
    # é melhor que esteja na ordem cronológica [x[k-d]...x[k]]
    x0_traj = copy(x0)

    if reverse
        reverse!(x0_traj)
    end
    
    for i in 1:passos
        x_atual = collect(x0_traj[end])
        x_atrasado = collect(x0_traj[end - d])
        x_novo = A * x_atual + BG * x_atrasado + E * r[i]

        push!(x0_traj, Tuple(x_novo))   
    end

    return x0_traj[d+1:end]
end