"""
    trajectory_segref(x0, A, E, r, passos) -> Vector{Float}

Simula a trajetória w.r.t

    x[k+1] = A*x(k) + Er(k)
"""
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

"""
    trajectory_segref_delay(x0, A, BG, E, r, passos, d; varying=false, reverse=true) -> Vector{Float}

Simula a trajetória w.r.t

    x[k+1] = A*x(k) + BG*x[k - d] + Er(k)

Admite atraso de tamanho d

"""

function trajectory_segref_delay(x0, A, BG, E, r, passos, d; varying=false, reverse=true)
    # Como o vetor vem na forma [x[k]...x[k-d]] para plotar a trajetória
    # é melhor que esteja na ordem cronológica [x[k-d]...x[k]]
    x0_traj = copy(x0)

    if reverse
        reverse!(x0_traj)
    end
    
    for i in 1:passos
        x_atual = x0_traj[end]
        x_atrasado = x0_traj[end - d]
        x_novo = A * x_atual + BG * x_atrasado + E * r[i]
        # precisamos resolver a tipagem na trajetoria
        push!(x0_traj, vec(x_novo))   
    end

    return x0_traj[d:end]
end

function expanded_state_trajectory(traj)
    # plotando o espaco aumentado do exemplo do artigo sem o v(k)
    plot(traj[:, 1:3], 
     label=["x1" "x2" "w"], 
     xlabel="k",  
     title="X[k]")
end