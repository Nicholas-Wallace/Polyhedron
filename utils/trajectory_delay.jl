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
function  trajectory_delay(x0, A, Ad, passos, d, umax, umin)
   x0_traj = copy(x0)
  
  function sat(x, umax, umin)
    if x > umax return umax
    elseif x < umin return umin
    else return x
    end
  end

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

