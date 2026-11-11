"""
  invariance_constraints!(model, system::System, X::Matrix; find::Bool=false)

Add constrainsts to ensure invariance for a polyhedron.

```
  x(k) \\in \\Omega => x(k+1) \\in \\Omega
```

We use Farka's Lemma to write constraints that need to be satisfacted.
"""
function invariance_constraints!(model, system::SLIT, X::Matrix; find::Bool=false)

  system.feedback == OL && error("realimentação OL não possui ganho a otimizar")

  if find

    lf = size(model[:F], 1)
    ones_f = ones(lf)

    n, m, p = checkdims(system.A, system.B, system.C)

    # Declaring problem variables
    @variable(model, 0 <= H[1:lf, 1:lf] <= 100) # Farka's aux matrix

    # The Gain matrix dimension depends on the feedback
    if system.feedback == SOF
      @variable(model, G[1:m, 1:p])
      A_cl = system.A + system.B*G*system.C # closed-loop
    else
      @variable(model, G[1:m, 1:n])
      A_cl = system.A + system.B*G # closed-loop
    end


    @constraint(model, contractiave_pi0, H*model[:F] == model[:F]*A_cl)
    @constraint(model, contractiave_pi1, H*ones_f <= ones_f*model[:lambda])

    # é preciso ainda que Rank(F) == n
    @variable(model, J[1:n, 1:lf])
    @constraint(model, pseudo_inv, J*model[:F] == I(n))

  else

    lx = size(X, 1)
    ones_x = ones(lx)

    n, m, p = checkdims(system.A, system.B, system.C)

    # Declaring problem variables
    @variable(model, 0 <= H[1:lx, 1:lx] <= 100) # Farka's aux matrix

    # The Gain matrix dimension depends on the feedback
    if system.feedback == SOF
      @variable(model, G[1:m, 1:p])
      A_cl = system.A + system.B*G*system.C # closed-loop
    else
      @variable(model, G[1:m, 1:n])
      A_cl = system.A + system.B*G # closed-loop
    end


    @constraint(model, contractiave_pi0, H*X == X*A_cl)
    @constraint(model, contractiave_pi1, H*ones_x <= ones_x*model[:lambda])   

  end
  return nothing
end


function invariance_constraints!(model, system::SystemDelayed, X::Matrix; find=true)

  # Restrições baseadas nas condições de invariancia delay-dependent

  system.feedback == OL && error("realimentação OL não possui ganho a otimizar")

  if find
    lf = size(model[:F], 1)
    ones_f = ones(lf)

    n, m, p = checkdims(system.A, system.Ad, system.C)

    # Declaring problem variable

    # aux variable to transformed model
    @variable(model, K[1:n, 1:n]) 
    
    # Farka's aux matrix
    @variable(model, 0 <= H[1:lf, 1:lf])
    @variable(model, 0 <= L[1:lf, 1:lf])
    @variable(model, 0 <= M[1:lf, 1:lf])
    @variable(model, 0 <= N[1:lf, 1:lf])

    # The Gain matrix dimension depends on the feedback
    if system.feedback == SOF
      @variable(model, G[1:m, 1:p])
      U = system.Ad*G*system.C # closed-loop
    else
      @variable(model, G[1:m, 1:n])
      U = system.Ad*G # closed-loop
    end


    @constraint(model, contractiave_pi_H, H*model[:F] == model[:F]*(system.A + K))
    @constraint(model, contractiave_pi_L, L*model[:F] == model[:F]*(U - K))
    @constraint(model, contractiave_pi_M, M*model[:F] == model[:F]*K*(system.A - I(n)))
    @constraint(model, contractiave_pi_N, N*model[:F] == model[:F]*K*U)
    
    @constraint(model, contractiave_pi, (H + L + system.d*(M+N))*ones_f <= ones_f*model[:lambda])

    # é preciso ainda que Rank(F) == n
    @variable(model, J[1:n, 1:lf])
    @constraint(model, pseudo_inv, J*model[:F] == I(n))

  else

    n, m, p = checkdims(system.A, system.Ad, system.C)
    lx = size(X, 1)
    ones_x = ones(lx)

    # Declaring problem variable

    # aux variable to transformed model
    @variable(model, K[1:n, 1:n]) 
    
    # Farka's aux matrix
    @variable(model, 0 <= H[1:lx, 1:lx])
    @variable(model, 0 <= L[1:lx, 1:lx])
    @variable(model, 0 <= M[1:lx, 1:lx])
    @variable(model, 0 <= N[1:lx, 1:lx])

    # The Gain matrix dimension depends on the feedback
    if system.feedback == SOF
      @variable(model, G[1:m, 1:p])
      U = system.Ad*G*system.C # closed-loop
    else
      @variable(model, G[1:m, 1:n])
      U = system.Ad*G # closed-loop
    end


    @constraint(model, contractiave_pi_H, H*X == X*(system.A + K))
    @constraint(model, contractiave_pi_L, L*X == X*(U - K))
    @constraint(model, contractiave_pi_M, M*X == X*K*(system.A - I(n)))
    @constraint(model, contractiave_pi_N, N*X == X*K*U)
    
    @constraint(model, contractiave_pi, (H + L + system.d*(M+N))*ones_x <= ones_x*model[:lambda])

  end

  return nothing
end

function invariance_constraints!(model, system::LTIDisturbanceDelayedInput, X::Matrix, D::Matrix; find::Bool=true)

  system.feedback == OL && error("realimentação OL não possui ganho a otimizar")

  if find
    lf = size(model[:F], 1)
    ld = size(D, 1)
    lr = size(system.E, 2) # ou ce
    ones_f = ones(lf)
    ones_d = ones(ld)

    n, m, p = checkdims(system.A, system.B, system.C)

    # Declaring problem variable

    # aux variable to transformed model
    @variable(model, K[1:n, 1:n]) 
    
    # Farka's aux matrix
    @variable(model, 0 <= H[1:lf, 1:lf])
    @variable(model, 0 <= L[1:lf, 1:lf])
    @variable(model, 0 <= M[1:lf, 1:lf])
    @variable(model, 0 <= N[1:lf, 1:lf])
    # aux variables to disturbance
    @variable(model, 0 <= P[1:ld, 1:ld])
    @variable(model, 0 <= Q[1:ld, 1:ld])

    # The Gain matrix dimension depends on the feedback
    if system.feedback == SOF
      @variable(model, G[1:m, 1:p])
      U = system.B*G*system.C # closed-loop
    else
      @variable(model, G[1:m, 1:n])
      U = system.B*G # closed-loop
    end


    @constraint(model, contractiave_pi_H, H*model[:F] == model[:F]*(system.A + K))
    @constraint(model, contractiave_pi_L, L*model[:F] == model[:F]*(U - K))
    @constraint(model, contractiave_pi_M, M*model[:F] == model[:F]*K*(system.A - I(n)))
    @constraint(model, contractiave_pi_N, N*model[:F] == model[:F]*K*U)
    @constraint(model, contractiave_pi_P, P*D == D*system.E)
    @constraint(model, contractiave_pi_Q, Q*D == -D*K*system.E)

    @constraint(model, contractiave_pi, (H + L + system.d*(M+N))*ones_f + (P+system.d*Q)*ones_d <= ones_f*model[:lambda])

    # é preciso ainda que Rank(F) == n
    @variable(model, J[1:n, 1:lf])
    @constraint(model, pseudo_inv, J*model[:F] == I(n))

  else

    n, m, p = checkdims(system.A, system.B, system.C)
    lx = size(X, 1)
    ld = size(D, 1)
    lr = size(system.E, 2) # ou ce
    ones_x = ones(lx)
    ones_d = ones(ld)

    # Declaring problem variable

    # aux variable to transformed model
    @variable(model, K[1:n, 1:n]) 
    
    # Farka's aux matrix
    @variable(model, 0 <= H[1:lx, 1:lx])
    @variable(model, 0 <= L[1:lx, 1:lx])
    @variable(model, 0 <= M[1:lx, 1:lx])
    @variable(model, 0 <= N[1:lx, 1:lx])
    # aux variables to disturbance
    @variable(model, 0 <= P[1:ld, 1:ld] <= 100)
    @variable(model, 0 <= Q[1:ld, 1:ld] <= 100)

    # The Gain matrix dimension depends on the feedback
    if system.feedback == SOF
      @variable(model, G[1:m, 1:p])
      U = system.B*G*system.C # closed-loop
    else
      @variable(model, G[1:m, 1:n])
      U = system.B*G # closed-loop
    end


    @constraint(model, contractiave_pi_H, H*X == X*(system.A + K))
    @constraint(model, contractiave_pi_L, L*X == X*(U - K))
    @constraint(model, contractiave_pi_M, M*X == X*K*(system.A - I(n)))
    @constraint(model, contractiave_pi_N, N*X == X*K*U)
    @constraint(model, contractiave_pi_P, P*D == D*system.E)
    @constraint(model, contractiave_pi_Q, Q*D == -D*K*system.E)

    @constraint(model, contractiave_pi, (H + L + system.d*(M+N))*ones_x + (P+system.d*Q)*ones_d <= ones_x*model[:lambda])
    
  end

  return nothing
end
