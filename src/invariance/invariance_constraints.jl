"""
  invariance_constraints(model, system::System)

Add constrainsts to ensure invariance for a polyhedrons.

```
  x(k) \in \Omega => x(k+1) \in \Omega
```

We use Farka's Lemma to write constraints that need to be satisfacted.
"""

function invariance_constraints!(model, system::System)
  return invariance_constraints(mode, system)

function invariance_constraints!(model, system::SLIT)

  lf = size(model[:F], 1) #use F or X? finding vs verificate
  ones_f = ones(lf)

  n, m, p = checkdims(system.A, system.B, system.C)

  # Declaring problem variables
  @variable(model, 0 <= H[1:lf, 1:lf] <= 100) # Farka's aux matrix

  # The Gain matrix dimension depends on the feedback
  if system.feedback == SOF
    @variable(model, G[1:m, 1:p]):
    A_cl = system.A + system.B*G*system.C # closed-loop
  else
    @variable(model, G[1:m, 1:n])
    A_cl = system.A + system.B*G # closed-loop
  end


  @constraint(model, contractiave_pi0, H*model[:F] == model[:F]*A_cl)
  @constraint(model, contractiave_pi1, H*ones_f <= ones_f*model[:lambda])
end
