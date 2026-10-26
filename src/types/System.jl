abstract type System end

@enum Feedback SOF SSF OL

"""
A Linear Time Invariant System defined in the space state as:

  x[k+1] = Ax[k] + Bu[k]
  y[k] = Cx[k]

Where u[k] depends on the feedback indication
  SOF -> u[k] = Gy[k] = GCx[k]
  SSF -> u[k] = Gx[k]
  OL  -> There is no feedback, it's a open-loop system 

"""
mutable struct SLIT <: System
  A::Matrix{Float64}
  B::Matrix{Float64}
  C::Matrix{Float64}
  feedback::Feedback
  G::Matrix{Float64}
end

"""
Dimensions of the feedback gain matrix G, where m = size(B, 2):

  SOF -> m x p (p = size(C, 1))
  SSF -> m x n (n = size(A, 1))
  OL  -> m x 0 (no feedback)
"""
function gdims(feedback::Feedback, A, B, C)
  m = size(B, 2)
  feedback == SOF && return (m, size(C, 1))
  feedback == SSF && return (m, size(A, 1))
  return (m, 0)
end

"""
    checkdims(A, B, C) -> (n, m, p)

Validates the state-space dimensions: `A` must be square (`n×n`), `B` must have
`n` rows (`n×m`) and `C` must have `n` columns (`p×n`). Throws `DimensionMismatch`
and otherwise returns `(n, m, p)`.
"""
function checkdims(A, B, C)
  n = size(A, 1)
  size(A, 2) == n ||
    throw(DimensionMismatch("A must be square (n×n), got $(size(A))"))
  size(B, 1) == n ||
    throw(DimensionMismatch("B must have $n rows (n×m), got $(size(B))"))
  size(C, 2) == n ||
    throw(DimensionMismatch("C must have $n columns (p×n), got $(size(C))"))
  return n, size(B, 2), size(C, 1)
end

"""
    checkdims(A, B, E, C) -> (n, m, p, q)

Validates the state-space dimensions of a system with a disturbance input:
`A` must be square (`n×n`), `B` and `E` must have `n` rows (`n×m` and `n×q`)
and `C` must have `n` columns (`p×n`). Throws `DimensionMismatch` and otherwise
returns `(n, m, p, q)`.
"""
function checkdims(A, B, E, C)
  n, m, p = checkdims(A, B, C)
  size(E, 1) == n ||
    throw(DimensionMismatch("E must have $n rows (n×q), got $(size(E))"))
  return n, m, p, size(E, 2)
end


function SLIT(A, B, C, feedback::Feedback)
  checkdims(A, B, C)
  SLIT(A, B, C, feedback, ones(gdims(feedback, A, B, C)...))
end

SLIT(A, B, C) = SLIT(A, B, C, OL)

"""
    set_feedback!(sys::SLIT, feedback::Feedback)

Updates the feedback of `sys` and resizes `G` to the matching dimensions.
Note that the previous value of `G` will be lost.
"""
function set_feedback!(sys::SLIT, feedback::Feedback)
  sys.feedback = feedback
  sys.G = ones(gdims(feedback, sys.A, sys.B, sys.C)...)
  return sys
end



"""
A Linear Time Invariant System defined in the space state.
The input u[k] is delayed "d" samples.

  x[k+1] = Ax[k] + Adu[k-d]
  y[k] = Cx[k]

Where u[k] depends on the feedback indication
  SOF -> u[k] = Gy[k] = GCx[k]
  SSF -> u[k] = Gx[k]
  OL  -> There is no feedback, it's a open-loop system 

"""
mutable struct SystemDelayed <: System
  A::Matrix{Float64}
  Ad::Matrix{Float64}
  C::Matrix{Float64}
  d::Int
  feedback::Feedback
  G::Matrix{Float64}
end

function SystemDelayed(A, Ad, C, d, feedback::Feedback)
  checkdims(A, Ad, C)
  SystemDelayed(A, Ad, C, d, feedback, ones(gdims(feedback, A, Ad, C)...))
end

SystemDelayed(A, Ad, C, d) = SystemDelayed(A, Ad, C, d, OL)

"""
    set_feedback!(sys::SystemDelayed, feedback::Feedback)

Updates the feedback of `sys` and resizes `G` to the matching dimensions.
"""
function set_feedback!(sys::SystemDelayed, feedback::Feedback)
  sys.feedback = feedback
  sys.G = ones(gdims(feedback, sys.A, sys.Ad, sys.C)...)
  return sys
end

"""
A Linear Time Invariant System defined in the space state, with a disturbance signal w(k).

  x(k+1) = Ax(k) + Bu(k) + Ew(k)
  y(k) = Cx(k)

Where u[k] depends on the feedback indication
  SOF -> u[k] = Gy[k] = GCx[k]
  SSF -> u[k] = Gx[k]
  OL  -> There is no feedback, it's a open-loop system 
"""

mutable struct LTIDisturbance <: System
  A::Matrix
  B::Matrix
  E::Matrix
  C::Matrix
  feedback::Feedback
  G::Matrix
end

function LTIDisturbance(A, B, E, C, feedback::Feedback)
  checkdims(A, B, E, C)
  LTIDisturbance(A, B, E, C, feedback, ones(gdims(feedback, A, B, C)...))
end

LTIDisturbance(A, B, E, C) = LTIDisturbance(A, B, E, C, OL)

"""
    set_feedback!(sys::LTIDisturbance, feedback::Feedback)

Updates the feedback of `sys` and resizes `G` to the matching dimensions.
Note that the previous value of `G` will be lost.
"""
function set_feedback!(sys::LTIDisturbance, feedback::Feedback)
  sys.feedback = feedback
  sys.G = ones(gdims(feedback, sys.A, sys.B, sys.C)...)
  return sys
end

"""
A Linear Time Invariant System defined in the space state,
with a disturbance signal w(k) and delayed input.

  x(k+1) = Ax(k) + Bu(k-d) + Ew(k)
  y(k) = Cx(k)

Where u[k] depends on the feedback indication
  SOF -> u[k] = Gy[k] = GCx[k]
  SSF -> u[k] = Gx[k]
  OL  -> There is no feedback, it's a open-loop system 
"""

mutable struct LTIDisturbanceDelayedInput <: System
  A::Matrix
  B::Matrix
  E::Matrix
  C::Matrix
  d::Int
  feedback::Feedback
  G::Matrix
end

function LTIDisturbanceDelayedInput(A, B, E, C, d, feedback::Feedback)
  checkdims(A, B, E, C)
  LTIDisturbanceDelayedInput(A, B, E, C, d, feedback, ones(gdims(feedback, A, B, C)...))
end

LTIDisturbanceDelayedInput(A, B, E, C, d) = LTIDisturbanceDelayedInput(A, B, E, C, d, OL)

"""
    set_feedback!(sys::LTIDisturbanceDelayedInput, feedback::Feedback)

Updates the feedback of `sys` and resizes `G` to the matching dimensions.
Note that the previous value of `G` will be lost.
"""
function set_feedback!(sys::LTIDisturbanceDelayedInput, feedback::Feedback)
  sys.feedback = feedback
  sys.G = ones(gdims(feedback, sys.A, sys.B, sys.C)...)
  return sys
end




