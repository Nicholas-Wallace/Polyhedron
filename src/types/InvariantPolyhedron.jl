"""
    InvariantPolyhedron{T}

Representa um poliedro p-invariante w.r.t. `system::System`:

    F x <= λ w

onde `x[k+1]` respeita a dinâmica do sistema.

Além dos dados algébricos (`F`, `w`, `λ`, `system`), guarda a representação
`Polyhedra`/`CDDLib` do poliedro, com sua H-representação (`hrep`) e
V-representação (`vrep`).
"""
struct InvariantPolyhedron{T<:Real}
  F::Matrix{T}
  w::Vector{T}
  λ::T
  system::System
  symetric::Bool
  poly::Polyhedra.Polyhedron{T}
  hrep::Polyhedra.HRepresentation{T}
  vrep::Polyhedra.VRepresentation{T}
end

"""
    InvariantPolyhedron(F, w, λ, system) -> InvariantPolyhedron

Constrói o poliedro p-invariante definido por `F x <= w`, com fator de
contração `λ`, computando a representação via `CDDLib` (H e V).
"""
function InvariantPolyhedron(F::AbstractMatrix, w::AbstractVector, λ::Real, system::System; symetric::Bool=false)
  T = Float64
  
  if symetric
    F = vcat(F, -F)
    w = vcat(w, w)
  end

  h = hrep(Matrix{T}(F), Vector{T}(w))
  p = polyhedron(h, CDDLib.Library())
  v = vrep(p)
  hr = hrep(p)

  return InvariantPolyhedron{T}(Matrix{T}(F), Vector{T}(w), T(λ), system, symetric, p, hr, v)
end
