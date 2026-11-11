module Polyhedron

using Reexport

# Public API reexports
@reexport using JuMP
@reexport using Plots

# Internal dependencies
using LinearAlgebra
using Base.Iterators: partition
using MathOptInterface
import HiGHS
import NEOSServer
import CDDLib
import Polyhedra
using Polyhedra: hrep, vrep, polyhedron, points, rays,
                 halfspaces, nhalfspaces, removehredundancy!

# --- types ------------------------------------------------------------------
include("types/Vertices.jl")
include("types/System.jl")
include("types/InvariantPolyhedron.jl")

# --- utilities --------------------------------------------------------------
include("utils/cond_iniciais_adm.jl")
include("utils/vet_eq_spc.jl")
include("utils/poly_projection.jl")
# include("utils/outros.jl")  # disabled: undefined `Poly`, leaks globals F/G

# --- trajectories -----------------------------------------------------------
include("trajectories/trajectory.jl")

# --- invariance -------------------------------------------------------------
include("invariance/state_constraints.jl")
include("invariance/input_constraints.jl")
include("invariance/invariance_constraints.jl")
include("invariance/is_pinvariant/is_pinvariant.jl")

# --- finding invariant polyhedra -------------------------------------------
include("invariance/finding_pinvariant/eq_spaced_vet_enl.jl")
include("invariance/finding_pinvariant/state_constraint_enl.jl")
include("invariance/finding_pinvariant/pinvariant_poly.jl")

# --- plots ------------------------------------------------------------------
include("plots/plot_poly.jl")
include("plots/plot_trajectory.jl")

# --- public API -------------------------------------------------------------
export System, Feedback, SOF, SSF, OL
export SLIT, SystemDelayed, LTIDisturbance, LTIDisturbanceDelayedInput
export set_feedback!, checkdims
export InvariantPolyhedron
export get_extVert_tuple, get_extreme_vertices
export poly_projection
export elimred, extended_F, extended_A, extended_A_Vector, allPossibleComb,
       admissable_initCond, calcular_v
export vet_eq_spc
export trajectory, trajectory_segref
export is_pinvariant
export invariance_constraints!, state_constraints!, input_constraints!
export pinvariant_poly
export eq_spaced_vet_enl!, state_constraint_enl!
export plot_poly, get_shape, plot_vertices, plot_trajectories,
       plot_expanded_state_trajectory, plot_v_state_trajectory, plot_along_k
export plot_trajectory, plot_trajectory!

end
