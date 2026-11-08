module Polyhedron

using Reexport

@reexport using LinearAlgebra
@reexport using DelimitedFiles
@reexport using JuMP
@reexport using HiGHS
@reexport using NEOSServer

@reexport using Polyhedra
@reexport using CDDLib
using Base.Iterators
@reexport using Plots
@reexport using MathOptInterface

include("types/Vertices.jl")
include("types/System.jl")
include("types/InvariantPolyhedron.jl")

include("utils/cond_iniciais_adm.jl")
include("utils/vet_eq_spc.jl")
include("utils/outros.jl")
include("trajectories/trajectory.jl")
include("utils/poly_projection.jl")


include("invariance/is_pinvariant/is_pinvariant.jl")

include("plots/plot_poly.jl")
include("plots/plot_trajectory.jl")

include("invariance/state_constraints.jl")
include("invariance/input_constraints.jl")
include("invariance/invariance_constraints.jl")

include("invariance/finding_pinvariant/eq_spaced_vet_enl.jl")
include("invariance/finding_pinvariant/state_constraint_enl.jl")
include("invariance/finding_pinvariant/finding_L_pinvariant.jl")
include("trajectories/trajectory_segref.jl")

include("pipelines/pipe_trajectory.jl")

include("utils/enl_simples.jl")

export pinvariant_poly
export invariance_constraints!, state_constraints!, input_constraints!
export eq_spaced_vet_enl!, state_constraint_enl!
export finding_L_pinvariant, finding_L_pinvariant_segref, finding_L_pinvariant_segref_delay, step1_saturation, step2_saturation, is_pinvariant_seg_ref
export step1_saturation_segref, step2_saturation_segref, step1_saturation_segref_alt, step2_saturation_segref_alt
export is_pinvariant
export plot_poly, get_shape, plot_vertices, plot_trajectories, plot_expanded_state_trajectory, plot_v_state_trajectory, plot_along_k
export plot_trajectory, plot_trajectory!
export Vertices, get_vertices, get_extreme_vertices, poly_projection, get_extVert_tuple
export SLIT, SystemDelayed, LTIDisturbance, LTIDisturbanceDelayedInput, Feedback, SOF, SSF, OL, set_feedback!, checkdims
export InvariantPolyhedron
export cond_iniciais_adm, mat_cond_iniciais_adm, elim_red, extended_F, extended_A, extended_A_Vector, allPossibleComb, admissable_initCond, calcular_v
export calcular_v

export pipe_trajectory

# falta outros.jl
export trajectory_delay, trajectory, vet_eq_spc, trajectory_segref, trajectory_delay_sat, trajectory_segref_delay

export enl_simples

end
