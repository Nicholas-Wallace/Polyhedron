module Polyhedron

using LinearAlgebra, DelimitedFiles, JuMP, HiGHS, NEOSServer, Polyhedra, CDDLib, Base.Iterators, Plots, MathOptInterface

include("types/Vertices.jl")
include("types/System.jl")
include("types/InvariantPolyhedron.jl")

include("utils/cond_iniciais_adm.jl")
include("utils/vet_eq_spc.jl")
include("utils/outros.jl")
include("trajectories/trajectory.jl")
include("utils/poly_projection.jl")


include("methods/is_pinvariant.jl")

include("plots/plot_poly.jl")
include("plots/plot_trajectory.jl")

include("finding_polyhedrons/finding_L_pinvariant.jl")
include("finding_polyhedrons/finding_L_pinvarint_segref_delay_sat.jl")
include("finding_polyhedrons/finding_L_pinvariant_segref.jl")
include("finding_polyhedrons/finding_L_pinvariant_segref_delay.jl")
include("trajectories/trajectory_segref.jl")

include("pipelines/pipe_trajectory.jl")

include("utils/enl_simples.jl")


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
