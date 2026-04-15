module Poly

using LinearAlgebra, DelimitedFiles, JuMP, HiGHS, NEOSServer, Polyhedra, CDDLib, Base.Iterators, Plots

include("types/Vertices.jl")

include("utils/cond_iniciais_adm.jl")
include("utils/vet_eq_spc.jl")
include("utils/outros.jl")
include("utils/trajectory.jl")
include("utils/trajectory_delay.jl")
include("utils/poly_projection.jl")
include("utils/calculate_admissible_references.jl")


include("methods/is_pinvariant.jl")
include("methods/is_pinvariant_delay.jl")


include("plots/plot_poly.jl")

include("finding_polyhedrons/finding_L_pinvariant.jl")
include("finding_polyhedrons/finding_L_pinvariant_delay.jl")
include("finding_polyhedrons/saturation.jl")
include("finding_polyhedrons/finding_L_pinvariant_segref.jl")
include("finding_polyhedrons/finding_L_pinvariant_segref_delay.jl")
include("finding_polyhedrons/ref_tracking_delay.jl")

include("trajectories/trajectory_segref.jl")


export finding_L_pinvariant, finding_L_pinvariant_delay, finding_L_pinvariant_segref, finding_L_pinvariant_segref_delay, finding_L_pinvariant_segref_delay_sim, step1_saturation, step2_saturation, is_pinvariant_seg_ref, finding_L_pinvariant_segref_delay2	
export is_pinvariant, is_pinvariant_delay
export plot_poly, get_shape
export Vertices, get_vertices, get_extreme_vertices, poly_projection, get_extVert_tuple
export cond_iniciais_adm, mat_cond_iniciais_adm, elim_red, extended_F, extended_A, extended_A_Vector, allPossibleComb, admissable_initCond, calculate_admissible_references
# falta outros.jl
export trajectory_delay, trajectory, vet_eq_spc, trajectory_segref, trajectory_delay_sat, trajectory_segref_delay


end