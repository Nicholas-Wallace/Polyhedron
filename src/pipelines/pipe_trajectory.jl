
function pipe_trajectory(F::Matrix{T}, A_exp::Matrix{T}, B_exp::Matrix{T}, G::Matrix{T}, d::Int, Ref::Float64; passos=300) where T
    w = ones(size(F, 1))
    nx = size(A_exp, 1)

    passos = 300
    r = ones(passos)*Ref # escolher a referência a se seguir 

    BG = B_exp * G

    init_cond_F, init_cond_w = Poly.admissable_initCond(A_exp, BG, F, d, w; fixed_d=true)
    ext_vertices_tuple = Poly.get_extVert_tuple(init_cond_F, init_cond_w, 2, nx)
    traj = Poly.trajectory_segref_delay(ext_vertices_tuple, A_exp, BG, E_exp, r, passos, d; reverse=false)

    Poly.plot_poly(T*F[:,1:2], T*ones(size(F, 1)))

    x1 = [p[1] for p in traj]
    x2 = [p[2] for p in traj]
    v = [p[3] for p in traj]
    w = [p[4] for p in traj]

    # Plot as connected trajectory with markers
    plot!(x1, x2,
        seriestype = :path,          # connects points with lines
        linewidth = 2,
        marker = :circle,
        markersize = 5,
        markercolor = :red,
        label = "d = 2",
        xlabel = "x1",
        ylabel = "x2")
end

