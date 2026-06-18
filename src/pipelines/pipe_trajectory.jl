
export pipe_trajectory

function pipe_trajectory(F::Matrix, A_exp::Matrix, B_exp::Matrix, E_exp::Matrix, G::Matrix, d::Int, Ref::Float64; passos=300)
    w = ones(size(F, 1))
    nx = size(A_exp, 1)

    T = Polyhedron.poly_projection(F[:,3:4])


    passos = 300
    r = ones(passos)*Ref # escolher a referência a se seguir 

    BG = B_exp * G

    init_cond_F, init_cond_w = Polyhedron.admissable_initCond(A_exp, BG, F, d, w; fixed_d=true)
    ext_vertices_tuple = Polyhedron.get_extVert_tuple(init_cond_F, init_cond_w, 2, nx)
    traj = Polyhedron.trajectory_segref_delay(ext_vertices_tuple, A_exp, BG, E_exp, r, passos, d; reverse=false)

    Polyhedron.plot_poly(T*F[:,1:2], T*ones(size(F, 1)))

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

