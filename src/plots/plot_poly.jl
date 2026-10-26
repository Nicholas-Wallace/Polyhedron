"""
    plot_poly(ip::InvariantPolyhedron) -> Plot

plota um poliedro p-invariante definido por `ip`. O poliedro deve estar em 2D.
"""
function plot_poly(ip::InvariantPolyhedron)
    size(ip.F, 2) == 2 ||
        error("plot_poly requires a 2D InvariantPolyhedron, got $(size(ip.F, 2)) columns")

    poly_F = get_shape(ip)

    # Plota a figura
    plt = plot(poly_F, fc=plot_color(:blue, 0.1), label=false, title="Polyhedron")
    xlabel!("x1")
    ylabel!("x2")

    return plt
end

"""
    plot_poly(ip1::InvariantPolyhedron, ip2::InvariantPolyhedron; label1, label2, legendpos) -> Plot

plota dois poliedros p-invariantes `ip1` e `ip2` em 2D, com um rótulo para cada um.
A posição da legenda pode ser definida por `legendpos` (padrão: `:bottomright`).
"""
function plot_poly(ip1::InvariantPolyhedron, ip2::InvariantPolyhedron;
                   label1="Poliedro Invariante", label2="Poliedro Auxiliar", title="Polyhedrons", legendpos=:bottomright)
    size(ip1.F, 2) == 2 && size(ip2.F, 2) == 2 ||
        error("plot_poly requires 2D InvariantPolyhedrons, got $(size(ip1.F, 2)) and $(size(ip2.F, 2)) columns")

    poly_F1 = get_shape(ip1)
    poly_F2 = get_shape(ip2)

    plt = plot(poly_F1, fc=plot_color(:blue, 0.1), label=label1, title=title,
               legend=legendpos, aspect_ratio=:equal)
    plot!(plt, poly_F2, fc=plot_color(:red, 0.1), label=label2)

    return plt
end

function get_shape(ip::InvariantPolyhedron)
    vertices = [Tuple(vet[1:2]) for vet in points(ip.vrep)]

    # Acha o valor do centróide do poliedro
    cx = sum(v[1] for v in vertices) / length(vertices)
    cy = sum(v[2] for v in vertices) / length(vertices)

    # Ordena com base no ângulo do vêrtice em relação ao centróide
    sort!(vertices, by = v -> atan(v[2] - cy, v[1] - cx))

    return Shape(vertices)
end

"""
    plot_poly(inq_matrix, w) -> Plot

plota um poliedro definido pela inequação

    inq_matrix <= w
"""
function plot_poly(inq_matrix, w)
    poly_F = get_shape(inq_matrix, w)
    
    # Plota a figura
    plt = plot(poly_F, fc=plot_color(:blue, 0.1), label=false, title="Polyhedron")
    xlabel!("x1")
    ylabel!("x2")
    
    return plt
end

function get_shape(inq_matrix, w)
    hrep_F = hrep(inq_matrix, w)
    P = vrep(polyhedron(hrep_F, CDDLib.Library()))

    vertices = [Tuple(vet[1:2]) for vet in points(P)]

    # Acha o valor do centróide do poliedro
    cx = sum(v[1] for v in vertices) / length(vertices)
    cy = sum(v[2] for v in vertices) / length(vertices)

    # Ordena com base no ângulo do vêrtice em relação ao centróide
    sort!(vertices, by = v -> atan(v[2] - cy, v[1] - cx))

    return Shape(vertices)
end

function plot_vertices(F, vertices)

    T = Polyhedron.poly_projection(F[:,3:4])
    plt = Polyhedron.plot_poly(T*F[:,1:2], T*ones(size(F, 1)))

    x1 = [v[1] for v in vertices]
    x2 = [v[2] for v in vertices]
    plot!(x1, x2, seriestype = :scatter, label = "Vertices", color = :red)
end

function plot_expanded_state_trajectory(trajs)
    plots_list = []
    
    for (i, traj) in enumerate(trajs)
        # Generate an individual subplot for the current trajectory
        p = plot([[p[1] for p in traj] [p[2] for p in traj] [p[4] for p in traj]],
                 label=["x1" "x2" "w"], 
                 xlabel="k",
                 ylims=(-15, 15),
                 st = :steppost,  
                 title="Traj $i")
        
        push!(plots_list, p)
    end
    
    # The `layout` parameter controls the grid shape (e.g., length(trajs) rows and 1 column)
    final_plot = plot(plots_list..., layout = (length(trajs), 1), size = (600, 300 * length(trajs)))
    
    return final_plot
end

function plot_along_k(trajs)
    plots_list = []
    
    for (i, traj) in enumerate(trajs)
        # Generate an individual subplot for the current trajectory
        p = plot([p for p in traj],
                 label="v", 
                 xlabel="k",
                 st = :steppost,  
                 title="Traj $i")
        
        push!(plots_list, p)
    end
    
    # The `layout` parameter controls the grid shape (e.g., length(trajs) rows and 1 column)
    final_plot = plot(plots_list..., layout = (length(trajs), 1), size = (600, 300 * length(trajs)))
    
    return final_plot
end

function plot_v_state_trajectory(trajs)
    plots_list = []
    
    for (i, traj) in enumerate(trajs)
        # Generate an individual subplot for the current trajectory
        p = plot([p[3] for p in traj],
                 label="v", 
                 xlabel="k",
                 st = :steppost,  
                 title="Traj $i")
        
        push!(plots_list, p)
    end
    
    # The `layout` parameter controls the grid shape (e.g., length(trajs) rows and 1 column)
    final_plot = plot(plots_list..., layout = (length(trajs), 1), size = (600, 300 * length(trajs)))
    
    return final_plot
end


function plot_trajectories(F, trajs)
    T = Polyhedron.poly_projection(F[:,3:4])
    plt = Polyhedron.plot_poly(T*F[:,1:2], T*ones(size(F, 1)))

    for (i, traj) in enumerate(trajs)
        x1 = [p[1] for p in traj]
        x2 = [p[2] for p in traj]
        plot!(x1, x2,
                 seriestype = :path,
                 linewidth = 2,
                 marker=:circle,
                 label = "Trajectory $i")
    end

    return plt
end
