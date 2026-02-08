#ANTES
hrep_L = hrep(cia_final, ones(242))
P = vrep(polyhedron(hrep_L, CDDLib.Library()))
vertices = collect(points(P))

vertices_tuple = Tuple{Float64, Float64}[]

for i in 1:2:22
  push!(vertices_tuple, (vertices[2][i], vertices[2][i+1]))
end

t_10 = trajectory_delay(vertices_tuple, A, Ad, 100, 10)

h_rep_F = hrep(vcat(F, -F), ones(4))
iter = points(vrep(polyhedron(h_rep_F, CDDLib.Library())))
vertices_F = []
for p in iter
  push!(vertices_F, Tuple(p))
end


# PLOTA O GRAFICO
plt = plot(Shape(vertices_F), fc = plot_color(:blue, 0.1), title = "Trajetórias", label = "F")

# PLOTA AS TRAJ

plt = plot(Shape(vertices_F), fc = plot_color(:blue, 0.1), title = "Trajetórias", label = "F")


x = [p[1] for p in t_1]
y = [p[2] for p in t_1]

# Plot as connected trajectory with markers
plot!(x, y,
     seriestype = :path,          # connects points with lines
     linewidth = 2,
     marker = :circle,
     markersize = 5,
     markercolor = :red,
     label = "d = 1",
     xlabel = "x1",
     ylabel = "x2")

x = [p[1] for p in t_10]
y = [p[2] for p in t_10]

# Plot as connected trajectory with markers
plot!(x, y,
     seriestype = :path,          # connects points with lines
     linewidth = 2,
     marker = :diamond,
     markersize = 5,
     markercolor = :purple,
     label = "d = 10")
