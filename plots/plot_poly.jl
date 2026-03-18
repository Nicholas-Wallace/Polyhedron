
# inq_matrix <= w
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