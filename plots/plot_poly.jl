
# inq_matrix <= w
function plot_poly(inq_matrix, w)
    hrep_L = hrep(inq_matrix, w)
    P = vrep(polyhedron(hrep_L, CDDLib.Library()))
    iter = points(P)
    vertices = []
    for vet in iter
    push!(vertices, Tuple(vet))
    end
    poly_F = Shape(vertices)
    plt = plot(poly_F, fc = plot_color(:blue, 0.1), label=false, title="poly")
    xlabel!("x1")
    ylabel!("x2")
    return plt

end
