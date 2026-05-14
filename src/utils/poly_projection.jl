function poly_projection(Gb)
    gb, u = size(Gb)

    # Construção das matrizes A1 e A2
    A1 = vcat(ones(1, gb), copy(transpose(Gb)))
    A2 = -Matrix{Float64}(I, gb, gb)

    # Construção dos vetores B1 e B2
    B1 = vcat(1.0, zeros(u))
    B2 = zeros(gb)

    # Concatenação das desigualdades: A1*x <= B1, -A1*x <= -B1, A2*x <= B2
    A = vcat(A1, -A1, A2)
    B = vcat(B1, -B1, B2)

    # Representação H e conversão para poliedro via CDD
    H = hrep(A, B)
    P = polyhedron(H, CDDLib.Library())

    # Extração das V-representações (pontos e raios extremos)
    pts = collect(points(P))
    rys = collect(rays(P))

    # Transforma os iteradores em matrizes onde cada linha é um vértice/raio
    # O uso do Vector() no raio extrai a direção numérica bruta do objeto Ray do Polyhedra.jl
    D = isempty(pts) ? zeros(0, gb) : vcat([v' for v in pts]...)
    R = isempty(rys) ? zeros(0, gb) : vcat([Vector(r)' for r in rys]...)

    # T = [D; R]
    T = vcat(D, R)

    return T
end