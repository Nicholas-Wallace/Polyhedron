function calver(G, w)
    vertex = Set{Tuple{Float64, Float64}}() 
    
    n_l = size(G, 1)
    
    # Loop through every unique pair of lines
    for i in 1 : n_l - 1
        for j in i + 1 : n_l # Start from i+1 to avoid self-intersection
            row1 = G[i, :]
            row2 = G[j, :]
            A = vcat(row1', row2')
            B = [w[i], w[j]]

            try
                solution = A \ B

                if satisfact_all(G, solution)
                    # Round the solution to a few decimal places to handle minor variations
                    # of the same point, then push to the set.
                    rounded_solution = Tuple(round.(solution, digits=8))
                    push!(vertex, rounded_solution)
                end

            catch e
                # This is expected for parallel lines, so no need to print unless debugging
                print(e)

            end
        end
    end

    # Convert the Set of vertices back to an Array
    P = Poly(G, reorderpoints(vertex))
    return P
end

# otimizar essa funcao (pense em pdroduto matricial)

function satisfact_all(G, x)

    tolerance = 1e-4 

    for v in eachrow(G)
        if v'*x > 1 + tolerance 
            return false
        end
    end
    return true
end

function reorderpoints(points)
    vet_in_order = Vector{Tuple{Float64, Float64}}() 
    # push!(vet_in_order, pop!(points))

    # print(vet_in_order)

    # finding the center of the Polyhedron

    pts = collect(points)

    cx = sum(p[1] for p in pts)/length(pts)
    cy = sum(p[2] for p in pts)/length(pts)

    vet_in_order = sort(pts, by = p -> atan(p[2] - cy, p[1] - cx)) # função anonima
    return vet_in_order
end 

