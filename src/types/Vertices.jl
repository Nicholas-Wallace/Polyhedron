function get_extVert_tuple(init_cond_F, init_cond_w, n)
    hrep_ic = hrep(init_cond_F, init_cond_w)
    #h = removehredundancy(hrep_ic, HiGHS.Optimizer)
    p = polyhedron(hrep_ic, CDDLib.Library())
    
    #v = removevredundancy(vrep(p), custom_highs)
    v = vrep(p)

    # seria bom saber de antemao o tamanho do vetor vertices
    vertices = Vector{Vector{Vector{Float64}}}()
    for pt in points(v)
        chunks = partition(collect(pt), n)
        push!(vertices, collect.(chunks)) # Adds each chunk as a Vector{Float64}
    end

    return vertices
end

function get_extreme_vertices(A, b, num_pontos=5)
    n_vars = size(A, 2)
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    @variable(model, x[1:n_vars])
    @constraint(model, A * x .<= b)
    
    lista_vertices = Vector{Vector{Float64}}()
    
    for i in 1:num_pontos
        c = randn(n_vars)
        @objective(model, Max, sum(c[j] * x[j] for j in 1:n_vars))
        
        optimize!(model)
        
        if termination_status(model) == MOI.OPTIMAL
            push!(lista_vertices, value.(x))
        end
    end
    return lista_vertices
end
