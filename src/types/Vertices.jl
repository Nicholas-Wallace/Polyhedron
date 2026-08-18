

# a cddlib retorna os vertices em um vetor, mas para as operações de plotar e calcular trajetórias 
# é melhor lidar com os vertices sendo um vetor de tuplas, em que cada tupla é uma cordenada no espaço de estados

struct Vertices{T <: Tuple} 
    points::Vector{T}
end

Base.getindex(v::Vertices, i::Int) = v.points[i]
Base.length(v::Vertices) = length(v.points)
Base.iterate(v::Vertices, state=1) = iterate(v.points, state)

# deve ser um construtor de Vertices que recebe um vetor com 
# vertices e reparte em chunks de tamanho n
function Vertices(vet ,n)
    chunks = Iterators.partition(vet, n)
    
    return [Tuple(chunk) for chunk in chunks if length(chunk) == n]
end

# ESSA FUNÇÃO RETORNA UMA TUPLA COM OS ESTADOS AGRUPADOS
# DE ALGUM VÉRTICE DO POLIEDRO DE CONDIÇÕES INICIAIS ADMISSÍVEIS
# "n" É A ORDEM DO SISTEMA E "i" É O ÍNDICE DO VÉRTICE A SER ESCOLHIDO

function calcular_v(F, E, R)
    # Pré-computa a matriz M = F * E
    M = F * E
    
    n_linhas = size(M, 1)
    n_cols_R = size(R, 2)  # dimensão do vetor r
    
    ones_r = ones(size(R, 1))
    
    v = zeros(n_linhas)
    
    # Cria UM único modelo e resolve repetidamente (mais eficiente que criar N modelos)
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    
    @variable(model, r[1:n_cols_R])
    @constraint(model, restricoes, R * r .<= ones_r)
    
    for i in 1:n_linhas
        # Define o objetivo como a i-ésima linha de M
        @objective(model, Max, sum(M[i, j] * r[j] for j in 1:n_cols_R))
        
        optimize!(model)
        
        if termination_status(model) == OPTIMAL
            v[i] = objective_value(model)
        else
            error("LP infeasible or unbounded at row $i")
        end
    end
    
    return v
end

function get_extVert_tuple(init_cond_F, init_cond_w, n)
    hrep_ic = hrep(init_cond_F, init_cond_w)
    #h = removehredundancy(hrep_ic, HiGHS.Optimizer)
    p = polyhedron(hrep_ic, CDDLib.Library(:exact))
    
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

function get_extVert_tuple_old(init_cond_F, init_cond_w, i, n)
    hrep_ic = hrep(init_cond_F, init_cond_w)
    P = vrep(polyhedron(hrep_ic, CDDLib.Library()))

    ext_vertices = collect(points(P))
    vertice_alvo = ext_vertices[i]

    ext_vertices_tuple = Tuple[]

    for i in 1:n:length(vertice_alvo)
        estado = vertice_alvo[i : i + n - 1]
        push!(ext_vertices_tuple, Tuple(estado))
    end

    return ext_vertices_tuple
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
