

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

function get_extVert_tuple(init_cond_F, init_cond_w, i, n)
    hrep_ic = hrep(init_cond, init_cond_w)
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