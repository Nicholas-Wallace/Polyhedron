module Polyhedron

using LinearAlgebra, DelimitedFiles, JuMP, NEOSServer


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
                # println("Rows $i and $j: singular matrix, skipping")
                print(e)

            end
        end
    end

    # Convert the Set of vertices back to an Array

    return reorderpoints(vertex)
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

# function satisfact_all(G, x)
#     tolerance = 1e-4 

#     return G*x > ones(size(G)[1]) + tolerance 
# end



function reorderpoints(points)
    vet_in_order = Vector{Tuple{Float64, Float64}}() 
    # push!(vet_in_order, pop!(points))

    # print(vet_in_order)

    # finding the center of the Polyhedron

    pts = collect(points)

    cx = sum(p[1] for p in pts)/length(pts)
    cy = sum(p[2] for p in pts)/length(pts)

    vet_in_order = sort(pts, by = p -> atan(p[2] - cy, p[1] - cx)) # função anonima
    
    # sum 

    # soma_x = 0
    # soma_y = 0
    # for p in points
    #   soma_x += p[1]
    #   soma_y += p[2]
    # end
  
    # soma_x /= length(points)
    # soma_y /= length(points)

    # c = (soma_x , soma_y)


    # for i in range(1, length(points) - 1)

    #     menor = 2*pi # upper bound
    #     next_pt = Tuple{Float64, Float64}
    #     for p in points
    #       a = [p[1] - c[1]; p[2] - c[2]]
    #       b = vet_in_order[i] #[vet_in_order[i][1], vet_in_order[i][2]]
    #       d = acos(dot(a, b)/(norm(a)*norm(b)))
          
    #       if d < menor
    #           menor = d
    #           next_pt = p
    #       end
    #     end

    #     pop!(points, next_pt)

    #     push!(vet_in_order, next_pt)

    # end

    # push!(vet_in_order, pop!(points))



    return vet_in_order
end 

#conjunto de vetores igualmente espaçados 
function vet_eq_spc(n)
    angle = 2*pi/n

    #matriz de rotação
    m_rot = [cos(angle) -sin(angle);
             sin(angle) cos(angle)]
    
    vet_array = Array{Vector{Float64}}(undef, n)
    
    vet_array[1] = [0, 1]

    for i in range(2, n)
        vet_array[i] = m_rot*vet_array[i-1] 
    end

    return vet_array

end

function trajectory(x0, G, passos)
    for i in range(2, passos)
        try
            x0 = hcat(x0, G*x0[:, i-1])    
        catch e
            print(e)
        end
        
    end

    return x0   
end

# The Linear Constrained Regulation Problem
# verifica se um certo poliedro X é invariante
function is_pinvariant(A, B, C, U, X, SOF = false) 
    
    model = Model() do
    return NEOSServer.Optimizer(; email = "wallace.lopes.162@ufrn.edu.br", solver = "Knitro")
    end

    #Parâmetros
    n = size(A, 1); #Ordem do sistema (Linhas de A)
    m = size(B, 2); #Número de Entradas (Colunas de B)
    p = size(C, 1); #Número de Saídas (Linhas de C)

    lx = size(X, 1); #Número de linhas de X
    lu = size(U, 1); #Número de linhas de U

    xi = ones(lx)
    phi = ones(lu);

    #Variáveis de decisão
    @variable(model, 0 <= H[1:lx,1:lx] <= 100)
    SOF ? @variable(model, K[1:m,1:p]) : @variable(model, F[1:m,1:n])
    @variable(model, 0 <= M[1:lu,1:lx] <= 100)
    @variable(model, 0 <= lambda <= 0.9999)

    #Objetivo

    func_obj = lambda

    @objective(model, Min, func_obj)

    # tentando encontrar um L, caso bilinear


    if SOF # realimentação de saida 
        @constraint(model, contractiave_pi0, H*X == X*(A + B*K*C) ) 
        @constraint(model, u_included0, M*X == U*K*C )


    else # realimentação de estado
        @constraint(model, contractiave_pi0, H*X == X*(A + B*F) ) 
        @constraint(model, u_included0, M*X == U*F )

    end
    
    @constraint(model, contractiave_pi1, H*xi <= xi*lambda )
    @constraint(model, u_included1, M*xi <= phi)

    #Variáveis encontradas pelo modelo
    H = value.(H)
    L = value.(L)
    F = value.(F)
    M = value.(M)
    lambda = value.(lambda)

    return lambda

end

function finding_L_pinvariant(A, B, C, U, X, SOF = false)
    
end

end
