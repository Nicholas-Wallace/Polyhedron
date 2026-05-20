"""
    trajectory(x0, A, passos) -> Vector{Tuple{Float64}}

    Gera e retorna a trajetoria dos estados de um sistema do tipo:

        x[k+1] = A*x[k]

    `x0` é o vetor de condição inicial e deve estar como Vector{Tuple{Float64}}
"""

function trajectory(x0, A, passos)
    for i in range(2, passos)
        try
            x0 = hcat(x0, A*x0[:, i-1])    
        catch e
            print(e)
        end
        
    end
    return x0   
end