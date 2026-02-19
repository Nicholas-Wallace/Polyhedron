
# cria um conjunto de vetores igualmente espaçados
# função auxiliar para o enl(.) dos problemas de otimização bilineares
function vet_eq_spc(n)
    angle = 2*pi/n

    #matriz de rotação
    m_rot = [cos(angle) -sin(angle);
             sin(angle) cos(angle)]
    
    vet_array = Array{Vector{Float64}}(undef, n)
    
    # por padrão iniciasse com um vetor unitário na direção do eixo y, e vai rotacionando ele  
    vet_array[1] = [0, 1]

    for i in range(2, n)
        vet_array[i] = m_rot*vet_array[i-1] 
    end
    return vet_array

end