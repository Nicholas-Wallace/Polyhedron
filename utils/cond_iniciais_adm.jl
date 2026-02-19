function cond_iniciais_adm(ext_F, ext_A, d)
    cond_iniciais_matrix = Matrix{Matrix{Matrix{Float64}}}(undef, d+1, 1)
    cond_iniciais_matrix[1] = ext_F
    for i in range(2, d+1)
        cond_iniciais_matrix[i] = ext_F*ext_A^(i-1)
    end

    

    return cond_iniciais_matrix
end

function mat_cond_iniciais_adm(A, Ad, F, d; symetric=true)
  ext_A = Polyhedron.extended_A(A, Ad, d)
  ext_F = Polyhedron.extended_F(F, d)
  cond_iniciais_adm = Polyhedron.cond_iniciais_adm(ext_F, ext_A, d)
  m, n = size(cond_iniciais_adm)#mxn
  p, q = size(cond_iniciais_adm[1])#pxq
  r, s = size(cond_iniciais_adm[1][1])#rxs

  cia_final = cond_iniciais_adm[1, 1][1, 1]
  temp = cond_iniciais_adm[1, 1][1, 1]

  for i in 1:m
    for j in 1:n
      for k in 1:p
        for l in 1:q
          if l == 1
            temp = cond_iniciais_adm[i, j][k, l]
          else
            temp = hcat(temp, cond_iniciais_adm[i, j][k, l])
          end
        end
        if i == j == k == 1
          cia_final = temp
        else
          cia_final = vcat(cia_final, temp)
        end
      end
    end
  end
  if symetric
    return vcat(cia_final, -cia_final)
  else
    return cia_final
  end
  
end

function elimred(G, ro)
    h = hrep(G, ro)
    p = polyhedron(h, CDDLib.Library())
    removehredundancy!(p)
    h_clean = hrep(p)

    n_h = nhalfspaces(h_clean)
    dim = size(G, 2)

    Gn = Matrix{Float64}(undef, n_h, dim)
    ron = Vector{Float64}(undef, n_h)

    for (i, half_space) in enumerate(halfspaces(h_clean))
        Gn[i, :] = half_space.a
        ron[i] = half_space.β
    end

    return Gn, ron
end

# Cria lista de matrizes F e forma a matriz diagonal extended_F
function extended_F(F, d)
    blocks = [F for _ in 1:(d+1)]
    # Concatena na diagonal (concatena tanto verticalmente quanto horizontalmente)
    # Isso gera a matriz diagonal, os zeros são colocados automaticamente
    return cat(blocks..., dims=(1, 2))
end

function extended_A(A, Ad, d; dm=d)
    n = size(A, 1)
    N = n * (dm + 1)
    id = (d-1) * n + 1 : d * n

    row1 = zeros(n, N-n)
    row1[:, id] = Ad
    row1 = hcat(A, row1)

    identity = Float64.(I(n*dm))
    bottom_zero = zeros(n*dm, n)

    return vcat(row1, hcat(identity, bottom_zero))
end

function extended_A_Vector(A, Ad, dm)
    n = size(A, 1)
    N = n * (dm + 1)
    x = [zeros(n, N) for _ in 1:dm]
    
    for i in 1:dm
        x[i] = extended_A(A, Ad, i; dm=dm)
    end

    return x
end

function allPossibleComb(j, dm)
    # Cria lista de ranges 1:dm
    ranges = [1:dm for _ in 1:j]

    # Retorna o produto cartesiano desses ranges numa matriz
    return collect(Iterators.product(ranges...))
end

function admissable_initCond(A, Ad, F, dm, w; symetric=false, fixed_d=false)
    ext_F = extended_F(F, dm)
    ext_w = repeat(w, dm + 1)

    init_cond_F = ext_F
    init_cond_w = ext_w 

    n = size(A, 1)
    N = n * (dm + 1)

    if fixed_d == false
        A_array = extended_A_Vector(A, Ad, dm)
        
        for i in 1:dm
            indexes = allPossibleComb(i, dm)
            for index in indexes
                product = Matrix{Float64}(I, N, N)
                for elem in index
                    product = A_array[elem] * product
                end
                init_cond_F = vcat(init_cond_F, ext_F * product)
                init_cond_w = vcat(init_cond_w, ext_w)
            end

            init_cond_F, init_cond_w = elimred(init_cond_F, init_cond_w)
        end
    else 
        ext_A = extended_A(A, Ad, dm)
        current_A_pow = ext_F

        for i in 1:dm
            current_A_pow *= ext_A
            init_cond_F = vcat(init_cond_F, current_A_pow)
            init_cond_w = vcat(init_cond_w, ext_w)
            
            init_cond_F, init_cond_w = elimred(init_cond_F, init_cond_w)
        end
    end
        
    if symetric
        return vcat(init_cond_F, -init_cond_F), vcat(init_cond_w, init_cond_w)
    end

    return init_cond_F, init_cond_w
end
