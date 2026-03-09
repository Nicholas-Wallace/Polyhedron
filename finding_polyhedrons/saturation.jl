# Passo 1 da saturação
function step1_saturation(A, B, X, f; lambda=0.99, symetric=true)
    model = Model() do
        return NEOSServer.Optimizer(; email = "wallace.lopes.162@ufrn.edu.br", solver = "Knitro")
    end

    # parâmetros

    n = size(A, 1) # ordem do sistema
    x = size(X, 1)
    w = ones(f)
    q = ones(x)

    @variable(model, 0 <= d <= 100)

    @variable(model, G[1:1, 1:n])
    @variable(model, K[1:n, 1:n])
    @variable(model, F[1:f, 1:n])
    @variable(model, Y[1:f, 1:n])
    @variable(model, J[1:n, 1:f])

    @constraint(model, Y == F*K)
    @constraint(model, J*F == I(n))

    @objective(model, Max, d)

    if symetric
        @variable(model, 0 <= H1[1:f, 1:f] <= 100)
        @variable(model, 0 <= H2[1:f, 1:f] <= 100)
        @variable(model, 0 <= L1[1:f, 1:f] <= 100)
        @variable(model, 0 <= L2[1:f, 1:f] <= 100)
        @variable(model, 0 <= M1[1:f, 1:f] <= 100)
        @variable(model, 0 <= M2[1:f, 1:f] <= 100)
        @variable(model, 0 <= N1[1:f, 1:f] <= 100)
        @variable(model, 0 <= N2[1:f, 1:f] <= 100)
        @variable(model, 0 <= R1[1:x, 1:f] <= 100)
        @variable(model, 0 <= R2[1:x, 1:f] <= 100)

        @constraint(model, (H1-H2)*F == F*(A + K))
        @constraint(model, (L1-L2)*F == F*(B*G - K))
        @constraint(model, (M1-M2)*F == -Y*(A - I(n)))
        @constraint(model, (N1-N2)*F == -Y*(B*G))
        @constraint(model, ((H1 + H2) + (L1 + L2) + d*((M1 + M2) + (N1 + N2)))*w .<= lambda * w)
        @constraint(model, (R1-R2)*F == X)
        @constraint(model, (R1+R2)*w .<= q)

        optimize!(model)
        
        d = value(d)
        F = value.(F)
        K = value.(K)
        H = value.(H1) - value.(H2)
        L = value.(L1) - value.(L2)
        M = value.(M1) - value.(M2)
        N = value.(N1) - value.(N2)
        R = value.(R1) - value.(R2)
        G = value.(G)

        print(termination_status(model))

        # retorna um dicionário com d e as matrizes
        result = Dict("d" => d, "F" => F, "G" => G, "K" => K, "H" => H, "L" => L, "M" => M, "N" => N, "R" => R)

        return result
    end
    
    @variable(model, 0 <= H[1:f, 1:f] <= 200)
    @variable(model, 0 <= L[1:f, 1:f] <= 200)
    @variable(model, 0 <= M[1:f, 1:f] <= 200)
    @variable(model, 0 <= N[1:f, 1:f] <= 200)
    @variable(model, 0 <= R[1:x, 1:f] <= 200)

    @constraint(model, H*F == F*(A + K))
    @constraint(model, L*F == F*(B*G - K))
    @constraint(model, M*F == -Y*(A - I(n)))
    @constraint(model, N*F == -Y*(B*G))
    @constraint(model, (H + L + d*(M + N))*w .<= lambda * w)
    @constraint(model, R*F == X)
    @constraint(model, R*w .<= q)

    optimize!(model)

    d = value(d)
    G = value.(G)
    F = value.(F)
    K = value.(K)
    H = value.(H)
    L = value.(L)
    M = value.(M)
    N = value.(N)
    R = value.(R)

    print(termination_status(model))

    result = Dict("d" => d, "F" => F, "G" => G, "K" => K, "H" => H, "L" => L, "M" => M, "N" => N, "R" => R)

    return result
end

# FALTA IMPLEMENTAR MUDANÇAS FEITAS NO PASSO 1 E MOVER O IF PRA CIMA DEPOIS DE DESCOBRIR 
# O QUE TEM DE ERRADO NO 1
function step2_saturation(A, B, X, G, umax, umin, d; f = 4, v=8, lambda=0.99, symetric=true)
    model = Model() do
        return NEOSServer.Optimizer(; email = "wallace.lopes.162@ufrn.edu.br", solver = "Knitro")
    end

    # parâmetros

    n = size(A, 1) # ordem do sistema
    x = size(X, 1)
    w = ones(f)
    q = ones(x)

    #vet = vet_eq_spc(v)
    #xt = ones(v)

    @variable(model, K[1:n, 1:n])
    @variable(model, F[1:f, 1:n])
    @variable(model, P[1:1, 1:n])
    @variable(model, igamma >= 1)
    

    # TESTAR SEM E DEPOIS COM ISSO
    @variable(model, J[1:n, 1:f])
    @constraint(model, J*F == I(n))

    #for i in 1:v
     #   @constraint(model, F * (gamma[i] * vet[i]) .<= w)
    #end

    @variable(model, 0 <= H[1:f, 1:f] <= 100)
    @variable(model, 0 <= L1[1:f, 1:f] <= 100)
    @variable(model, 0 <= L2[1:f, 1:f] <= 100)
    @variable(model, 0 <= M[1:f, 1:f] <= 100)
    @variable(model, 0 <= N1[1:f, 1:f] <= 100)
    @variable(model, 0 <= N2[1:f, 1:f] <= 100)
    @variable(model, 0 <= R[1:x, 1:f] <= 100)
    @variable(model, 0 <= Q1[1:1, 1:f] <= 100)
    @variable(model, 0 <= Q2[1:1, 1:f] <= 100)

    # talvez ajude na otimização criar uma variavel auxiliar FK
    @variable(model, FK[1:f, 1:n])
    @constraint(model, FK == F*K)

    @variable(model, BG[1:n, 1:n])
    @constraint(model, BG == B*G)


    @variable(model, BP[1:n, 1:n])
    @constraint(model, BP == B*P)

    @variable(model, KA[1:n, 1:n])
    @constraint(model, KA == K*(A- I(n)))

    @variable(model, KG[1:n, 1:n])
    @constraint(model, KG == K*BG)

    @variable(model, KP[1:n, 1:n])
    @constraint(model, KP == K*BP)

    @constraint(model, H*F == F*(A + K))
    @constraint(model, L1*F == F*(BG - K))
    @constraint(model, L2*F == F*(BP - K))
    @constraint(model, M*F == -F*KA)
    @constraint(model, N1*F == -F*KG)
    @constraint(model, N2*F == -F*KP)

    @constraint(model, (H + L1 + d*(M + N1))*w .<= lambda * w)
    @constraint(model, (H + L1 + d*(M + N2))*w .<= lambda * w)
    @constraint(model, (H + L2 + d*(M + N1))*w .<= lambda * w)
    @constraint(model, (H + L2 + d*(M + N2))*w .<= lambda * w)

    @constraint(model, Q1*F == P)
    @constraint(model, Q2*F == -P)

    @constraint(model, Q1*w .<= umax)
    @constraint(model, Q2*w .<= umin)

    @constraint(model, R*F == X)
    @constraint(model, R*w .<= q*igamma)

    pond = 0.00000001;

    func_obj = igamma

    @objective(model, Min, func_obj)
    
    optimize!(model)

    igamma = value.(igamma)
    F = value.(F)
    G = value.(G)
    K = value.(K)
    H = value.(H)
    L1 = value.(L1)
    L2 = value.(L2)

    print(termination_status(model))

    result = Dict("F" => F, "K" => K, "H" => H, "L1" => L1, "L2" => L2, "G" => G)

    return result


    # if symetric
    #     @variable(model, 0 <= Hp[1:f, 1:f])
    #     @variable(model, 0 <= Hm[1:f, 1:f])
    #     @variable(model, 0 <= L1p[1:f, 1:f])
    #     @variable(model, 0 <= L1m[1:f, 1:f])
    #     @variable(model, 0 <= L2p[1:f, 1:f])
    #     @variable(model, 0 <= L2m[1:f, 1:f])
    #     @variable(model, 0 <= Mp[1:f, 1:f])
    #     @variable(model, 0 <= Mm[1:f, 1:f])
    #     @variable(model, 0 <= N1p[1:f, 1:f])
    #     @variable(model, 0 <= N1m[1:f, 1:f])
    #     @variable(model, 0 <= N2p[1:f, 1:f])
    #     @variable(model, 0 <= N2m[1:f, 1:f])
    #     @variable(model, 0 <= Rp[1:x, 1:f])
    #     @variable(model, 0 <= Rm[1:x, 1:f])
    #     @variable(model, 0 <= Qp[1:1, 1:f])
    #     @variable(model, 0 <= Qm[1:1, 1:f])
    
    #     @constraint(model, (Hp - Hm)*F == F*(A + K))
    #     @constraint(model, (L1p - L1m)*F == F*(B*G - K))
    #     @constraint(model, (L2p - L2m)*F == F*(B*P - K))
    #     @constraint(model, (Mp - Mm)*F == -Y*(A - I(n)))
    #     @constraint(model, (N1p - N1m)*F == -Y*(B*G))
    #     @constraint(model, (N2p - N2m)*F == -Y*(B*P))
    #     @constraint(model, ((Hp + Hm) + (L1p + L1m) + d*((Mp + Mm) + (N1p + N1m)))*w .<= lambda * w)
    #     @constraint(model, ((Hp + Hm) + (L1p + L1m) + d*((Mp + Mm) + (N2p + N2m)))*w .<= lambda * w)
    #     @constraint(model, ((Hp + Hm) + (L2p + L2m) + d*((Mp + Mm) + (N1p + N1m)))*w .<= lambda * w)
    #     @constraint(model, ((Hp + Hm) + (L2p + L2m) + d*((Mp + Mm) + (N2p + N2m)))*w .<= lambda * w)

    #     @constraint(model, (Rp - Rm)*F == X)
    #     @constraint(model, (Rp + Rm)*w .<= q)

    #     @constraint(model, (Qp - Qm)*F == P)
    #     @constraint(model, (Qp + Qm)*w .<= umax)

    #     optimize!(model)
        
    #     F = value(F)
    #     K = value.(K)
    #     H = value.(Hp) - value.(Hm)
    #     L1 = value.(L1p) - value.(L1m)
    #     L2 = value.(L2p) - value.(L2m)
    #     print(termination_status(model))

    #     # retorna um dicionário com matrizes F, K, H e L
    #     result = Dict("F" => F, "K" => K, "H" => H, "L1" => L1, "L2" => L2, "G" => G)

    #     return result
    # end
end