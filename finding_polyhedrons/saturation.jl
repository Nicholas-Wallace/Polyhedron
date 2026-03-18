function step1_saturation(A, B, X, f; lambda=0.99, pond=0.5, symetric=true)
    model = Model() do
        return NEOSServer.Optimizer(; email = "wallace.lopes.162@ufrn.edu.br", solver = "Knitro")
    end

    # --- Opções do Solver (Knitro) ---
    set_optimizer_attribute(model, "nlp_algorithm", 1)
    set_optimizer_attribute(model, "feastol_abs", 1e-8)
    set_optimizer_attribute(model, "opttol_abs", 1e-8)
    set_optimizer_attribute(model, "scale", 1)
    set_optimizer_attribute(model, "xtol_iters", 500)
    set_optimizer_attribute(model, "ms_maxsolves", 1000)
    set_optimizer_attribute(model, "outlev", 2)
    set_optimizer_attribute(model, "numthreads", 4)
    set_optimizer_attribute(model, "honorbnds", 1)
    set_optimizer_attribute(model, "ms_enable", 1)
    set_optimizer_attribute(model, "ms_seed", 0)

    # parâmetros
    n = size(A, 1) # ordem do sistema 
    m = size(B, 2) # entradas 
    x = size(X, 1) # restrições de estado 
    w = ones(f)
    q = ones(x)

    @variable(model, 1 <= d <= 100) 
    @variable(model, 0.001 <= igama <= 1000)

    @variable(model, -100 <= G[1:m, 1:n] <= 100)
    @variable(model, -1000 <= K[1:n, 1:n] <= 1000)
    @variable(model, -50 <= F[1:f, 1:n] <= 50)
    @variable(model, J[1:n, 1:f])
    
    @variable(model, FK[1:f, 1:n])
    @variable(model, BG[1:n, 1:n])
    @variable(model, KA[1:n, 1:n])
    @variable(model, KG[1:n, 1:n])

    @constraint(model, FK .== F*K)
    @constraint(model, BG .== B*G)
    @constraint(model, KA .== K*(A - Matrix(I, n, n)))
    @constraint(model, KG .== K*BG)
    @constraint(model, J*F .== Matrix(I, n, n))

    @objective(model, Min, pond * igama + (1 - pond) / d)

    # --- Inicialização Global (Warm Start) ---
    set_start_value(d, 10.0)
    set_start_value(igama, 10.0)
    set_start_value.(G, 0.1)
    
    for i in 1:n, j in 1:n
        set_start_value(K[i,j], i == j ? 0.1 : 0.0)
    end

    for i in 1:f, j in 1:n
        val = i <= x ? X[i,j] : (i % 2 == 0 ? 0.01 : -0.01)
        set_start_value(F[i,j], val)
    end

    set_start_value.(J, 0.0)

    # Pré-cálculo para inicializar as auxiliares e ajudar o Jacobiano inicial
    K_init = [i == j ? 0.1 : 0.0 for i in 1:n, j in 1:n]
    F_init = [i <= x ? X[i,j] : (i % 2 == 0 ? 0.01 : -0.01) for i in 1:f, j in 1:n]
    BG_init = B * fill(0.1, m, n)
    
    set_start_value.(FK, F_init * K_init)
    set_start_value.(BG, BG_init)
    set_start_value.(KA, K_init * (A - Matrix(I, n, n)))
    set_start_value.(KG, K_init * BG_init)

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
        
        @variable(model, 0 <= jota1[1:f, 1:x] <= 100)
        @variable(model, 0 <= jota2[1:f, 1:x] <= 100)

        @variable(model, dM1[1:f, 1:f] >= 0)
        @variable(model, dM2[1:f, 1:f] >= 0)
        @variable(model, dN1[1:f, 1:f] >= 0)
        @variable(model, dN2[1:f, 1:f] >= 0)

        # --- Inicialização Bloco Simétrico ---
        for i in 1:f, j in 1:f
            set_start_value(H1[i,j], i == j ? 0.1 : 0.0)
            set_start_value(L1[i,j], i == j ? 0.1 : 0.0)
            set_start_value(M1[i,j], i == j ? 0.01 : 0.0)
            set_start_value(N1[i,j], i == j ? 0.01 : 0.0)
            
            set_start_value(dM1[i,j], i == j ? 10.0 * 0.01 : 0.0)
            set_start_value(dN1[i,j], i == j ? 10.0 * 0.01 : 0.0)

            # As parcelas negativas iniciam zeradas
            set_start_value(H2[i,j], 0.0)
            set_start_value(L2[i,j], 0.0)
            set_start_value(M2[i,j], 0.0)
            set_start_value(N2[i,j], 0.0)
            set_start_value(dM2[i,j], 0.0)
            set_start_value(dN2[i,j], 0.0)
        end

        set_start_value.(R1, 0.01)
        set_start_value.(R2, 0.0)
        set_start_value.(jota1, 0.1)
        set_start_value.(jota2, 0.0)

        @constraint(model, dM1 .== d * M1)
        @constraint(model, dM2 .== d * M2)
        @constraint(model, dN1 .== d * N1)
        @constraint(model, dN2 .== d * N2)

        @constraint(model, (H1-H2)*F .== F*A + FK)
        @constraint(model, (L1-L2)*F .== F*(BG) - FK)
        @constraint(model, (M1-M2)*F .== -F*KA)
        @constraint(model, (N1-N2)*F .== -F*KG)
        
        @constraint(model, ((H1 + H2) + (L1 + L2) + (dM1 + dM2) + (dN1 + dN2))*w .<= lambda * w)
        
        @constraint(model, (R1-R2)*F .== X)
        @constraint(model, (R1+R2)*w .<= q)

        @constraint(model, (jota1-jota2)*X .== F)
        @constraint(model, (jota1+jota2)*q .<= igama * w)

        optimize!(model)
        
        d_val = value(d)
        igama_val = value(igama)
        F_val = value.(F)
        K_val = value.(K)
        H_val = value.(H1) - value.(H2)
        L_val = value.(L1) - value.(L2)
        M_val = value.(M1) - value.(M2)
        N_val = value.(N1) - value.(N2)
        R_val = value.(R1) - value.(R2)
        jota_val = value.(jota1) - value.(jota2)
        G_val = value.(G)

        print(termination_status(model))

        result = Dict("d" => d_val, "igama" => igama_val, "F" => F_val, "G" => G_val, "K" => K_val, "H" => H_val, "L" => L_val, "M" => M_val, "N" => N_val, "R" => R_val, "jota" => jota_val)

        return result
    end
    
    @variable(model, 0 <= H[1:f, 1:f] <= 200)
    @variable(model, 0 <= L[1:f, 1:f] <= 200)
    @variable(model, 0 <= M[1:f, 1:f] <= 200)
    @variable(model, 0 <= N[1:f, 1:f] <= 200)
    @variable(model, 0 <= R[1:x, 1:f] <= 200)
    @variable(model, 0 <= jota[1:f, 1:x] <= 200)

    @variable(model, dM[1:f, 1:f] >= 0)
    @variable(model, dN[1:f, 1:f] >= 0)

    # --- Inicialização Bloco Não-Simétrico ---
    for i in 1:f, j in 1:f
        set_start_value(H[i,j], i == j ? 0.1 : 0.0)
        set_start_value(L[i,j], i == j ? 0.1 : 0.0)
        set_start_value(M[i,j], i == j ? 0.01 : 0.0)
        set_start_value(N[i,j], i == j ? 0.01 : 0.0)
        
        set_start_value(dM[i,j], i == j ? 10.0 * 0.01 : 0.0)
        set_start_value(dN[i,j], i == j ? 10.0 * 0.01 : 0.0)
    end
    set_start_value.(R, 0.01)
    set_start_value.(jota, 0.1)

    @constraint(model, dM .== d * M)
    @constraint(model, dN .== d * N)

    @constraint(model, H*F .== F*A + FK)
    @constraint(model, L*F .== F*(BG) - FK)
    @constraint(model, M*F .== -F*KA)
    @constraint(model, N*F .== -F*KG)
    
    @constraint(model, (H + L + dM + dN)*w .<= lambda * w)
    
    @constraint(model, R*F .== X)
    @constraint(model, R*w .<= q)

    @constraint(model, jota*X .== F)
    @constraint(model, jota*q .<= igama * w)

    optimize!(model)

    d_val = value(d)
    igama_val = value(igama)
    G_val = value.(G)
    F_val = value.(F)
    K_val = value.(K)
    H_val = value.(H)
    L_val = value.(L)
    M_val = value.(M)
    N_val = value.(N)
    R_val = value.(R)
    jota_val = value.(jota)

    print(termination_status(model))

    result = Dict("d" => d_val, "igama" => igama_val, "F" => F_val, "G" => G_val, "K" => K_val, "H" => H_val, "L" => L_val, "M" => M_val, "N" => N_val, "R" => R_val, "jota" => jota_val)

    return result
end

# FALTA IMPLEMENTAR MUDANÇAS FEITAS NO PASSO 1 E MOVER O IF PRA CIMA DEPOIS DE DESCOBRIR 
# O QUE TEM DE ERRADO NO 1
function step2_saturation(A, B, X, G, umax, umin, d; f = 4, lambda=0.99, symetric=true)
    model = Model() do
        return NEOSServer.Optimizer(; email = "wallace.lopes.162@ufrn.edu.br", solver = "Knitro")
    end

    set_optimizer_attribute(model, "nlp_algorithm", 1)
    set_optimizer_attribute(model, "feastol_abs", 1e-8)
    set_optimizer_attribute(model, "opttol_abs", 1e-8)
    set_optimizer_attribute(model, "scale", 1)
    set_optimizer_attribute(model, "xtol_iters", 500)
    set_optimizer_attribute(model, "ms_maxsolves", 10000)
    set_optimizer_attribute(model, "outlev", 1)
    set_optimizer_attribute(model, "honorbnds", 1)
    set_optimizer_attribute(model, "ms_enable", 1)
    set_optimizer_attribute(model, "ms_seed", 0)

    # parâmetros
    n = size(A, 1) # ordem do sistema
    x = size(X, 1)
    m = size(B, 2)
    w = ones(f)
    q = ones(x)

    @variable(model, -100 <= K[1:n, 1:n] <= 100)
    @variable(model, -100 <= F[1:f, 1:n] <= 100)
    @variable(model, -100 <= P[1:m, 1:n] <= 100)

    @variable(model, 0.001 <= igama <= 1000)
    
    @variable(model, J[1:n, 1:f])
    
    @variable(model, 0 <= jota[1:f, 1:x] <= 100)

    @variable(model, 0 <= Q1[1:m, 1:f] <= 100)
    @variable(model, 0 <= Q2[1:m, 1:f] <= 100)
    @variable(model, 0 <= R[1:x, 1:f] <= 100)

    if symetric
        @variable(model, 0 <= Hp[1:f, 1:f] <= 100);  @variable(model, 0 <= Hm[1:f, 1:f] <= 100)
        @variable(model, 0 <= L1p[1:f, 1:f] <= 100); @variable(model, 0 <= L1m[1:f, 1:f] <= 100)
        @variable(model, 0 <= L2p[1:f, 1:f] <= 100); @variable(model, 0 <= L2m[1:f, 1:f] <= 100)
        @variable(model, 0 <= Mp[1:f, 1:f] <= 100);  @variable(model, 0 <= Mm[1:f, 1:f] <= 100)
        @variable(model, 0 <= N1p[1:f, 1:f] <= 100); @variable(model, 0 <= N1m[1:f, 1:f] <= 100)
        @variable(model, 0 <= N2p[1:f, 1:f] <= 100); @variable(model, 0 <= N2m[1:f, 1:f] <= 100)

        H = Hp - Hm; L1 = L1p - L1m; L2 = L2p - L2m; 
        M = Mp - Mm; N1 = N1p - N1m; N2 = N2p - N2m
    else
        @variable(model, 0 <= H[1:f, 1:f] <= 100)
        @variable(model, 0 <= L1[1:f, 1:f] <= 100)
        @variable(model, 0 <= L2[1:f, 1:f] <= 100)
        @variable(model, 0 <= M[1:f, 1:f] <= 100)
        @variable(model, 0 <= N1[1:f, 1:f] <= 100)
        @variable(model, 0 <= N2[1:f, 1:f] <= 100)
    end

    @variable(model, -100 <= BG[1:n, 1:n] <= 100)
    @variable(model, -100 <= BP[1:n, 1:n] <= 100)
    @variable(model, -100 <= KA[1:n, 1:n] <= 100)
    @variable(model, -100 <= KG[1:n, 1:n] <= 100)
    @variable(model, -100 <= KP[1:n, 1:n] <= 100)
    @variable(model, -100 <= FK[1:f, 1:n] <= 100)

    set_start_value(igama, 10.0)

    K_init = Matrix(0.1 * I(n))
    P_init = fill(0.1, m, n)
    F_init = zeros(f, n)

    for i in 1:m, j in 1:n
        set_start_value(P[i,j], P_init[i,j])
    end

    for i in 1:n, j in 1:n
        set_start_value(K[i,j], K_init[i,j])
    end

    for i in 1:n, j in 1:f
        set_start_value(J[i,j], 0.0)
    end

    for i in 1:f, j in 1:n
        if i <= x
            F_init[i,j] = X[i,j]
            set_start_value(F[i,j], F_init[i,j])
        else
            val = i % 2 == 0 ? 0.01 : -0.01
            F_init[i,j] = val
            set_start_value(F[i,j], val)
        end
    end

    for i in 1:x, j in 1:f
        set_start_value(R[i,j], 0.01) 
        set_start_value(jota[j,i], 0.1)
    end

    for i in 1:f, j in 1:f
        diag = (i == j)
        if symetric
            set_start_value(Hp[i,j], diag ? 0.1 : 0.0);  set_start_value(Hm[i,j], 0.0)
            set_start_value(L1p[i,j], diag ? 0.1 : 0.0); set_start_value(L1m[i,j], 0.0)
            set_start_value(L2p[i,j], diag ? 0.1 : 0.0); set_start_value(L2m[i,j], 0.0)
            set_start_value(Mp[i,j], diag ? 0.01 : 0.0); set_start_value(Mm[i,j], 0.0)
            set_start_value(N1p[i,j], diag ? 0.01 : 0.0); set_start_value(N1m[i,j], 0.0)
            set_start_value(N2p[i,j], diag ? 0.01 : 0.0); set_start_value(N2m[i,j], 0.0)
        else
            set_start_value(H[i,j], diag ? 0.1 : 0.0);   set_start_value(L1[i,j], diag ? 0.1 : 0.0)
            set_start_value(L2[i,j], diag ? 0.1 : 0.0);  set_start_value(M[i,j], diag ? 0.01 : 0.0)
            set_start_value(N1[i,j], diag ? 0.01 : 0.0); set_start_value(N2[i,j], diag ? 0.01 : 0.0)
        end
    end

    for i in 1:m, j in 1:f
        set_start_value(Q1[i,j], umax[i]/2)
        set_start_value(Q2[i,j], -umin[i]/2)
    end

    for i in 1:n, j in 1:n
        set_start_value(BG[i,j], BG_init[i,j]); set_start_value(BP[i,j], BP_init[i,j])
        set_start_value(KA[i,j], KA_init[i,j]); set_start_value(KG[i,j], KG_init[i,j])
        set_start_value(KP[i,j], KP_init[i,j])
    end
    for i in 1:f, j in 1:n
        set_start_value(FK[i,j], FK_init[i,j])
    end

    BG_init = B * G; BP_init = B * P_init
    KA_init = K_init * (A - I); KG_init = K_init * BG_init
    KP_init = K_init * BP_init; FK_init = F_init * K_init

    @constraint(model, H*F == F*(A + K))
    @constraint(model, L1*F == F*(B*G - K))
    @constraint(model, L2*F == F*(B*P - K))
    @constraint(model, M*F == -FK*(A - I(n)))
    @constraint(model, N1*F == -FK*(B*G))
    @constraint(model, N2*F == -FK*(B*P))

    @constraint(model, (H + L1 + d*(M + N1))*w .<= lambda * w)
    @constraint(model, (H + L1 + d*(M + N2))*w .<= lambda * w)
    @constraint(model, (H + L2 + d*(M + N1))*w .<= lambda * w)
    @constraint(model, (H + L2 + d*(M + N2))*w .<= lambda * w)

    @constraint(model, Q1*F == P)
    @constraint(model, Q2*F == -P)

    @constraint(model, Q1*w .<= umax)
    @constraint(model, Q2*w .<= umin)

    @constraint(model, R*F == X)
    @constraint(model, R*w .<= q)

    @objective(model, Max, sum(gamma)/v)
    
    optimize!(model)

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