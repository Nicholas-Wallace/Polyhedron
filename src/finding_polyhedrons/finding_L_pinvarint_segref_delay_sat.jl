function step1_saturation_segref(A, B, E, S, R, V; f=12, lambda=0.99, pond=0.5, symetric=true)
    model = Model() do
        return NEOSServer.Optimizer(; email = "wallace.lopes.162@ufrn.edu.br", solver = "Knitro")
    end

    # --- Opções do Solver ---
    set_optimizer_attribute(model, "nlp_algorithm", 1)
    set_optimizer_attribute(model, "feastol_abs", 1e-8)
    set_optimizer_attribute(model, "opttol_abs", 1e-8)
    set_optimizer_attribute(model, "outlev", 2)
    set_optimizer_attribute(model, "ms_enable", 1)
    set_optimizer_attribute(model, "ms_maxsolves", 100)

    # Parâmetros
    n = size(A, 1)       # Ordem do sistema
    m = size(B, 2)       # Entradas de controle
    s = size(S, 1)       # Linhas de restrição de estado
    r = size(R, 1)       # Linhas de restrição de referência
    v = size(V, 1)       # Linhas de restrição de controle

    ones_f = ones(f)
    ones_s = ones(s)
    ones_r = ones(r)
    ones_v = ones(v)

    # Variáveis Principais
    @variable(model, 1 <= d <= 100) 
    @variable(model, 0.001 <= igamma <= 1000)

    @variable(model, -100 <= G[1:m, 1:n] <= 100)
    @variable(model, -1000 <= K[1:n, 1:n] <= 1000)
    @variable(model, -50 <= F[1:f, 1:n] <= 50)
    @variable(model, J[1:n, 1:f])

    if symetric
        @variable(model, 0 <= H1[1:f, 1:f] <= 100); @variable(model, 0 <= H2[1:f, 1:f] <= 100)
        @variable(model, 0 <= L1[1:f, 1:f] <= 100); @variable(model, 0 <= L2[1:f, 1:f] <= 100)
        @variable(model, 0 <= M1[1:f, 1:f] <= 100); @variable(model, 0 <= M2[1:f, 1:f] <= 100)
        @variable(model, 0 <= N1[1:f, 1:f] <= 100); @variable(model, 0 <= N2[1:f, 1:f] <= 100)
        @variable(model, 0 <= R1[1:s, 1:f] <= 100); @variable(model, 0 <= R2[1:s, 1:f] <= 100)
        @variable(model, 0 <= jota1[1:f, 1:s] <= 100); @variable(model, 0 <= jota2[1:f, 1:s] <= 100)
        
        @variable(model, 0 <= U1[1:v, 1:f] <= 100)
        @variable(model, 0 <= U2[1:v, 1:f] <= 100)

        @variable(model, 0 <= Pref1[1:f, 1:r] <= 100)
        @variable(model, 0 <= Pref2[1:f, 1:r] <= 100)

        @variable(model, dM1[1:f, 1:f] >= 0); @variable(model, dM2[1:f, 1:f] >= 0)
        @variable(model, dN1[1:f, 1:f] >= 0); @variable(model, dN2[1:f, 1:f] >= 0)

        H = H1 - H2; L = L1 - L2; M = M1 - M2; N = N1 - N2
        Pref = Pref1 - Pref2; dM = dM1 - dM2; dN = dN1 - dN2
    else
        @variable(model, 0 <= H[1:f, 1:f] <= 200);   @variable(model, 0 <= L[1:f, 1:f] <= 200)
        @variable(model, 0 <= M[1:f, 1:f] <= 200);   @variable(model, 0 <= N[1:f, 1:f] <= 200)
        @variable(model, 0 <= R[1:s, 1:f] <= 200);   @variable(model, 0 <= jota[1:f, 1:s] <= 200)
        @variable(model, 0 <= U[1:v, 1:f] <= 200)
        @variable(model, dM[1:f, 1:f] >= 0);         @variable(model, dN[1:f, 1:f] >= 0)
        @variable(model, 0 <= Pref[1:f, 1:r] <= 200)
    end

    @variable(model, FK[1:f, 1:n])
    @variable(model, BG[1:n, 1:n])
    @variable(model, KA[1:n, 1:n])
    @variable(model, KG[1:n, 1:n])

    # Restrições de Definição Bilinear
    @constraint(model, FK .== F*K)
    @constraint(model, BG .== B*G)
    @constraint(model, KA .== K*(A - Matrix(I, n, n)))
    @constraint(model, KG .== K*BG)
    @constraint(model, J*F .== Matrix(I, n, n))

    # Restrições do Lema de Farkas (Igualdades)
    @constraint(model, dM .== d * M)
    @constraint(model, dN .== d * N)
    @constraint(model, H*F .== F*A + FK)
    @constraint(model, L*F .== F*(BG) - FK)
    @constraint(model, M*F .== -F*KA)
    @constraint(model, N*F .== -F*KG)
    @constraint(model, Pref * R .== F * E)

    # Restrições de Inclusão e Desigualdades Principais
    if symetric
        @constraint(model, ((H1 + H2) + (L1 + L2) + (dM1 + dM2) + (dN1 + dN2))*ones_f + (Pref1+Pref2)*ones_r .<= lambda * ones_f)
        
        @constraint(model, (R1-R2)*F .== S)
        @constraint(model, (R1+R2)*ones_f .<= ones_s)
        
        @constraint(model, (U1-U2)*F .== V*G)
        @constraint(model, (U1+U2)*ones_f .<= ones_v)

        @constraint(model, (jota1-jota2)*S .== F)
        @constraint(model, (jota1+jota2)*ones_s .<= igamma * ones_f)
    else
        @constraint(model, (H + L + dM + dN)*ones_f + Pref*ones_r .<= lambda * ones_f)
        
        @constraint(model, R*F .== S)
        @constraint(model, R*ones_f .<= ones_s)

        @constraint(model, U*F .== V*G)
        @constraint(model, U*ones_f .<= ones_v)

        @constraint(model, jota*S .== F)
        @constraint(model, jota*ones_s .<= igamma * ones_f)
    end

    @objective(model, Min, pond * igamma + (1 - pond) / d)

    optimize!(model)

    return Dict("d" => value(d), "igamma" => value(igamma), "F" => value.(F), "G" => value.(G), "Pref" => value.(Pref))
end

function step2_saturation_segref(A, B, E, S, R, V, G, umax, umin, d; f=12, lambda=0.99, symetric=true)
    model = Model() do
        return NEOSServer.Optimizer(; email = "wallace.lopes.162@ufrn.edu.br", solver = "Knitro")
    end

    set_optimizer_attribute(model, "nlp_algorithm", 1)
    set_optimizer_attribute(model, "feastol_abs", 1e-8)
    set_optimizer_attribute(model, "opttol_abs", 1e-8)
    set_optimizer_attribute(model, "outlev", 1)
    set_optimizer_attribute(model, "ms_enable", 1)

    # Parâmetros
    n = size(A, 1)
    s = size(S, 1)
    m = size(B, 2)
    r = size(R, 1)
    v = size(V, 1)
    
    ones_f = ones(f)
    ones_s = ones(s)
    ones_r = ones(r)
    ones_v = ones(v)

    @variable(model, -100 <= K[1:n, 1:n] <= 100)
    @variable(model, -100 <= F[1:f, 1:n] <= 100)
    @variable(model, -100 <= Psat[1:m, 1:n] <= 100) # Rebatizada para evitar colisão com a referência

    @variable(model, 0.001 <= igamma <= 1000)
    @variable(model, J[1:n, 1:f])

    @variable(model, 0 <= Q1[1:m, 1:f] <= 100)
    @variable(model, 0 <= Q2[1:m, 1:f] <= 100)

    if symetric
        @variable(model, 0 <= Hp[1:f, 1:f] <= 100);  @variable(model, 0 <= Hm[1:f, 1:f] <= 100)
        @variable(model, 0 <= L1p[1:f, 1:f] <= 100); @variable(model, 0 <= L1m[1:f, 1:f] <= 100)
        @variable(model, 0 <= L2p[1:f, 1:f] <= 100); @variable(model, 0 <= L2m[1:f, 1:f] <= 100)
        @variable(model, 0 <= Mp[1:f, 1:f] <= 100);  @variable(model, 0 <= Mm[1:f, 1:f] <= 100)
        @variable(model, 0 <= N1p[1:f, 1:f] <= 100); @variable(model, 0 <= N1m[1:f, 1:f] <= 100)
        @variable(model, 0 <= N2p[1:f, 1:f] <= 100); @variable(model, 0 <= N2m[1:f, 1:f] <= 100)
        @variable(model, 0 <= T1[1:s, 1:f] <= 100);  @variable(model, 0 <= T2[1:s, 1:f] <= 100)
        @variable(model, 0 <= Pref1[1:f, 1:r] <= 100); @variable(model, 0 <= Pref2[1:f, 1:r] <= 100)
        
        @variable(model, 0 <= U1[1:v, 1:f] <= 100)
        @variable(model, 0 <= U2[1:v, 1:f] <= 100)

        @variable(model, 0 <= jota1[1:f, 1:s] <= 100); @variable(model, 0 <= jota2[1:f, 1:s] <= 100)

        H = Hp - Hm; L1 = L1p - L1m; L2 = L2p - L2m 
        M = Mp - Mm; N1 = N1p - N1m; N2 = N2p - N2m
        Pref = Pref1 - Pref2; T = T1 - T2; U = U1 - U2
        jota = jota1 - jota2
    else
        @variable(model, 0 <= H[1:f, 1:f] <= 200)
        @variable(model, 0 <= L1[1:f, 1:f] <= 200)
        @variable(model, 0 <= L2[1:f, 1:f] <= 200)
        @variable(model, 0 <= M[1:f, 1:f] <= 200)
        @variable(model, 0 <= N1[1:f, 1:f] <= 200)
        @variable(model, 0 <= N2[1:f, 1:f] <= 200)
        @variable(model, 0 <= T[1:s, 1:f] <= 200)

        @variable(model, 0 <= Pref[1:f, 1:r] <= 200)
        @variable(model, 0 <= U[1:v, 1:f] <= 200)
        @variable(model, 0 <= jota[1:f, 1:s] <= 200)
    end

    @variable(model, -100 <= BG[1:n, 1:n] <= 100)
    @variable(model, -100 <= BPsat[1:n, 1:n] <= 100)
    @variable(model, -100 <= KA[1:n, 1:n] <= 100)
    @variable(model, -100 <= KG[1:n, 1:n] <= 100)
    @variable(model, -100 <= KPsat[1:n, 1:n] <= 100)
    @variable(model, -100 <= FK[1:f, 1:n] <= 100)

    @constraint(model, BG .== B*G)
    @constraint(model, BPsat .== B*Psat)
    @constraint(model, KA .== K*(A - I))
    @constraint(model, KG .== K*BG)
    @constraint(model, KPsat .== K*BPsat)
    @constraint(model, FK .== F*K)
    @constraint(model, Pref * R .== F * E) # Inclusão da referência

    @constraint(model, H*F .== F*A + FK)
    @constraint(model, L1*F .== F*BG - FK)
    @constraint(model, L2*F .== F*BPsat - FK)
    @constraint(model, M*F .== -F*KA)
    @constraint(model, N1*F .== -F*KG)
    @constraint(model, N2*F .== -F*KPsat)

    @constraint(model, J * F .== I(n))

    # Saturação de Entrada
    @constraint(model, Q1*F .== Psat)
    @constraint(model, Q2*F .== -Psat)
    @constraint(model, Q1*ones_f .<= umax)
    @constraint(model, Q2*ones_f .<= -umin)

    # Restrições de Estado S e Controle V
    @constraint(model, T*F .== S)

    # Restrição no controle
    @constraint(model, U*F .== V*G)

    # Restrição do igamma
    @constraint(model, jota*S .== F)
    
    if symetric
        @constraint(model, ((Hp + Hm) + (L1p + L1m) + d*((Mp + Mm) + (N1p + N1m)))*ones_f + (Pref1+Pref2)*ones_r .<= lambda * ones_f)
        @constraint(model, ((Hp + Hm) + (L1p + L1m) + d*((Mp + Mm) + (N2p + N2m)))*ones_f + (Pref1+Pref2)*ones_r .<= lambda * ones_f)
        @constraint(model, ((Hp + Hm) + (L2p + L2m) + d*((Mp + Mm) + (N1p + N1m)))*ones_f + (Pref1+Pref2)*ones_r .<= lambda * ones_f)
        @constraint(model, ((Hp + Hm) + (L2p + L2m) + d*((Mp + Mm) + (N2p + N2m)))*ones_f + (Pref1+Pref2)*ones_r .<= lambda * ones_f)

        @constraint(model, (U1+U2)*ones_f .<= ones_v)
        @constraint(model, (T1+T2)*ones_f .<= ones_s)
        @constraint(model, (jota1+jota2)*ones_s .<= igamma*ones_f)
    else
        @constraint(model, (H + L1 + d*(M + N1))*ones_f + Pref*ones_r .<= lambda * ones_f)
        @constraint(model, (H + L1 + d*(M + N2))*ones_f + Pref*ones_r .<= lambda * ones_f)
        @constraint(model, (H + L2 + d*(M + N1))*ones_f + Pref*ones_r .<= lambda * ones_f)
        @constraint(model, (H + L2 + d*(M + N2))*ones_f + Pref*ones_r .<= lambda * ones_f)

        @constraint(model, U*ones_f .<= ones_v)
        @constraint(model, T*ones_f .<= ones_s)
        @constraint(model, jota*ones_s .<= igamma*ones_f)
    end

    @objective(model, Min, igamma)
    optimize!(model)

    println(termination_status(model))

    return Dict("F" => value.(F), "K" => value.(K), "Psat" => value.(Psat), "G" => G, "Pref" => value.(Pref), "igamma" => value(igamma))    
end