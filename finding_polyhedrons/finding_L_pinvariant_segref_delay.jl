
#############################################
### essa função é referente a um sistema: ###
###     x[k+1] = Ax(k) BGx(k-d) + Er(k)   ### 
#############################################


###        Farkas para o poliedro ser invaiante         ###
# H*F == F(A+K)                                         ###
# L*F == F(BG - K)                                      ###
# M*F == -FK(A -I)                                      ###
# N*F == -FKBG                                          ###
# P*R == F*E                                            ###    
# (H + L + d(M+N))*ones_f + P*ones_r .<= lambda*ones_f  ###
#                                                       ###
###  Para estar contido em S (poliedro das restricoes)  ###
#                                                       ###
# T*F == S                                              ###
# T*ones_f .<= ones_f                                   ###
###########################################################

function finding_L_pinvariant_segref_delay(A, B, E, S, R, d; lambda=0.99, lf=6, symmetric=false) 
    
    model = Model() do
        return NEOSServer.Optimizer(; email = "wallace.lopes.162@ufrn.edu.br", solver = "Knitro")
    end

    #0  Auto (Knitro escolhe sozinho).
    #1: Active Set 
    #2: Interior Point / CG
    #3: Interior Point / Direct
    #4: SQP (Programação Quadrática Sequencial)

    set_optimizer_attribute(model, "nlp_algorithm", 1)
    set_optimizer_attribute(model, "feastol_abs", 1e-8)
    set_optimizer_attribute(model, "opttol_abs", 1e-8)
    set_optimizer_attribute(model, "scale", 1)
    set_optimizer_attribute(model, "xtol_iters", 500)
    set_optimizer_attribute(model, "ms_maxsolves", 30)
    set_optimizer_attribute(model, "outlev", 2)
    set_optimizer_attribute(model, "numthreads", 4)
    set_optimizer_attribute(model, "honorbnds", 1)
    set_optimizer_attribute(model, "ms_enable", 1)
    set_optimizer_attribute(model, "ms_seed", 0)

    #Parâmetros
    n = size(A, 1); #Ordem do sistema (Linhas de A)
    m = size(B, 2); #Número de Entradas (Colunas de B)
    r = size(R, 1); #Número de linhas das restrições em r[k]
    s = size(S, 1); #Número de linhas da matriz das restrições de X[k]

    ones_f = ones(lf)
    ones_r = ones(r)
    ones_s = ones(s)

    @variable(model, K[1:n, 1:n]) #Variavel auxiliar trans model
    @variable(model, G[1:m, 1:n]) #Ganho do controlador 
    @variable(model, -500 <= F[1:lf,1:n] <= 500) #Poliedro que queremos achar
    @variable(model, J[1:n, 1:lf]) #Pseudo inversa
    @variable(model, 0.001 <= igamma <= 1)

    # Auxiliary variables for bilinear products
    @variable(model, BG[1:n, 1:n])  # B*G
    @variable(model, FK[1:lf, 1:n])  # F*K
    @variable(model, KA[1:n, 1:n])  # K*(A-I)
    @variable(model, KBG[1:n, 1:n])  # K*B*G

    # Constraints for common auxiliary variables
    @constraint(model, BG .== B*G)
    @constraint(model, FK .== F*K)
    @constraint(model, KA .== K*(A - I(n)))
    @constraint(model, KBG .== K*BG)
    @constraint(model, J*F .== I(n))

    # --- Warm Start Initialization ---
    set_start_value(igamma, 10.0)
    set_start_value.(G, 0.1)

    for i in 1:n, j in 1:n
        set_start_value(K[i,j], i == j ? 0.1 : 0.0)
    end

    for i in 1:lf, j in 1:n
        set_start_value(F[i,j], i <= s ? S[i,j] : (i % 2 == 0 ? 0.01 : -0.01))
    end

    set_start_value.(J, 0.0)
    
    # Pre-computed initializations
    G_init = fill(0.1, m, n)
    K_init = [i == j ? 0.1 : 0.0 for i in 1:n, j in 1:n]
    F_init = [i <= s ? S[i,j] : (i % 2 == 0 ? 0.01 : -0.01) for i in 1:lf, j in 1:n]
    BG_init = B * G_init
    FK_init = F_init * K_init
    KA_init = K_init * (A - I(n))
    KBG_init = K_init * BG_init
    
    set_start_value.(BG, BG_init)
    set_start_value.(FK, FK_init)
    set_start_value.(KA, KA_init)
    set_start_value.(KBG, KBG_init)
    
    if symmetric
        # Symmetric formulation: variables split into positive and negative parts
        @variable(model, 0 <= Hp[1:lf,1:lf] <= 300)
        @variable(model, 0 <= Hm[1:lf,1:lf] <= 300)
        @variable(model, 0 <= Lp[1:lf,1:lf] <= 300)
        @variable(model, 0 <= Lm[1:lf,1:lf] <= 300)
        @variable(model, 0 <= Mp[1:lf,1:lf] <= 300)
        @variable(model, 0 <= Mm[1:lf,1:lf] <= 300)
        @variable(model, 0 <= Np[1:lf,1:lf] <= 300)
        @variable(model, 0 <= Nm[1:lf,1:lf] <= 300)
        @variable(model, 0 <= Pp[1:lf,1:r] <= 300)
        @variable(model, 0 <= Pm[1:lf,1:r] <= 300)
        @variable(model, 0 <= Tp[1:s, 1:lf] <= 300)
        @variable(model, 0 <= Tm[1:s, 1:lf] <= 300)
        @variable(model, 0 <= Zp[1:lf, 1:s] <= 300)
        @variable(model, 0 <= Zm[1:lf, 1:s] <= 300)
        
        # Auxiliary variables for d*Mp, d*Mm, etc.
        @variable(model, dMp[1:lf,1:lf])
        @variable(model, dMm[1:lf,1:lf])
        @variable(model, dNp[1:lf,1:lf])
        @variable(model, dNm[1:lf,1:lf])
        @constraint(model, dMp .== d * Mp)
        @constraint(model, dMm .== d * Mm)
        @constraint(model, dNp .== d * Np)
        @constraint(model, dNm .== d * Nm)
        H = Hp - Hm
        L = Lp - Lm
        M = Mp - Mm
        N = Np - Nm
        P = Pp - Pm
        T = Tp - Tm
        Z = Zp - Zm

        # Initialize symmetric pairs
        for i in 1:lf, j in 1:lf
            set_start_value(Hp[i,j], i == j ? 0.1 : 0.0)
            set_start_value(Lp[i,j], i == j ? 0.1 : 0.0)
            set_start_value(Mp[i,j], i == j ? 0.01 : 0.0)
            set_start_value(Np[i,j], i == j ? 0.01 : 0.0)
            set_start_value(Hm[i,j], 0.0)
            set_start_value(Lm[i,j], 0.0)
            set_start_value(Mm[i,j], 0.0)
            set_start_value(Nm[i,j], 0.0)
        end
        for i in 1:lf, j in 1:r
            set_start_value(Pp[i,j], 0.01)
            set_start_value(Pm[i,j], 0.0)
        end
        for i in 1:s, j in 1:lf
            set_start_value(Tp[i,j], 0.01)
            set_start_value(Tm[i,j], 0.0)
        end
        for i in 1:lf, j in 1:s
            set_start_value(Zp[i,j], 0.1)
            set_start_value(Zm[i,j], 0.0)
        end

        # Start values for dMp, dMm, dNp, dNm
        Mp_init = [i == j ? 0.01 : 0.0 for i in 1:lf, j in 1:lf]
        Mm_init = zeros(lf, lf)
        Np_init = [i == j ? 0.01 : 0.0 for i in 1:lf, j in 1:lf]
        Nm_init = zeros(lf, lf)
        set_start_value.(dMp, d * Mp_init)
        set_start_value.(dMm, d * Mm_init)
        set_start_value.(dNp, d * Np_init)
        set_start_value.(dNm, d * Nm_init)

        # Auxiliary products for symmetric formulation
        @variable(model, HFp[1:lf, 1:n])
        @variable(model, HFm[1:lf, 1:n])
        @variable(model, LFp[1:lf, 1:n])
        @variable(model, LFm[1:lf, 1:n])
        @variable(model, MFp[1:lf, 1:n])
        @variable(model, MFm[1:lf, 1:n])
        @variable(model, NFp[1:lf, 1:n])
        @variable(model, NFm[1:lf, 1:n])
        @variable(model, PRp[1:lf, 1:r])
        @variable(model, PRm[1:lf, 1:r])
        @variable(model, TFp[1:s, 1:n])
        @variable(model, TFm[1:s, 1:n])
        @variable(model, ZSp[1:lf, 1:n])
        @variable(model, ZSm[1:lf, 1:n])

        @constraint(model, HFp .== Hp*F)
        @constraint(model, HFm .== Hm*F)
        @constraint(model, LFp .== Lp*F)
        @constraint(model, LFm .== Lm*F)
        @constraint(model, MFp .== Mp*F)
        @constraint(model, MFm .== Mm*F)
        @constraint(model, NFp .== Np*F)
        @constraint(model, NFm .== Nm*F)
        @constraint(model, PRp .== Pp*R)
        @constraint(model, PRm .== Pm*R)
        @constraint(model, TFp .== Tp*F)
        @constraint(model, TFm .== Tm*F)
        @constraint(model, ZSp .== Zp*S)
        @constraint(model, ZSm .== Zm*S)

        @constraint(model, HFp .- HFm .== F*(A+K))
        @constraint(model, LFp .- LFm .== F*(BG - K))
        @constraint(model, MFp .- MFm .== -F*KA)
        @constraint(model, NFp .- NFm .== -F*KBG)
        @constraint(model, PRp .- PRm .== F*E)
        @constraint(model, TFp .- TFm .== S)
        @constraint(model, ZSp .- ZSm .== F)

        @constraint(model, ((Hp + Hm) + (Lp + Lm) + ((dMp + dMm) + (dNp + dNm)))*ones_f + (Pp + Pm)*ones_r .<= ones_f*lambda)
        @constraint(model, (Tp + Tm)*ones_f .<= ones_s)
        @constraint(model, (Zp + Zm)*ones_s .<= ones_f*igamma)
        
    else
        # Standard formulation: non-negative variables only
        @variable(model, 0 <= H[1:lf,1:lf] <= 300)
        @variable(model, 0 <= L[1:lf,1:lf] <= 300) #Farkas para A
        @variable(model, 0 <= M[1:lf,1:lf] <= 300) #Farkas para A
        @variable(model, 0 <= N[1:lf,1:lf] <= 300) #Farkas para A
        @variable(model, 0 <= P[1:lf,1:r] <= 300) #Farkas para E
        @variable(model, 0 <= T[1:s, 1:lf] <= 300) #Pra fazer a inclusão de L em S
        @variable(model, 0 <= Z[1:lf, 1:s] <= 300) #Pra fazer a restrição do igamma

        # Auxiliary variables for d*M and d*N
        @variable(model, dM[1:lf,1:lf])
        @variable(model, dN[1:lf,1:lf])
        @constraint(model, dM .== d * M)
        @constraint(model, dN .== d * N)
        for i in 1:lf, j in 1:lf
            set_start_value(H[i,j], i == j ? 0.1 : 0.0)
            set_start_value(L[i,j], i == j ? 0.1 : 0.0)
            set_start_value(M[i,j], i == j ? 0.01 : 0.0)
            set_start_value(N[i,j], i == j ? 0.01 : 0.0)
        end
        for i in 1:lf, j in 1:r
            set_start_value(P[i,j], 0.01)
        end
        for i in 1:s, j in 1:lf
            set_start_value(T[i,j], 0.01)
        end
        for i in 1:lf, j in 1:s
            set_start_value(Z[i,j], 0.1)
        end

        # Start values for dM and dN
        M_init = [i == j ? 0.01 : 0.0 for i in 1:lf, j in 1:lf]
        N_init = [i == j ? 0.01 : 0.0 for i in 1:lf, j in 1:lf]
        set_start_value.(dM, d * M_init)
        set_start_value.(dN, d * N_init)

        # Auxiliary variables for standard formulation
        @variable(model, HF[1:lf, 1:n])
        @variable(model, LF[1:lf, 1:n])
        @variable(model, MF[1:lf, 1:n])
        @variable(model, NF[1:lf, 1:n])
        @variable(model, PR[1:lf, 1:r])
        @variable(model, TF[1:s, 1:n])
        @variable(model, ZS[1:lf, 1:n])

        @constraint(model, HF .== H*F)
        @constraint(model, LF .== L*F)
        @constraint(model, MF .== M*F)
        @constraint(model, NF .== N*F)
        @constraint(model, PR .== P*R)
        @constraint(model, TF .== T*F)
        @constraint(model, ZS .== Z*S)
        
        @constraint(model, HF .== F*(A+K))
        @constraint(model, LF .== F*(BG - K))
        @constraint(model, MF .== -F*KA)
        @constraint(model, NF .== -F*KBG)
        @constraint(model, PR .== F*E)
        @constraint(model, TF .== S)
        @constraint(model, ZS .== F)

        @constraint(model, (H + L + (dM + dN))*ones_f + P*ones_r .<= ones_f*lambda)
        @constraint(model, T*ones_f .<= ones_s)
        @constraint(model, Z*ones_s .<= ones_f*igamma)
    end

    @objective(model, Min, igamma)

    optimize!(model)

    F = value.(F)
    G = value.(G)

    result = Dict("F" => F, "G" => G)

    print(termination_status(model))

    return result
end
