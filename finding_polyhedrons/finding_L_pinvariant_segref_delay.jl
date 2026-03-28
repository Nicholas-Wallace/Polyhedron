
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
        
        H = Hp - Hm
        L = Lp - Lm
        M = Mp - Mm
        N = Np - Nm
        P = Pp - Pm
        T = Tp - Tm
        Z = Zp - Zm

        # Farkas invariance constraint with absolute values (symmetric)
        @constraint(model, ((Hp + Hm) + (Lp + Lm) + d*((Mp + Mm) + (Np + Nm)))*ones_f + (Pp + Pm)*ones_r .<= ones_f*lambda)
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
        
        # Farkas invariance constraint
        @constraint(model, (H + L + d*(M + N))*ones_f + P*ones_r .<= ones_f*lambda)
        @constraint(model, T*ones_f .<= ones_s)
        @constraint(model, Z*ones_s .<= ones_f*igamma)
    end
    
    @variable(model, J[1:n, 1:lf]) #Pseudo inversa
    @variable(model, igamma >= 1)

    @objective(model, Min, igamma)

    @constraint(model, H*F == F*(A+K)) 
    @constraint(model, L*F == F*(B*G - K))
    @constraint(model, M*F == -F*K*(A-I(n)))
    @constraint(model, N*F == -F*K*B*G)
    @constraint(model, P*R == F*E)
    @constraint(model, T*F == S)
    @constraint(model, J*F == I(n))

    # Restrição do igamma
    @constraint(model, Z*S == F)

    #0  Auto (Knitro escolhe sozinho).
    #1: Active Set 
    #2: Interior Point / CG
    #3: Interior Point / Direct
    #4: SQP (Programação Quadrática Sequencial)

    set_optimizer_attribute(model, "outlev", 2)
    set_optimizer_attribute(model, "algorithm", 1)
    set_optimizer_attribute(model, "maxtime", 3600.0)
    set_optimizer_attribute(model, "maxit", 10000)
    set_optimizer_attribute(model, "ms_enable", 1)
    set_optimizer_attribute(model, "ms_maxsolves", 100)

    optimize!(model)

    F = value.(F)
    G = value.(G)

    result = Dict("F" => F, "G" => G)

    print(termination_status(model))

    return result
end
