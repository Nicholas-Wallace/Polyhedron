
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

function finding_L_pinvariant_segref_delay(A, B, E, S, R, d; lambda=0.99, time=10, lf=10) 
    
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
    @variable(model, 0 <= H[1:lf,1:lf] <= 300)
    @variable(model, 0 <= L[1:lf,1:lf] <= 300) #Farkas para A
    @variable(model, 0 <= M[1:lf,1:lf] <= 300) #Farkas para A
    @variable(model, 0 <= N[1:lf,1:lf] <= 300) #Farkas para A
    @variable(model, 0 <= P[1:lf,1:r] <= 300) #Farkas para E
    @variable(model, 0 <= T[1:s, 1:lf] <= 300) #Pra fazer a inclusão de L em S
    @variable(model, J[1:n, 1:lf]) #Pseudo inversa

    @variable(model, 0.001 <= igamma <= 1000)
    @variable(model, 0 <= Z[1:lf, 1:s] <= 300)

    @objective(model, Min, igamma)

    @constraint(model, H*F == F*(A+K)) 
    @constraint(model, L*F == F*(B*G - K))
    @constraint(model, M*F == -F*K*(A-I(n)))
    @constraint(model, N*F == -F*K*B*G)
    @constraint(model, P*R == F*E)
    @constraint(model, T*F == S)
    @constraint(model, J*F == I(n))

    # Restrições do igamma
    @constraint(model, Z*S == F)
    @constraint(model, Z*ones_s .<= ones_f*igamma)
    
    @constraint(model, (H + L + d*(M + N))*ones_f + P*ones_r .<= ones_f*lambda)
    @constraint(model, T*ones_f .<= ones_s)

    #0  Auto (Knitro escolhe sozinho).
    #1: Active Set 
    #2: Interior Point / CG
    #3: Interior Point / Direct
    #4: SQP (Programação Quadrática Sequencial)

    set_optimizer_attribute(model, "outlev", 2)
    set_optimizer_attribute(model, "algorithm", 2)
    set_optimizer_attribute(model, "maxtime", time*60.0)
    set_optimizer_attribute(model, "maxit", 10000)
    set_optimizer_attribute(model, "ms_enable", 1)
    set_optimizer_attribute(model, "ms_maxsolves", 100)

    optimize!(model)

    F = value.(F)
    G = value.(G)
    K = value.(K)
    H = value.(H)
    L = value.(L)
    M = value.(M)
    N = value.(N)
    P = value.(P)
    T = value.(T)
    Z = value.(Z)
    J = value.(J)

    result = Dict(
    "F" => F,
    "G" => G,
    "K" => K,
    "H" => H,
    "L" => L, 
    "M" => M,
    "N" => N,
    "P" => P, 
    "T" => T,
    "Z" => Z
    )

    # Salvar matrizes em arquivo
    open("results_pinvariant_segref_delay.txt", "w") do io
        println(io, "=== Resultados da otimização ===")
        println(io, "Status: $(termination_status(model))\n")
        
        println(io, "=== Matriz F (Poliedro) ===")
        writedlm(io, F)
        println(io, "")
        
        println(io, "=== Matriz G (Ganho do controlador) ===")
        writedlm(io, G)
        println(io, "")
        
        println(io, "=== Matriz K (Variável auxiliar) ===")
        writedlm(io, K)
        println(io, "")
        
        println(io, "=== Matriz H (Farkas) ===")
        writedlm(io, H)
        println(io, "")
        
        println(io, "=== Matriz L (Farkas) ===")
        writedlm(io, L)
        println(io, "")
        
        println(io, "=== Matriz M (Farkas) ===")
        writedlm(io, M)
        println(io, "")
        
        println(io, "=== Matriz N (Farkas) ===")
        writedlm(io, N)
        println(io, "")
        
        println(io, "=== Matriz P (Farkas para E) ===")
        writedlm(io, P)
        println(io, "")
        
        println(io, "=== Matriz T (Inclusão em S) ===")
        writedlm(io, T)
        println(io, "")
        
        println(io, "=== Matriz Z ===")
        writedlm(io, Z)
        println(io, "")
    end

    print(termination_status(model))

    return result
end

function finding_L_pinvariant_segref_delay2(A, B, S, d; a=0.5, time=10, lf=10) 
    
    model = Model() do
        return NEOSServer.Optimizer(; email = "wallace.lopes.162@ufrn.edu.br", solver = "Knitro")
    end

    # Parâmetros
    n = size(A, 1) # Ordem do sistema (Linhas de A)
    m = size(B, 2) # Número de Entradas (Colunas de B)
    s = size(S, 1) # Número de linhas da matriz das restrições de X[k]

    ones_f = ones(lf)
    ones_s = ones(s)

    # Variáveis
    @variable(model, K[1:n, 1:n]) 
    @variable(model, G[1:m, 1:n]) 
    @variable(model, -500 <= F[1:lf, 1:n] <= 500) 
    @variable(model, 0 <= H[1:lf, 1:lf] <= 300)
    @variable(model, 0 <= L[1:lf, 1:lf] <= 300) 
    @variable(model, 0 <= M[1:lf, 1:lf] <= 300) 
    @variable(model, 0 <= N[1:lf, 1:lf] <= 300) 
    @variable(model, 0 <= T[1:s, 1:lf] <= 300) 
    @variable(model, J[1:n, 1:lf]) 

    @variable(model, 0.001 <= igamma <= 1000)
    @variable(model, 0.001 <= lambda <= 0.999) # Lambda agora é variável de decisão
    @variable(model, 0 <= Z[1:lf, 1:s] <= 300)

    # 1. & 2. Nova Função Objetivo (min a * igama + (1-a) * lambda)
    @objective(model, Min, a * igamma + (1 - a) * lambda)

    # Restrições de Igualdade (Lema de Farkas sem P/E/R)
    @constraint(model, H*F .== F*(A+K)) 
    @constraint(model, L*F .== F*(B*G - K))
    @constraint(model, M*F .== -F*K*(A-I(n)))
    @constraint(model, N*F .== -F*K*B*G)
    @constraint(model, T*F .== S)
    @constraint(model, J*F .== I(n))

    # Restrições de igamma (Volume do Poliedro)
    @constraint(model, Z*S .== F)
    @constraint(model, Z*ones_s .<= ones_f*igamma)
    
    # 3. Nova Desigualdade Principal
    @constraint(model, (H + L + d*(M + N))*ones_f .<= ones_f*lambda)
    @constraint(model, T*ones_f .<= ones_s)

    # Configurações do Solver Knitro
    set_optimizer_attribute(model, "outlev", 2)
    set_optimizer_attribute(model, "algorithm", 2)
    set_optimizer_attribute(model, "maxtime", time*60)
    set_optimizer_attribute(model, "maxit", 100000)
    set_optimizer_attribute(model, "ms_enable", 1)
    set_optimizer_attribute(model, "ms_maxsolves", 100)

    optimize!(model)

    F_val = value.(F)
    G_val = value.(G)
    H_val = value.(H)
    L_val = value.(L)
    M_val = value.(M)
    N_val = value.(N)
    T_val = value.(T)
    Z_val = value.(Z)
    K_val = value.(K)
    lambda_val = value(lambda)
    igamma_val = value(igamma)

    result = Dict("F" => F_val, "G" => G_val, "K" => K_val, "H" => H_val, 
                  "L" => L_val, "M" => M_val, "N" => N_val, "T" => T_val, 
                  "Z" => Z_val, "lambda" => lambda_val, "igamma" => igamma_val)

    # Salvar matrizes em arquivo
    open("results_pinvariant_segref_delay.txt", "w") do io
        println(io, "=== Resultados da otimização ===")
        println(io, "Status: $(termination_status(model))\n")
        
        println(io, "=== Matriz F (Poliedro) ===")
        writedlm(io, F)
        println(io, "")
        
        println(io, "=== Matriz G (Ganho do controlador) ===")
        writedlm(io, G)
        println(io, "")
        
        println(io, "=== Matriz K (Variável auxiliar) ===")
        writedlm(io, K)
        println(io, "")
        
        println(io, "=== Matriz H (Farkas) ===")
        writedlm(io, H)
        println(io, "")
        
        println(io, "=== Matriz L (Farkas) ===")
        writedlm(io, L)
        println(io, "")
        
        println(io, "=== Matriz M (Farkas) ===")
        writedlm(io, M)
        println(io, "")
        
        println(io, "=== Matriz N (Farkas) ===")
        writedlm(io, N)
        println(io, "")
        
        println(io, "=== Matriz T (Inclusão em S) ===")
        writedlm(io, T)
        println(io, "")
        
        println(io, "=== Matriz Z ===")
        writedlm(io, Z)
        println(io, "")

        println(io, "=== lambda ===")
        writedlm(io, lambda)
        println(io, "")

        println(io, "=== igamma ===")
        writedlm(io, igamma)
        println(io, "")
    end

    print(termination_status(model))

    return result
end

function finding_L_pinvariant_segref_delay_sim(A, B, E, S, R, d; lambda=0.99, time=10, lf=10) 
    
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
    @variable(model, 0 <= Hm[1:lf,1:lf] <= 300)
    @variable(model, 0 <= Hp[1:lf,1:lf] <= 300)
    @variable(model, 0 <= Lm[1:lf,1:lf] <= 300) #Farkas para A
    @variable(model, 0 <= Lp[1:lf,1:lf] <= 300)
    @variable(model, 0 <= Mm[1:lf,1:lf] <= 300) #Farkas para A
    @variable(model, 0 <= Mp[1:lf,1:lf] <= 300)
    @variable(model, 0 <= Nm[1:lf,1:lf] <= 300) #Farkas para A
    @variable(model, 0 <= Np[1:lf,1:lf] <= 300)
    @variable(model, 0 <= Pm[1:lf,1:r] <= 300) #Farkas para E
    @variable(model, 0 <= Pp[1:lf,1:r] <= 300)
    @variable(model, 0 <= Tp[1:s, 1:lf] <= 300) #Pra fazer a inclusão de L em S
    @variable(model, 0 <= Tm[1:s, 1:lf] <= 300)
    @variable(model, 0 <= Zp[1:lf, 1:s] <= 300)
    @variable(model, 0 <= Zm[1:lf, 1:s] <= 300)
    @variable(model, J[1:n, 1:lf]) #Pseudo inversa

    @variable(model, igamma >= 1)

    @objective(model, Min, igamma)

    @constraint(model, (Hp - Hm)*F == F*(A+K)) 
    @constraint(model, (Lp - Lm)*F == F*(B*G - K))
    @constraint(model, (Mp - Mm)*F == -F*K*(A-I(n)))
    @constraint(model, (Np - Nm)*F == -F*K*B*G)
    @constraint(model, (Pp - Pm)*R == F*E)
    @constraint(model, (Tp - Tm)*F == S)
    @constraint(model, J*F == I(n))

    # Restrições do igamma
    @constraint(model, (Zp - Zm)*S == F)
    @constraint(model, (Zp + Zm)*ones_s .<= ones_f*igamma)
    
    @constraint(model, ((Hp + Hm) + (Lp + Lm) + d*((Mp + Mm) + (Np + Nm)))*ones_f + (Pp + Pm)*ones_r .<= ones_f*lambda)
    @constraint(model, (Tp + Tm)*ones_f .<= ones_s)

    #0  Auto (Knitro escolhe sozinho).
    #1: Active Set 
    #2: Interior Point / CG
    #3: Interior Point / Direct
    #4: SQP (Programação Quadrática Sequencial)

    set_optimizer_attribute(model, "outlev", 2)
    set_optimizer_attribute(model, "algorithm", 1)
    set_optimizer_attribute(model, "maxtime", time*60)
    set_optimizer_attribute(model, "maxit", 10000)
    set_optimizer_attribute(model, "ms_enable", 1)
    set_optimizer_attribute(model, "ms_maxsolves", 100)

    optimize!(model)

    F = value.(F)
    G = value.(G)
    Hp = value.(Hp)
    Hm = value.(Hm)
    Lp = value.(Lp)
    Lm = value.(Lm)
    Mp = value.(Mp)
    Mm = value.(Mm)
    Np = value.(Np)
    Nm = value.(Nm)
    Pp = value.(Pp)
    Pm = value.(Pm)
    Tp = value.(Tp)
    Tm = value.(Tm)
    Zp = value.(Zp)
    Zm = value.(Zm)
    K = value.(K)
    J = value.(J)

    result = Dict("F" => F, "G" => G, "K" => K, "J" => J, 
                  "H_diff" => Hp - Hm, "L_diff" => Lp - Lm, "M_diff" => Mp - Mm, "N_diff" => Np - Nm, 
                  "P_diff" => Pp - Pm, "T_diff" => Tp - Tm, "Z_diff" => Zp - Zm,
                  "H_sum" => Hp + Hm, "L_sum" => Lp + Lm, "M_sum" => Mp + Mm, "N_sum" => Np + Nm,
                  "P_sum" => Pp + Pm, "T_sum" => Tp + Tm, "Z_sum" => Zp + Zm,
                  "Hp" => Hp, "Hm" => Hm, "Lp" => Lp, "Lm" => Lm, 
                  "Mp" => Mp, "Mm" => Mm, "Np" => Np, "Nm" => Nm,
                  "Pp" => Pp, "Pm" => Pm, "Tp" => Tp, "Tm" => Tm, 
                  "Zp" => Zp, "Zm" => Zm)

    # Salvar matrizes em arquivo (valores simétricos como a diferença)
    open("results_pinvariant_segref_delay_sim.txt", "w") do io
        println(io, "=== Resultados da otimização (versão simétrica) ===")
        println(io, "Status: $(termination_status(model))\n")
        
        println(io, "=== Matriz F (Poliedro) ===")
        writedlm(io, F)
        println(io, "")
        
        println(io, "=== Matriz G (Ganho do controlador) ===")
        writedlm(io, G)
        println(io, "")
        
        println(io, "=== Matriz K (Variável auxiliar) ===")
        writedlm(io, K)
        println(io, "")
        
        println(io, "=== Matriz J (Pseudo inversa) ===")
        writedlm(io, J)
        println(io, "")
        
        println(io, "=== Matriz H (Farkas - diferença Hp - Hm) ===")
        writedlm(io, Hp - Hm)
        println(io, "")
        
        println(io, "=== Matriz L (Farkas - diferença Lp - Lm) ===")
        writedlm(io, Lp - Lm)
        println(io, "")
        
        println(io, "=== Matriz M (Farkas - diferença Mp - Mm) ===")
        writedlm(io, Mp - Mm)
        println(io, "")
        
        println(io, "=== Matriz N (Farkas - diferença Np - Nm) ===")
        writedlm(io, Np - Nm)
        println(io, "")
        
        println(io, "=== Matriz P (Farkas para E - diferença Pp - Pm) ===")
        writedlm(io, Pp - Pm)
        println(io, "")
        
        println(io, "=== Matriz T (Inclusão em S - diferença Tp - Tm) ===")
        writedlm(io, Tp - Tm)
        println(io, "")
        
        println(io, "=== Matriz Z (diferença Zp - Zm) ===")
        writedlm(io, Zp - Zm)
        println(io, "")

        println(io, "=== Valores de H em valor absoluto (Hp + Hm) ===")
        writedlm(io, Hp + Hm)
        println(io, "")
        
        println(io, "=== Valores de L em valor absoluto (Lp + Lm) ===")
        writedlm(io, Lp + Lm)
        println(io, "")
        
        println(io, "=== Valores de M em valor absoluto (Mp + Mm) ===")
        writedlm(io, Mp + Mm)
        println(io, "")
        
        println(io, "=== Valores de N em valor absoluto (Np + Nm) ===")
        writedlm(io, Np + Nm)
        println(io, "")
        
        println(io, "=== Valores de P em valor absoluto (Pp + Pm) ===")
        writedlm(io, Pp + Pm)
        println(io, "")
        
        println(io, "=== Valores de T em valor absoluto (Tp + Tm) ===")
        writedlm(io, Tp + Tm)
        println(io, "")
        
        println(io, "=== Valores de Z em valor absoluto (Zp + Zm) ===")
        writedlm(io, Zp + Zm)
        println(io, "")
    end

    print(termination_status(model))

    return result
end

function finding_L_pinvariant_segref_delay_sim2(A, B, E, S, R, V, d; lambda=0.99, time=10, lf=10) 
    
    model = Model() do
        return NEOSServer.Optimizer(; email = "wallace.lopes.162@ufrn.edu.br", solver = "Knitro")
    end

    #Parâmetros
    n = size(A, 1); #Ordem do sistema (Linhas de A)
    m = size(B, 2); #Número de Entradas (Colunas de B)
    r = size(R, 1); #Número de linhas das restrições em r[k]
    s = size(S, 1); #Número de linhas da matriz das restrições de X[k]
    v = size(V, 1); 

    ones_f = ones(lf)
    ones_r = ones(r)
    ones_s = ones(s)
    ones_v = ones(v)


    @variable(model, K[1:n, 1:n]) #Variavel auxiliar trans model
    @variable(model, G[1:m, 1:n]) #Ganho do controlador 
    @variable(model, -500 <= F[1:lf,1:n] <= 500) #Poliedro que queremos achar
    @variable(model, 0 <= Hm[1:lf,1:lf] <= 300)
    @variable(model, 0 <= Hp[1:lf,1:lf] <= 300)
    @variable(model, 0 <= Lm[1:lf,1:lf] <= 300) #Farkas para A
    @variable(model, 0 <= Lp[1:lf,1:lf] <= 300)
    @variable(model, 0 <= Mm[1:lf,1:lf] <= 300) #Farkas para A
    @variable(model, 0 <= Mp[1:lf,1:lf] <= 300)
    @variable(model, 0 <= Nm[1:lf,1:lf] <= 300) #Farkas para A
    @variable(model, 0 <= Np[1:lf,1:lf] <= 300)
    @variable(model, 0 <= Pm[1:lf,1:r] <= 300) #Farkas para E
    @variable(model, 0 <= Pp[1:lf,1:r] <= 300)
    @variable(model, 0 <= Tp[1:s, 1:lf] <= 300) #Pra fazer a inclusão de L em S
    @variable(model, 0 <= Tm[1:s, 1:lf] <= 300)
    @variable(model, 0 <= Zp[1:lf, 1:s] <= 300)
    @variable(model, 0 <= Zm[1:lf, 1:s] <= 300)
    @variable(model, 0 <= Up[1:v, 1:lf] <= 300)
    @variable(model, 0 <= Um[1:v, 1:lf] <= 300)
    @variable(model, J[1:n, 1:lf]) #Pseudo inversa

    @variable(model, igamma >= 1)

    @objective(model, Min, igamma)

    @constraint(model, (Hp - Hm)*F == F*(A+K)) 
    @constraint(model, (Lp - Lm)*F == F*(B*G - K))
    @constraint(model, (Mp - Mm)*F == -F*K*(A-I(n)))
    @constraint(model, (Np - Nm)*F == -F*K*B*G)
    @constraint(model, (Pp - Pm)*R == F*E)
    @constraint(model, (Tp - Tm)*F == S)
    @constraint(model, J*F == I(n))

    # Restrições do igamma
    @constraint(model, (Zp - Zm)*S == F)
    @constraint(model, (Zp + Zm)*ones_s .<= ones_f*igamma)

    # Inclusão do poliedro dos estados no poliedro das restrições do u[k]
    @constraint(model, (Up - Um)*F == V*G)
    @constraint(model, (Up + Um)*ones_f .<= ones_v)
    
    @constraint(model, ((Hp + Hm) + (Lp + Lm) + d*((Mp + Mm) + (Np + Nm)))*ones_f + (Pp + Pm)*ones_r .<= ones_f*lambda)
    @constraint(model, (Tp + Tm)*ones_f .<= ones_s)

    #0  Auto (Knitro escolhe sozinho).
    #1: Active Set 
    #2: Interior Point / CG
    #3: Interior Point / Direct
    #4: SQP (Programação Quadrática Sequencial)

    set_optimizer_attribute(model, "outlev", 2)
    set_optimizer_attribute(model, "algorithm", 1)
    set_optimizer_attribute(model, "maxtime", time*60)
    set_optimizer_attribute(model, "maxit", 10000)
    set_optimizer_attribute(model, "ms_enable", 1)
    set_optimizer_attribute(model, "ms_maxsolves", 100)

    optimize!(model)

    F = value.(F)
    G = value.(G)
    Hp = value.(Hp)
    Hm = value.(Hm)
    Lp = value.(Lp)
    Lm = value.(Lm)
    Mp = value.(Mp)
    Mm = value.(Mm)
    Np = value.(Np)
    Nm = value.(Nm)
    Pp = value.(Pp)
    Pm = value.(Pm)
    Tp = value.(Tp)
    Tm = value.(Tm)
    Zp = value.(Zp)
    Zm = value.(Zm)
    Up = value.(Up)
    Um = value.(Um)
    K = value.(K)
    J = value.(J)

    result = Dict("F" => F, "G" => G, "K" => K, "J" => J, 
                  "H_diff" => Hp - Hm, "L_diff" => Lp - Lm, "M_diff" => Mp - Mm, "N_diff" => Np - Nm, 
                  "P_diff" => Pp - Pm, "T_diff" => Tp - Tm, "Z_diff" => Zp - Zm,
                  "H_sum" => Hp + Hm, "L_sum" => Lp + Lm, "M_sum" => Mp + Mm, "N_sum" => Np + Nm,
                  "P_sum" => Pp + Pm, "T_sum" => Tp + Tm, "Z_sum" => Zp + Zm,
                  "U_diff" => Up - Um, "U_sum" => Up + Um,
                  "Hp" => Hp, "Hm" => Hm, "Lp" => Lp, "Lm" => Lm, 
                  "Mp" => Mp, "Mm" => Mm, "Np" => Np, "Nm" => Nm,
                  "Pp" => Pp, "Pm" => Pm, "Tp" => Tp, "Tm" => Tm, 
                  "Zp" => Zp, "Zm" => Zm, "Up" => Up, "Um" => Um)

    # Salvar matrizes em arquivo (valores simétricos como a diferença)
    open("results_pinvariant_segref_delay_sim.txt", "w") do io
        println(io, "=== Resultados da otimização (versão simétrica) ===")
        println(io, "Status: $(termination_status(model))\n")
        
        println(io, "=== Matriz F (Poliedro) ===")
        writedlm(io, F)
        println(io, "")
        
        println(io, "=== Matriz G (Ganho do controlador) ===")
        writedlm(io, G)
        println(io, "")
        
        println(io, "=== Matriz K (Variável auxiliar) ===")
        writedlm(io, K)
        println(io, "")
        
        println(io, "=== Matriz J (Pseudo inversa) ===")
        writedlm(io, J)
        println(io, "")
        
        println(io, "=== Matriz H (Farkas - diferença Hp - Hm) ===")
        writedlm(io, Hp - Hm)
        println(io, "")
        
        println(io, "=== Matriz L (Farkas - diferença Lp - Lm) ===")
        writedlm(io, Lp - Lm)
        println(io, "")
        
        println(io, "=== Matriz M (Farkas - diferença Mp - Mm) ===")
        writedlm(io, Mp - Mm)
        println(io, "")
        
        println(io, "=== Matriz N (Farkas - diferença Np - Nm) ===")
        writedlm(io, Np - Nm)
        println(io, "")
        
        println(io, "=== Matriz P (Farkas para E - diferença Pp - Pm) ===")
        writedlm(io, Pp - Pm)
        println(io, "")
        
        println(io, "=== Matriz T (Inclusão em S - diferença Tp - Tm) ===")
        writedlm(io, Tp - Tm)
        println(io, "")
        
        println(io, "=== Matriz Z (diferença Zp - Zm) ===")
        writedlm(io, Zp - Zm)
        println(io, "")

        println(io, "=== Matriz U (diferença Up - Um) ===")
        writedlm(io, Up - Um)
        println(io, "")

        println(io, "=== Valores de H em valor absoluto (Hp + Hm) ===")
        writedlm(io, Hp + Hm)
        println(io, "")
        
        println(io, "=== Valores de L em valor absoluto (Lp + Lm) ===")
        writedlm(io, Lp + Lm)
        println(io, "")
        
        println(io, "=== Valores de M em valor absoluto (Mp + Mm) ===")
        writedlm(io, Mp + Mm)
        println(io, "")
        
        println(io, "=== Valores de N em valor absoluto (Np + Nm) ===")
        writedlm(io, Np + Nm)
        println(io, "")
        
        println(io, "=== Valores de P em valor absoluto (Pp + Pm) ===")
        writedlm(io, Pp + Pm)
        println(io, "")
        
        println(io, "=== Valores de T em valor absoluto (Tp + Tm) ===")
        writedlm(io, Tp + Tm)
        println(io, "")
        
        println(io, "=== Valores de Z em valor absoluto (Zp + Zm) ===")
        writedlm(io, Zp + Zm)
        println(io, "")

        println(io, "=== Valores de U em valor absoluto (Up + Um) ===")
        writedlm(io, Up + Um)
        println(io, "")
    end

    print(termination_status(model))

    return result
end
