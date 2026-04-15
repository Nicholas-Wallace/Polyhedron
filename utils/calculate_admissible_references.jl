function calculate_admissible_references(F, A, E)
    lf, n = size(F)
    Ksi = zeros(lf)

    # 4. Acha o vetor Ksi
    # Para cada linha i de F, resolve um PL: max (F_i * A * x) sujeito a Fx <= 1
    for i in 1:lf
        # O objetivo é F_i * A * x. F_i é a linha i de F.
        obj_coefs = F[i, :]' * A
        
        model = Model(HiGHS.Optimizer)
        set_silent(model)
        
        @variable(model, x[1:n])
        @constraint(model, F * x .<= ones(lf))
        @objective(model, Max, sum(obj_coefs[j] * x[j] for j in 1:n))
        
        optimize!(model)
        
        if termination_status(model) == MOI.OPTIMAL
            Ksi[i] = objective_value(model)
        else
            println("Aviso: PL não convergiu para a linha $i de Ksi.")
            Ksi[i] = NaN
        end
    end

    # 5. Calcula o poliedro FEr <= 0.999*ones() - Ksi
    # O resultado são as matrizes que definem o poliedro H-rep para r: H_r * r <= k_r
    H_r = F * E
    k_r = 0.999 .* ones(lf) - Ksi

    return Ksi, H_r, k_r
end