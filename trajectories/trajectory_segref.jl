#############################################
### essa função é referente a um sistema: ###
###     x[k+1] = A*x(k) + Er(k)     ### 
#############################################

function trajectory_segref(x0, A, E,passos)
    for i in range(2, passos)
        try
            x0 = hcat(x0, A*x0[:, i-1] + E*r[i-1])    
        catch e
            print(e)
        end
        
    end
    return x0   
end