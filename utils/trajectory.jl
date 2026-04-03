function trajectory(x0, A, passos)
    for i in range(2, passos)
        try
            x0 = hcat(x0, A*x0[:, i-1])    
        catch e
            print(e)
        end
        
    end
    return x0   
end