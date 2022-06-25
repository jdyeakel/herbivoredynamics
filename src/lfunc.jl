function lfunc(x,A,massvec)
    #Calculate Likelihood of proposed alpha, beta, gamma
    alpha = x[1];
    beta = x[2];
    gamma = x[3];
    
    numprey = size(A)[1];
    numpred = size(A)[2];
    massvec_pred = copy(massvec);
    massvec_prey = copy(massvec);
    
    L = Array{Float64}(undef,1);

    #instead of looping over ALL SPECIES, loop only across consumers
    numpreyperpred = vec(sum(A,dims=1));
    consumers = findall(!iszero,numpreyperpred);
    
    let cumLij = 0.0
        for j=consumers
        # for i=1:numpred
            for i=10:numprey
                # if i != j
                    #NOTE: PREY/PRED mass ratio
                    aij = copy(A[i,j]);
                    mi = massvec_prey[i]; #j is prey
                    mj = massvec_pred[j]; #i is predator
                    
                    if mi > 50 && mj > 50
                        Lij = 0.0;
                        
                        pij = exp(alpha + beta*log(mi/mj) + gamma*log(mi/mj)^2)/(1+exp(alpha + beta*log(mi/mj) + gamma*log(mi/mj)^2));
                        
                        if aij == 1 #ERROR: 7/15/20 THIS WAS A ZERO HA
                            Lij = copy(pij);
                        else
                            Lij = copy(1 - pij);
                        end
                        
                        cumLij += log(Lij);
                    end
                # end
            end
        end
        L[1] = -cumLij
    end
    return L[1]
end