function walkeranalysis()
    ##############################
    # HERBIVORE POSTCRANIAL DATA
    ##############################

    pcdata = CSV.read("$(homedir())/Dropbox/Postdoc/2018_eigenvec/pc_bones.csv",header=true);
    sp = pcdata[:common_name];
    genus = pcdata[:Genus];
    species = pcdata[:Species];
    _array = Array{Char}(length(genus)); _array[1:length(genus)] = ' ';
    gensp = mapslices(join, [Array{String}(genus) _array Array{String}(species)], 2);
    nsp = size(pcdata)[1];
    nmeas = size(pcdata)[2];

    #Archive body mass measure
    bodymass = Array(pcdata[:mass_max_kg]);
    #take out body mass data info
    remove = [find(x->x==:mass_max_kg,names(pcdata));find(x->x==:mass_min_kg,names(pcdata))]
    keep = deleteat!(collect(1:nmeas),remove);
    pcdata = pcdata[:,keep];
    nmeas = size(pcdata)[2];
    pcdatatr = copy(pcdata[:,6:nmeas]);
    #Delete measurements that are missing for ALL species
    nmeas = size(pcdatatr)[2];
    todelete = Array{Int64}(0);
    for i=1:nmeas
        if all(ismissing.(pcdatatr[:,i]))
            push!(todelete,i);
        end
    end
    pcdatatr = pcdatatr[:,setdiff(collect(1:nmeas),todelete)];
    meas = Array{String}(names(pcdatatr));
    nmeas = size(pcdatatr)[2];

    # pcdatatr = Array(pcdatatr) ./ Array(pcdatatr[:femur_length]);

    #scale to 'mass'
    #Femur length
    fl = copy(Array(pcdatatr[:femur_length]));
    for i=1:nmeas
        pcdatatr[:,i] = Array(pcdatatr[:,i]) ./ (bodymass.^(0.25));
    end

    #Eliminate measures that are empty for all species
    nomeas = Array{Int64}(nsp);
    for i=1:nsp
        nomeas[i] = length(find(ismissing,Array(pcdatatr[i,:])));
    end
    bydataamt = sortperm(nomeas);

    PC = Array{Float64}(nsp,nsp);
    # measmeans = mean(pcdatatr[!isnothing(pcdatatr)],1);
    #Build similarity matrix
    for i = 0:(nsp^2 - 1)
        a = mod(i,nsp) + 1;
        b = Int64(floor(i/nsp)) + 1;
        if a == b
            PC[a,b] = 0.0;
            continue
        end
        ct = 0;
        ctones = 0;
        for j = 1:nmeas
            if !ismissing(pcdatatr[a,j]) && !ismissing(pcdatatr[b,j])
                ct += log(minimum([pcdatatr[a,j],pcdatatr[b,j]])/maximum([pcdatatr[a,j],pcdatatr[b,j]]));
                ctones += 1;
            end
        end
        ctscaled = exp(ct/ctones);
        PC[a,b] = Float64(ctscaled); #/Float64(ctones);s
    end

    S = laplacian(PC,10);
    ev = eigs(S; nev=10,which=:SR);
    eval = ev[1];
    evecs = ev[2];
    
    return(
    pcdata,
    pcdatatr,
    gensp,
    eval,
    evecs
    )
    
end