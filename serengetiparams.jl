using(DataFrames)
using(CSV)
using(RCall)
using(LinearAlgebra)
using(Distributions)
using(Arpack)
using(Optim)
using(Distributed)
include("$(homedir())/Dropbox/PostDoc/2021_TaranWebs/src/laplacian.jl")
include("$(homedir())/Dropbox/PostDoc/2021_TaranWebs/src/eigencluster.jl")
# @everywhere include("$(homedir())/Dropbox/PostDoc/2021_TaranWebs/src/walkeranalysis.jl")
include("$(homedir())/Dropbox/PostDoc/2021_TaranWebs/src/eigendistance.jl")


include("$(homedir())/Dropbox/PostDoc/2021_TaranWebs/src/logitlikelihood.jl")
include("$(homedir())/Dropbox/PostDoc/2021_TaranWebs/src/lfunc.jl")
include("$(homedir())/Dropbox/PostDoc/2021_TaranWebs/src/fc.jl")


pcdata = CSV.read("$(homedir())/Dropbox/PostDoc/2018_eigenvec/pc_bones.csv",header=true,DataFrame);
sp = pcdata[!,:common_name];
genus = pcdata[!,:Genus];
species = pcdata[!,:Species];
_array = Array{Char}(undef,length(genus)); _array[1:length(genus)] .= ' ';
gensp = mapslices(join, [Array{String}(genus) _array Array{String}(species)], dims=2);
nsp = size(pcdata)[1];
nmeas = size(pcdata)[2];

#Archive body mass measure
bodymass = Array(pcdata[!,:mass_max_kg]);
#take out body mass data info
remove = [findall(x->x==:mass_max_kg,names(pcdata));findall(x->x==:mass_min_kg,names(pcdata))]
keep = deleteat!(collect(1:nmeas),remove);
pcdata = pcdata[:,keep];
nmeas = size(pcdata)[2];
pcdatatr = copy(pcdata[:,6:nmeas]);
#Delete measurements that are missing for ALL species
nmeas = size(pcdatatr)[2];
todelete = Array{Int64}(undef,0);
for i=1:nmeas
    if all(ismissing.(pcdatatr[:,i]))
        push!(todelete,i);
    end
end
pcdatatr = pcdatatr[:,setdiff(collect(1:nmeas),todelete)];
meas = Array{String}(String.(names(pcdatatr)));
nmeas = size(pcdatatr)[2];




baskspcodes = CSV.read("$(homedir())/Dropbox/PostDoc/2018_eigenvec/baskervilledata/Table_S1.csv",header=true,DataFrame);
baskcode = Array(baskspcodes[!,:code]);
baskgensp = Array(baskspcodes[!,:species]);


basktrophic = CSV.read("$(homedir())/Dropbox/PostDoc/2018_eigenvec/baskervilledata/Table_S2.csv",header=true,DataFrame);
preds = basktrophic[!,:pred_code];
preys = basktrophic[!,:prey_code];
npreds = length(preds);

trophic = Array{String}(undef,length(preds),2);
predsunique = unique(preds);
preysunique = unique(preys);
for i = 1:length(predsunique)
    id = findall(x->x==predsunique[i],preds);
    codetosp = findall(x->x==predsunique[i],baskcode);
    trophic[id,1] = repeat(baskgensp[codetosp],outer=length(id));
end
for i = 1:length(preysunique)
    id = findall(x->x==preysunique[i],preys);
    codetosp = findall(x->x==preysunique[i],baskcode);
    trophic[id,2] = repeat(baskgensp[codetosp],outer=length(id));
end
trophiclevel = [repeat([2],inner=9);repeat([1],inner=(32-9))]
#Include only consumers in the serengeti system
spnames = unique(trophic[:,1]);
SerengetiAdj = zeros(Int64,32,9);
for i=1:9
    predloc = findall(x->x==spnames[i],trophic[:,1])
    preyid = trophic[predloc,2];
    preyloc = Array{Int64}(undef,0);
    for j=1:length(preyid)
        sp_position = findall(x->x==preyid[j],spnames)[1];
        push!(preyloc,sp_position)
    end
    SerengetiAdj[preyloc,i] .= 1;
end


#Mass vector
massvec = Array{Float64}(undef,32);
for i=1:32
    spid = spnames[i];
    masspos = findall(x->x==spid,vec(gensp));
    if length(masspos) > 0
        massvec[i] = bodymass[masspos[1]];
    else
        massvec[i] = 0;
    end
end

# Missing masses
# "Canis aureus"
# "Damaliscus korrigum"
# "Heterohyrax brucei"
# "Papio anubis"
# "Pedetes capensis"
# "Procavia capensis"
# "Rhabdomys pumilio"
missingmass = [10,108,2.4,20,3.5,3.8,0.0397];
massvec[findall(iszero,massvec)] .= missingmass;


#Assess logit parameters


## Serengeti
x0 = [0.0,0.0,0.0];
results_ser = optimize(x->lfunc(x,SerengetiAdj,massvec),x0,NelderMead());
results_ser.minimizer
xmax = results_ser.minimizer;
fcorr_bg, Apredict_ser = fc(xmax,SerengetiAdj,massvec)

predtopreyratio = exp(xmax[2]/(2*xmax[3]))



#HAYWARD ANALYSIS
prefdata = CSV.read("$(homedir())/Dropbox/PostDoc/2021_TaranWebs/data/data_preypreference.csv",header=true,DataFrame);
#create a bodymass vector

#Replace some names with closely related species
prefdata[!,:prey][findall(x->x=="Tragelaphus strepsiceros",prefdata[!,:prey])] .= "Tragelaphus eurycerus"
prefdata[!,:prey][findall(x->x=="Tragelaphus angasi",prefdata[!,:prey])] .= "Tragelaphus eurycerus"
prefdata[!,:prey][findall(x->x=="Kobus kob",prefdata[!,:prey])] .= "Kobus ellipsiprymnus"
prefdata[!,:prey][findall(x->x=="Oryx gazelle",prefdata[!,:prey])] .= "Oryx beisa"
prefdata[!,:prey][findall(x->x=="Equus burchelli",prefdata[!,:prey])] .= "Equus quagga"


spnameshay = [unique(prefdata[!,:pred]);unique(prefdata[!,:prey])]
massvechay = Array{Float64}(undef,length(spnameshay));

for i=1:length(spnameshay)
    spid = spnameshay[i];
    masspos = findall(x->x==spid,vec(gensp));
    if length(masspos) > 0
        massvechay[i] = bodymass[masspos[1]];
    else
        massvechay[i] = 0;
    end
end
spnameshay[findall(iszero,massvechay)]
# "Gazella thompsoni"
# "Antidorcas marsupialis"
# "Gazella granti"
# "Damaliscus dorcas phillipsi"
# "Giraffe camelopardis"
missingmasshay = [20,35,65,150,1045];
massvechay[findall(iszero,massvechay)] .= missingmasshay;

#Caculcate mean preferred mass per predator
preds = unique(prefdata[!,:pred]);
prefmassmean = Array{Float64}(undef,length(preds));
predmass = Array{Float64}(undef,length(preds));
maxpreymass = Array{Float64}(undef,length(preds));
for i=1:length(preds)
    pref = prefdata[!,:actual][findall(x->x==preds[i],prefdata[!,:pred])];
    preyi = prefdata[!,:prey][findall(x->x==preds[i],prefdata[!,:pred])]
    predmass[i] = massvechay[findall(x->x==preds[i],spnameshay)[1]];
    #Average non-zero entries
    nonzeromean = mean(pref[findall(!iszero,pref)]);
    pref[findall(iszero,pref)].=nonzeromean;
    #create normalized pref index
    prefnorm = pref./ sum(pref) ; #
   
    preymass = Array{Float64}(undef,length(pref));
    for j=1:length(pref)
        preymass[j] = massvechay[findall(x->x==preyi[j],spnameshay)[1]];
    end
    maxpreymass[i] = maximum(preymass);
    prefmassmean[i] = sum(preymass .* prefnorm);
end

R"""
plot($predmass,$prefmassmean,pch=16,col='black',log='xy',xlim=c(1,300),ylim=c(1,1000),xlab='Predator mass (kg)',ylab='Mean preferred prey mass (kg)')
lines(seq(0,1000),seq(0,1000))
"""
R"""
model = lm(log($prefmassmean) ~ log($predmass))
summary(model)
"""
R"""
plot(log($predmass),log($prefmassmean),pch=16,col='black',xlim=c(1,6),ylim=c(1,6),xlab='Predator mass (kg)',ylab='Mean preferred prey mass (kg)')
abline(model)
lines(seq(0,1000),seq(0,1000),lty=2)
"""

R"""
model2 = lm(log($maxpreymass) ~ log($predmass))
summary(model2)
"""
R"""
plot(log($predmass),log($maxpreymass),pch=16,col='black',xlim=c(1,11),ylim=c(1,11),xlab='Predator mass (kg)',ylab='Max prey mass (kg)')
abline(model2)
lines(seq(0,1000),seq(0,1000),lty=2)
"""



#Caculcate mean preferred mass per predator
preds = unique(prefdata[!,:pred]);
#SIMULATION FIT
function genpredprey()
preyinds = Array{Float64}(undef,0);
predinds = Array{Float64}(undef,0);
for i=1:length(preds)
    pref = prefdata[!,:actual][findall(x->x==preds[i],prefdata[!,:pred])];
    #Average non-zero entries
    nonzeromean = mean(pref[findall(!iszero,pref)]);
    pref[findall(iszero,pref)].=nonzeromean;
    pref = round.((pref ./ sum(pref))*100);
    preyind = sum(pref);

    meanpredmass = massvechay[findall(x->x==preds[i],spnameshay)[1]];
    predmassSD = 0.1*meanpredmass;
    predbodysizedist = Normal(meanpredmass,predmassSD);
    predinds_draw = rand(predbodysizedist,Int64(preyind));
    predinds = [predinds; predinds_draw];
    preyi = prefdata[!,:prey][findall(x->x==preds[i],prefdata[!,:pred])];
    
    for j=1:length(preyi)
        #draw body masses
        meanmass = massvechay[findall(x->x==preyi[j],spnameshay)[1]];
        massSD = 0.1*meanmass;
        preybodysizedist = Normal(meanmass,massSD);
        numbers = Int64(pref[j]);
        preyinds_draw = rand(preybodysizedist,numbers);
        preyinds = [preyinds; preyinds_draw];
    end
end
return(predinds,preyinds)
end
predinds,preyinds = genpredprey();
R"""
model3 = lm(log($preyinds) ~ log($predinds))
summary(model3)
CI = confint(model3,level=0.99)
int99 = CI[[3]];
slope99 = CI[[4]];
"""
namespace = string(homedir(),"/Dropbox/PostDoc/2021_TaranWebs/figures/fig_predpreyratioanalysis.pdf")
R"""
pdf($namespace,height=5,width=6)
plot(log($predinds),log($preyinds),pch=16,col='black',xlim=c(1,6),ylim=c(1,6),xlab='Predator mass (kg)',ylab='Mean preferred prey mass (kg)')
abline(model3,col='blue')
lines(seq(0,1000),seq(0,1000),lty=2)
lines(seq(1,10),int99+slope99*seq(1,10),col='red')
dev.off()
"""



#FULL HAYWARD DATASET
haywardfulldata = CSV.read("$(homedir())/Dropbox/PostDoc/2021_TaranWebs/data/data_hayward_all.csv",header=true,DataFrame);
# haypredmass = haywardfulldata[!,:Predbodymaskg];
# haypreymass = haywardfulldata[!,:Preybodymasskg];
# haypercent = haywardfulldata[!,:PercentOfKills];

# "Panthera leo"
# "Crocuta crocuta"
# "Panthera pardus"
# "Cuon alpinus"
# "Lycaon pictus"
# "Acinonyx jubatus"
# "Panthera tigris"
groupsize = [4.1, 2.36, 1, 1, 4, 1, 1]

#Caculcate mean preferred mass per predator
preds = unique(haywardfulldata[!,:Predator]);
# preds = preds[[1,2,3,4,6]]
#SIMULATION FIT
function genpredprey()
preyinds = Array{Float64}(undef,0);
predinds = Array{Float64}(undef,0);
for i=1:length(preds)
    pref = haywardfulldata[!,:PercentOfKills][findall(x->x==preds[i],haywardfulldata[!,:Predator])];
    pref[findall(isnan,pref)].=0.;
    #Average non-zero entries
    # nonzeromean = mean(pref[findall(!iszero,pref)]);
    # pref[findall(iszero,pref)].=nonzeromean;
    # pref = round.((pref)*10000);
    pref = round.((pref ./ sum(pref))*10000);
    preyind = sum(pref);

    meanpredmass = mean(haywardfulldata[!,:Predbodymasskg][findall(x->x==preds[i],haywardfulldata[!,:Predator])]);
    # meanpredmass *= groupsize[i];
    predmassSD = 0.2*meanpredmass;
    predbodysizedist = Normal(meanpredmass,predmassSD);
    predinds_draw = rand(predbodysizedist,Int64(preyind));
    predinds = [predinds; predinds_draw];
    preyi = haywardfulldata[!,:Prey][findall(x->x==preds[i],haywardfulldata[!,:Predator])];
    
    for j=1:length(preyi)
        if pref[j] > 0.
            #draw body masses
            #Preybodymasskg
            #Preybodymasskg34adultfemalemass
            meanmass = haywardfulldata[!,:Preybodymasskg34adultfemalemass][findall(x->x==preds[i],haywardfulldata[!,:Predator])][j];
            massSD = 0.2*meanmass;
            preybodysizedist = Normal(meanmass,massSD);
            numbers = Int64(pref[j]);
            preyinds_draw = rand(preybodysizedist,numbers);
            preyinds = [preyinds; preyinds_draw];
        end
    end
end
return(predinds,preyinds)
end
predinds,preyinds = genpredprey();
R"""
model3 = lm(log($preyinds) ~ log($predinds))
summary(model3)
CI = confint(model3,level=0.99)
int99 = CI[[3]];
slope99 = CI[[4]];
"""
namespace = string(homedir(),"/Dropbox/PostDoc/2021_TaranWebs/figures/fig_predpreyratioanalysis_alldata34.pdf")
R"""
pdf($namespace,height=5,width=6)
plot(log($predinds),log($preyinds),pch=16,col='black',xlim=c(0.1,8),ylim=c(0.1,8),xlab='Predator mass (kg)',ylab='Mean preferred prey mass (kg)',cex=0.5)
abline(model3,col='blue')
lines(seq(0,1000),seq(0,1000),lty=2)
lines(seq(1,10),int99+slope99*seq(1,10),col='red')
dev.off()
"""

sizebins = 100;
maxprey = log10(maximum(preyinds));
minprey = log10(100); #Run this only for prey = 100 KG to max KG
# minprey = log10(minimum(preyinds));
stepsize = (maxprey-minprey)/sizebins;
preysizeclass = collect(minprey:stepsize:maxprey);
meanpredsize = Array{Float64}(undef,(length(preysizeclass)-1));
for i=2:length(preysizeclass)
    preysizeinds_pos = findall(x->((x>10^preysizeclass[i-1]) && (x < 10^preysizeclass[i])), preyinds);
    if length(preysizeinds_pos) > 0
        meanpredsize[i-1] = mean(predinds[preysizeinds_pos]);
    else 
        meanpredsize[i-1] =NaN;
    end
end
filledspots = findall(!isnan,meanpredsize)
scatterplot(10 .^preysizeclass[filledspots],meanpredsize[filledspots])
R"""
    model4 = lm(log($(meanpredsize[filledspots])) ~ log($(10 .^preysizeclass[filledspots])))
    summary(model4)
    fitintercept = model4[[1]][[1]];
    fitslope = model4[[1]][[2]];
    CI = confint(model4,level=0.95)
    intlow = CI[[1]];
    slopelow = CI[[2]];
    inthigh = CI[[3]];
    slopehigh = CI[[4]];
"""
namespace = string(homedir(),"/Dropbox/PostDoc/2021_TaranWebs/figures/fig_meanpredsize34.pdf")
R"""
pdf($namespace,height=5,width=6)
plot(log($(10 .^preysizeclass[filledspots])),log($(meanpredsize[filledspots])),pch=16,col='black',xlim=c(4.5,8.5),ylim=c(4.5,6),xlab='Prey mass (kg)',ylab='Mean predator mass (kg)',cex=0.5)
lines(seq(0,3000),seq(0,3000),lty=2)
abline(model4,col='blue')
lines(seq(1,10),inthigh+slopehigh*seq(1,10),col='red')
dev.off()
"""
fitintercept = @rget fitintercept;
fitslope = @rget fitslope;
fitinterceptlow = @rget intlow;
fitslopelow = @rget slopelow;
fitintercepthigh = @rget inthigh;
fitslopehigh = @rget slopehigh;


#Export for the mathematica notebook
fit_table = DataFrame([fitintercept fitinterceptlow fitintercepthigh; fitslope fitslopelow fitslopehigh]);
rename!(fit_table,[:Fit,:FitLow,:FitHigh])
CSV.write("$(homedir())/Dropbox/PostDoc/2021_TaranWebs/data/ppmr_fit_table.csv",fit_table; header=true);

## CALCULATE THE OPPOSITE ACROSS PREDATOR CLASSES




# CONVERT TO GRAMS TO MATCH MATHEMATICA notebook
fitslopegrams = copy(fitslope);
fitinterceptgrams = exp(fitintercept)*((1/1000)^fitslope)*1000;

namespace = string(homedir(),"/Dropbox/PostDoc/2021_TaranWebs/figures/fig_meanpredsize34_log.pdf")
R"""
pdf($namespace,height=5,width=6)
plot($((10 .^preysizeclass[filledspots]).*1000),$(meanpredsize[filledspots].*1000),pch=16,col='black',xlim=c(exp(4.5)*1000,exp(10)*1000),ylim=c(exp(4.5)*1000,exp(7)*1000),xlab='Prey mass (g)',ylab='Expected predator mass (g)',cex=0.5,log='xy')
lines(seq(0,3000000),seq(0,3000000),lty=2)
lines($(10 .^preysizeclass[filledspots])*1000, ($fitinterceptgrams)*($((10 .^preysizeclass[filledspots]).*1000))^$fitslopegrams,col='blue')
chi1 = -0.97
chi2 = 1.5
lines($(collect(100:20000))*1000, ($fitinterceptgrams*(1+chi1))*($(collect(100:20000))*1000)^($fitslopegrams*(1+chi2)),col='red')
# chi1b = 1.37
# chi2b = 0.28
# lines($(collect(100:20000)*1000), ($fitinterceptgrams*(1+chi1b))*($(collect(100:20000))*1000)^($fitslopegrams*(1+chi2b)),col='green')
dev.off()
"""
@rget chi1; @rget chi2
fitinterceptgrams*(1+chi1)
fitslopegrams*(1+chi2)