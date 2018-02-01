###Selection Test####
###Based on dNdSloc - estimation of mutation rates per gene, basing neutral mutation rate on synonymous mutations###
###Alberto Vicens###

genemuts = data.frame(gene_name = sapply(RefCDS, function(x) x$gene_name), n_syn=NA, n_mis=NA, n_non=NA, n_spl=NA, exp_syn=NA, exp_mis=NA, exp_non=NA, exp_spl=NA)
genemuts[,2:5] = t(sapply(RefCDS, function(x) colSums(x$N)))
mutrates = sapply(substmodel[,1], function(x) prod(parmle[base::strsplit(x,split="\\*")[[1]]])) # Expected rate per available site
genemuts[,6:9] = t(sapply(RefCDS, function(x) colSums(x$L*mutrates)))
numrates = length(mutrates)

y = apply(genemuts[,-1],2, as.numeric)

# a. Neutral model: wmis==1, wnon==1, wspl==1
mrfold = sum(y[,1:4])/sum(y[,5:8]) # Correction factor of "t" based on the obs/exp ratio of "neutral" mutations under the model
ll0 = sum(dpois(N, lambda = L*mutrates*mrfold*t(array(c(1,1,1,1),dim=c(4,numrates))))) # lik null model

# b. Missense model: wmis==1, free wnon, free wspl
mrfold = max(1e-10, sum(y[,c(1,2)])/sum(y[,c(5,6)])) # Correction factor of "t" based on the obs/exp ratio of "neutral" mutations under the model
wfree = y[,3:4]/y[,7:8]/mrfold; wfree[y[,3:4]==0] = 0
llmis = sum(dpois(N, lambda= L*mutrates*mrfold*t(array(c(1,1,wfree),dim=c(4,numrates))))) # lik free wmis

# c. free wmis, wnon and wspl
mrfold = max(1e-10, y[,1]/y[,5]) # Correction factor of "t"
w = y[,2:4]/y[,6:8]/mrfold; w[y[,2:4]==0] = 0 # MLE of dN/dS based on the local rate (using syn muts as neutral)
llall = sum(dpois(x=N, lambda= L*mutrates*mrfold*t(array(c(1,w),dim=c(4,numrates))))) # lik free wmis, wnon, wspl

#Likelihood Ratio Test and p-value
lrt = pchisq(2*(llall-c(llmis,ll0)),df=c(1,3))
names(lrt) <- c("lall-lmis", "lall-l0")
p = 1-lrt
names(p) <- names(lrt)
lrtout = as.data.frame(cbind(lrt,p)))