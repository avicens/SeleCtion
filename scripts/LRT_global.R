###Selection Test####
###Based on  global (whole genome) parameters###
###Alberto Vicens###

mutrates = sapply(substmodel[,1], function(x) prod(parmle[base::strsplit(x,split="\\*")[[1]]])) # Expected rate per available site
numrates = length(mutrates)

# a. Neutral model: wmis==1, wnon==1, wspl==1
mrfold = sum(N)/sum(L) # Correction factor of "t" based on the obs/exp ratio of "neutral" mutations under the model
ll0 = sum(dpois(N, lambda = L*mutrates*mrfold*t(array(c(1,1,1,1),dim=c(4,numrates))))) # lik null model

# b. Missense model: wmis==1, free wnon, free wspl
mrfold = max(1e-10, sum(N[,c(1,2)])/sum(L[,c(1,2)])) # Correction factor of "t" based on the obs/exp ratio of "neutral" mutations under the model
wfree = N[,3:4]/L[,3:4]/mrfold; wfree[y[,3:4]==0] = 0
llmis = sum(dpois(N, lambda= L*mutrates*mrfold*t(array(c(1,1,wfree),dim=c(4,numrates))))) # lik free wmis

# c. free wmis, wnon and wspl
mrfold = sum(N[,1])/sum(L[,1]) # Correction factor of "t"
w = N[,2:4]/L[,2:4]/mrfold; w[N[,2:4]==0] = 0 # MLE of dN/dS based on the local rate (using syn muts as neutral)
llall = sum(dpois(x=N, lambda= L*mutrates*mrfold*t(array(c(1,w),dim=c(4,numrates))))) # lik free wmis, wnon, wspl
w[w>1e4] = 1e4

#Likelihood Ratio Test and p-value
lrt = pchisq(2*(llall-c(llmis,ll0)),df=c(1,3))
names(lrt) <- c("lall-lmis", "lall-l0")
p = 1-lrt
names(p) <- names(lrt)