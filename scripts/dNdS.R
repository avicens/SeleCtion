library("GenomicRanges")
library("dndscv")

###Estimation of pN/pS for somatic mutations (based on dNdScv with default parameters)

##Function dNdS
dnds = function(mutations, refdb = "hg19", sm = "192r_3w") {

#1. Environment
  
  # [Input] Reference database
  if (refdb == "hg19") {
    data("refcds_hg19", package="dndscv")
  } else {
    load(refdb)
  }
  
  # [Input] Substitution model (The user can also input a custom substitution model as a matrix)
  if (length(sm)==1) {
    data(list=sprintf("submod_%s",sm), package="dndscv")
  } else {
    substmodel = sm
  }

   # Expanding the reference sequences [for faster access]
  for (j in 1:length(RefCDS)) {
    RefCDS[[j]]$seq_cds = base::strsplit(as.character(RefCDS[[j]]$seq_cds), split="")[[1]]
    RefCDS[[j]]$seq_cds1up = base::strsplit(as.character(RefCDS[[j]]$seq_cds1up), split="")[[1]]
    RefCDS[[j]]$seq_cds1down = base::strsplit(as.character(RefCDS[[j]]$seq_cds1down), split="")[[1]]
    if (!is.null(RefCDS[[j]]$seq_splice)) {
      RefCDS[[j]]$seq_splice = base::strsplit(as.character(RefCDS[[j]]$seq_splice), split="")[[1]]
      RefCDS[[j]]$seq_splice1up = base::strsplit(as.character(RefCDS[[j]]$seq_splice1up), split="")[[1]]
      RefCDS[[j]]$seq_splice1down = base::strsplit(as.character(RefCDS[[j]]$seq_splice1down), split="")[[1]]
    }
  }

#2. Mutation annotation

mutations[,c(1,2,4,5)] = lapply(mutations[,c(1,2,4,5)], as.character)  
colnames(mutations) = c("sampleID","chr","pos","ref","mut")
nt = c("A","C","G","T")
  
  #Trinucleotide substitutions
  trinucs = paste(rep(nt,each=16,times=1),rep(nt,each=4,times=4),rep(nt,each=1,times=16), sep="")
  trinucinds = setNames(1:64, trinucs)
  trinucsubs = NULL
  for (j in 1:length(trinucs)) {
    trinucsubs = c(trinucsubs, paste(trinucs[j], paste(substr(trinucs[j],1,1), setdiff(nt,substr(trinucs[j],2,2)), substr(trinucs[j],3,3), sep=""), sep=">"))
  }
  trinucsubsind = setNames(1:192, trinucsubs)
  
  #Retrieving genomic ranges 
  ind = setNames(1:length(RefCDS), sapply(RefCDS,function(x) x$gene_name))
  gr_genes_ind = ind[gr_genes$names] #GRanges object with gene ranges
  
  # Warning about possible unannotated dinucleotide substitutions
  if (any(diff(mutations$pos)==1)) {
    warning("Mutations observed in contiguous sites within a sample. Please annotate or remove dinucleotide or complex substitutions for best results.")
  }
  
  # Mapping mutations to genes
  gr_muts = GRanges(mutations$chr, IRanges(mutations$pos,mutations$pos))
  ol = as.matrix(findOverlaps(gr_muts, gr_genes, type="any", select="all"))
  mutations = mutations[ol[,1],] # Duplicating subs if they hit more than one gene
  mutations$geneind = gr_genes_ind[ol[,2]]
  mutations$gene = sapply(RefCDS,function(x) x$gene_name)[mutations$geneind]
  
  # Additional annotation of substitutions
  snv = (mutations$ref %in% nt & mutations$mut %in% nt)
  indels = mutations[!snv,]
  mutations = mutations[snv,]
  mutations$ref_cod = mutations$ref
  mutations$mut_cod = mutations$mut
  compnt = setNames(rev(nt), nt)
  
  mutations$strand = sapply(RefCDS,function(x) x$strand)[mutations$geneind]
  isminus = (mutations$strand==-1)
  mutations$ref_cod[isminus] = compnt[mutations$ref[isminus]]
  mutations$mut_cod[isminus] = compnt[mutations$mut[isminus]]
  
  # Initialising the matrices of observed mutations (N)
  for (j in 1:length(RefCDS)) {
    RefCDS[[j]]$N = array(0, dim=c(192,4)) 
  }
  
  # Subfunction: obtaining the codon positions of a coding mutation given the exon intervals
  
  chr2cds = function(pos,cds_int,strand) {
    if (strand==1) {
      return(which(pos==unlist(apply(cds_int, 1, function(x) x[1]:x[2]))))
    } else if (strand==-1) {
      return(which(pos==rev(unlist(apply(cds_int, 1, function(x) x[1]:x[2])))))
    }
  }
  
  # Annotating the functional impact of each substitution and populating the N matrices
  
  ref3_cod = mut3_cod = wrong_ref = aachange = ntchange = impact = array(NA, nrow(mutations))
  
  for (j in 1:nrow(mutations)) {
    
    geneind = mutations$geneind[j]
    pos = mutations$pos[j]
    if (any(pos == RefCDS[[geneind]]$intervals_splice)) { # Essential splice-site substitution
      
      impact[j] = "Essential_Splice"; impind = 4
      pos_ind = (pos==RefCDS[[geneind]]$intervals_splice)
      cdsnt = RefCDS[[geneind]]$seq_splice[pos_ind]
      ref3_cod[j] = sprintf("%s%s%s", RefCDS[[geneind]]$seq_splice1up[pos_ind], RefCDS[[geneind]]$seq_splice[pos_ind], RefCDS[[geneind]]$seq_splice1down[pos_ind])
      mut3_cod[j] = sprintf("%s%s%s", RefCDS[[geneind]]$seq_splice1up[pos_ind], mutations$mut_cod[j], RefCDS[[geneind]]$seq_splice1down[pos_ind])
      aachange[j] = ntchange[j] = "-"
      
    } else { # Coding substitution
      
      pos_ind = chr2cds(pos, RefCDS[[geneind]]$intervals_cds, RefCDS[[geneind]]$strand)
      cdsnt = RefCDS[[geneind]]$seq_cds[pos_ind]
      ref3_cod[j] = sprintf("%s%s%s", RefCDS[[geneind]]$seq_cds1up[pos_ind], RefCDS[[geneind]]$seq_cds[pos_ind], RefCDS[[geneind]]$seq_cds1down[pos_ind])
      mut3_cod[j] = sprintf("%s%s%s", RefCDS[[geneind]]$seq_cds1up[pos_ind], mutations$mut_cod[j], RefCDS[[geneind]]$seq_cds1down[pos_ind])
      codon_pos = c(ceiling(pos_ind/3)*3-2, ceiling(pos_ind/3)*3-1, ceiling(pos_ind/3)*3)
      old_codon = as.character(as.vector(RefCDS[[geneind]]$seq_cds[codon_pos]))
      pos_in_codon = pos_ind-(ceiling(pos_ind/3)-1)*3
      new_codon = old_codon; new_codon[pos_in_codon] = mutations$mut_cod[j]
      old_aa = seqinr::translate(old_codon)
      new_aa = seqinr::translate(new_codon)
      aachange[j] = sprintf('%s%0.0f%s',old_aa,ceiling(pos_ind/3),new_aa)
      ntchange[j] = sprintf('%s%0.0f%s',mutations$ref_cod[j],pos_ind,mutations$mut_cod[j])
      
      # Annotating the impact of the mutation
      if (new_aa == old_aa){ 
        impact[j] = "Synonymous"; impind = 1
      } else if (new_aa == "*"){
        impact[j] = "Nonsense"; impind = 3
      } else if (old_aa != "*"){
        impact[j] = "Missense"; impind = 2
      } else if (old_aa=="*") {
        impact[j] = "Stop_loss"; impind = NA
      }
    }
    
    if (mutations$ref_cod[j] != as.character(cdsnt)) { # Incorrect base annotation in the input mutation file (the mutation will be excluded with a warning)
      wrong_ref[j] = 1
    } else if (!is.na(impind)) { # Correct base annotation in the input mutation file
      trisub = trinucsubsind[ paste(ref3_cod[j], mut3_cod[j], sep=">") ]
      RefCDS[[geneind]]$N[trisub,impind] = RefCDS[[geneind]]$N[trisub,impind] + 1 # Adding the mutation to the N matrices
    }
    
    if (round(j/1e4)==(j/1e4)) { message(sprintf('    %0.3g %%...', round(j/nrow(mutations),2)*100)) }
  }
  
  mutations$ref3_cod = ref3_cod
  mutations$mut3_cod = mut3_cod
  mutations$aachange = aachange
  mutations$ntchange = ntchange
  mutations$impact = impact
  mutations$pid = sapply(RefCDS,function(x) x$protein_id)[mutations$geneind]
  
  if (any(!is.na(wrong_ref))) {
    stop(sprintf('%0.0f mutations have a wrong reference base, please correct and rerun.',sum(!is.na(wrong_ref)))) # This can be made into a mere warning and the rest of the code will work
    wrong_refbase = mutations[!is.na(wrong_ref),]
    mutations = mutations[is.na(wrong_ref),]
  }
  
  if (any(nrow(indels))) { # If there are indels we concatenate the tables of subs and indels
    indels = cbind(indels, data.frame(ref_cod=".", mut_cod=".", strand=".", ref3_cod=".", mut3_cod=".", aachange=".", ntchange=".", impact="no-SNV", pid=sapply(RefCDS,function(x) x$protein_id)[indels$geneind]))
    annot = rbind(mutations, indels)
  } else {
    annot = mutations
  }
  annot = annot[order(annot$sampleID, annot$chr, annot$pos),]
  
#3. Estimation of the global rates

Lall = array(sapply(RefCDS, function(x) x$L), dim=c(192,4,length(RefCDS)))
Nall = array(sapply(RefCDS, function(x) x$N), dim=c(192,4,length(RefCDS)))
L = apply(Lall, c(1,2), sum)
N = apply(Nall, c(1,2), sum)
  
# Subfunction: fitting substitution model

fit_substmodel = function(N, L, substmodel) {
  
  l = c(L); n = c(N); r = c(substmodel)
  n = n[l!=0]; r = r[l!=0]; l = l[l!=0] #Discard rate paramters where number of expected mutations is 0
  
  params = unique(base::strsplit(x=paste(r,collapse="*"), split="\\*")[[1]])
  indmat = as.data.frame(array(0, dim=c(length(r),length(params))))
  colnames(indmat) = params
  for (j in 1:length(r)) {
    indmat[j, base::strsplit(r[j], split="\\*")[[1]]] = 1
  }
  
  model = glm(formula = n ~ offset(log(l)) + . -1, data=indmat, family=poisson(link=log))
  mle = exp(coefficients(model)) # Maximum-likelihood estimates for the rate params
  ci = exp(confint.default(model)) # Confidence intervals
  par = data.frame(name=gsub("\`","",rownames(ci)), mle=mle[rownames(ci)], cilow=ci[,1], cihigh=ci[,2])
  return(list(par=par, model=model))
}

# Fitting all mutation rates and the 3 global selection parameters

poissout = fit_substmodel(N, L, substmodel) # Original substitution model
par = poissout$par
poissmodel = poissout$model
parmle =  setNames(par[,2], par[,1])
mle_submodel = par
rownames(mle_submodel) = NULL

# Fitting models with 1 and 2 global selection parameters

s1 = gsub("wmis","wall",gsub("wnon","wall",gsub("wspl","wall",substmodel)))
par1 = fit_substmodel(N, L, s1)$par # Substitution model with 1 selection parameter
s2 = gsub("wnon","wtru",gsub("wspl","wtru",substmodel))
par2 = fit_substmodel(N, L, s2)$par # Substitution model with 2 selection parameter
globaldnds = rbind(par, par1, par2)[c("wmis","wnon","wspl","wtru","wall"),]

#4. Likelihood Ratio Tests (based on dNdSloc)

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
p = 1-lrt
lrtout<-data.frame(cbind(lrt,p),row.names = c("lall-lmis", "lall-l0"))

pnpsout = list(globaldnds = globaldnds, annotmuts = annot, mle_submodel = mle_submodel, poissmodel = poissmodel, lrt = lrtout)

}
