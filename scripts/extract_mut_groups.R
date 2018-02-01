###global dNdS estimation in single cell###
##Alberto Vicens##

#########Inputs#########
#data frames with genotypes matrix (gt) and mutations (mut) for all SNVs to be analyzed
gt<-read.table("Dropbox/phylogenomics_lab_dbx/sc_selection/W55/W55.all.gt",header=TRUE)
mutations<-read.table("/home/user/Dropbox/phylogenomics_lab_dbx/sc_selection/W55/W55.all.mutlist",header = T)


###Function extract_mut_groups
#Group the mutations by cells on VCF data from Single Cell

extract_mut_groups <- function (gt,mut)  {
  
  
  mut_cells<-apply(gt[,3:ncol(gt)], 1,function(x) {length(which(x=="0/1" | x =="1/1"))})
  gt2<-cbind(gt,mut_cells) #Append it as a new column
  
  mut_groups.list<-list()
  
  for (i in unique(mut_cells)) {
    positions <- gt2[gt2$mut_cells==i,1:2]
    matched<-which(paste(mutations$chr,"_",mutations$pos) %in% paste(positions$CHROM,"_",positions$POS))
    mut<-mutations[matched,]
    mut.list[[i]]<-mut
    return(mut.list)
    }
  
}

#libraries
library("seqinr")
library("Biostrings")
library("MASS")
library("GenomicRanges")
library("dndscv")
