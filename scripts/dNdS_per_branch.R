###Collect mutations in each branch of a phylogenetic tree

#Open Genotype file (Colum names of genotypes must coincide with tip labels of the tree)
gt<-read.table("Dropbox/phylogenomics_lab_dbx/sc/W55/W55.all.gt",header=TRUE)

#Open mutations file
muts<-read.table("/home/user/Dropbox/phylogenomics_lab_dbx/sc/W55/W55.all.mutlist",header = T)

#Load phylogenetic tree
library(ape)

sctree<-read.tree(sctree)

#Extract clade and tip labels of this clade
cl70<-extract.clade(sctree,70)
tip70<-cl70$tip.label

#Retrieve mutations carried by the cells within the clade

