#Open Genotype file
gt<-read.table("Dropbox/phylogenomics_lab_dbx/sc_selection/W55/W55.all.GT.FORMAT",header=TRUE)

#For each mutation, count the number of cells carrying such mutations
mut_cells<-apply(gt[,3:ncol(gt)], 1,function(x) {length(which(x=="0/1" | x =="1/1"))})
gt2<-cbind(gt,mut_cells) #Append it as a new column

plot(table(mut_cells),xlab="Number of cells", ylab="Number of mutations", main= "Wang 2014 - ERBC - Single Cell")

#Open mutations file
mutations<-read.table("/home/user/Dropbox/phylogenomics_lab_dbx/sc_selection/W55/W55.all.mutlist",header = T)

#Extract mutations for different groups
for (i in unique(mut_cells)) {
  positions <- gt2[gt2$mut_cells==i,1:2]
  matched<-which(paste(mutations$chr,"_",mutations$pos) %in% paste(positions$CHROM,"_",positions$POS))
  write.table(mutations[matched,],file=paste("Dropbox/phylogenomics_lab_dbx/sc_selection/W55/muts/W55",i,"muts",sep="_"),quote =  F, row.names = F)
}
  

  -------------------------------------------------------------
  #Genotype file (Removing genotypes with <5 reads)
  gt_recode<-read.table("Dropbox/phylogenomics_lab_dbx/sc_selection/W55/W55.all.recode.GT.FORMAT",header=TRUE)
  mut_cells_recode<-apply(gt_recode[,3:ncol(gt)], 1,function(x) {length(which(x=="0/1" | x =="1/1"))})
  gt2_rec<-cbind(gt_recode,mut_cells_recode) #Append it as a new column
  
   #Extract mutations for different groups
  mut_cells_rec<-mut_cells_recode[mut_cells_recode!=0]
  for (i in unique(mut_cells_rec)) {
    positions <- gt2_rec[gt2_rec$mut_cells_recode==i,1:2]
    matched<-which(paste(mutations$chr,"_",mutations$pos) %in% paste(positions$CHROM,"_",positions$POS))
    write.table(mutations[matched,],file=paste("Dropbox/phylogenomics_lab_dbx/sc_selection/W55/muts2/W55",i,"muts",sep="_"),quote =  F, row.names = F)
  }
  