library(dndscv)

# [Input] Reference database
data("refcds_hg19", package="dndscv")

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

#Retrieving genomic ranges 
ind = setNames(1:length(RefCDS), sapply(RefCDS,function(x) x$gene_name))
gr_genes_ind = ind[gr_genes$names] #GRanges object with gene ranges

#Coordinates file
mutations<-read.table("Dropbox/phylogenomics_lab_dbx/sc_selection/W55/infer_ancestral_states/W55.all.coding.mut.coordinates")
colnames(mutations) = c("chr","pos")

# Mapping mutations to genes
gr_muts = GRanges(mutations$chr, IRanges(mutations$pos,mutations$pos))
ol = as.matrix(findOverlaps(gr_muts, gr_genes, type="any", select="first"))
mutations$geneind = gr_genes_ind[ol[,1]]
mutations$gene = sapply(RefCDS,function(x) x$gene_name)[mutations$geneind]

#Export a file of gene names
gn<-mutations$gene
write(gn, "Dropbox/phylogenomics_lab_dbx/sc_selection/W55/infer_ancestral_states/W55.all.coding.genenames.txt")
