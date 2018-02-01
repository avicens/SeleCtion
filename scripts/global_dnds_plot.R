library(ggplot2)

##Global dN/dS estimates
load(file ="phylogenomics_lab/sc_selection/data/W55_dnds.RData")
load(file ="phylogenomics_lab/sc_selection/data/X25_dnds.RData")

global_dnds<-rbind(w55_dNdS$globaldnds,x25_dNdS$globaldnds)
dataset<-rep(c("W55","X25"),each=5)
global_dnds <-cbind(dataset,global_dnds)
subset_global_dnds<-global_dnds[(global_dnds$name!= "wnon" &  global_dnds$name!= "wspl"),] #filter out nonsense and splicing mutations due to Inf estimates

ggplot(subset_global_dnds,aes(x=name,y=mle, colour=dataset)) + geom_errorbar(aes(ymin=cilow,ymax=cihigh),width=.2, 
                              position=position_dodge(.9)) + geom_point(position=position_dodge(0.9),size=2) + xlab("Mutation type") + ylab("Global dN/dS MLE [95% CI]") + ggtitle("Global dN/dS estimates - SCS data") +  scale_x_discrete(limit = c("wmis", "wtru","wall"),labels = c("Missense","Truncating", "All"))

##Clonal and Subclonal dN/dS estimates

#W55$globaldnds)
load(file ="phylogenomics_lab/sc_selection/data/W55_cl_dnds.RData")
load(file ="phylogenomics_lab/sc_selection/data/W55_subcl_dnds.RData")

subtype=rep(c("Clonal","Subclonal"),each=5)

w55_dnds=rbind(w55_cl_dNdS$globaldnds,w55_subcl_dNdS)
w55_dnds=cbind(subtype,w55_dnds)

#X25
load(file ="phylogenomics_lab/sc_selection/data/X25_cl_dnds.RData")
load(file ="phylogenomics_lab/sc_selection/data/W55_subcl_dnds.RData")
x25_dnds=rbind(x25_cl_dNdS$globaldnds,x25_subcl_dNdS$globaldnds)
x25_dnds=cbind(subtype,x25_dnds)

#Merge estimates form SC datasets
dataset<-rep(c("W55","X25"),each=10)

sc_global_dnds<-rbind(w55_dnds,x25_dnds)
sc_global_dnds<-cbind(dataset,sc_global_dnds)
sc_mis_dnds<-sc_global_dnds[(sc_global_dnds$name== "wmis"),]

