# Compute "expresssion" - the overlap of reads and a set of annotations by means of genomeIntervals
library(genomeIntervals);

args <- commandArgs(TRUE)

ncRNAs.file = args[1]
reads.file  = args[2]
output.file = args[3]
rundir = args[4]

reads <- read.table(reads.file, stringsAsFactors=F)
reads.ivals <- new("Genome_intervals_stranded", as.matrix(reads[,2:3]), closed=T, annotation=data.frame(seq_name=reads[,1], inter_base=F, strand=factor(reads[,6], levels=c("+", "-"))))

ncRNAs <- read.table(ncRNAs.file, stringsAsFactors=F)
ncRNAs.ivals <- new("Genome_intervals_stranded", as.matrix(ncRNAs[,2:3]), closed=T, annotation=data.frame(seq_name=ncRNAs[,1], inter_base=F, strand=factor(ncRNAs[,6], levels=c("+", "-"))))

n <- interval_overlap(reads.ivals, ncRNAs.ivals)
n_idx <- which(unlist(lapply(n, length)) > 0)

write.table(reads[n_idx,], output.file,row.names=F,col.names=F,append=T,quote=F,sep="\t")
