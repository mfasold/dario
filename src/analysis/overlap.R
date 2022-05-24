# Compute "expresssion" - the overlap of reads and a set of annotations by means of genomeIntervals
library(genomeIntervals);

args <- commandArgs(TRUE)
ncRNAs.annotations = args[1]
exons.annotations = args[2]
introns.annotations = args[3]
upload.reads = args[4]
ncRNAs.output = args[5]
unknown.output = args[6]
overlap.info = args[7]
rundir = args[8]

reads <- read.table(upload.reads, stringsAsFactors=F)
reads.ivals <- new("Genome_intervals_stranded", as.matrix(reads[,2:3]), closed=T, annotation=data.frame(seq_name=reads[,1], inter_base=F, strand=factor(reads[,6], levels=c("+", "-"))))

exons <- read.table(exons.annotations, stringsAsFactors=F)
ivals <- new("Genome_intervals_stranded", as.matrix(exons[,2:3]), closed=T, annotation=data.frame(seq_name=exons[,1], inter_base=F, strand=factor(exons[,6], levels=c("+", "-")), type="exon"))

introns <- read.table(introns.annotations, stringsAsFactors=F)
ivals <- c(ivals, new("Genome_intervals_stranded", as.matrix(introns[,2:3]), closed=T, annotation=data.frame(seq_name=introns[,1], inter_base=F, strand=factor(introns[,6], levels=c("+", "-")), type="intron")))

ncRNAs <- read.table(ncRNAs.annotations, stringsAsFactors=F)
ivals <- c(ivals, new("Genome_intervals_stranded", as.matrix(ncRNAs[,2:3]), closed=T, annotation=data.frame(seq_name=ncRNAs[,1], inter_base=F, strand=factor(ncRNAs[,6], levels=c("+", "-")), type=ncRNAs[,7])))

n <- interval_overlap(reads.ivals, ivals)
n_idx <- which(unlist(lapply(n, length)) > 0)

ncRNAs.labels <- unique(ncRNAs[,7])
ncRNAs.cnt <- rep(0, length(ncRNAs.labels))
ncRNAs.exp <- rep(0, length(ncRNAs.labels))
ncRNA.cnt <- 0
ncRNA.exp <- 0
exon.cnt <- 0
exon.exp <- 0
intron.cnt <- 0
intron.exp <- 0

ncRNA.out <- NULL
exon.out <- NULL
unknown.out <- NULL

tmp <- unlist(sapply(n_idx, function(i){
	types <- ivals@annotation$type[n[[i]]]
	if (!all(grepl("exon", types) | grepl("intron", types))){
		ncRNA.out <<- c(ncRNA.out, i)
		ncRNA.cnt <<- ncRNA.cnt + reads[i,7]
		ncRNA.exp <<- ncRNA.exp + reads[i,5]
		for (j in 1:length(ncRNAs.labels)){
			if (any(grepl(ncRNAs.labels[j], types))){
				ncRNAs.cnt[j] <<- ncRNAs.cnt[j] + reads[i,7]
				ncRNAs.exp[j] <<- ncRNAs.exp[j] + reads[i,5]
			}
		}
	}
	else {
		if (any(grepl("intron", types))){
			unknown.out <<- c(unknown.out, i)
			intron.cnt <<- intron.cnt + reads[i,7]
			intron.exp <<- intron.exp + reads[i,5]	
		}
		else {
			exon.out <<- c(exon.out, i)
			exon.cnt <<- exon.cnt + reads[i,7]
			exon.exp <<- exon.exp + reads[i,5]
		}
	}
}))


write.table(reads[ncRNA.out,], ncRNAs.output,row.names=F,col.names=F,append=T,quote=F,sep="\t")
write.table(reads[unknown.out,], unknown.output,row.names=F,col.names=F,append=T,quote=F,sep="\t")

n_idx <- which(unlist(lapply(n, length)) == 0)
write.table(reads[n_idx,], unknown.output,row.names=F,col.names=F,append=T,sep="\t")
intergenic.cnt <- length(n_idx)
intergenic.exp <- sum(reads[n_idx,5])

sink(overlap.info)
cat("ncRNAs",ncRNA.cnt,ncRNA.exp,prettyNum(ncRNA.cnt, big.mark = ",", decimal.mark="."),prettyNum(ncRNA.exp, big.mark = ",", decimal.mark="."), sep=":")
cat("\n")
cat("exons",exon.cnt,exon.exp,prettyNum(exon.cnt, big.mark = ",", decimal.mark="."),prettyNum(exon.exp, big.mark = ",", decimal.mark="."), sep=":")
cat("\n")
cat("introns",intron.cnt,intron.exp,prettyNum(intron.cnt, big.mark = ",", decimal.mark="."),prettyNum(intron.exp, big.mark = ",", decimal.mark="."), sep=":")
cat("\n")
for (j in 1:length(ncRNAs.labels)){
	cat(ncRNAs.labels[j],ncRNAs.cnt[j],ncRNAs.exp[j],prettyNum(ncRNAs.cnt[j], big.mark = ",", decimal.mark="."),prettyNum(ncRNAs.exp[j], big.mark = ",", decimal.mark="."), sep=":")
	cat("\n")
}
cat("intergenic",intergenic.cnt,intergenic.exp,prettyNum(intergenic.cnt, big.mark = ",", decimal.mark="."),prettyNum(intergenic.exp, big.mark = ",", decimal.mark="."), sep=":")
sink()
