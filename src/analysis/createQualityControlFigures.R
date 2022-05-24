#! /usr/bin/Rscript --no-environ
##
## @todo
##  * define colGroups for read files -> speed up
## 

args <- commandArgs(TRUE)

# Set image parameters
imageSizes = list( smallImage = list(x = 420, y = 300, res = 100, ext = "small"),
                   bigImage   = list(x = 1000, y = 800, res = 150, ext = "big"))

# Graphical settings
library(ggplot2)
library(RColorBrewer)
library(colorspace)

getColors <- function(n, variant=1) {
  schema = c("#586B7A", "#BEAF85", "#BE9585", "#48535C", "#1D394F", "#8E866E", "#7B672B", "#8E776E", "#7B422B", "#820F00", "#024250", "#497800")

  if (n <= length(schema)) {
    return(schema[1:n])
  } else {
  return(rainbow_hcl(n,c=60,l=75))
  }
}

# personal shortcuts
do <- function() dev.off()
l <- `length`
spaste <- function(...) paste(..., sep = "/")
dpaste <- function(...) paste(..., sep = ".")
npaste <- function(...) paste(..., sep = "")


# Input and output directories
dir = Sys.glob(args[1])
outdir = dir

## ------------------------------------------------------------------------------------------
## Read Length histogram 
## ------------------------------------------------------------------------------------------
dl = read.table(spaste(dir, "length.out"), colClasses=c("character", "numeric"))
dl = read.table(spaste(dir, "length.out") ,as.is=T)
colnames(dl) = c("ReadLength", "Frequency")

for (imgSize in imageSizes) {
  png(spaste(outdir, dpaste("ReadLengthHistogram", imgSize$ext, "png")), imgSize$x, imgSize$y, res=imgSize$res)
  p = ggplot(dl, aes(x=ReadLength, y=Frequency)) + geom_bar(stat="identity")
  print(p)
  do()
}

## ------------------------------------------------------------------------------------------
## Multiple Mappings histogram
## ------------------------------------------------------------------------------------------
dl = read.table(spaste(dir, "multipleMappings.out"), colClasses=c("numeric", "numeric"))
colnames(dl) = c("id","NumberOfMappings")

for (imgSize in imageSizes) {
  png(spaste(outdir, dpaste("MultipleMappingsHistogram", imgSize$ext, "png")), imgSize$x, imgSize$y, res=imgSize$res)
  p = qplot(dl$id, dl$NumberOfMappings, size = I(2)) + ylab("Frequency") + xlab("# of Mappings per Read") + scale_x_log10() + scale_y_log10() 
  print(p)
  do()
}


## ------------------------------------------------------------------------------------------
## PIE CHARTS
## ------------------------------------------------------------------------------------------
myTheme <- theme(
                 panel.grid  = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 axis.line = element_blank(),
                 axis.text.y = element_blank(),
                 axis.text.x = element_blank(),
                 axis.ticks = element_blank()
                 )

dl = read.table(spaste(dir, "reads.info"), row.names=1, sep=":")

normalizedReads.idx = 2
is.not.ncrna = rownames(dl) %in% c("ncRNAs", "exons" , "introns", "intergenic")
slices = dl[is.not.ncrna, normalizedReads.idx]
labels = rownames(dl)[is.not.ncrna]
percentages <- round(slices/sum(slices)*100, digits = 1)
labels = npaste(labels, " (", percentages, "%)")
mycols = getColors(length(labels)) 

for (imgSize in imageSizes) {
  png(spaste(outdir, dpaste("ReadMappingsOverview1", imgSize$ext, "png")), imgSize$x, imgSize$y, res=imgSize$res)
  df = data.frame(Region = labels, value = slices) 
  p = ggplot(df, aes(x = "", y = value, fill = Region))  + myTheme +
      geom_bar(width = 1, stat = "identity") +  scale_fill_manual(values = mycols) +
      coord_polar("y", start=0) + ylab(NULL) + xlab(NULL) 
  print(p)
  do()
}

# Remove intergenic
dl = dl[!is.not.ncrna,]
slices = dl[,normalizedReads.idx]
labels = rownames(dl)
percentages <- round(slices/sum(slices)*100, digits = 1)
labels = npaste(labels, " (", percentages, "%)")
mycols = getColors(length(labels)) 

for (imgSize in imageSizes) {
  png(spaste(outdir, dpaste("ReadMappingsOverview2", imgSize$ext, "png")), imgSize$x, imgSize$y, res=imgSize$res)
  df = data.frame(ncRNA = labels, value = slices) 
  p = ggplot(df, aes(x = "", y = value, fill = ncRNA)) + myTheme +
    geom_bar(width = 1, stat = "identity") +  scale_fill_manual(values = mycols) +
      coord_polar("y", start=0) + ylab(NULL) + xlab(NULL)
  print(p)
  do()
}


## ------------------------------------------------------------------------------------------
## Classifier PIE
## ------------------------------------------------------------------------------------------
if(file.exists(spaste(dir, "predictions.expression.bed"))) {
  p = try({ # add try clause: if fails, do not yield error, just do not create graphics
    
    ## Parse file
    lines = readLines(spaste(dir, "my.modelstat"))
    idx = which( lines == "=== Stratified cross-validation ===") ## line index of cross-validation results

    ## get number of Correctly/Incorrectly Classified Instances
    ## by splitting line in file at >=2 consecutive spaces
    validate.correct = as.numeric(strsplit(lines[idx+2], "[[:blank:]]{2,}", perl=TRUE)[[1]][2]) 
    validate.uncorrect = as.numeric(strsplit(lines[idx+3], "[[:blank:]]{2,}", perl=TRUE)[[1]][2])

    slices = c(validate.correct, validate.uncorrect)
    labels = c("Correctly Classified", "Incorrectly Classified")
    percentages <- round(slices/sum(slices)*100, digits = 1)
    labels = npaste(labels, " (", percentages, "%)")
    mycols = getColors(length(labels)) 

    for (imgSize in imageSizes) {
      png(spaste(outdir, dpaste("ClassificationOverview", imgSize$ext, "png")), imgSize$x, imgSize$y, res=imgSize$res)
      df = data.frame(Classification = labels, value = slices) 
      p = ggplot(df, aes(x = "", y = value, fill = Classification)) + myTheme +
        geom_bar(width = 1, stat = "identity") +  scale_fill_manual(values = mycols) +
          coord_polar("y", start=0) +  ylab(NULL) + xlab(NULL) 
      print(p)
      do()
    }
  })
  if(inherits(p, "try-error")) {
    cat("Could not create prediction statistics. \n")
  }
}
