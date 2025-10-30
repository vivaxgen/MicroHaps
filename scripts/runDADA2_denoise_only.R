#!/bin/R env

library(dada2)
library(limma)
library(argparse)
library(data.table)

# Custom filtering, denoising parameters (if not default) can be provided as a separate config file?

parser <- ArgumentParser()
parser$add_argument("-p", "--path_to_meta", help="Path to input meta file listing fastqs (required)")
parser$add_argument("-d", "--dir", help="Working directory path for writing all dada2 output files")
parser$add_argument("-o", "--output_filename", help="output tab-separated filename (required)")
parser$add_argument("-s", "--save_run", help="save Run as R workspace image")
parser$add_argument("-ee", "--maxEE", type="integer",
                    help="Maximum expected errors for filtering forward and reverse read")
parser$add_argument("-tR", "--trimRight", type="integer",
                    help="Length for trimming from right for both forward and reverse reads")
parser$add_argument("-mL", "--minLen", type="integer",
                    help="minimum length required for reads on both end. Shorter reads are discarded")
parser$add_argument("-tQ", "--truncQ", type="integer",
                    help="truncate reads to first occurence of truncQ. All filtered reads have quality >= truncQ")
parser$add_argument("-mC", "--max_consist", type="integer",
                    help="Maximum cycles for error model until consistency. If no convergence, error values at max_consist cycle are used")
parser$add_argument("-wA", "--omega_a", type="double",
                    help="P-value threshold in sample inference for forming a new partition")
parser$add_argument("--bimera", action='store_true', help="Optionally output list of sequences identified as bimeras")
parser$add_argument("--threads", type="integer", default=1, help="Number of threads to use for parallel processing")
args <- parser$parse_args()

# Universal parameters
work_dir <- args$dir
path_to_meta <- args$path_to_meta
if (file.exists(path_to_meta)) {
  metafile <- fread(path_to_meta, sep = "\t", header=FALSE)
  sample.names <- metafile$V1
  fnMergeds <- metafile$V2
} else {
  stop(paste("metafile",path_to_meta,"not found!"))
}

# obtain/initialize Parameters
# (Universal) Parameters
randomize=TRUE
selfConsist=TRUE
filter = TRUE

# DADA2 and Filtering parameters
maxEE <- args$maxEE
trimRight <- args$trimRight
minLen <- args$minLen
truncQ <- args$truncQ
max_consist <- args$max_consist
omega_a <- args$omega_a

#Output parameters
if (dirname(args$output_filename) != ".") {
  output_filename <- args$output_filename
  } else {
    output_filename <- paste0(work_dir,"/",args$output_filename)
  }

#Datatable to summarize parmeters
parameter_df <- data.frame(maxEE=maxEE,
		trimRight=trimRight,
  		minLen=minLen,
		truncQ=truncQ,
		max_consist=max_consist,
		randomize=randomize,
		selfConsist=selfConsist,
		OMEGA_A=omega_a)

print(parameter_df)

if (length(fnMergeds) == 0) {
	stop("fastq files incomplete or not found")
}

# Plot Quality profiles before filering
png(paste0(work_dir,"/qualityMerged.png"), height = 800, width = 700)
try(print(plotQualityProfile(fnMergeds[1:2])), silent = TRUE)
dev.off()

filtMergeds <- file.path(work_dir, "filtered", paste0(sample.names, "_filt_Merged.fastq.gz"))
names(filtMergeds) <- sample.names

# Filter read
if (filter == TRUE) {
	print("filtering samples...")
	out <- filterAndTrim(fnMergeds, filtMergeds,
            maxN=0, maxEE=maxEE, trimRight=trimRight, truncQ=truncQ, minLen=minLen,
            rm.phix=TRUE, compress=TRUE, multithread=args$threads, verbose=TRUE)
	print("filtering done!")
} else {
	print("skipping filter except mandatory removal of N's... ")
	out <- filterAndTrim(fnMergeds, filtMergeds, truncQ=c(0), maxN=0, rm.phix=TRUE,
            compress=TRUE, multithread=args$threads, verbose=TRUE)
}

# Report and Correct for samples with zero reads after filter
zeros <- row.names(out)[out[,2] == 0]
write.table(zeros, paste0(work_dir,"/zeroReadSamples.txt"), sep = "\t", quote = FALSE)
filtMergeds <- filtMergeds[out[,2] != 0]
sample.names <- sample.names[out[,2] != 0]

# Update Out table
out <- out[(out[,2] != 0),]

#Compute the error model
print("starting error model learning for merged reads...")
errM <- learnErrors(filtMergeds, multithread=args$threads, verbose=2, randomize=randomize, MAX_CONSIST=max_consist)

#Plot the Errors
png(paste0(work_dir,"/errM.png"), height = 800, width = 700)
try(print(plotErrors(errM, nominalQ=TRUE)), silent = TRUE)
dev.off()

#DeReplicate the reads
derepMs <- derepFastq(filtMergeds, verbose = TRUE)
names(derepMs) <- sample.names

#Run core DADA2 algorithm
print("starting dada2 for merged reads...")
dadaMs <- dada(derepMs, err=errM, selfConsist=selfConsist, multithread=args$threads, verbose=TRUE, OMEGA_A=omega_a)

#Generate sequence table
print("generating sequence table...")
seqtab <- makeSequenceTable(dadaMs)
print("Number of sequences in table")
print(dim(seqtab))
# Inspect distribution of sequence lengths
print(table(nchar(getSequences(seqtab))))

#Remove Chimeras
if(args$bimera) {
  print("identifying bimeric sequences...")
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=args$threads, verbose=TRUE)
  print("Number of non-bimeric sequences:")
  print(dim(seqtab.nochim)[2])
  print("Percentage of reads which are non-bimeric:")
  print(sum(seqtab.nochim)/sum(seqtab))
  bimeras <- !(colnames(seqtab) %in% colnames(seqtab.nochim))
  write.table(data.frame(sequence = colnames(seqtab), bimera = bimeras), file=paste0(work_dir,"/ASVBimeras.txt"),
    quote=FALSE, sep="\t", row.names=FALSE)
} else {
  print("skipping Bimera identification..")
  seqtab.nochim <- seqtab
}

# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaMs, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedM", "nonchim")
rownames(track) <- sample.names

# sink summary from stdout to a file
sink(paste0(work_dir,"/reads_summary.txt"))
print(track)
#close sink
sink()


#Show the barplot of length distribution
png(paste0(work_dir,"/sequences_barplot.png"), height = 800, width = 700)
print(barplot(table(nchar(getSequences(seqtab)))))
dev.off()

#Generate output: sequence table to a tsv
write.table(seqtab, file=output_filename, quote = FALSE, sep = "\t")

# Save Run as R workspace image (Optional)
if (is.null(args$save_run)||args$save_run == '') {
    print("--save_run not found or empty. skip saving Rdata image to a file")
  } else {
    rm(fnMergeds,filtMergeds,derepMs,errM)
    save.image(paste0(work_dir,"/",args$save_run))
  }