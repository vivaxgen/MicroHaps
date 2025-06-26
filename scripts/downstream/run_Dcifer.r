#!/usr/bin/env Rscript

library("argparse")
library(dcifer)
parser <- ArgumentParser(description = "Running DCIFER analysis on microhaplotype data")
parser$add_argument(
  "--infile",
  type = "character",
  required = TRUE,
  help = "Path to the filtered microhaplotype data from prepare_data (tsv)"
)
parser$add_argument(
  "--outdir",
  type = "character",
  required = TRUE,
  help = "Output directory for DCIFER results"
)

args <- parser$parse_args()
if (!dir.exists(args$outdir)) {
  dir.create(args$outdir, recursive = TRUE)
}

## Load all functions by running everything from line 5-92
# use case: alleles in dsmp are not present in the population allele 
pad.afreq <- function(afreq, dsmp) {
  # assumes the sample allele frequencies are more diverse
  smpafreq <- dsmp[[1]]
  
  # normalise number of locus in population
  locus.diff <- setdiff(names(smpafreq), names(afreq))
  afreq <- c(afreq, setNames(vector("list", length(locus.diff)), locus.diff))
  afreq <- afreq[order(names(afreq))]
  
  # normalise number of allele in population
  for (i in seq_along(smpafreq)) {
    allele.diff <- setdiff(names(smpafreq[[i]]), names(afreq[[i]]))
    
    # set the frequencies for the missing alleles
    small <- min(afreq[[i]]) / length(allele.diff)
    smalls <- setNames(rep(small, length(allele.diff)), allele.diff)
    
    afreq[[i]] <- c(afreq[[i]], smalls)
    afreq[[i]] <- afreq[[i]][order(names(afreq[[i]]))]
    
    # normalise frequencies to 1
    afreq[[i]] <- afreq[[i]] / sum(afreq[[i]])
  }
  
  return(afreq)
}


get.m1.estimate <- function(dres0) {
  dmat <- dres0[, , "estimate"]
  m1.estimate <- as.data.frame(as.table(dmat))
  m1.estimate <- m1.estimate[!is.na(m1.estimate[["Freq"]]), ]
  names(m1.estimate) <- c("sample_id2", "sample_id1", "M1")
  
  return(m1.estimate)
}


analyse.significantly.related.samples <-
  function(dsmp, coi, afreq, dres0, alpha = 0.05) {
    isig <- which(dres0[, , "p_value"] <= alpha, arr.ind = TRUE)[, 2:1]
    
    # if no significant relatedness found
    if (nrow(isig) == 0) {
      return(NULL)
    }
    
    sig2 <- vector("list", nrow(isig))
    for (i in 1:nrow(isig)) {
      sig2[[i]] <- ibdEstM(dsmp[isig[i, ]], coi[isig[i, ]], afreq, equalr = TRUE)
    }
    M2 <- sapply(sig2, length)
    rtotal2 <- sapply(sig2, sum)
    
    samples <- names(dsmp)
    sig <- data.frame(
      sample_id1 = samples[isig[, 1]],
      sample_id2 = samples[isig[, 2]],
      M = M2,
      rtotal = rtotal2
    )
    
    return(sig)
  }


calculate.overall.relatedness.estimate <- function(m1.estimate, sig, coi) {
  if (is.null(sig)) {
    mall.estimate <- cbind(m1.estimate, M = NA, rtotal = NA)
  } else {
    mall.estimate <- merge(m1.estimate, sig, all = TRUE)
  }
  
  coi.df <-
    data.frame(
      sample_id = sapply(strsplit(names(coi), ".", fixed = TRUE), "[[", 1),
      coi = coi,
      row.names = NULL
    )
  
  mall.estimate <-
    merge(mall.estimate, coi.df, by.x = "sample_id2", by.y = "sample_id")
  mall.estimate <-
    merge(
      mall.estimate,
      coi.df,
      by.x = "sample_id1",
      by.y = "sample_id",
      suffixes = c("2", "1")
    )
  
  mall.estimate[["scaled_r"]] <-
    mall.estimate[["rtotal"]] /
    (pmin(mall.estimate[["coi1"]], mall.estimate[["coi2"]]))
  
  mall.estimate[["relatedness"]] <- mall.estimate[["scaled_r"]]
  mall.estimate[is.na(mall.estimate[["relatedness"]]), "relatedness"] <-
    mall.estimate[is.na(mall.estimate[["relatedness"]]), "M1"]
  
  return(mall.estimate)
}


## Run the commands below step-by-step
## Change FIXME lines to the appropriate parameters

# FIXME: change to path of filtered microhaplotype data
sfile <- args$infile

dlong <- read.delim(sfile)
dsmp <- formatDat(dlong, svar = "sample_id", lvar = "locus", avar = "allele")
coi   <- getCOI(dsmp)
afreq <- calcAfreq(dsmp, coi)


## only one pair of strains between two infections can be related
dres0 <- ibdDat(dsmp, coi, afreq)

m1.estimate <- get.m1.estimate(dres0)


## allow multiple pairs of strains to be related between two infections
sig <- analyse.significantly.related.samples(dsmp, coi, afreq, dres0)


## combine results from restrained and unrestrained M
mall.estimate <- calculate.overall.relatedness.estimate(m1.estimate, sig, coi)


# FIXME: change to path of between-infection relatedness estimate
mall.estimate.file <- file.path(args$outdir, "mhap_between_relatedness_estimate.tsv")

write.table(
  mall.estimate,
  mall.estimate.file,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)
