#!/usr/bin/env Rscript

library("argparse")
parser <- ArgumentParser(description = "Running MOIRE analysis on microhaplotype data")
parser$add_argument(
  "--infile",
  type = "character",
  required = TRUE,
  help = "Path to the filtered microhaplotype data from prepare_data (tsv)"
)
parser$add_argument(
  "--burnin",
  type = "integer",
  default = 200,
  help = "Number of burn-in iterations (default: 200)"
)
parser$add_argument(
  "--samples_per_chain",
  type = "integer",
  default = 10000,
  help = "Number of samples per chain (default: 10000)"
)
parser$add_argument(
  "--pt_chains",
  type = "integer",
  default = 40,
  help = "Number of parallel tempering chains (default: 40)"
)
parser$add_argument(
  "--pt_num_threads",
  type = "integer",
  default = 10,
  help = "Number of threads for parallel tempering (default: 10)"
)
parser$add_argument(
  "--outdir",
  type = "character",
  required = TRUE,
  help = "Output directory for MCMC results"
)

args <- parser$parse_args()
if (!dir.exists(args$outdir)) {
  dir.create(args$outdir, recursive = TRUE)
}


plot.posterior.burnin.trace <- function(mcmc_results) {
  for (i in seq_along(mcmc_results[["chains"]])) {
    chain <- mcmc_results[["chains"]][[i]]
    
    plot(
      chain[["posterior_burnin"]],
      type = "l",
      main = paste0("Chain ", i),
      xlab = "Iteration",
      ylab = "Posterior (burn-in)"
    )
  }
}

plot.posterior.sample.trace <- function(mcmc_results) {
  for (i in seq_along(mcmc_results[["chains"]])) {
    chain <- mcmc_results[["chains"]][[i]]
    
    plot(
      chain[["posterior_sample"]],
      type = "l",
      main = paste0("Chain ", i),
      xlab = "Iteration",
      ylab = "Posterior (sample)"
    )
  }
}

plot.posterior.sample.density <- function(mcmc_results) {
  for (i in seq_along(mcmc_results[["chains"]])) {
    chain <- mcmc_results[["chains"]][[i]]
    
    plot(
      density(chain[["posterior_sample"]], bw = "SJ"),
      main = paste0("Chain ", i),
      ylab = "Posterior (sample) density"
    )
  }
}

plot.posterior.sample.acf <- function(mcmc_results) {
  for (i in seq_along(mcmc_results[["chains"]])) {
    chain <- mcmc_results[["chains"]][[i]]
    
    acf <-
      acf(
        chain[["posterior_sample"]],
        lag.max = length(chain[["posterior_sample"]]) / 5,
        plot = FALSE
      )
    
    plot(
      acf,
      main = paste0("Chain ", i),
      ylab = "Posterior (sample) ACF"
    )
  }
}


## Change FIXME lines to the appropriate parameters

## If you already have the mcmc_results saved as RDS,
## you don't need to re-do run_mcmc (use readRDS() instead)

# FIXME: change to path of filtered microhaplotype data
sfile <- args$infile

df <- read.delim(sfile)
data <- moire::load_long_form_data(df)


# these default parameters can be adjusted accordingly:
# burn-in = 200 [burnin], samples per chain = 10,000 [samples_per_chain]
# parallel tempering chains = 40 [pt_chains]
# adjust the number of threads [pt_num_threads] accordingly
## You might want to start with 1,000 samples per chain before committing to a larger value
burnin <- args$burnin
pt_num_threads <- args$pt_num_threads
pt_chains <- args$pt_chains
samples_per_chain <- args$samples_per_chain

mcmc_results <-
  moire::run_mcmc(
    data = data,
    is_missing = data[["is_missing"]],
    burnin = burnin,
    samples_per_chain = samples_per_chain,
    pt_chains = pt_chains,
    pt_num_threads = pt_num_threads
 )

# FIXME: change to path of MOIRE MCMC 
file <- file.path(args$outdir, "mhap_MOIRE.rds")

saveRDS(mcmc_results, file = file)


## produce MCMC diagnostics for assessment, see:
## https://www.statlect.com/fundamentals-of-statistics/Markov-Chain-Monte-Carlo-diagnostics
nchain <- length(mcmc_results[["chains"]])
m <- ceiling(sqrt(nchain))

# FIXME: change to path of posterior (burn-in) trace plot (adjusted)
posterior_burnin_trace_plot_file <- file.path(args$outdir, "mhap_MOIRE_posterior_burnin_trace.pdf")
pdf(file = posterior_burnin_trace_plot_file, width = 7, height = 7)

par(mfrow = c(m, m))
plot.posterior.burnin.trace(mcmc_results)
dev.off()

# FIXME: change to path of posterior (sample) trace plot (adjusted)
posterior_sample_trace_plot_file <- file.path(args$outdir, "mhap_MOIRE_posterior_sample_trace.pdf")
pdf(file = posterior_sample_trace_plot_file, width = 7, height = 7)

par(mfrow = c(m, m))
plot.posterior.sample.trace(mcmc_results)
dev.off()

# FIXME: change to path of posterior (sample) ACF plot (adjusted)
posterior_sample_acf_plot_file <- file.path(args$outdir, "mhap_MOIRE_posterior_sample_ACF.pdf")
pdf(file = posterior_sample_acf_plot_file, width = 7, height = 7)

par(mfrow = c(m, m))
plot.posterior.sample.acf(mcmc_results)
dev.off()

## STOP and determine if the MCMC parameters need to be re-adjusted
## https://www.statlect.com/fundamentals-of-statistics/Markov-Chain-Monte-Carlo-diagnostics
## section: Trace plots & Autocorrelation function (ACF) plots 


# FIXME: change to path of COI summary
coi_summary_file <- file.path(args$outdir, "mhap_COI_summary.tsv")

coi_summary <- moire::summarize_coi(mcmc_results)
write.table(
  coi_summary,
  coi_summary_file,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)


# FIXME: change to path of polyclonal status
polyclonal_status_file <- file.path(args$outdir, "mhap_polyclonal_status.tsv")

polyclonal_status <-
  data.frame(
    sample_id = coi_summary[["sample_id"]],
    is_polyclonal = coi_summary[["prob_polyclonal"]] >= 0.5
  )

write.table(
  polyclonal_status,
  polyclonal_status_file,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)


# FIXME: change to path of effective COI summary
effective_coi_summary_file <- file.path(args$outdir, "mhap_effective_COI_summary.tsv")

effective_coi_summary <- moire::summarize_effective_coi(mcmc_results)
write.table(
  effective_coi_summary,
  effective_coi_summary_file,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)
