#
# R parquet benchmark
#
library(arrow)
library(readr)
library(microbenchmark)

# load input data
df <- read_csv(snakemake@input[[1]], col_types = cols())

# benchmark write performance
write_times <- microbenchmark(
  write_parquet(df, snakemake@output[['data']]),
  times = snakemake@config$benchmark$num_times
)
min_write_time <- min(write_times$time) / 1e9

# benchmark read performance
read_times <- microbenchmark(
  read_parquet(snakemake@output[['data']]),
  times = snakemake@config$benchmark$num_times
)
min_read_time <- min(read_times$time) / 1e9

# get file size
fsize <- file.size(snakemake@output[['data']]) / 1e6

# save timings and size info
entry <- c(snakemake@wildcards[['dataset']], 'r', 'r-arrow', 'parquet', 
           'none', min_read_time, min_write_time, fsize)

entry <- paste0(entry, collapse = ', ')

fp <- file(snakemake@output[['timings']], 'w')
write(entry, file = fp)
close(fp)
