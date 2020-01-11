#
# R readr csv benchmark
#
library(readr)
library(microbenchmark)

# benchmark read performance
read_times <- microbenchmark(
  read_csv(snakemake@input[[1]], col_types = cols()),
  times = snakemake@config$benchmark$num_times
)
min_read_time <- min(read_times$time) / 1e9

# benchmark write performance
df <- read_csv(snakemake@input[[1]], col_types = cols())

write_times <- microbenchmark(
  # exclude NA placeholders to be more consistent with pandas
  write_csv(df, snakemake@output[['data']], na = ""),
  times = snakemake@config$benchmark$num_times
)
min_write_time <- min(write_times$time) / 1e9

# get file size
fsize <- file.size(snakemake@output[['data']]) / 1e6

if (endsWith(snakemake@output[[1]], '.gz')) {
  compression <- 'gzip'
} else {
  compression <- 'none'
}

# save timings and size info
entry <- c(snakemake@wildcards[['dataset']], 'r', 'readr', 'read::read_csv', 
           compression, min_read_time, min_write_time, fsize)

entry <- paste0(entry, collapse = ', ')

fp <- file(snakemake@output[['timings']], 'w')
write(entry, file = fp)
close(fp)
