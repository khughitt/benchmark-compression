R/Python Data Compression Benchmark Results Summary
================
2020-01-11

``` r
library(tidyverse)
library(knitr)

opts_chunk$set(fig.width = 5.625, fig.height = 5.625, dpi = 192)
```

# Overview

``` r
res <- read_csv('/data/benchmark/compression/timings/all_timings.csv', col_types = cols())

# add a more useful "Method" column and reorder columns
res <- res %>%
  mutate(Method = sprintf("%s / %s (%s)", Library, `File Format`, Compression)) %>%
  select(Dataset, Method, everything())

# convert "Method" field to a factor with the desired order
method_levels <- res %>%
  arrange(`File Format`, desc(Compression), Library) %>%
  pull(Method) %>%
  unique()

res$Method <- factor(res$Method, levels = method_levels)
  
# finally, let's add a fields representing the average compression ratio, and the
# compression ratio to average i/o ratio
data_sizes <- res %>%
  filter(Library == 'pandas' & `File Format` == 'csv' & Compression == 'none') %>%
  select(dataset=Dataset, size=`File Size (MB)`)

res$natural_size <- data_sizes$size[match(res$Dataset, data_sizes$dataset)]

res <- res %>%
  mutate(`Compression Ratio` = natural_size / `File Size (MB)`) %>%
  mutate(`Compression to I/O Ratio` = `Compression Ratio` / (0.5 * (`Read Time (Secs)` + `Write Time (Secs)`)))

res %>%
  select(-Library, -Compression) %>%
  kable(digits = 2)
```

| Dataset   | Method                       | Language | File Format | Read Time (Secs) | Write Time (Secs) | File Size (MB) | natural\_size | Compression Ratio | Compression to I/O Ratio |
| :-------- | :--------------------------- | :------- | :---------- | ---------------: | ----------------: | -------------: | ------------: | ----------------: | -----------------------: |
| itraq     | pandas / csv (none)          | python   | csv         |             0.13 |              1.17 |          14.16 |         14.16 |              1.00 |                     1.54 |
| golub1999 | pandas / csv (none)          | python   | csv         |             0.04 |              0.31 |           2.04 |          2.04 |              1.00 |                     5.72 |
| pancan    | pandas / csv (none)          | python   | csv         |             9.43 |             20.22 |         226.34 |        226.34 |              1.00 |                     0.07 |
| yeast     | pandas / csv (none)          | python   | csv         |             0.08 |              0.72 |          10.07 |         10.07 |              1.00 |                     2.51 |
| itraq     | pandas / csv (gzip)          | python   | csv         |             0.19 |              3.28 |           5.81 |         14.16 |              2.44 |                     1.40 |
| golub1999 | pandas / csv (gzip)          | python   | csv         |             0.05 |              1.11 |           0.74 |          2.04 |              2.74 |                     4.73 |
| pancan    | pandas / csv (gzip)          | python   | csv         |            10.73 |             39.70 |          76.81 |        226.34 |              2.95 |                     0.12 |
| yeast     | pandas / csv (gzip)          | python   | csv         |             0.13 |              1.29 |           4.80 |         10.07 |              2.10 |                     2.95 |
| itraq     | pyarrow / feather (none)     | python   | feather     |             0.01 |              0.05 |           9.30 |         14.16 |              1.52 |                    57.76 |
| golub1999 | pyarrow / feather (none)     | python   | feather     |             0.01 |              0.01 |           3.96 |          2.04 |              0.52 |                    43.78 |
| pancan    | pyarrow / feather (none)     | python   | feather     |             0.21 |              1.44 |         133.13 |        226.34 |              1.70 |                     2.06 |
| yeast     | pyarrow / feather (none)     | python   | feather     |             0.00 |              0.01 |           4.54 |         10.07 |              2.22 |                   331.29 |
| itraq     | pyarrow / parquet (lz4)      | python   | parquet     |             0.02 |              0.06 |           3.50 |         14.16 |              4.05 |                   109.34 |
| golub1999 | pyarrow / parquet (lz4)      | python   | parquet     |             0.02 |              0.03 |           1.07 |          2.04 |              1.91 |                    79.81 |
| pancan    | pyarrow / parquet (lz4)      | python   | parquet     |             2.23 |              2.66 |         145.48 |        226.34 |              1.56 |                     0.64 |
| yeast     | pyarrow / parquet (lz4)      | python   | parquet     |             0.01 |              0.03 |           2.47 |         10.07 |              4.09 |                   201.30 |
| itraq     | fastparquet / parquet (lz4)  | python   | parquet     |             0.06 |              0.10 |           5.31 |         14.16 |              2.67 |                    35.23 |
| golub1999 | fastparquet / parquet (lz4)  | python   | parquet     |             0.03 |              0.09 |           1.53 |          2.04 |              1.33 |                    21.62 |
| pancan    | fastparquet / parquet (lz4)  | python   | parquet     |             6.01 |              8.68 |         121.36 |        226.34 |              1.87 |                     0.25 |
| yeast     | fastparquet / parquet (lz4)  | python   | parquet     |             0.02 |              0.05 |           3.01 |         10.07 |              3.35 |                    96.44 |
| itraq     | pyarrow / parquet (snappy)   | python   | parquet     |             0.02 |              0.06 |           3.49 |         14.16 |              4.06 |                   110.81 |
| golub1999 | pyarrow / parquet (snappy)   | python   | parquet     |             0.02 |              0.03 |           1.07 |          2.04 |              1.91 |                    76.29 |
| pancan    | pyarrow / parquet (snappy)   | python   | parquet     |             2.27 |              2.67 |         145.09 |        226.34 |              1.56 |                     0.63 |
| yeast     | pyarrow / parquet (snappy)   | python   | parquet     |             0.01 |              0.03 |           2.61 |         10.07 |              3.86 |                   200.58 |
| itraq     | pyarrow / parquet (zstd)     | python   | parquet     |             0.01 |              0.06 |           3.30 |         14.16 |              4.29 |                   109.02 |
| golub1999 | pyarrow / parquet (zstd)     | python   | parquet     |             0.02 |              0.04 |           0.85 |          2.04 |              2.41 |                    93.23 |
| pancan    | pyarrow / parquet (zstd)     | python   | parquet     |             2.23 |              3.34 |         137.31 |        226.34 |              1.65 |                     0.59 |
| yeast     | pyarrow / parquet (zstd)     | python   | parquet     |             0.01 |              0.04 |           2.34 |         10.07 |              4.30 |                   173.31 |
| itraq     | fastparquet / parquet (zstd) | python   | parquet     |             0.06 |              0.11 |           3.29 |         14.16 |              4.30 |                    48.51 |
| golub1999 | fastparquet / parquet (zstd) | python   | parquet     |             0.03 |              0.10 |           0.83 |          2.04 |              2.45 |                    36.27 |
| pancan    | fastparquet / parquet (zstd) | python   | parquet     |             5.81 |              9.43 |         112.67 |        226.34 |              2.01 |                     0.26 |
| yeast     | fastparquet / parquet (zstd) | python   | parquet     |             0.02 |              0.07 |           2.15 |         10.07 |              4.69 |                    98.69 |
| itraq     | readr / csv (none)           | r        | csv         |             0.57 |              0.20 |          12.42 |         14.16 |              1.14 |                     2.97 |
| golub1999 | readr / csv (none)           | r        | csv         |             0.06 |              0.09 |           2.04 |          2.04 |              1.00 |                    13.40 |
| pancan    | readr / csv (none)           | r        | csv         |            30.49 |              4.69 |         201.73 |        226.34 |              1.12 |                     0.06 |
| yeast     | readr / csv (none)           | r        | csv         |             0.37 |              0.11 |          10.03 |         10.07 |              1.00 |                     4.16 |
| itraq     | readr / csv (gzip)           | r        | csv         |             0.57 |              1.07 |           5.46 |         14.16 |              2.59 |                     3.16 |
| golub1999 | readr / csv (gzip)           | r        | csv         |             0.06 |              0.26 |           0.76 |          2.04 |              2.69 |                    16.35 |
| pancan    | readr / csv (gzip)           | r        | csv         |            30.28 |             15.51 |          72.07 |        226.34 |              3.14 |                     0.14 |
| yeast     | readr / csv (gzip)           | r        | csv         |             0.37 |              0.68 |           4.79 |         10.07 |              2.10 |                     4.00 |
| itraq     | r-feather / feather (none)   | r        | feather     |             0.01 |              0.01 |           9.30 |         14.16 |              1.52 |                   149.98 |
| golub1999 | r-feather / feather (none)   | r        | feather     |             0.01 |              0.01 |           3.96 |          2.04 |              0.52 |                    64.67 |
| pancan    | r-feather / feather (none)   | r        | feather     |             0.31 |              0.16 |         133.13 |        226.34 |              1.70 |                     7.27 |
| yeast     | r-feather / feather (none)   | r        | feather     |             0.00 |              0.01 |           4.54 |         10.07 |              2.22 |                   448.41 |
| itraq     | r-arrow / parquet (none)     | r        | parquet     |             0.01 |              0.05 |           3.93 |         14.16 |              3.60 |                   124.33 |
| golub1999 | r-arrow / parquet (none)     | r        | parquet     |             0.01 |              0.03 |           1.55 |          2.04 |              1.32 |                    66.69 |
| pancan    | r-arrow / parquet (none)     | r        | parquet     |             1.51 |              1.46 |         137.93 |        226.34 |              1.64 |                     1.10 |
| yeast     | r-arrow / parquet (none)     | r        | parquet     |             0.00 |              0.03 |           2.72 |         10.07 |              3.70 |                   239.57 |

``` r
# include original dataset sizes in plot sub-headers
res$Dataset <- sprintf("%s (%0.1f MB)", res$Dataset, res$natural_size)
```

# Read Time (seconds)

``` r
ggplot(res, aes(x = Method, y = `Read Time (Secs)`, fill = `File Format`)) +
  geom_bar(stat = 'identity') + 
  facet_wrap(~Dataset, scales = 'free_y') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Read Time (Lower is Better)") +
  ylab("Time (seconds)")
```

![](summarize_results_files/figure-gfm/read_time-1.png)<!-- -->

# Write Time (seconds)

``` r
ggplot(res, aes(x = Method, y = `Write Time (Secs)`, fill = `File Format`)) +
  geom_bar(stat = 'identity') + 
  facet_wrap(~Dataset, scales = 'free_y') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Write Time (Lower is Better)") +
  ylab("Time (seconds)")
```

![](summarize_results_files/figure-gfm/write_time-1.png)<!-- -->

# Filesize (MB)

``` r
ggplot(res, aes(x = Method, y = `File Size (MB)`, fill = `File Format`)) +
  geom_bar(stat = 'identity') + 
  facet_wrap(~Dataset, scales = 'free_y') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("File Size (Lower is Better)") +
  ylab("Size (MB)")
```

![](summarize_results_files/figure-gfm/file_size-1.png)<!-- -->

# Compression to Average I/O Speed Ratio

``` r
ggplot(res, aes(x = Method, y = `Compression to I/O Ratio`, fill = `File Format`)) +
  geom_bar(stat = 'identity') + 
  facet_wrap(~Dataset, scales = 'free_y') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Compression to Average I/O Speed Ratio (Higher is Better)") +
  ylab("Compression to Average I/O Speed Ratio")
```

![](summarize_results_files/figure-gfm/compression_to_speed_ratio-1.png)<!-- -->
