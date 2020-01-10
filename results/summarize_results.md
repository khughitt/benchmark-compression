R/Python Data Compression Benchmark Results Summary
================
2020-01-09

``` r
library(tidyverse)
library(knitr)

opts_chunk$set(fig.width = 5.626, fig.height = 4.25, dpi = 192)
```

# Overview

``` r
res <- read_csv('/data/benchmark/compression/timings/all_timings.csv', col_types = cols())

# drop the "Language" and "Method" columns for now since everything is just using pandas..
res <- res %>%
  select(-Language, -Method)

# add a more useful "Method" column and reorder columns
res <- res %>%
  mutate(Method = sprintf("%s (%s)", Library, Compression)) %>%
  select(Dataset, Method, everything()) %>%
  select(-Library, -Compression)

# finally, let's add a fields representing the average compression ratio, and the
# compression ratio to average i/o ratio
data_sizes <- res %>%
  filter(Method == 'pandas (none)') %>%
  select(dataset=Dataset, size=`File Size (MB)`)

res$natural_size <- data_sizes$size[match(res$Dataset, data_sizes$dataset)]

res <- res %>%
  mutate(`Compression Ratio` = natural_size / `File Size (MB)`) %>%
  mutate(`Compression to I/O Ratio` = `Compression Ratio` / (0.5 * (`Read Time (Secs)` + `Write Time (Secs)`)))

knitr::kable(res, digits = 2)
```

| Dataset   | Method             | Read Time (Secs) | Write Time (Secs) | File Size (MB) | natural\_size | Compression Ratio | Compression to I/O Ratio |
| :-------- | :----------------- | ---------------: | ----------------: | -------------: | ------------: | ----------------: | -----------------------: |
| itraq     | pandas (none)      |             0.12 |              1.13 |          14.23 |         14.23 |              1.00 |                     1.60 |
| golub1999 | pandas (none)      |             0.04 |              0.31 |           2.07 |          2.07 |              1.00 |                     5.77 |
| pancan    | pandas (none)      |             9.48 |             20.71 |         226.34 |        226.34 |              1.00 |                     0.07 |
| yeast     | pandas (none)      |             0.08 |              0.71 |          10.10 |         10.10 |              1.00 |                     2.52 |
| itraq     | pandas (gzip)      |             0.19 |              3.27 |           5.84 |         14.23 |              2.43 |                     1.41 |
| golub1999 | pandas (gzip)      |             0.05 |              1.08 |           0.76 |          2.07 |              2.72 |                     4.79 |
| pancan    | pandas (gzip)      |            10.52 |             39.12 |          76.82 |        226.34 |              2.95 |                     0.12 |
| yeast     | pandas (gzip)      |             0.13 |              1.27 |           4.82 |         10.10 |              2.10 |                     2.99 |
| itraq     | feather (none)     |             0.01 |              0.04 |           9.30 |         14.23 |              1.53 |                    67.60 |
| golub1999 | feather (none)     |             0.01 |              0.01 |           3.96 |          2.07 |              0.52 |                    47.52 |
| pancan    | feather (none)     |             0.22 |              1.42 |         133.13 |        226.34 |              1.70 |                     2.08 |
| yeast     | feather (none)     |             0.00 |              0.01 |           4.54 |         10.10 |              2.22 |                   349.88 |
| itraq     | pyarrow (lz4)      |             0.01 |              0.06 |           3.50 |         14.23 |              4.07 |                   113.56 |
| golub1999 | pyarrow (lz4)      |             0.02 |              0.03 |           1.07 |          2.07 |              1.94 |                    82.04 |
| pancan    | pyarrow (lz4)      |             2.17 |              2.62 |         145.48 |        226.34 |              1.56 |                     0.65 |
| yeast     | pyarrow (lz4)      |             0.01 |              0.03 |           2.47 |         10.10 |              4.10 |                   206.91 |
| itraq     | fastparquet (lz4)  |             0.06 |              0.09 |           5.31 |         14.23 |              2.68 |                    36.10 |
| golub1999 | fastparquet (lz4)  |             0.03 |              0.09 |           1.53 |          2.07 |              1.35 |                    22.47 |
| pancan    | fastparquet (lz4)  |             5.33 |              8.66 |         121.36 |        226.34 |              1.87 |                     0.27 |
| yeast     | fastparquet (lz4)  |             0.02 |              0.05 |           3.01 |         10.10 |              3.36 |                    95.87 |
| itraq     | pyarrow (snappy)   |             0.01 |              0.06 |           3.49 |         14.23 |              4.07 |                   115.28 |
| golub1999 | pyarrow (snappy)   |             0.02 |              0.03 |           1.07 |          2.07 |              1.94 |                    82.39 |
| pancan    | pyarrow (snappy)   |             2.26 |              2.64 |         145.09 |        226.34 |              1.56 |                     0.64 |
| yeast     | pyarrow (snappy)   |             0.01 |              0.03 |           2.61 |         10.10 |              3.87 |                   202.20 |
| itraq     | pyarrow (zstd)     |             0.01 |              0.06 |           3.30 |         14.23 |              4.31 |                   110.34 |
| golub1999 | pyarrow (zstd)     |             0.02 |              0.03 |           0.85 |          2.07 |              2.45 |                    96.35 |
| pancan    | pyarrow (zstd)     |             2.18 |              3.33 |         137.31 |        226.34 |              1.65 |                     0.60 |
| yeast     | pyarrow (zstd)     |             0.01 |              0.04 |           2.34 |         10.10 |              4.32 |                   174.61 |
| itraq     | fastparquet (zstd) |             0.06 |              0.11 |           3.29 |         14.23 |              4.32 |                    49.66 |
| golub1999 | fastparquet (zstd) |             0.03 |              0.10 |           0.83 |          2.07 |              2.50 |                    36.87 |
| pancan    | fastparquet (zstd) |             5.60 |              9.16 |         112.67 |        226.34 |              2.01 |                     0.27 |
| yeast     | fastparquet (zstd) |             0.02 |              0.07 |           2.15 |         10.10 |              4.70 |                   105.27 |

# Read Time (seconds)

``` r
ggplot(res, aes(x = Method, y = `Read Time (Secs)`, fill = Method)) +
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
ggplot(res, aes(x = Method, y = `Write Time (Secs)`, fill = Method)) +
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
ggplot(res, aes(x = Method, y = `File Size (MB)`, fill = Method)) +
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
ggplot(res, aes(x = Method, y = `Compression to I/O Ratio`, fill = Method)) +
  geom_bar(stat = 'identity') + 
  facet_wrap(~Dataset, scales = 'free_y') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Compression to Average I/O Speed Ratio (Higher is Better)") +
  ylab("Compression to Average I/O Speed Ratio")
```

![](summarize_results_files/figure-gfm/compression_to_speed_ratio-1.png)<!-- -->
