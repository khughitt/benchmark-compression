"""
CSV compression benchmarck
KH Jan 9, 2020
"""
import os
import pandas as pd
import timeit

# load test data
df = pd.read_csv(snakemake.input[0])

# benchmark write performance
num_times = int(snakemake.config['benchmark']['num_times'])
write_times = timeit.repeat("df.to_csv('{}')".format(snakemake.output['data']), 
                            globals=globals(), number=1, repeat=num_times)

# benchmark read performance
read_times = timeit.repeat("pd.read_csv('{}')".format(snakemake.output['data']), 
                            globals=globals(), number=1, repeat=num_times)

# get output file size
fsize = os.path.getsize(snakemake.output['data']) / 1e6

# determine output label
compression = 'gzip' if snakemake.output[0].endswith('.gz') else 'none'

with open(snakemake.output['timings'], 'w') as fp:
    entry = [snakemake.wildcards['dataset'], 'python', 'pandas', 'pandas.read_csv', compression, 
             str(min(read_times)), str(min(write_times)), str(fsize)]
    fp.write(", ".join(entry) + "\n")
