"""
Parquet compression benchmarck
KH Jan 9, 2020
"""
import os
import pandas as pd
import timeit

# load test data
df = pd.read_csv(snakemake.input[0])

# parameters
engine = snakemake.params['engine']
compression = snakemake.params['compression']

# benchmark write performance
num_times = int(snakemake.config['benchmark']['num_times'])

write_cmd = "df.to_parquet('{}', engine='{}', compression='{}')".format(snakemake.output['data'], engine, compression)
write_times = timeit.repeat(write_cmd, globals=globals(), number=1, repeat=num_times)

# benchmark read performance
read_cmd = "pd.read_parquet('{}', engine='{}')".format(snakemake.output['data'], engine, compression)
read_times = timeit.repeat(read_cmd, globals=globals(), number=1, repeat=num_times)

# get output file size
fsize = os.path.getsize(snakemake.output['data']) / 1e6

with open(snakemake.output['timings'], 'w') as fp:
    entry = [snakemake.wildcards['dataset'],
            'python', 
             engine,
             'pandas.read_parquet',
             compression,
             str(min(read_times)), 
             str(min(write_times)), 
             str(fsize)]
    fp.write(", ".join(entry) + "\n")
