"""
Compressed data performance benchmark
KH Jan 09, 2020
"""
import os
import shutil
import numpy as np
import pandas as pd
from zipfile import ZipFile 

out_dir = config['output_dir']

rule combine_timings:
    input:
        expand(os.path.join(out_dir, 'timings', '{dataset}', 'pandas_csv_uncompressed.csv'), dataset=config['datasets']),
        expand(os.path.join(out_dir, 'timings', '{dataset}', 'pandas_csv_gzip.csv'), dataset=config['datasets']),
        expand(os.path.join(out_dir, 'timings', '{dataset}', 'pandas_feather.csv'), dataset=config['datasets']),
        expand(os.path.join(out_dir, 'timings', '{dataset}', 'pandas_parquet_pyarrow_lz4.csv'), dataset=config['datasets']),
        expand(os.path.join(out_dir, 'timings', '{dataset}', 'pandas_parquet_fastparquet_lz4.csv'), dataset=config['datasets']),
        expand(os.path.join(out_dir, 'timings', '{dataset}', 'pandas_parquet_pyarrow_snappy.csv'), dataset=config['datasets']),
        expand(os.path.join(out_dir, 'timings', '{dataset}', 'pandas_parquet_pyarrow_zstd.csv'), dataset=config['datasets']),
        expand(os.path.join(out_dir, 'timings', '{dataset}', 'pandas_parquet_fastparquet_zstd.csv'), dataset=config['datasets'])
    output: 
        os.path.join(out_dir, 'timings', 'all_timings.csv')
    run:
        dfs = []

        for infile in input:
            dfs.append(pd.read_csv(infile, header=None))

        dat = pd.concat(dfs)

        dat.columns = ['Dataset', 'Language', 'Library', 'Method', 'Compression', 'Read Time (Secs)', 'Write Time (Secs)', 'File Size (MB)']
        dat.to_csv(output[0], index=False)

rule benchmark_parquet_pyarrow_lz4:
    input: os.path.join(out_dir, 'datasets', '{dataset}.csv'),
    output:
        data=os.path.join(out_dir, 'output', '{dataset}', 'pyarrow_lz4.parquet'),
        timings=os.path.join(out_dir, 'timings', '{dataset}', 'pandas_parquet_pyarrow_lz4.csv')
    params:
        engine='pyarrow',
        compression='lz4'
    threads: config['benchmark']['num_threads']
    script: 'src/benchmark_parquet.py'

rule benchmark_parquet_fastparquet_lz4:
    input: os.path.join(out_dir, 'datasets', '{dataset}.csv'),
    output:
        data=os.path.join(out_dir, 'output', '{dataset}', 'fastparquet_lz4.parquet'),
        timings=os.path.join(out_dir, 'timings', '{dataset}', 'pandas_parquet_fastparquet_lz4.csv')
    params:
        engine='fastparquet',
        compression='lz4'
    threads: config['benchmark']['num_threads']
    script: 'src/benchmark_parquet.py'

rule benchmark_parquet_pyarrow_snappy:
    input: os.path.join(out_dir, 'datasets', '{dataset}.csv'),
    output:
        data=os.path.join(out_dir, 'output', '{dataset}', 'pyarrow_snappy.parquet'),
        timings=os.path.join(out_dir, 'timings', '{dataset}', 'pandas_parquet_pyarrow_snappy.csv')
    params:
        engine='pyarrow',
        compression='snappy'
    threads: config['benchmark']['num_threads']
    script: 'src/benchmark_parquet.py'

rule benchmark_parquet_pyarrow_zstd:
    input: os.path.join(out_dir, 'datasets', '{dataset}.csv'),
    output:
        data=os.path.join(out_dir, 'output', '{dataset}', 'pyarrow_zstd.parquet'),
        timings=os.path.join(out_dir, 'timings', '{dataset}', 'pandas_parquet_pyarrow_zstd.csv')
    params:
        engine='pyarrow',
        compression='zstd'
    threads: config['benchmark']['num_threads']
    script: 'src/benchmark_parquet.py'

rule benchmark_parquet_fastparquet_zstd:
    input: os.path.join(out_dir, 'datasets', '{dataset}.csv'),
    output:
        data=os.path.join(out_dir, 'output', '{dataset}', 'fastparquet_zstd.parquet'),
        timings=os.path.join(out_dir, 'timings', '{dataset}', 'pandas_parquet_fastparquet_zstd.csv')
    params:
        engine='fastparquet',
        compression='zstd'
    threads: config['benchmark']['num_threads']
    script: 'src/benchmark_parquet.py'

rule benchmark_feather:
    input: os.path.join(out_dir, 'datasets', '{dataset}.csv'),
    output:
        data=os.path.join(out_dir, 'output', '{dataset}', 'dat.feather'),
        timings=os.path.join(out_dir, 'timings', '{dataset}', 'pandas_feather.csv')
    threads: config['benchmark']['num_threads']
    script: 'src/benchmark_feather.py'

rule benchmark_pandas_uncompressed_csv:
    input: os.path.join(out_dir, 'datasets', '{dataset}.csv'),
    output:
        data=os.path.join(out_dir, 'output', '{dataset}', 'pandas_csv_uncompressed.csv'),
        timings=os.path.join(out_dir, 'timings', '{dataset}', 'pandas_csv_uncompressed.csv')
    threads: config['benchmark']['num_threads']
    script: 'src/benchmark_pandas_csv.py'

rule benchmark_pandas_compressed_csv:
    input: os.path.join(out_dir, 'datasets', '{dataset}.csv'),
    output:
        data=os.path.join(out_dir, 'output', '{dataset}', 'pandas_csv_gzip.csv.gz'),
        timings=os.path.join(out_dir, 'timings', '{dataset}', 'pandas_csv_gzip.csv')
    threads: config['benchmark']['num_threads']
    script: 'src/benchmark_pandas_csv.py'

rule download_test_data:
    run:
        import kaggle
        from zipfile import ZipFile

        kaggle.api.authenticate()

        # retrieve datasets
        for dataset in config['datasets'].keys():
            # check if file already exists, and if so, skip
            outfile = os.path.join(out_dir, 'datasets', dataset + '.csv')

            if os.path.exists(outfile):
                continue

            # download zipfile for specified dataset/filename
            kaggle.api.dataset_download_file(config['datasets'][dataset]['id'],
                                             config['datasets'][dataset]['file'],
                                             path='/tmp')

            # unzip dataset
            zip_file = os.path.join('/tmp', config['datasets'][dataset]['file'] + '.zip')

            with ZipFile(zip_file, 'r') as zobj:
                zobj.extractall('/tmp')

            # move to brenchmark data directory and rename
            shutil.move(zip_file.replace('.zip', ''), outfile)
