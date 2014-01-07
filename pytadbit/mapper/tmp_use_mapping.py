"""
07 Jan 2014


"""
from hiclib import mapping
from mapping import iterative_mapping_gem


PATH   = '/home/fransua/Tools/hiclib/examples/tutorial/'
INFILE = PATH + 'data/datatestshort.txt'

mapping.iterative_mapping_gem(
    gem_index_path       = PATH + 'data/hg19_AXYM.gem',
    fastq_path           = INFILE,
    out_sam_path         = PATH + 'data/SRR027056_llg1.txt',
    min_seq_len          = 25,
    len_step             = 5,
    seq_start            = 0,
    seq_end              = 75,
    nthreads             = 8,  # on intel corei7 CPUs 4 threads are as fast as
                               # 8, but leave some room for you other applications
    # max_reads_per_chunk = 10000000,  #optional, on low-memory machines
    temp_dir             = PATH + 'data/tmp',  # optional, keep temporary files here
    # bash_reader         = 'bin/sratoolkit.2.1.16-ubuntu32/bin/fastq-dump -Z'
    bash_reader          = 'cat',
    single_end           = True
    )



iterative_mapping_gem(
    gem_index_path       = PATH + 'data/hg19_AXYM.gem',
    fastq_path           = INFILE,
    out_sam_path         = PATH + 'data/SRR027056_lltdbg1.txt',
    min_seq_len          = 25,
    len_step             = 5,
    seq_start            = 0,
    seq_end              = 75,
    nthreads             = 8,  # on intel corei7 CPUs 4 threads are as fast as
                               # 8, but leave some room for you other applications
    # max_reads_per_chunk = 10000000,  #optional, on low-memory machines
    temp_dir             = PATH + 'data/tmp',  # optional, keep temporary files here
    # bash_reader         = 'bin/sratoolkit.2.1.16-ubuntu32/bin/fastq-dump -Z'
    single_end           = True
    )

