"""
12 Dec 2013


"""
from os.path import abspath, expanduser, join as joinpath, split as splitpath
from os      import remove
import tempfile
import gem
import pysam
import gzip


def trimming(raw_seq_len, seq_start, min_seq_len):
    return seq_start, raw_seq_len - seq_start - min_seq_len


def iterative_mapping_gem(gem_index_path, fastq_path, out_sam_path,
                          min_seq_len, len_step, **kwargs):
    """
    :param fastq_path: 152 bases first 76 from one end, next 76 from the other
       end. Both to be read from left to right.
    """
    gem_index_path      = abspath(expanduser(gem_index_path))
    fastq_path          = abspath(expanduser(fastq_path))
    out_sam_path        = abspath(expanduser(out_sam_path))
    single_end          = kwargs.get('single_end'          , False)
    seq_start           = kwargs.get('seq_start'           , 0)
    seq_end             = kwargs.get('seq_end'             , None)
    nthreads            = kwargs.get('nthreads'            , 4)
    max_edit_distance   = kwargs.get('max_edit_distance'   , 0.04)
    mismatches          = kwargs.get('mismatches'          , 0.04)
    nthreads            = kwargs.get('nthreads'            , 4)
    temp_dir = abspath(expanduser(
        kwargs.get('temp_dir', tempfile.gettempdir())))
    fastqh = open(fastq_path)
    raw_seq_len = int(fastqh.next().strip().split('length=')[1])
    fastqh.close()

    # end position according to sequence in the file
    local_seq_end = min(raw_seq_len, seq_end) if seq_end else raw_seq_len

    # we grow seq to map up to 75 bases, no more
    if min_seq_len > local_seq_end - seq_start:
        return

    # define what we trim
    trim_5, trim_3 = trimming(raw_seq_len, seq_start, min_seq_len)

    # output
    local_out_sam = out_sam_path + '.' + str(min_seq_len)
    # input
    inputf = gem.files.open(fastq_path)

    # trimming
    trimmed = gem.filter.run_filter(inputf, 
                                    ['--hard-trim', '%d,%d' % (trim_5, trim_3)],
                                    threads=nthreads, paired=not single_end)
    # mapping
    mapped = gem.mapper(trimmed, gem_index_path, min_decoded_strata=0,
                        max_decoded_matches=2, unique_mapping=True,
                        max_edit_distance=max_edit_distance,
                        mismatches=mismatches,
                        output='/tmp/test.map',
                        threads=nthreads)

    # convert to sam
    sam = gem.gem2sam(mapped, index=gem_index_path, output=local_out_sam,
                      threads=nthreads, single_end=single_end)

    # Check if the next iteration is required.
    if len_step <= 0:
        return
    if min_seq_len + len_step > local_seq_end - seq_start:
        return

    # Recursively go to the next iteration.
    unmapped_fastq_path = joinpath(
        temp_dir, splitpath(fastq_path)[1] + '.%d' % min_seq_len)
    _filter_unmapped_fastq(fastq_path, local_out_sam, unmapped_fastq_path)

    iterative_mapping_gem(gem_index_path, unmapped_fastq_path,
                          out_sam_path,
                          min_seq_len=min_seq_len + len_step,
                          len_step=len_step, **kwargs)

    remove(unmapped_fastq_path)

    

def _filter_unmapped_fastq(in_fastq, in_sam, nonunique_fastq):
    '''Read raw sequences from **in_fastq** and alignments from
    **in_sam** and save the non-uniquely aligned and unmapped sequences
    to **unique_sam**.
    '''
    samfile = pysam.Samfile(in_sam)

    nonunique_ids = set()
    for read in samfile:
        tags_dict = dict(read.tags)
        read_id = read.qname
        # If exists, the option 'XS' contains the score of the second
        # best alignment. Therefore, its presence means non-unique alignment.
        if 'XS' in tags_dict or read.is_unmapped or (
            'NM' in tags_dict and int(tags_dict['NM']) > 1):
            nonunique_ids.add(read_id)

    _filter_fastq(nonunique_ids, in_fastq, nonunique_fastq)


def _filter_fastq(ids, in_fastq, out_fastq):
    '''Filter FASTQ sequences by their IDs.

    Read entries from **in_fastq** and store in **out_fastq** only those
    the whose ID are in **ids**.
    '''
    out_file = open(out_fastq, 'w')
    in_file = _gzopen(in_fastq)
    while True:
        line = in_file.readline()
        if not line:
            break

        if not line.startswith('@'):
            raise Exception(
                '{0} does not comply with the FASTQ standards.'.format(in_fastq))

        fastq_entry = [line, in_file.readline(),
                       in_file.readline(), in_file.readline()]
        read_id = line.split()[0][1:]
        if read_id.endswith('/1') or read_id.endswith('/2'):
            read_id = read_id[:-2]
        if read_id in ids:
            out_file.writelines(fastq_entry)


def _gzopen(path):
    if path.endswith('.gz'):
        return gzip.open(path)
    else:
        return open(path)

