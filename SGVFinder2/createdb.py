import argparse
from collections import defaultdict
from glob import glob
from os.path import join, basename, splitext
from Bio import SeqIO
import logging

from pandas import to_pickle

from .helpers.ICRAUtils import _open_gz_indif, _set_logging_handlers

log_ = logging.getLogger()

def create_db_from_reps(input_path, out_prefix):
    dlen_dict = {}
    lengths_dict = {}
    dests_dict = {}

    log_.info('Starting...')

    with open(out_prefix + '.fasta', 'wt') as o_fa_h:
        for f in sum((glob(join(input_path, glb)) for glb in ('*.fasta', '*.fasta.gz', '*.fa', '*.fa.gz')), list()):
            destid = splitext(basename(f.replace('.gz', '')))[0].replace(' ', '_')
            assert destid not in dlen_dict
            dlen_dict[destid] = 0

            log_.info(f'Processing {destid}...')

            with _open_gz_indif(f) as in_fa_spec:
                for r in SeqIO.parse(in_fa_spec, 'fasta'):
                    r.description = r.description.replace(' ', '_')
                    if not r.description.startswith(destid + '.'):
                        if r.description.startswith(destid):
                            r.id = r.description.replace(destid, destid + '.')
                        else:
                            r.id = f'{destid}.{r.description}'
                    else:
                        r.id = r.description
                    r.description = ''
                    dpid = r.id#description
                    assert dpid not in lengths_dict and dpid not in lengths_dict
                    dests_dict[dpid] = (destid, dlen_dict[destid])
                    lengths_dict[dpid] = len(r)
                    dlen_dict[destid] += lengths_dict[dpid]
                    SeqIO.write(r, o_fa_h, 'fasta')
    to_pickle(dlen_dict, out_prefix + '.dlen')
    to_pickle(lengths_dict, out_prefix + '.lengths')
    to_pickle(dests_dict, out_prefix + '.dests')
    log_.info(f"Done. Please run bowtie2-build {out_prefix+'.fasta'} {out_prefix}")


if __name__ == '__main__':
    _set_logging_handlers()
    parser = argparse.ArgumentParser()
    parser.add_argument('input_path', help='A folder containing one fasta file per genome. Files need to end with '
                        'either fasta, fasta.gz, fa, or fa.gz')
    parser.add_argument('output_prefix',
                        help='prefix for the created database. This is what you\'ll later supply to ICRA or SGVF.')
    args = parser.parse_args()
    create_db_from_reps(args.input_path, args.output_prefix)
