import logging
import ujson
import gzip
import os
import numpy as np
from pathlib import Path
from os.path import join, basename, exists
from pandas import read_pickle
from datetime import datetime
from collections import defaultdict
from .helpers.Bowtie2WrapperSlim import do_pair_simple, MapPreset
from .helpers.ICRAUtils import timeit, _open_gz_indif
from .helpers.ReadContainer import ReadContainer, _get_flank_size
from .helpers.sam2pmp import sam2pmp, SourceReadSAM

log_ = logging.getLogger('ICRA')

GENOMES = 'genomes'
CAMI = 'cami'

def get_ujson_splitting_point(delta):
    no_aln = 0
    lens = []
    for idx in range(len(delta)):
        no_aln += len(delta[idx][1])
        lens.append([idx, len(delta[idx][1])])
    print(f'no of alignments is {str(no_aln)}')

    counter = 0
    for index, rd_cnts in lens:
        counter += rd_cnts
        if counter >= no_aln / 2:
            return index


def single_file(
        fq1, fq2, 
        outfol=os.getcwd(),
        max_mismatch=8,
        consider_lengths=False,
        epsilon=1e-6,
        max_iterations=100,
        min_bins=10,
        max_bins=100,
        min_reads=100, 
        dense_region_coverage=60, 
        length_minimum=1e5,
        length_maximum=2e7, 
        dbpath=None, 
        use_theta=False, 
        threads=1, 
        senspreset=MapPreset.VERY_SENSITIVE,
        report_alns=20, 
        max_ins=750
    ):
    if not os.path.exists(outfol):
        print('Creating directory %s...'%outfol)
        Path(outfol).mkdir(parents=True, exist_ok=True)

    if fq2 is None:
        print('Running ICRA on single-end read!')
        print('Single-',fq1)
    else:
        print('Running ICRA on paired-end reads!')
        print('Forward-',fq1)
        print('Reverse-',fq2)
    print('--------------------------------------------\n')

    log_.debug('Loading database...')
    length_db_f, indexf, dest_dictf, dlen_db_f = dbpath + '.lengths', dbpath, dbpath + '.dests', dbpath + '.dlen'
    outpref = join(outfol,
                   basename(fq1).replace('_1.fastq.gz', '').replace('_1.fastq', '').replace('.fastq.gz', '').replace(
                       '.fastq', ''))
    log_.debug('Loaded. Mapping file...')

    #AYA added Feb 2022
    if not exists(outpref + '.bam'):
        do_pair_simple(fq1, fq2, outpref, indexf, senspreset, report_alns, max_ins, threads)
    else:
        print('%s already exists, skipping mapping step...'%(outpref + '.bam'))

    #do_pair_simple(fq1, fq2, outpref, indexf, senspreset, report_alns, max_ins, threads)
    log_.debug('Mapped. Converting to pmp')
    if not exists(outpref + '.pmp'):
        sam2pmp(outpref + '.bam', outpref + '.pmp', full=True)
    #AYA DEBUG - i commented it out to save bam
    #_tryrm(outpref + '.bam')
    log_.debug('PMP ready. Running ICRA...')
    #if not exists('delta_new_for_debug.pkl'):
    if not exists(outpref + '_1.jsdel'):

        read_container, pi, theta1, average_read_length, lengthdb = \
            _initialize(fq1, fq2, outpref + '.pmp', dlen_db_f, max_mismatch, consider_lengths, length_minimum,
                        length_maximum,
                        min_reads, min_bins, max_bins, dense_region_coverage, dest_dictf, use_theta)
        if len(pi) == 0:
            delta = {}
        else:
            delta, pi = _runIterative(pi, theta1, read_container, min_bins, max_bins,
                                      min_reads, average_read_length, max_mismatch,
                                      lengthdb, dense_region_coverage, consider_lengths,
                                      epsilon, max_iterations, use_theta)
        delta_new = [("", mappings) for i, (read_id, mappings) in enumerate(delta)]
        #else:
            # delta_new = read_pickle('delta_new_for_debug.pkl')
            # average_read_length = 75
        print('this is the len of delta:')
        print(len(delta_new))

        with gzip.open(outpref + '.jspi', 'wt') as of:
            ujson.dump(pi, of)
        #half = get_ujson_splitting_point(delta_new)
        half = int(len(delta_new)/2)
        twenieth = int(len(delta_new)/20)
        # if len(delta_new) > 50: #50000000:

        #     lst1 = delta_new[0:half] # first half
        #     print(1, len(lst1))
        #     with gzip.open(outpref + '_1.jsdel', 'wt') as of:
        #         ujson.dump([average_read_length, lst1], of)

        #     start = half
        #     for num in range(2,12):
        #         if num == 11 :
        #             lst = delta_new[start:]
        #             print(num, len(lst))
        #         else:
        #             lst = delta_new[start: start+twenieth]
        #             print(num, len(lst))
        #         with gzip.open(outpref + f'_{num}.jsdel', 'wt') as of:
        #             ujson.dump([average_read_length, lst], of)
        #         start += twenieth

        # else:
        with gzip.open(outpref + '.jsdel', 'wt') as of:
            ujson.dump([average_read_length, delta_new], of)
        return outpref + '.jspi', outpref + '.jsdel'
    else:
        return outpref + '.jspi', outpref + '.jsdel'


@timeit(log_, "%s: Parameter init complete. Time: %s", logging.INFO)
def _initialize(fq1, fq2, pmpf, dlen_db_f, max_mismatch, consider_lengths,
                length_minimum, length_maximum, min_reads, min_bins,
                max_bins, dense_region_coverage, dest_dictf, use_theta):
    dest_dict = read_pickle(dest_dictf) if dest_dictf is not None else None
    average_read_length = _getaveragereadlength(fq1)
    lengthdb = read_pickle(dlen_db_f)
    read_container = ReadContainer()
    pe = fq2 is not None
    allowed_dests = {k for k, v in lengthdb.items() if length_minimum <= v <= length_maximum}
    for sr in _load_from_file(pmpf):
        read_container.add_pmp_sr(sr, pe, allowed_dests, max_mismatch, None, dest_dict)
    del dest_dict, allowed_dests
    pi = read_container.get_pi()

    read_container.remove_dests([k for k, val in pi.items() if val <= min_reads])
    pi = {k: val for k, val in pi.items() if val > min_reads}

    read_container.init_cov_dict(pi, min_bins, max_bins, min_reads, average_read_length, lengthdb)
    # note: this will create a fully updated coverage dict in the process
    # (in the original implementation this was explicit and not implicit)
    pi = read_container.get_dense_region_coverage_pi(dense_region_coverage)

    read_container.remove_dests([k for k, val in pi.items() if val <= min_reads])
    pi = {k: v for k, v in pi.items() if v > min_reads}

    if consider_lengths:
        flank_sizes = _get_flank_sizes(pi, lengthdb, min_bins, max_bins, average_read_length)
        pi = {k: (v / (lengthdb[k] - (2 * flank_sizes[k]))) for k, v in pi.items()}
    pisum = sum(pi.values())
    pi = {k: v / float(pisum) for k, v in pi.items()}

    theta1 = read_container.get_theta() if use_theta else None

    return read_container, pi, theta1, average_read_length, lengthdb


def _runIterative(pi, theta1, read_container, min_bins, max_bins, min_reads, average_read_length, max_mismatch,
                  lengthdb, dense_region_coverage, consider_lengths, epsilon, max_iterations, use_theta):
    starttime = datetime.now()
    pi_dict = {}
    i = 1

    while True:
        pi_dict[i] = pi
        delta = read_container.do_delta(pi, theta1, 1. / ((len(pi) ** 2) * (10 ** max_mismatch)), use_theta)
        read_container.recalc_cov_bins()
        prevPi = pi
        # reminder: implicitly creates an updated covdic
        pi = read_container.get_dense_region_coverage_pi(dense_region_coverage, delta)
        read_container.remove_dests([k for k, val in pi.items() if val < min_reads])
        pi = {k: v for k, v in pi.items() if v >= min_reads}

        if len(pi) == 0:
            log_.info("No adequately covered strains found")
            return {}, {}

        theta1 = read_container.get_theta() if use_theta else None

        if consider_lengths:
            flank_sizes = _get_flank_sizes(pi, lengthdb, min_bins, max_bins, average_read_length)
            pi = {k: (v / (lengthdb[k] - (2 * flank_sizes[k]))) for k, v in pi.items()}
        pisum = sum(pi.values())
        pi = {k: v / float(pisum) for k, v in pi.items()}

        prevPi = {k: v for k, v in prevPi.items() if k in pi}
        dPi = _LogDistDict(pi, prevPi)

        i += 1
        log_.info("Iteration {} - Time: {}, dPi = {:.2e}, nPi = {}".format(i, datetime.now() - starttime, dPi,
                                                                           len(pi)))

        if dPi < epsilon or i > max_iterations:
            break
    log_.info("Final result - Time: {}".format(datetime.now() - starttime))
    return read_container.get_full_delta(pi, theta1, 1. / ((len(pi) ** 2) * (10 ** max_mismatch)), use_theta), pi


def _logdist(u, v):
    x = np.abs(np.log10(np.asarray(u)) - np.log10(np.asarray(v)))
    return np.percentile(x, 90)


def _LogDistDict(u, v):
    x = []
    y = []
    if not (set(u.keys()).intersection(set(v.keys())) == set(u.keys()) == set(v.keys())):
        raise KeyError("Keys not overlapping between the two dictionaries")
    for key in u:
        if u[key] > 1e-5 or v[key] > 1e-5:
            x.append(u[key])
            y.append(v[key])
    return _logdist(x, y)


def _getaveragereadlength(fq):
    with _open_gz_indif(fq) as inf:
        sum_bps = 0
        for i, read in enumerate(inf):
            if i % 4 == 1:
                sum_bps += len(read) - 1  # has a /n at the end
            if i >= 40000: break
    arl = int(np.round(sum_bps / 10000))
    log_.info("Average read length for {} set to {}".format(basename(fq), arl))
    return arl


def _get_len_dct(length_db_f, read_dest_dict=None):
    lengthdb = read_pickle(length_db_f)
    if read_dest_dict is not None:
        newlengthdb = defaultdict(int)
        for l in lengthdb:
            newlengthdb[read_dest_dict[l][0]] += lengthdb[l]
        return dict(newlengthdb)
    else:
        return lengthdb


def _get_flank_sizes(pi, lengthdb, min_bins, max_bins, average_read_length):
    numbins = {k: max(min_bins, min(max_bins, int(v / 100))) + 2 for k, v in pi.items()}
    return {k: _get_flank_size(lengthdb[k], numbins[k], average_read_length) for k, v in pi.items()}


def _load_from_file(fpath):
    with gzip.open(fpath) as fin:
        for r in _load_iterable_fromdesc(fin):
            yield r


def _load_iterable_fromdesc(desc):
    try:
        while True:
            yield SourceReadSAM.from_ser(ujson.loads(desc.readline()))
    except ValueError as ve:
        if ve.args[0] == 'No JSON object could be decoded' or \
                ve.args[0] == 'Expected object or value':
            return
        else:
            raise ve
    except EOFError:
        return
