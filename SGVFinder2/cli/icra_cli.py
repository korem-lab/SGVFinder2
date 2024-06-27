import logging
import argparse
from pandas import to_pickle
from ..icra import single_file
from ..svfinder import get_sample_map
from ..helpers.Bowtie2WrapperSlim import MapPreset, do_pair_simple
from ..helpers.sam2pmp import sam2pmp
from os.path import join, basename, exists

parser = argparse.ArgumentParser()

parser.add_argument('--outfol', help='Where to write results')
parser.add_argument('--bamfol', help='Path to generate or to the exists bam file (defaults to outfol)')
parser.add_argument('--pmpfol', help='Path to generate or to the exists pmp file (defaults to outfol)')
parser.add_argument('--fq1', help='Path to first fastq file')
parser.add_argument('--fq2', help='Path to second fastq file (defaults to None)', default=None)
parser.add_argument('--db', help='Database path (including db prefix)', default=None)
parser.add_argument('--max_mismatch', help='How many mismatch are considered acceptable', default=8, type=int)
parser.add_argument('--consider_lengths', help='Should genome lengths be considered when calculating abundances', action='store_false')
parser.add_argument('--epsilon', help='The stop criteria. Epsilon is the euclidean distance between the internal vectors of element abundances under which ICRA stops.', default=1e-6, type=float)
parser.add_argument('--max_iterations', help='An upper limit to the number of iterations ICRA will run', default=100, type=int)
parser.add_argument('--min_bins', help='The minimum number of bins per genomic element.', default=10, type=int)
parser.add_argument('--max_bins', help='The maximum number of bins per genomic element', default=100, type=int)
parser.add_argument('--min_reads',help='Minimal number of reads mapped to a genomic element for it to be considered present in the sample.', default=100, type=int)
parser.add_argument('--dense_region_coverage', help='The percentage of the genome examined for coverage purposes.', default=60, type=int)
parser.add_argument('--length_minimum', help='Minimal genome length considered,', default=1e5, type=float)
parser.add_argument('--length_maximum', help='Maximal genome length considered.', default=2e7, type=float)
parser.add_argument('--use_theta',help='Theta is an extra algorithmic components that causes ICRA to converge faster',action='store_false')
parser.add_argument('--debug', help='Logs debug output', action='store_true')
parser.add_argument('--threads', help='Number of threads to use', default=1,type=float)
parser.add_argument('--sensitivity', help='Bowtie2 sensitivity (defaults to very-sensitive)', default=MapPreset.VERY_SENSITIVE)
parser.add_argument('--report_alignments', help='Bowtie2 alignment mode, specify `all` for -a or an integer for -k', default=20)
parser.add_argument('--max_ins', help='Bowtie2 -X parameter, the maximum fragment length for valid paired-end alignments', default=750)
parser.add_argument('--x_coverage', help = '`get_sample_map` param- the desired coverage across the genome in units of 100bp reads. This parameter is used to determine bin size: bin_size = rate_param/x_coverage (Default = 0.01)', type=float, default = 0.01)
parser.add_argument('--rate_param', help = '`get_sample_map` param- the lower limit for the median number of reads per genomic bin. Genomes with coverage lower than rate_param will be discarded from the analysis (default = 10)', type = int, default = 10)

parser.add_argument('--generate_bam', action='store_true', help='Only generate BAM files from FASTQ files using bowtie2 and samtools.')
parser.add_argument('--bam_to_pmp', action='store_true', help='Only convert BAM files to PMP files. Use this after generating BAM files with bowtie2 and samtools (`generate_bam` command).')

args = parser.parse_args()
logging.basicConfig(level=(logging.DEBUG if args.debug else logging.INFO),
                    format='%(asctime)s-%(levelname)s: %(message)s')
logging.debug(args)

def run():
    fq1= args.fq1
    fq2 = args.fq2
    if args.bamfol is None:
        args.bamfol = args.outfol
    if args.pmpfol is None:
        args.pmpfol = args.outfol
    outfol = args.outfol
    bamfol = args.bamfol
    pmpfol = args.pmpfol
    file_basename = basename(fq1).replace('_1.fastq.gz', '').replace('_1.fastq', '').replace('.fastq.gz', '').replace(
            '.fastq', '').replace('_1.fq.gz', '').replace('_1.fq', '').replace('.fq.gz', '').replace('.fq', '')
    bampref = join(bamfol, file_basename) 
    pmppref = join(pmpfol, file_basename) 
    if args.generate_bam:
        print('Running ICRA `generate_bam` command...')
        do_pair_simple(fq1, fq2, bampref, 
                       dbpath = args.db, 
                       senspreset=args.sensitivity, 
                       report_alns=args.report_alignments, 
                       max_ins=args.max_ins, 
                       threads=args.threads)
        print('Mapped. BAM ready')
    elif args.bam_to_pmp:
        print('Running ICRA `bam_to_pmp` command...')
        if exists(bampref + '.bam'):
            print('%s exists, converting to pmp...'%(bampref + '.bam'))
            sam2pmp(bampref + '.bam', pmppref + '.pmp', full=True)
            print('Finished running bam_to_pmp, saving results to %s'%args.pmpfol)
        else:
            print('NO BAM files found in %s. Please run ICRA `generate_bam` command first!'%(bamfol))
    else:
        print('Running ICRA `single_file` command...')
        jspi_file, jsdel_file = single_file(
            fq1, fq2, bamfol, pmpfol, outfol,
            max_mismatch=args.max_mismatch,
            consider_lengths=args.consider_lengths,
            epsilon=args.epsilon,
            max_iterations=args.max_iterations,
            min_bins=args.min_bins,
            max_bins=args.max_bins,
            min_reads=args.min_reads,
            dense_region_coverage=args.dense_region_coverage,
            length_minimum=args.length_minimum,
            length_maximum=args.length_maximum,
            dbpath=args.db,
            use_theta=args.use_theta,
            threads=args.threads,
            senspreset=args.sensitivity,
            report_alns=args.report_alignments,
            max_ins=args.max_ins
        )
        print('Finished running ICRA, saving results to %s'%args.outfol)
        print('Running `get_sample_map` on %s, output will be saved to %s'%(jsdel_file,jsdel_file.replace('.jsdel','.smp')))
        sample_map = get_sample_map(jsdel_file,args.db+'.dlen',args.x_coverage, args.rate_param) 
        to_pickle(sample_map,jsdel_file.replace('.jsdel','.smp'))