import pandas as pd
import logging
import argparse
import sys
from ..svfinder import work_on_collection, get_sample_map

BETAPRIME = 'betaprime'
CHISQ = 'ncx2'

parser = argparse.ArgumentParser()
parser.add_argument('--debug', help='Logs debug output', action='store_true')

sub_parsers = parser.add_subparsers(dest='command')

parser_woc = sub_parsers.add_parser('work_on_collection')
parser_gsm = sub_parsers.add_parser('get_sample_map')

# WORK_ON_COLLECTION commands
parser_woc.add_argument('--samp_to_map_dir', help='Path to folder with `.smp` files')
parser_woc.add_argument('--output_dsgv', help='Output path for the deletion-sgv dataframe. By default a pickled pandas dataframe')
parser_woc.add_argument('--output_vsgv', help='Output path for the variable-sgv dataframe. By default a pickled pandas dataframe')
parser_woc.add_argument('--max_spacing', help='Max spacing param (default 10)', type=float, default=10)
parser_woc.add_argument('--x_coverage', help='The desired coverage across the genome in units of 100bp reads. This parameter is used to determine bin size: bin_size = rate_param/x_coverage (Default = 0.1)', type=float, default=0.1)
parser_woc.add_argument('--rate_param', help='The lower limit for the median number of reads per genomic bin. Genomes with coverage lower than rate_param will be discarded from the analysis (Default = 10)', type=int, default=10)
parser_woc.add_argument('--min_samp_cutoff', help='Minimum number of samples in which a microbe exists with sufficient coverage to be considered in the analysis (Default=2)',type=int, default=2)
parser_woc.add_argument('--dels_detect_thresh', help='Determines the minimum and maximum ratio of samples for which a bin is considered a deletion-SGV. (Default=0.25, setting the ratio at 0.25-0.75)',type=float, default=0.25)
parser_woc.add_argument('--real_del_thresh', help='Threshold above which a bin is considered deleted for all individuals (Default=0.95)', type=float, default=0.95)
parser_woc.add_argument('--vsgv_dissim_thresh',
                    help='Maximal correlation dissimilarity for concatenation and clustering of variable SGV '
                            'bins. Correlation dissimilarity is defined as calculated as 1-((rho(u,v)+1)/2), '
                            'where rho is the Spearman correlation and u, v are the bin vectors being compared ('
                            'Default=0.125)',
                    type=float, default=0.125)
parser_woc.add_argument('--dels_cooc_thresh',
                    help='Maximal cooccurrence dissimilarity for concatenation and clustering of deletion SGV '
                            'bins. Coocurrence dissimilarity is defined as the proportion of samples which are in '
                            'disagreement on the deletion-state of the two bins being compared (wherein one bin is '
                            'deleted and one is retained for the same sample) out of all samples that harbor the '
                            'microbe (Default=0.25)',
                    type=float, default=0.25)
parser_woc.add_argument('--vsgv_clip_quantile',
                    help='Determines clipping performed on the distribution of bin values prior to fitting a '
                            'distribution for the detection of variable-SGVs (Default=0.02 corresponding to clipping '
                            'outside the 2nd to 98th percentiles)',
                    type=float, default=0.02)
parser_woc.add_argument('--vsgv_fit_interval',
                    help='Significance cutoff for the fitted distribution above which a bin is considered '
                            'variable (Default=0.95)',
                    type=float, default=0.95)
parser_woc.add_argument('--vsgv_fit_method',
                    help='Determines the distribution being fit on bin values (either a Beta-prime or Chi-square '
                            'distribution; Default=betaprime)',
                    type=str, default=BETAPRIME, choices=[BETAPRIME, CHISQ])
parser_woc.add_argument('--vsgv_dense_perc',
                    help='The percent of the data that is considered when standardizing the bin values of a '
                            'microbe in a sample. The algorithm chooses the densest part of the data. If a '
                            'percentage p is selected, the algorithm calculates a subset x of the vector of bins '
                            'such that max(x)-min(x) is minimal and |x| = p*length(bins). The mean and standard '
                            'deviation of this vector are calcuated and used to standardize the bin vector ('
                            'Default=85)',
                    type=float, default=85)
parser_woc.add_argument('--browser_path',
                    help='Optional; a path for the html output of SGV-Browser (Default=None, resulting in no output)',
                    type=str, default=None)
parser_woc.add_argument('--taxonomy_path',
                    help='Optional; taxonomy database path (must end in `.taxonomy.df`)',
                    type=str, default=None)
parser_woc.add_argument('--genepos_path',
                    help='Optional; gene position path (must end in `.genepos.df`)',
                    type=str, default=None)
parser_woc.add_argument('--bac_frames_path',
                    help='Optional; directory to save bacterial dataframes',
                    type=str, default=None)

# GET_SAMPLE_MAP commands
parser_gsm.add_argument('--delta_file', help = 'The ICRA .jsdel file to process')
parser_gsm.add_argument('--db_path', help = 'Path to database created by `create_db_from_reps`')
parser_gsm.add_argument('--x_coverage', help = 'The desired coverage across the genome in units of 100bp reads. This parameter is used to determine bin size: bin_size = rate_param/x_coverage (Default = 0.01)', type=float, default = 0.01)
parser_gsm.add_argument('--rate_param', help = 'The lower limit for the median number of reads per genomic bin. Genomes with coverage lower than rate_param will be discarded from the analysis (default = 10)', type = int, default = 10)
    
args = parser.parse_args()
logging.basicConfig(level=(logging.DEBUG if args.debug else logging.INFO),
                    format='%(asctime)s-%(levelname)s: %(message)s')
logging.debug(args)

def run():
	if args.command == 'get_sample_map':
		print('Running `get_sample_map')
		sample_map = get_sample_map(args.delta_file,args.db_path+'.dlen',args.x_coverage,args.rate_param) 
		print('Finished! Writing processed .smp file to %s...'%args.delta_file.replace('.jsdel','.smp'))
		pd.to_pickle(sample_map,args.delta_file.replace('.jsdel','.smp'))
	elif args.command == 'work_on_collection':
		print('Running SGVFinder `work_on_collection`...')
		vsgv, dsgv = work_on_collection(
			samp_to_map=args.samp_to_map_dir,
			max_spacing=args.max_spacing,
			min_samp_cutoff=args.min_samp_cutoff, 
			delsdetectthresh=args.dels_detect_thresh, 
			real_del_thresh=args.real_del_thresh, 
			dels_cooc_thresh=args.dels_cooc_thresh,
			vsgv_dissim_thresh=args.vsgv_dissim_thresh, 
			vsgv_clip_quantile=args.vsgv_clip_quantile, 
			vsgv_fit_interval=args.vsgv_fit_interval, 
			vsgv_fit_method=args.vsgv_fit_method,
			x_coverage=args.x_coverage,
			rate_param=args.rate_param, 
			vsgv_dense_perc=args.vsgv_dense_perc, 
			browser_path=args.browser_path, 
			taxonomypath=args.taxonomy_path,
			genepospath=args.genepos_path,
			frames_path=args.bac_frames_path
		)
		if args.csv_output:
				vsgv.to_csv(args.output_vsgv)
				dsgv.to_csv(args.output_dsgv)
		else:
				vsgv.to_pickle(args.output_vsgv)
				dsgv.to_pickle(args.output_dsgv)
	else:
		print('Unrecognized command `%s`'%args.command)
		print('Available commands are [`get_sample_map`,`work_on_collection`]')
		sys.exit(0)
       