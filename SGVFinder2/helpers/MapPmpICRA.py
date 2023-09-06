from tools.Bowtie2WrapperSlim import MapPreset, dose, dope
from tools.sam2pmp import sam2pmp
from ICRA2 import single_file, OCEAN_GENES
from os.path import dirname, basename, exists
from ICRA_BAM import add_ICRA_probs
from Utils import tryrm, mkdirifnotexists, Write, Load
from os import stat
import logging
log_ = logging.getLogger(__name__)
del logging
import os
    
default_map_param = dict(preset=MapPreset.SENSITIVE, report_alns=1, minins=0, maxins=500, 
                         no_mixed=False, no_discordant=False, dovetail=False, no_contain=False, 
                         no_overlap=False)
default_icra_param_genes = dict(max_mismatch=8, consider_lengths=True, epsilon=1e-6, \
                               max_iterations=100, min_bins=4, max_bins=100, min_reads=10, 
                               dense_region_coverage=60, length_minimum=300, \
                               length_maximum=2e5, use_theta=False, average_read_length=None, 
                               force_save_delta=True)

def do_single_mappmp(fq1, fq2, output_prefix, index_file, threads, map_param_dict):
    base_fq = basename(fq1)
    met_fname = output_prefix+'.met'
    if exists(met_fname) and stat(met_fname).st_size > 0:
        log_.info('{} already mapped'.format(base_fq))
    else:
        log_.info('Running mapping on {}'.format(base_fq))
        pe = fq2 is not None
        if pe: 
            dope(fq1, fq2, output_prefix, index_file, bam=True, local=False, 
                 metfile=met_fname, threads=threads, **default_map_param)
        else:
            map_param_se = {k:v for k,v in map_param_dict.items() if k in ['preset','report_alns']}
            dose(fq1, output_prefix, index_file, bam=True, local=False, 
                 metfile=met_fname, threads=threads, **map_param_se)
        if os.stat(met_fname).st_size > 0:
            log_.info('Mapping completed successfully')
        else:
            raise RuntimeError('Mapping of {} unsuccessful'.format(fq1))
    if exists(output_prefix + '.pmpdone') and Load(output_prefix + '.pmpdone'):
        log_.info('{} already converted to pmp'.format(output_prefix + '.bam'))
    else:
        log_.info('Converting BAM to PMP')
        sam2pmp(output_prefix+'.bam', output_prefix+'.pmp',full=True)
        Write(output_prefix + '.pmpdone', True)
        log_.info('Conversion completed')


def do_single_icrabam(fq1, fq2, output_prefix, index_file, icra_usage,  
              icra_param_dict, remove_unmapped=True, remove_not_delta=False, delta_thresh=0.999,
              delete_pmp=False, delete_old_bam=True):
    jsdel_f = output_prefix + '.jsdel'
    if exists(jsdel_f) and stat(jsdel_f).st_size > 0:
        log_.info('ICRA was run successfully')
    else:
        log_.info('Running ICRA...')
        single_file(fq1, fq2, outfol=dirname(output_prefix), usage=icra_usage, 
                    pmpf=output_prefix+'.pmp', sam_based=True, **icra_param_dict)
    
    if exists(output_prefix + '.icrabamdone'):
        log_.info('ICRA update performed... nothing to do here')
    else:
        add_ICRA_probs(output_prefix+'.jsdel', output_prefix+'.bam', output_prefix+'.icra.bam', 
                       remove_unmapped=remove_unmapped, remove_not_in_delta=remove_not_delta, 
                       delta_thresh=delta_thresh)
        Write(output_prefix + '.icrabamdone', True)
        log_.info('ICRA update completed')
    
    log_.info('Cleaning up...')
    if delete_pmp:
        tryrm(output_prefix+'.pmp')            
    if delete_old_bam:
        tryrm(output_prefix+'.bam')