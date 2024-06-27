from .icra import single_file
from .svfinder import get_sample_map, work_on_collection
from .createdb import create_db_from_reps
from distutils.spawn import find_executable
from .helpers.Bowtie2WrapperSlim import do_pair_simple
from .helpers.sam2pmp import sam2pmp

assert find_executable('bowtie2') is not None, 'Cannot find `bowtie2` installation'
assert find_executable('samtools') is not None, 'Cannot find `samtools` installation'