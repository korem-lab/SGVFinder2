import logging
from helpers.sam2pmp import sam2pmp
import sys

args = sys.argv
outpref = args[1]
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s-%(levelname)s: %(message)s')
logging.debug(args)
log_ = logging.getLogger('ICRA')
log_.debug("{} Converting to pmp".format(outpref))

sam2pmp(outpref + '.bam', outpref + '.pmp', full=True)

log_.debug('PMP ready')