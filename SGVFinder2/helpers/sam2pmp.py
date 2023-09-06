import sys
sys.path.insert(0,'cy_ext/')

from collections import namedtuple
import pysam
from itertools import groupby
from operator import attrgetter
import logging
import gzip
import ujson
import sam2pmp_helper as s2ph
from .ICRAUtils import timeit

log_ = logging.getLogger('sam2pmp')
del logging


@timeit(log_)
def sam2pmp(samfile, pmpfile, full=True):
    unmapped = 0
    save = pysam.set_verbosity(0)
    f = pysam.AlignmentFile(samfile)
    pysam.set_verbosity(save)
    with gzip.open(pmpfile, 'wt') as of:
        for _, grp in groupby(f, attrgetter('query_name')):
            alngrp = list(grp)
            if _readunmapped(alngrp):
                unmapped += 1
                continue
            curr = SourceReadSAM.from_alngroup(alngrp, full)
            if len(curr) > 0:
                curr.sort()
                ujson.dump(curr.to_ser(), of)
                of.write('\n')


def _get_map_quality(aln):
    alp = aln.get_aligned_pairs()
    ref = aln.get_reference_sequence()
    quals = aln.query_qualities
    return s2ph.calc_quality(alp, quals, ref)  # @UndefinedVariable


def _readunmapped(alngrp):
    return len(alngrp) == 1 and alngrp[0].is_unmapped \
           or len(alngrp) == 2 and alngrp[0].is_unmapped and alngrp[1].is_unmapped


DestMap_se = namedtuple('DestMap_se', ['dest_id', 'mismatches', 'strand', 'pos', 'cigar',
                                       'end', 'mismatch_pos', 'alnscore', 'map_quality'])
DestMap_pe = namedtuple('DestMap_pe', ['dest_id', 'mismatches1', 'mismatches2', 'strand1',
                                       'strand2', 'pos1', 'pos2', 'cigar1', 'cigar2',
                                       'mismatch_pos1', 'mismatch_pos2', 'alnscore1',
                                       'alnscore2', 'map_quality1', 'map_quality2'])
DestMap_se_part = namedtuple('DestMap_se_part', ['dest_id', 'mismatches', 'strand', 'pos', 'end',
                                                 'map_quality'])
DestMap_pe_part = namedtuple('DestMap_pe_part', ['dest_id', 'mismatches1', 'mismatches2', 'strand1',
                                                 'strand2', 'pos1', 'pos2', 'map_quality1',
                                                 'map_quality2'])


def _get_nm_md_as(aln):
    # NM - edit distance (mismatches); MD - mismatch positions string; AS - alignment score
    retdct = dict(aln.tags)
    return retdct['NM'], retdct['MD'], retdct['AS']


def _get_strand(aln):
    return '-' if aln.is_reverse else '+'


def _get_pe_end(aln):
    return 0 if aln.is_read1 else 1


def _se_from_aln(aln, unique, full):
    qual = 1 if unique else _get_map_quality(aln)
    refname = aln.reference_name.split()[0]
    if refname == 'SAMN02603102.565034.SAMN02603102.CP001357':
        refname = 'SAMN02603102.565034.SAMN02603102.CP001357\t<annotation:_organism="Serpulina_phage_VSH-1"_taxon="58620"/>'
    if full:
        nm, md, asc = _get_nm_md_as(aln)
        return DestMap_se(refname, nm, _get_strand(aln), aln.pos,
                          aln.cigarstring, _get_pe_end(aln), md, asc, qual)
    return DestMap_se_part(refname, dict(aln.tags)['NM'], _get_strand(aln),
                           aln.pos, _get_pe_end(aln), qual)


def _pe_from_aln(aln1, aln2, unique, full):
    qual1 = 1 if unique else _get_map_quality(aln1)
    qual2 = 1 if unique else _get_map_quality(aln2)
    if full:
        nm1, md1, asc1 = _get_nm_md_as(aln1)
        nm2, md2, asc2 = _get_nm_md_as(aln2)
        return DestMap_pe(aln1.reference_name.split()[0], nm1, nm2,
                          _get_strand(aln1), _get_strand(aln2), aln1.pos, aln2.pos,
                          aln1.cigarstring, aln2.cigarstring, md1, md2, asc1, asc2,
                          qual1, qual2)
    return DestMap_pe_part(aln1.reference_name.split()[0], dict(aln1.tags)['NM'],
                           dict(aln2.tags)['NM'], _get_strand(aln1), _get_strand(aln2),
                           aln1.pos, aln2.pos, qual1, qual2)


class SourceReadSAM(object):
    def __init__(self, rid, seq, qual, full=True):
        self.rid = rid
        self.seq = [str(s) for s in seq]
        self.quals = qual
        self.full = full
        self.se_maps = []
        self.pe_maps = []

    def to_ser(self):
        if self.full:
            return (self.rid, self.seq, self.quals, [tuple(x) for x in self.se_maps],
                    [tuple(y) for y in self.pe_maps])
        else:
            return (self.rid, self.seq, self._letterqual(), [tuple(x) for x in self.se_maps],
                    [tuple(y) for y in self.pe_maps])

    def _letterqual(self):
        return [''.join([chr(i + 33) for i in q]) for q in self.quals]

    @staticmethod
    def _numberqual(letterqual):
        return [[ord(c) - 33 for c in st] for st in letterqual]

    @classmethod
    def from_alngroup(cls, alngrp, full):
        fst = alngrp[0]
        res = cls(fst.query_name, [fst.query_sequence], [fst.query_qualities.tolist()], full)
        read1 = None
        seen_read1 = fst.is_read1
        seen_read2 = fst.is_read2
        unique = (len(alngrp) == 1) \
                 or ((len(alngrp) == 2) and (alngrp[0].is_read1 ^ alngrp[1].is_read1))
        for aln in alngrp:
            if not seen_read1 and aln.is_read1:
                res.seq.insert(0, aln.query_sequence)
                res.quals.insert(0, aln.query_qualities.tolist())
                seen_read1 = True
            if not seen_read2 and aln.is_read2:
                res.seq.append(aln.query_sequence)
                res.quals.append(aln.query_qualities.tolist())
                seen_read2 = True
            if aln.is_proper_pair:
                if aln.is_read1:
                    read1 = aln
                elif aln.is_read2:
                    assert (read1.reference_name == aln.reference_name) \
                           and (read1.is_reverse ^ aln.is_reverse)
                    res.pe_maps.append(_pe_from_aln(read1, aln, unique, full))
                    read1 = None
                else:
                    raise RuntimeError('huh?')
            else:
                if not aln.is_unmapped:
                    res.se_maps.append(_se_from_aln(aln, unique, full))
        return res

    @staticmethod
    def from_ser(ser):
        if type(ser[2][0]) == str:
            res = SourceReadSAM(ser[0], ser[1], SourceReadSAM._numberqual(ser[2]), True)
            res.se_maps = [DestMap_se_part(*x) for x in ser[3]]
            res.pe_maps = [DestMap_pe_part(*y) for y in ser[4]]
        else:
            res = SourceReadSAM(ser[0], ser[1], ser[2], False)
            res.se_maps = [DestMap_se(*x) for x in ser[3]]
            res.pe_maps = [DestMap_pe(*y) for y in ser[4]]
        return res

    def __len__(self):
        return len(self.se_maps) + len(self.pe_maps)

    def sort(self):
        self.se_maps.sort(key=attrgetter('map_quality'), reverse=True)
        self.pe_maps.sort(key=lambda m: m.map_quality1 + m.map_quality2, reverse=True)
