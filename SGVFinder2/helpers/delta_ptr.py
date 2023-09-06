from lmfit import minimize, Parameters
from lmfit import __version__ as lmfitv
from numpy import asarray, sqrt, add
from numpy.random import permutation
from pandas import DataFrame, Series, concat, \
    read_pickle, to_pickle
import numpy as np
from os.path import splitext, split, join, basename
from glob import glob
import logging
import ujson
import gzip

NEWLMFIT = True
__version__ = '2.0.0'
log_ = logging.getLogger('PTRC')


class genomelsqfit(object):
    def __init__(self, genomelength):
        self._params = Parameters()
        self._params.add('genomelength', value=genomelength, vary=False)
        self._params.add('troughx', value=genomelength / 2, min=0, max=genomelength)
        self._params.add('troughy', value=10, min=0)
        self._params.add('peakx', value=0, min=0, max=genomelength)
        self._params.add('peaky', value=10, min=0)
        self.genomelength = genomelength
        self._pvalue = None

    @property
    def troughx(self):
        return self._params['troughx'].value

    @property
    def troughy(self):
        return self._params['troughy'].value

    @property
    def peakx(self):
        return self._params['peakx'].value

    @property
    def peaky(self):
        return self._params['peaky'].value

    def bilinfunc(self, x):
        return _bilinfunc(x, self.genomelength, self.troughx, self.troughy, self.peakx, self.peaky)

    def dofit(self, x, data):
        self._lastx = x
        self._lastdata = data
        self._params['peaky'].value = data.mean() + data.std()
        self._params['troughy'].value = data.mean() - data.std()
        res = minimize(_getresiduals, self._params, args=(x, data))
        if NEWLMFIT:
            self._params = res.params
        if self.peaky < self.troughy:
            self._params['peakx'], self._params['troughx'] = self._params['troughx'], self._params['peakx']
            self._params['peaky'], self._params['troughy'] = self._params['troughy'], self._params['peaky']
        return res

    def determination(self):
        yhat = self.bilinfunc(self._lastx)
        ybar = np.average(self._lastdata)
        ssreg = sum((yhat - ybar) ** 2)
        sstot = sum((self._lastdata - ybar) ** 2)
        if sstot == 0:
            return np.nan
        return ssreg / sstot

    def pValue(self, num_permutations=100):
        if not self._pvalue:
            num_permuted_wins = 0
            selfdet = self.determination()
            if selfdet == np.nan:
                return 1
            for _ in range(num_permutations):
                res = genomelsqfit(self.genomelength)
                res.dofit(self._lastx, permutation(self._lastdata))
                resdet = res.determination()
                if resdet == np.nan or resdet >= selfdet:
                    num_permuted_wins += 1
            self._pvalue = 1 if num_permutations == 0 else (max(1, num_permuted_wins) / float(num_permutations))
        return self._pvalue


class SampleMap(object):
    def __init__(self, init_res=100):
        self.bacid_maps = {}
        self.res = init_res
        self.modbacid_maps = None
        self.modres = None

    def adddelta(self, destid, pos1, delta, destlength):
        if destid not in self.bacid_maps:
            self.bacid_maps[destid] = \
                Series(np.zeros((destlength // self.res) + \
                                (1 if destlength % self.res != 0 else 0)))
        self.bacid_maps[destid][pos1 // self.res] += delta

    @staticmethod
    def FromDeltaFile(filename, destid_length):
        smap = SampleMap()
        with gzip.open(filename, 'rt') as inf:
            deltl = ujson.load(inf)[1]
        for rid, mapngs in deltl:
            for dest_id, pos1, pos2, koef, q_koef in mapngs:
                if dest_id in destid_length:
                    smap.adddelta(dest_id, max(pos1, 0), koef, destid_length[dest_id])
        return smap

    def sumevery(self, x, sumwindow, slide, covper):
        x = concat((x[-sumwindow // 2:], x, x[:sumwindow // 2])).rolling(window=sumwindow,
                                                                       min_periods=int(covper*sumwindow),
                                                                       center=True).sum()
        x = x[sumwindow // 2:-sumwindow // 2]
        # x = rolling_sum(concat((x[-sumwindow / 2:], x, x[:sumwindow / 2])),
        #                 window=sumwindow, min_periods=int(covper * sumwindow),
        #                 center=True)[sumwindow / 2:-sumwindow / 2]
        # Fix for floating point inaccuracy
        x[(x < 0) & (x > -1e-10)] = 0
        return x[list(range(0, x.shape[0], slide))]

    def medianevery(self, x, window, slide, covper):
        x = concat((x[-window // 2:], x, x[:window // 2])).rolling(window=window, min_periods=int(covper * window),
                                                                 center=True).median()
        x = x[window // 2:-window // 2]
        # x = rolling_median(concat((x[-window / 2:], x, x[:window / 2])),
        #                    window=window, min_periods=int(covper * window),
        #                    center=True)[window / 2:-window / 2]
        return x[list(range(0, x.shape[0], slide))]

    def mff(self, x, medianfoldfilter):
        x = x.copy()
        med = x.median()
        mn = med / sqrt(medianfoldfilter)
        mx = med * sqrt(medianfoldfilter)
        x[(x < mn) | (x > mx)] = np.nan
        return x

    def ManipulateData(self, sumwindow=1, sumslide=1, medianwindow=1, medianslide=1, medianfoldfilter=None,
                       minimalPercentageOfCoverage=0.6, minimalbincount=10, mediancoveragefilter=0):
        self.modbacid_maps = {}
        self.modres = self.res * (sumslide if sumwindow != 1 else 1) * \
                      (medianslide if medianwindow != 1 else 1)
        for ky, bacvec in self.bacid_maps.items():
            bacmap = bacvec.copy()
            bacmap = self.sumevery(bacmap, sumwindow, sumslide, minimalPercentageOfCoverage)
            bacmap = self.mff(bacmap, medianfoldfilter)
            bacmap = self.medianevery(bacmap, medianwindow, medianslide, minimalPercentageOfCoverage)
            bacmap.index = bacmap.index * self.res
            cbins = (bacmap > 0).sum()
            if bacmap.shape[0] > 0 and \
                    cbins / float(bacmap.shape[0]) >= minimalPercentageOfCoverage \
                    and cbins >= minimalbincount \
                    and bacmap.median() >= mediancoveragefilter:
                bacmap = bacmap.dropna()
                self.modbacid_maps[ky] = bacmap


def _bilinfunc(x, genomelength, troughx, troughy, peakx, peaky):
    assert all((x >= 0) & (x < genomelength))
    a = (troughy - peaky) / (troughx - peakx)
    levx = np.array(
        [troughx if b else peakx for b in (troughx < peakx and x < peakx) | (peakx < troughx and x >= troughx)])
    levy = np.array([troughy if b else peaky for b in \
                     (troughx < peakx and x < peakx) | (peakx < troughx and x >= troughx)])
    moda = np.array([1 if b else -1 for b in \
                     (troughx < peakx and ((troughx < x) & (x < peakx))) | (
                             peakx < troughx and ((peakx < x) & (x < troughx)))])
    return np.multiply(moda * a, x) + (levy - np.multiply(moda * a, levx))


def _troughlimfunct(peakx, troughx, genomelength):
    ret = troughx
    fact = -1 if troughx < peakx else 1
    mn = 0.45 * genomelength
    mx = 0.55 * genomelength
    dif = abs(peakx - troughx)
    if dif > mx:
        ret = peakx + fact * mx
    elif dif < mn:
        ret = peakx + fact * mn
    return ret % genomelength


def _getresiduals(params, x, data):
    params['troughx'].value = _troughlimfunct(params['peakx'].value,
                                              params['troughx'].value,
                                              params['genomelength'].value)
    return data - _bilinfunc(x, params['genomelength'].value,
                             params['troughx'].value, params['troughy'].value,
                             params['peakx'].value, params['peaky'].value)


def _get_sample_predictions(smap, threshold, pvcutoff, destid_length):
    predictions = {}
    for bac, ser in smap.modbacid_maps.items():
        ser = ser.dropna()
        ftr = genomelsqfit(destid_length[bac])
        try:
            res = ftr.dofit(ser.index.values, np.log2(ser.values))
        except:
            continue

        if ((threshold and (2 ** ftr.peaky) / (2 ** ftr.troughy) >= threshold)
            or not threshold) and ftr.pValue(1000) <= pvcutoff:
            predictions[bac] = (ftr.troughx, ftr.peakx, ftr.troughy, ftr.peaky, ftr.pValue())
    return predictions


def _merge_preds(dcts):
    retval = {}
    for curdic in dcts:
        for k, v in curdic.items():
            retval.setdefault(k, []).append(v)
    return retval


def _circularmedianbasedpredictioncensus(pks, trs, lenx):
    bstrng = lenx * 2
    bstpk = None
    bsttr = None
    for t in pks:
        tmpx = np.mod(np.array(pks) - t, lenx)
        rng = tmpx.max() - tmpx.min()
        if rng < bstrng:
            bstrng = rng
            bstpk = (np.median(tmpx) + t) % lenx
            bsttr = (np.median(np.mod(np.array(trs) - t, lenx)) + t) % lenx
    return bsttr, bstpk


def _circular_median_census_on_preds(predsraw, destid_lengths, minpredcount):
    preds = {}
    for ky in predsraw.keys():
        if len(predsraw[ky]) >= minpredcount:
            trs, pks, _, _, _ = zip(*predsraw[ky])
            preds[ky] = _circularmedianbasedpredictioncensus(pks, trs, destid_lengths[ky])
    return preds


def _gety(ser, offset, res):
    ser = ser.dropna()
    offset = int(offset)
    if offset in ser:
        result = ser[offset]
        if result != 0:
            return result
    vicinity = []
    if (offset + res) in ser:
        vicinity.append(ser[offset + res])
    if (offset - res) in ser:
        vicinity.append(ser[offset - res])
    if len(vicinity) == 2:
        return np.average(vicinity)
    elif len(vicinity) == 0:
        return None
    if (offset + 2 * res) in ser:
        vicinity.append(ser[offset + 2 * res])
    if (offset - 2 * res) in ser:
        vicinity.append(ser[offset - 2 * res])
    if len(vicinity) > 1:
        return np.average(vicinity)
    else:
        return None


def _calculateptr(smap, preds, tr_thresh=0):
    p2t = {}
    for bac, ser in smap.modbacid_maps.items():
        if bac in preds:
            # extra safe for no apparent reason
            tr = np.round(np.round((preds[bac][0] / smap.modres)) * smap.modres)
            pk = np.round(np.round((preds[bac][1] / smap.modres)) * smap.modres)
            pk_y = _gety(ser, pk, smap.modres)
            if pk_y:
                pk_y = np.round(pk_y, 6)
            tr_y = _gety(ser, tr, smap.modres)
            if tr_y:
                tr_y = np.round(tr_y, 6)
            if tr_y and pk_y and tr_y != 0 and tr_y >= tr_thresh:
                p2t[bac] = (pk_y if pk_y != 0 else 1) / tr_y
    return p2t


def _ptrframefromdct(dct):
    res = DataFrame(dct)
    try:
        res.drop('unclassified', inplace=True)
    except:
        pass
    return res.T


def coverage_analysis_fromEM(deltaf, outfol, cov_thresh, destid_length):
    smap = SampleMap.FromDeltaFile(deltaf, destid_length)
    smap.ManipulateData(100, 1, 10000, 100, 8, 0.6, 10, cov_thresh)
    preds = _get_sample_predictions(smap, 1.1, 0.05, destid_length)
    to_pickle((smap, preds), join(outfol, splitext(split(deltaf)[1])[0] + '.ptr'))


def calculate_ptrs_strain(inputfolder, outpath, destid_length, minpredcount=3,
                          tr_thresh=5, opreds=None, csv_output=False):
    predses = _merge_preds([read_pickle(f)[1] \
                            for f in glob(join(inputfolder, '*.ptr'))])
    preds = _circular_median_census_on_preds(predses, destid_length, minpredcount)
    if opreds is not None:
        to_pickle(preds, opreds)
    res = _ptrframefromdct({splitext(basename(f))[0]: _calculateptr(read_pickle(f)[0], preds, tr_thresh)
                            for f in glob(join(inputfolder, '*.ptr'))})
    if csv_output:
        res.to_csv(outpath)
    else:
        res.to_pickle(outpath)
