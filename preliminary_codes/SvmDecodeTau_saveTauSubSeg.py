import matplotlib as mpl

mpl.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import itertools
import h5py
import numpy as np
import scipy.stats as scstat
import random
import pandas as pd
from Pos_file_processing import *
# import caiman.base.rois
import scipy
from scipy.io import loadmat
import pickle
import time
from sklearn import linear_model
from sklearn.preprocessing import StandardScaler
from random import randint
import random

expt = 'CircleSquare'

pathAll = '/Volumes/GoogleDrive/My Drive/Fenton_Lab/CellPress_CircleSquare/MNS_ANALYSIS/Data/' + expt + '/'

TimeBinSize_AVC = 10 * 1000
TimeBinSize_TAU = 1000

AV_plot_lim = [-0.3, 0.3]  # [-0, 0.15]
# AV_plot_lim_mean = [-0, 0.15]
AV_plot_lim_mean = [-0, 0.05]  # for z
tauPlotLimMean = [-0.05, 0.3]  # [-0.005, 0.025]
tauPlotLim = [-0.05, 1]
tauNSeg = nEqualSegs = 21
scatterSubSamp = 0.5

# tauplotLim = [-0.2, 1]
tauplotLim = [-0.5, 1]  # for Z

AvPlotSize = 2 * 5 * 60 * 1000 / TimeBinSize_AVC
# AvPlotSize = 2*5*60*1000/TimeBinSize_TAU # for Z score

threshMult = 0.02
thresh = threshMult * (TimeBinSize_TAU / 1000)

tauSubSeg = 60  # length (in dt [curr 1dt =1s] ) of smaller segments to calculate Tau multiple times

animalList = ['M29'] #, 'M20', 'M29' 'M34', 'M35', 'M39']

strTimeScale = str(int(TimeBinSize_AVC)) + 'msAV_' + str(int(TimeBinSize_TAU)) + 'msKT'

nAni = len(animalList)

plotEveryRecCorr = False
day_thresh = 3
daySepMin = 0

doPTI = False

print(doPTI)

strParam = str(int(daySepMin)) + 'daysep_' + str(int(day_thresh)) + 'dayThresh_' + str(TimeBinSize_TAU) + 'ms'

perfAll = []

for iAni, animal in enumerate(animalList):

    if doPTI:
        with open(pathAll + 'analysisFiles/' + 'Zlinear_List_' +
                  str(int(TimeBinSize_AVC)) + 'msAV_' + str(int(TimeBinSize_TAU)) + 'msKT'
                                                                                    '_CircleSquare_' + animal + '.file',
                  'rb') as f:
            outputComb = pickle.load(f)
    else:
        with open(pathAll + 'analysisFiles/' + 'TempCorr_List_' +
                  str(int(TimeBinSize_AVC)) + 'msAV_' + str(int(TimeBinSize_TAU)) + 'msKT'
                                                                                    '_CircleSquare_' + animal + '.file',
                  'rb') as f:
            outputComb = pickle.load(f)

    if doPTI:
        tauVecSingleList = outputComb[0][7][0]
        sConvTauList = outputComb[0][4][0]
        # [dayRealList, _ , _, sConvAvcList, [zConvTauList, sConvtauList], isSeparatedList,
        #   tauPairsSingleList, [zTauVecSingleList, tauVecSingleList], trainingList] = outputComb[0]
    else:
        tauVecSingleList = outputComb[0][7]
        sConvTauList = outputComb[0][4]
        # [dayRealList, _ , _, sConvAvcList, sConvtauList, isSeparatedList,
        #  tauPairsSingleList, tauVecSingleList, trainingList] = outputComb[0]
    dayCountList = outputComb[0][1]
    sessList = outputComb[0][2]
    isSeparatedList = outputComb[0][5]
    tauPairsSingleList = outputComb[0][6]

    outputComb = []

    tAllFilteredSave = []
    sessIndicAllSave = []
    tauPairsFilteredSave = []
    daySave = []

    # [trackList, _, aMaxLocList] = outputComb[1]
    # [atnSplitList, iDayList] = outputComb[2]

    # filter Tau Pairs
    # create filter for significance
    # isSignificant = [[(tPair[0][1]<0.05 and tPair[1][1]<0.05) for tPair in zip(*tauVecComb)] for tauVecComb in tauVec]

    perfAni = []

    for iDay, day in enumerate(range(1, 10)):

        print(animal)
        print(iDay)

        if day not in dayCountList: continue

        daySave.append(day)

        # find index of non-HMC recording from the 'day' defined above
        recIdx = np.where([d == day and 'HMC' not in s for d, s in zip(dayCountList, sessList)])[0]

        # get vector of correlations
        tauVec = [[t[0] for t in tauVecSingleList[idx]] for idx in recIdx]
        tauPairs = [tauPairsSingleList[idx] for idx in recIdx]
        sess = [sessList[idx] for idx in recIdx]
        sConvTau = [sConvTauList[idx] for idx in recIdx]

        ###  FOR EIGENVALUES  ###
        # filter tau values
        # filter for distanced cells
        isSepDay = isSeparatedList[recIdx[0]]
        # create filter for Nan's
        isAnyNan = np.any([np.isnan(t) for t in tauVec], 0)
        # combine filters
        filterVec = [not anyNan for isSep, anyNan in zip(isSepDay, isAnyNan)]
        # apply filter and extract tau corr
        tauVecFiltered = [[t for t, f in zip(tV, filterVec) if f] for tV in tauVec]
        tauPairsFiltered = [[t for t, f in zip(tV, filterVec) if f] for tV in tauPairs]

        ###  FOR DECODING  ###
        tAll = []
        sessIndicAll = []

        # compute tau values for rec segments of 100s
        for sDay, sessDay in zip(sConvTau, sess):
            for dt in range(len(sDay[0]) // tauSubSeg):
                s = [sCell[0 + tauSubSeg * dt:tauSubSeg + tauSubSeg * dt] for sCell in sDay]
                s_all = []
                for s_tmp in s:
                    #print(s_tmp)
                    random.shuffle(s_tmp)
                    s_all.append(s_tmp)
                    #print(s_tmp)

                tAll.append([scstat.kendalltau(c1, c2)[0] for c1, c2 in itertools.combinations(s_all, 2)])
                if 'CYL' in sessDay:
                    sessIndicAll.append(0)
                else:
                    sessIndicAll.append(1)

        tauPairs = [(x1, x2) for x1, x2 in itertools.combinations(range(len(sConvTau[0])), 2)]

        # filter tau values from these segments
        nanSegFilt = np.any(np.isnan(tAll), 0)

        filterSegVec = [anyNan or not isSep for isSep, anyNan in zip(isSepDay, nanSegFilt)]

        tAllFiltered = np.delete(tAll, np.where(filterSegVec), axis=1)
        tauPairsFiltered = np.delete(tauPairs, np.where(filterSegVec), axis=0)

        tAllFilteredSave.append(tAllFiltered)
        sessIndicAllSave.append(sessIndicAll)
        tauPairsFilteredSave.append(tauPairsFiltered)

        # check for anticorrelated pairs
        # isAntiCorr =

        # filter for rate of activity
        # rateFilter = [[(iR[0][0] > thresh and iR[0][1] > thresh and iR[1][0] > thresh and iR[1][1] > thresh)
        #                for iR in zip(*ratePair)]
        #               for ratePair in ratePairAllComb]

        tAllFiltered_tmp = np.take(tAllFiltered, np.random.permutation(tAllFiltered.shape[0]), axis=0)
        tAllFiltered_tmp2 = np.take(tAllFiltered_tmp, np.random.permutation(tAllFiltered_tmp.shape[1]), axis=1)

        # normalization by Std of each cell pair
        # tauList_norm = tauList / np.std(tauList,axis=0)

        # apply the in-built normalization for SVM decoding
        scaler = StandardScaler()
        scaler.fit(tAllFiltered_tmp2)
        tAllFiltered_norm = scaler.transform(tAllFiltered_tmp2)

        random_data = np.random.rand(2000,50)
        scaler = StandardScaler()
        scaler.fit(random_data)
        random_data_norm = scaler.transform(random_data)

        random_sess = np.round(np.random.rand(2000))

        if len(tAllFiltered[0]) > 5:
             perf, coef,roc, _ = SVM_decode(random_data_norm, random_sess, do_cv=False, cv_rep=10, cv_param=[1e-10, 10000, 1])
             perfAni.append(perf)
             print(perf)

             sortAbsCoef = np.argsort(np.abs(coef))

        perfAll.append(perfAni)
        print(perfAni)

    print('saving...')
    listsToSave = [tAllFilteredSave, tauPairsFilteredSave, sessIndicAllSave, [daySave, dayCountList]]

    if doPTI:
        with open(pathAll + 'analysisFiles_Simon/' +
                  'TauSvm_PTI_' + str(int(tauSubSeg)) + 'perSeg_' +
                  str(int(TimeBinSize_TAU)) + 'msKT' +
                  '_' + animal + '.file', 'wb') as f:
            pickle.dump(listsToSave, f, pickle.HIGHEST_PROTOCOL)
    else:
        with open(pathAll + 'analysisFiles_Simon/' +
                  'TauSvm_' + str(int(tauSubSeg)) + 'perSeg_' +
                  str(int(TimeBinSize_TAU)) + 'msKT' +
                  '_' + animal + '.file', 'wb') as f:
            pickle.dump(listsToSave, f, pickle.HIGHEST_PROTOCOL)
    print('saved!')

print(perfAll)

plt.show(block=True)
