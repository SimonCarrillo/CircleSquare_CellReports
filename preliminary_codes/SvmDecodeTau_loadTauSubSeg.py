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
import caiman.base.rois
import scipy
from scipy.io import loadmat
import pickle
import time
from sklearn import linear_model

expt = 'CircleSquare'

pathAll = '/Users/eliott/Dropbox/dDocuments/NYU/Fenton/MNS_ANALYSIS/Data/'+ expt +'/'
# pathFigures = '/Users/eliott/Dropbox/dDocuments/NYU/Fenton/MNS_ANALYSIS/Figures/' + expt + '/savedFigs/'

saveFig = False

BinNb = 12  # for loading place maps
TimeBinSize_AVC = 1 * 1000
TimeBinSize_TAU = 1000

cvRep = 100

AV_plot_lim = [-0.3, 0.3] #[-0, 0.15]
# AV_plot_lim_mean = [-0, 0.15]
AV_plot_lim_mean = [-0, 0.05] # for z
tauPlotLimMean = [-0.05, 0.3]#[-0.005, 0.025]
tauPlotLim = [-0.05, 1]
tauNSeg = nEqualSegs = 21
scatterSubSamp = 0.5

# tauplotLim = [-0.2, 1]
tauplotLim = [-0.5, 1] # for Z

AvPlotSize = 2*5*60*1000/TimeBinSize_AVC
# AvPlotSize = 2*5*60*1000/TimeBinSize_TAU # for Z score

threshMult = 0.02
thresh = threshMult * (TimeBinSize_TAU / 1000)

tauSubSeg = 60 # length (in dt [curr 1dt =1s] ) of smaller segments to calculate Tau multiple times

animalList = ['M19', 'M35', 'M20', 'M29', 'M34', 'M39']

strTimeScale = str(int(TimeBinSize_AVC)) + 'msAV_' + str(int(TimeBinSize_TAU)) + 'msKT'
# pdffilepath = pathFigure + 'loadTemporal_thresh' + str(100*threshMult) + '_' + strTimeScale + '.pdf'

nAni = len(animalList)

plotEveryRecCorr = False
day_thresh = 3
daySepMin = 0

doPTI = True
doPTICells = True

if doPTI: strType = 'PTI_'
else: strType = 'rate_'

if doPTI:
    pathFigure = '/Users/eliott/Dropbox/dDocuments/NYU/Fenton/MNS_ANALYSIS/Figures/' + expt + '/savedFigs/z/' # for z
else:
    pathFigure = '/Users/eliott/Dropbox/dDocuments/NYU/Fenton/MNS_ANALYSIS/Figures/' + expt + '/savedFigs/'

print(doPTI)
print(doPTICells)
print(pathFigure)

# bins to cut tau Vec
if doPTICells:
    tauBins = [-0.2, -0.1, -0.06, -0.04, -0.03, -0.025, -0.02, -0.0175, -0.015, -0.0125, -0.01, -0.005,
               0.005, 0.01, 0.02, 0.03, 0.04, 0.06, 0.1, 0.2, 0.35]
else:
    tauBins = [-0.1, -0.06, -0.04, -0.03, -0.025, -0.02, -0.0175, -0.015, -0.0125, -0.01, -0.005,
               0.005, 0.01, 0.02, 0.03, 0.04, 0.06, 0.1, 0.2]
tauPlot = [tauBins[0] - (tauBins[1] - tauBins[0]) / 2]
tauPlot.extend([b + d / 2 for b, d in zip(tauBins, np.diff(tauBins))])
tauPlot.extend([tauBins[-1] - (tauBins[-2] - tauBins[-1]) / 2])

strParam = str(int(daySepMin)) +'daysep_'+ str(int(day_thresh)) +'dayThresh_' +  str(TimeBinSize_TAU) +'ms'

perfAll = []
perfSegAll = []
percCellSegAll = []
distCellDiff = []
distTauDiff = []
ratioPc = []

tauSvmMeanAll = [[], [], []]

aTauSvmMean = plt.subplots(nAni, 9)[1]
aTauSvmMean2 = plt.subplots(nAni, 9)[1]
aTauSvm = plt.subplots(nAni, 9)[1]

for iAni, animal in enumerate(animalList):

    print(animal)

    # load Tau/PTI-Tau computed on short segments of recording
    if doPTI:
        with open(pathAll + 'analysisFiles/' +
                  'TauSvm_PTI_' + str(int(tauSubSeg)) +'perSeg_' +
                  str(int(TimeBinSize_TAU)) + 'msKT' +
                  '_' + animal + '.file', 'rb') as f:
            savedList = pickle.load(f)
    else:
        with open(pathAll + 'analysisFiles/' +
                  'TauSvm_' + str(int(tauSubSeg)) +'perSeg_' +
                  str(int(TimeBinSize_TAU)) + 'msKT' +
                  '_' + animal + '.file', 'rb') as f:
            savedList = pickle.load(f)

    tAllFilteredSave, tauPairsFilteredSave, sessIndicAllSave, daySaveTmp = savedList
    #sessIndicAllSave: CYL = 0, RCT = 1
    daySave, dayCountList = daySaveTmp

    if doPTICells:
        with open(pathAll + 'analysisFiles/' + 'Zlinear_List_' +
                  str(int(TimeBinSize_AVC)) + 'msAV_' + str(int(TimeBinSize_TAU)) + 'msKT'
                  '_CircleSquare_' + animal +'.file', 'rb') as f:
            outputComb = pickle.load(f)
    else:
        with open(pathAll + 'analysisFiles/' + 'TempCorr_List_' +
                  str(int(TimeBinSize_AVC)) + 'msAV_' + str(int(TimeBinSize_TAU)) + 'msKT'
                  '_CircleSquare_' + animal +'.file', 'rb') as f:
            outputComb = pickle.load(f)
    if doPTICells:
        # PTItauVecSingleList = outputComb[0][7][0]
        tauVecSingleList = outputComb[0][7][0]
        sConvTauList = outputComb[0][4][1]
        # [dayRealList, _ , _, sConvAvcList, [zConvTauList, sConvTauList], isSeparatedList,
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

    # load PC analysis
    loadpath = pathAll + 'analysisFiles/' + 'MapCorr_' + str(int(BinNb)) + '_Bins_' + animal + '.file'
    with open(loadpath, 'rb') as f:
        output = pickle.load(f)
    dayRealListMap = output[1][0]
    dayCountListMap = output[1][1]
    sessListMap = output[1][2]
    PcIdxList = output[1][4]
    output = []

    # [trackList, _, aMaxLocList] = outputComb[1]
    # [atnSplitList, iDayList] = outputComb[2]

    # filter Tau Pairs
    #create filter for significance
    # isSignificant = [[(tPair[0][1]<0.05 and tPair[1][1]<0.05) for tPair in zip(*tauVecComb)] for tauVecComb in tauVec]

    fSvm, aSvm = plt.subplots(1, 9)
    fSvmNotOrdered, aSvmNotOrdered = plt.subplots(1, 9)
    fTau, aTau = plt.subplots(3, 2)
    fRatio, aRatio = plt.subplots(9,3)

    percList10 = [[x, x + 10] for x in [90, 80, 70, 60, 50, 40, 30, 20, 10, 0]]
    # percList5 = [[x, x + 5] for x in [95, 85, 75, 65, 55, 40, 30, 20, 10, 0]]

    distTau = [[] for _ in percList10]
    distCell = [[] for _ in percList10]
    randDistCell = []
    perfSeg = [[np.nan for _ in percList10] for d in range(9)]
    percCellSeg = [[np.nan for _ in percList10] for d in range(9)]
    perfAni = [np.nan for d in range(9)]
    ratioPcAni = []
    daySaveRatio = []

    for iDay, tauList, pairList, sessIndic in zip(daySave, tAllFilteredSave, tauPairsFilteredSave, sessIndicAllSave):
        print(iDay)
        daySaveRatio.append(iDay)
        if len(tauList[0])>25 and len(tauList)>10 :
            perf, coef, _ = SVM_decode(tauList, np.array(sessIndic), do_cv=True, cv_rep=cvRep, cv_param=[1e-10, 10000, 0])
            perfAni[iDay-1] = perf
            aSvm[iDay-1].imshow([np.transpose(tauList)[x] for x in np.argsort(coef)], aspect='auto',
                              vmin=-0.2, vmax=0.5)
            zPlot = aSvmNotOrdered[iDay-1].imshow(np.transpose(tauList), aspect='auto', vmin=-0.2, vmax=0.5)
            plt.colorbar(zPlot, ax=aSvmNotOrdered[iDay-1])

            print('day' +  str(iDay))
            print( str(len(tauList[0])) + 'cell pair correlations')

            aTauSvmMean[iAni][iDay-1].hist2d(np.mean(tauList[np.where([s == 0 for s in sessIndic])[0]], 0), coef/np.std(coef))#, s=1, color='b')
            # aTauSvmMean[iAni][iDay-1].hist2d(np.mean(tauList[np.where([s == 0 for s in sessIndic])[0]], 0)/np.std(np.mean(tauList[np.where([s == 0 for s in sessIndic])[0]], 0))
            #                                  , coef)#, s=1, color='r')

            tauSvmMeanAll[0].extend(np.mean(tauList[np.where([s == 0 for s in sessIndic])[0]], 0))
            tauSvmMeanAll[1].extend(np.mean(tauList[np.where([s == 1 for s in sessIndic])[0]], 0))
            tauSvmMeanAll[2].extend(coef/np.std(coef))

            aTauSvmMean2[iAni][iDay - 1].hist2d(tauSvmMeanAll[0], tauSvmMeanAll[2], [200,200])

            for s, t in zip(sessIndic, tauList):
                if s > 0:
                    aTauSvm[iAni][iDay-1].scatter(t, coef, s=1, color='b')
                else:
                    aTauSvm[iAni][iDay-1].scatter(t, coef, s=1, color='r')

            # sortAbsCoef = np.argsort(np.abs(coef))
            isSegAbsWeight = [[np.percentile(np.abs(coef), p[0]) < np.abs(x) <= np.percentile(np.abs(coef), p[1])
                            for x in coef] for p in percList10]
            # sortCoef = np.argsort(coef)
            isSegWeight = [[np.percentile(coef, p[0]) < x <= np.percentile(coef, p[1])
                            for x in coef] for p in percList10]
            # isSegWeight = [[np.percentile(np.abs(coef), p) < np.abs(x)
            #                 for x in coef] for p in [99, 95, 90, 80, 50]]
            totCells = len(np.unique(pairList))

            # do random draw to compute a random dist of distCell
            randDistCellTmp = []
            for nRepeat in range(100):
                nPairs = np.sum(isSegAbsWeight[0])
                randIdx = np.random.permutation(np.arange(len(pairList)))[:nPairs]
                randPairs = [pairList[x] for x in randIdx]
                randDistCellTmp.extend([np.sum(randPairs == x) / (totCells-1) for x in np.unique(randPairs)])
            randDistCell.append(randDistCellTmp)

            for iSeg, isWeight in enumerate(isSegAbsWeight):
                tAllFiltered20 = np.delete(tauList, np.where([not x for x in isWeight]), axis=1)
                tauPairsFiltered20 = np.delete(pairList, np.where([not x for x in isWeight]), axis=0)

                tauPairsFiltered20NonAbs = np.delete(pairList, np.where([not x for x in isSegWeight[iSeg]]), axis=0)

                perf, coef, _ = SVM_decode(tAllFiltered20, np.array(sessIndic), do_cv=True, cv_rep=cvRep,
                                           cv_param=[1e-10, 10000, 0])

                perfSeg[iDay-1][iSeg] = perf
                percCellSeg[iDay-1][iSeg] = len(np.unique(tauPairsFiltered20))/totCells

                distTau[iSeg].extend([t_ for t in tAllFiltered20 for t_ in t])
                distCell[iSeg].extend([np.sum(tauPairsFiltered20 == x) / (totCells-1) for x in np.unique(tauPairsFiltered20)])

                if iSeg == 0:
                    threshDistCell = np.mean(randDistCellTmp)+1.5*np.std(randDistCellTmp)
                    impttCells = [x for x in np.unique(tauPairsFiltered20)
                                  if np.sum(tauPairsFiltered20 == x) / totCells > threshDistCell]
                    impttCellsHigh = [x for x in np.unique(tauPairsFiltered20NonAbs)
                                  if np.sum(tauPairsFiltered20NonAbs == x) / totCells > threshDistCell]
                if iSeg == len(isSegAbsWeight)-1:
                    impttCellsLow = [x for x in np.unique(tauPairsFiltered20NonAbs)
                                      if np.sum(tauPairsFiltered20NonAbs == x) / totCells > threshDistCell]


            # PC ANALYSIS
            # find index of non-HMC recordings from the 'day' defined above to extract PC-categories
            recIdxMap = np.where([d == iDay and 'HMC' not in s for d, s in zip(dayCountListMap, sessListMap)])[0]
            PcIdx = np.any([PcIdxList[rIdx] for rIdx in recIdxMap], 0)

            ratioPcAni.append(np.mean([PcIdx[x] for x in impttCells]) / np.mean(PcIdx))

            # find index of non-HMC recording from the 'day' defined above to extract tau corr
            recIdx = np.where([d == iDay and 'HMC' not in s for d, s in zip(dayCountList, sessList)])[0]
            for rIdx in recIdx:
                # get vector of correlations
                tauVec = [t[0] for t in tauVecSingleList[rIdx]]
                tauPairs = tauPairsSingleList[rIdx]
                sess = sessList[rIdx]
                sConvTau = sConvTauList[rIdx]

                # filter tau values
                # filter for distanced cells
                isSepDay = isSeparatedList[rIdx]
                # create filter for Nan's
                isAnyNan = [np.isnan(t) for t in tauVec]
                # # filter out cells that are rarely active
                meanRate = [np.mean([s > 0 for s in c]) for c in sConvTau]
                rate = [meanRate[p[0]] > thresh and meanRate[p[1]] > thresh for p in tauPairs]

                # create filter for spatial separation and significance and non nan correlation
                tauFilter = [sep and r and not iNan for sep, iNan, r in zip(isSepDay, isAnyNan, rate)]
                meanRate = rate = isAnyNan = isSeparated = []

                # filter tau values and pair numbers
                tauVec = [t for t,f in zip(tauVec, tauFilter) if f]
                idxTauPairs = [idx for idx,f in zip(tauPairs, tauFilter) if f]

                keptCells = np.unique([x for pair in idxTauPairs for x in pair])
                # svmCoef = [cCv[iCell] for iCell in keptCells]
                nKeptCells = len(keptCells)

                nBins = len(tauBins) + 1
                binList = [*range(nBins)]
                # sort pairs in different segments
                equalSegIdx = np.digitize(tauVec, tauBins)  # a vector containing the bin number for each pair

                # for each segment, list the pairs it contains
                equalSegPairList = [[idxTauPairs[iPair] for iPair in np.where(equalSegIdx == iBin)[0]]
                                    for iBin in binList]

                # aBin.plot(tauPlot, [len(x) for x in equalSegPairList], '.')

                # evaluate how many times each cell is in each segment
                segPairsCountTmp = [[np.sum([np.any([p[0] == i, p[1] == i]) for p in pairsSeg]) for i in keptCells]
                                    for pairsSeg in equalSegPairList]

                # normalize for the number of pairs a cell is in
                segPairsCount = [x / np.sum(segPairsCountTmp, 0) for x in segPairsCountTmp]

                if 'CYL' in sess:
                    # impttCellsNonAbs = impttCellsLow
                    linCol = 'b'
                else:
                    # impttCellsNonAbs = impttCellsHigh
                    linCol = 'r'
                ratio = [np.mean([s[np.where(keptCells == x)][0] for x in impttCellsLow if
                                  np.size(np.where(keptCells == x)) > 0]) / np.mean(s) for s in segPairsCount]
                aRatio[iDay - 1][0].plot(tauPlot, [1 for _ in tauPlot], 'k--')
                aRatio[iDay - 1][0].plot(tauPlot,ratio, '.-', color=linCol)

                ratio = [np.mean([s[np.where(keptCells == x)][0] for x in impttCellsHigh if
                                  np.size(np.where(keptCells == x)) > 0]) / np.mean(s) for s in segPairsCount]
                aRatio[iDay - 1][1].plot(tauPlot, [1 for _ in tauPlot], 'k--')
                aRatio[iDay - 1][1].plot(tauPlot,ratio, '.-', color=linCol)

                ratio = [np.mean([s[np.where(keptCells == x)][0] for x in impttCells if
                                  np.size(np.where(keptCells == x)) > 0]) / np.mean(s) for s in segPairsCount]
                aRatio[iDay - 1][2].plot(tauPlot, [1 for _ in tauPlot], 'k--')
                aRatio[iDay - 1][2].plot(tauPlot,ratio, '.-', color=linCol)

        else:
            ratioPcAni.append(np.nan)
            print('add nan')

    h, u = np.histogram(distCell[0], np.linspace(0, 0.3, 21), density=True)
    uxCell = u[:-1]+np.diff(u)[0]/2 # go from bin to x loc of bin center
    aTau[0][1].plot(uxCell,h)
    aTau[0][1].errorbar(uxCell, np.mean([np.histogram(d, np.linspace(0, 0.3, 21), density=True)[0] for d in distCell[1:]], 0),
                        np.std([np.histogram(d, np.linspace(0, 0.3, 21), density=True)[0] for d in distCell[1:]], 0))
    for d in randDistCell:
        aTau[0][1].hist(d, np.linspace(0, 0.3, 21), density=True, color=[0.7,0.7,0.7], histtype='step')

    aTau[0][1].errorbar(uxCell, np.mean([np.histogram(d, np.linspace(0, 0.3, 21), density=True)[0] for d in randDistCell], 0),
                        np.std([np.histogram(d, np.linspace(0, 0.3, 21), density=True)[0] for d in randDistCell], 0),
                        color=[0.5,0.5,0.5])
    distCellDiff.append(h-np.mean([np.histogram(d, np.linspace(0, 0.3, 21), density=True)[0] for d in distCell[1:]], 0))

    # for iSeg in range(len(percList)):
    #     aTau[0][2].hist(distCell[iSeg],np.linspace(0,0.3,21), histtype='step', density=True)
    aTau[0][0].hist(distTau[0], np.linspace(-1, 1, 61), histtype='step', density=True)
    aTau[0][0].hist([t_ for tauList in tAllFilteredSave for t in tauList for t_ in t], np.linspace(-1, 1, 61), histtype='step', density=True)
    aTau[1][0].imshow(perfSeg, vmin=0.5, vmax=1)
    aTau[1][1].imshow(percCellSeg, vmin=0.8, vmax=1)

    aTau[2][0].boxplot([[p_ for p_ in p if np.isfinite(p_)] for p in zip(*perfSeg)])
    aTau[2][1].boxplot([[p_ for p_ in p if np.isfinite(p_)] for p in zip(*percCellSeg)])

    distTauDiff.append(np.histogram(distTau[0], np.linspace(-1, 1, 61), density=True)[0]-
                       np.histogram([t_ for tauList in tAllFilteredSave for t in tauList for t_ in t],
                                    np.linspace(-1, 1, 61), density=True)[0])


    perfAll.append(perfAni)
    perfSegAll.append(perfSeg)
    percCellSegAll.append(percCellSeg)

    # reordering ratioPcAni for slicing per week
    ratioPcAniOrdered = []
    for iDayFind in range(1,10):
        idxDayRatio = np.where([d==iDayFind for d in daySaveRatio])[0]
        if len(idxDayRatio)>0:
            ratioPcAniOrdered.append(ratioPcAni[idxDayRatio[0]])
        else:
            ratioPcAniOrdered.append(np.nan)

    ratioPc.append(ratioPcAniOrdered)

    print()
    print(perfAni)
    print(ratioPcAni)
    print(ratioPcAniOrdered)

    if saveFig:
        fTau.savefig(pathFigures + strType + 'decodeTau_' + animal + '.eps', format='eps')
        fSvm.savefig(pathFigures + strType + 'decodeTauOrderedPairs_' + animal + '.eps', dpi=240, format='eps')
        fSvmNotOrdered.savefig(pathFigures + strType + 'decodeTauPairs_' + animal + '.eps', dpi=240, format='eps')

print(ratioPc)
fPcRatio, aPcRatio  = plt.subplots(1,1)

ratioPcBoxPlot = [[r for rAni in ratioPc for r in rAni[0 + iWk*3:3 + iWk*3] if np.isfinite(r)] for iWk in [0, 1, 2]]
tagRatioPc = [[[animalList[iAni], iWk, iDay] for iAni, rAni in enumerate(ratioPc) for iDay, r in enumerate(rAni[0 + iWk*3:3 + iWk*3]) if np.isfinite(r)] for iWk in [0, 1, 2]]
aPcRatio.boxplot(ratioPcBoxPlot, positions = [2, 9, 16], notch=False)

if saveFig: fPcRatio.savefig(pathFigure + 'SvmTau_ratioPc.eps', format='eps')

# export data for statistics
df = pd.DataFrame()
dataStr = ['ratioPC_PtiSvm']  # create index string
df = df.append(pd.DataFrame(data=np.concatenate(ratioPcBoxPlot), columns=dataStr).transpose())
dataStr = [s + '_tag_ratioCorr' + s for s in ['ani', 'week', 'day']]  # create index string
df = df.append(pd.DataFrame(data=np.concatenate(tagRatioPc), columns=dataStr).transpose())
if saveFig: df.to_csv(pathFigure + 'dataSvmTau_ratioPc.csv')

f,a = plt.subplots(2,3)

uxTau = np.linspace(-1, 1, 61)[:-1]+np.diff(np.linspace(-1, 1, 61))[0]/2

for iAni in range(nAni):
    a[0][0].plot(uxTau, distTauDiff[iAni])
    a[0][1].plot(uxCell, distCellDiff[iAni])
    a[0][2].plot(perfAll[iAni], '.-')


decodePerfData = [[p for pAni in perfAll for p in pAni[0 + iWk:3 + iWk] if np.isfinite(p)] for iWk in [0, 3, 6]]
decodePerfTag = [[[animalList[iAni], iWk/3, iDay] for iAni, pAni in enumerate(perfAll) for iDay, p in enumerate(pAni[0 + iWk:3 + iWk]) if np.isfinite(p)] for iWk in [0, 3, 6]]

dfPerf = pd.DataFrame()
dataStr = ['perf']  # create index string
dfPerf = dfPerf.append(pd.DataFrame(data=np.concatenate(decodePerfData), columns=dataStr).transpose())
dataStr = [s + '_tag_perf' + s for s in ['ani', 'week', 'day']]  # create index string
dfPerf = dfPerf.append(pd.DataFrame(data=np.concatenate(decodePerfTag), columns=dataStr).transpose())
dfPerf.to_csv(pathFigure + 'dataSvmDecodeTau_perf.csv')

a[0][2].errorbar([*range(9)], np.nanmean(perfAll,0), np.nanstd(perfAll,0), color='k')
a[1][0].boxplot(decodePerfData, positions = [2, 9, 16], notch=False)

a[1][1].boxplot([[pDay[iSeg] for pAni in perfSegAll for pDay in pAni if np.isfinite(pDay[iSeg])] for iSeg in range(10)],
             positions=[*range(10)], notch=False)

a[1][2].boxplot([[pDay[iSeg] for pAni in percCellSegAll for pDay in pAni if np.isfinite(pDay[iSeg])] for iSeg in range(10)],
             positions=[*range(10)], notch=False)

a[1][0].set_ylim([-.05,1.05])
a[1][1].set_ylim([-.05,1.05])
a[1][2].set_ylim([-.05,1.05])

if saveFig:
    f.savefig(pathFigure+'decodeTauSummary.eps', format='eps')

f,a = plt.subplots(1,2)
for iPlot in range(2):
    h, x, y, zPlot = a[iPlot].hist2d(tauSvmMeanAll[iPlot],tauSvmMeanAll[2],
                                     bins=[np.linspace(-1,1,101),np.linspace(-6,6,101)],
                                     normed=True, vmax=3)
# zPlot1 = a[0].hist2d(tauSvmMeanAll[0],tauSvmMeanAll[2],bins=[np.linspace(-1,1,101),np.linspace(-6,6,101)], density=True)
# zPlot2 = a[1].hist2d(tauSvmMeanAll[1],tauSvmMeanAll[2],bins=[np.linspace(-1,1,101),np.linspace(-6,6,101)], density=True)

    regr = linear_model.LinearRegression()
    regr.fit(np.reshape(tauSvmMeanAll[iPlot], [-1, 1]),tauSvmMeanAll[2])
    x1 = np.min(tauSvmMeanAll[iPlot])
    x2 = np.max(tauSvmMeanAll[iPlot])
    a[iPlot].plot([x1, x2], [x1* regr.coef_ + regr.intercept_, x2* regr.coef_ + regr.intercept_], 'b--')
    a[iPlot].set_ylim([-4,4])
    a[iPlot].set_xlim([-0.2,0.4])
    r,p = scstat.pearsonr(tauSvmMeanAll[iPlot], tauSvmMeanAll[2])
    a[iPlot].set_title('c = {:.3f} p = {:.3f}'.format(r,p))

    plt.colorbar(zPlot, ax=a[iPlot])

if saveFig:
    f.savefig(pathFigure + strType + 'TauWeightsHeatmap.eps', format='eps')

print(perfAll)

plt.show(block=True)