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
from sklearn import svm
import time

expt = 'CircleSquare'
frame = 'Arena'
saveOutput = True

BinNb = 12  # for place maps
TimeBinSize = 1 * 1000
TimeBinSize_AVC = 1 * 1000
TimeBinSize_TAU = 1000
min_cell_nb = 15
nPercReshape = 10
x_day = np.array([1, 2, 3, 8, 9, 10, 15, 16, 17])

doPerc = False

cvRep = 100
maxIter = 2000000
verbose = 0
cvParams = [1e-10, maxIter, verbose]
threshMult = 2
thresh = 0.01* threshMult*(TimeBinSize_TAU/1000)

pathAll = '/Users/eliott/Dropbox/dDocuments/NYU/Fenton/MNS_ANALYSIS/Data/' + expt + '/'
ExcelPath = pathAll + 'FilenameMatch.xlsx'

if doPerc:
    save_path = pathAll + 'analysisFiles/' + 'CS_Decoder_SVM_' + str(
        int(TimeBinSize)) + 'msBins_' + str(
        int(TimeBinSize_TAU)) + 'msKT_perc_thresh' + str(int(threshMult))+ '.file'
    tauBins = []
else:
    save_path = pathAll + 'analysisFiles/' + 'CS_Decoder_SVM_' + str(
        int(TimeBinSize)) + 'msBins_' + str(
        int(TimeBinSize_TAU)) + 'msKT_seg_thresh' + str(int(threshMult))+ '.file'

    # bins to cut tau Vec
    # nEqualSegs = 45
    # tauBins = np.linspace(-0.2, 0.7, nEqualSegs)
    # tauBins = [-0.1, -0.06, -0.04, -0.03, -0.025, -0.02, -0.0175, -0.015, -0.0125, -0.01, -0.005,
    #            0.005, 0.01, 0.02, 0.03, 0.04, 0.06, 0.1, 0.2]
    tauBins = [-0.05, -0.025, -0.015,-0.005, 0.005, 0.04, 0.1, 0.2]
    tauPlot = [tauBins[0] - (tauBins[1]-tauBins[0])/2]
    tauPlot.extend([b+d/2 for b,d in zip(tauBins,np.diff(tauBins))])
    tauPlot.extend([tauBins[-1] - (tauBins[-2]-tauBins[-1])/2])
    print(tauPlot)

if not saveOutput:
    print('not saving')

animalList = ['M19', 'M20', 'M29', 'M34', 'M35', 'M39']
nAni = len(animalList)

corrTauWeightAll = []
SegWeightRatioAll = []

zScoreRankedAni = []

aTauWeight = plt.subplots(nAni, 4)[1]
aTauWeight2d = plt.subplots(nAni,6)[1]
aTauSvmScat = plt.subplots(nAni,2)[1]
aZScoreRanked = plt.subplots(nAni, 2)[1]

for iAni, animal in enumerate(animalList):

    aBin = plt.subplots(1, 1)[1]
    fRankedZscore, aRankedZscore = plt.subplots(9,4)
    aRankedZscore2 = plt.subplots(9,4)[1]

    print(animal)

    # process excel file
    xlread = pd.ExcelFile(ExcelPath)
    xlAll = pd.read_excel(xlread, animal)
    xlAll = xlAll[pd.notna(xlAll.videoname)]

    Days = np.unique(list(xlAll.day))
    print(Days)

    # load PC analysis
    loadpath = pathAll + 'analysisFiles/' + 'MapCorr_' + str(int(BinNb)) + '_Bins_' + animal + '.file'
    with open(loadpath, 'rb') as f:
        output = pickle.load(f)
    dayRealListMap = output[1][0]
    sessListMap = output[1][2]
    PcIdxList = output[1][4]

    #load temporal analysis
    with open(pathAll + 'analysisFiles/TempCorr_List_'+
               str(int(TimeBinSize_AVC)) + 'msAV_' + str(int(TimeBinSize_TAU)) + 'msKT_'
              + expt + '_' + animal + '.file', 'rb') as f:
        outputTemp = pickle.load(f)
    dayRealListTemp = outputTemp[0][0]
    sessListTemp = outputTemp[0][2]
    sConvTauList = outputTemp[0][4]
    isSeparatedList = outputTemp[0][5]
    tauPairsSingleList = outputTemp[0][6]
    tauVecSingleList = outputTemp[0][7]

    # pre-allocate
    corrTauWeight = [[], []]
    SegWeightRatio = []
    tauWeight_toPlot = [[[], []] for _ in range(6)]
    zScoreRankedDay = []

    # create figures
    fTauDay, aTauDay = plt.subplots(3, 3)

    for i_day, day in enumerate(Days):

        # load data
        xlOut, fileNames, physio = loadPhysiology(day, animal, xlAll, pathAll, expt, frame)
        VideoName, DayCount, Sessions, doAnalysis, Training, _ = xlOut
        MnsTSFileName = fileNames[0]
        sOutput = physio[0]
        rec_idx = physio[2]
        physio = None

        sess_indic = []
        s_all = []

        # extract data for the day
        for i_vid, vid in enumerate(VideoName):

            if doAnalysis[i_vid] > 0:

                # extract physiology from longer file
                s_tmp = np.transpose(sOutput[rec_idx[i_vid][0]:rec_idx[i_vid][1]] / np.std(sOutput, 0))
                # s_tmp = np.transpose((S_all[rec_idx[i_vid][0]:rec_idx[i_vid][1]]-np.mean(S_all, 0)) / np.std(S_all, 0))
                # s_tmp = np.transpose(S_all[rec_idx[i_vid][0]:rec_idx[i_vid][1]])

                # load video time file and adjust for delay to start MNS
                MnsFrameCount, MnsTS = read_timestamps(MnsTSFileName[i_vid])

                # correction for lost frames during tiff conversion
                if Sessions[i_vid][:3] == 'WIR' or Sessions[i_vid][:3] == 'PTC' or len(s_tmp[0])>3333:
                    MnsTS_corrected_for_tiff = MnsTS[1:1000] + MnsTS[10001:10001 + 3 * ((len(MnsTS) - 10001) // 3)]
                else:
                    MnsTS_corrected_for_tiff = MnsTS[1:1 + 3 * ((len(MnsTS) - 1) // 3)]

                # correction for delay between scope and camera
                MnsTS_delay_adjusted = [ts + 100 for ts in MnsTS_corrected_for_tiff]

                # correction for re sampling
                MnsTS_resampld = np.min(np.reshape(MnsTS_delay_adjusted, [-1, 3]), 1)  # used for binning later

                # binning for calcium trace (w/o tacking tracker into account, see older code for that)
                # AVC
                maxTimeBin = np.ceil(np.max(MnsTS_resampld)) + 1 + TimeBinSize
                minTimeBin = np.floor(np.min(MnsTS_resampld)) - 1
                TimeBins = np.arange(minTimeBin, maxTimeBin, TimeBinSize)
                MnsMatchTimeBins = np.digitize(MnsTS_resampld, TimeBins)
                sBinned = [
                    [np.mean([s[ii] for ii in np.where(MnsMatchTimeBins == i + 1)[0]]) for i in range(len(TimeBins) - 1)]
                    for s in s_tmp]

                if 'HMC' in Sessions[i_vid][:3]:
                    print('hmc')
                    # sess_indic.append(np.zeros_like(S_binned_avc[0]))
                elif 'CYL' in Sessions[i_vid][:3]:
                    sess_indic.append(np.ones_like(sBinned[0]))
                    s_all.append(sBinned)
                elif 'RCT' in Sessions[i_vid][:3]:
                    sess_indic.append(np.zeros_like(sBinned[0]))
                    s_all.append(sBinned)

        if len(sess_indic)>2:

            #concatenate data of the day
            sess_concat_tmp = np.concatenate(sess_indic, 0)
            S_concat_tmp = np.transpose(np.concatenate(s_all, 1))

            S_concat = S_concat_tmp[np.any(~np.isnan(S_concat_tmp), 1)]
            sess_concat = sess_concat_tmp[np.any(~np.isnan(S_concat_tmp), 1)]

            # self training recording with CV
            pCv, cCv, _ = SVM_decode(S_concat, sess_concat, do_cv=True, cv_rep=cvRep, cv_param=[1e-10, 100000, 0])

            #EVALUATE TAU vs Weights
            #find day in temporal recording
            idxDayTemp = np.where([(d == day) and ('HMC' not in sess) for d, sess in zip(dayRealListTemp, sessListTemp)])[0]
            idxDayMap = np.where([(d == day) and ('HMC' not in sess) for d, sess in zip(dayRealListMap, sessListMap)])[0]
            SegWeightRatioDay = []
            corrTauWeightDay = [[], []]

            if i_day>6:
                f, a = plt.subplots(1, 2)
                isPcList = [np.any(x) for x in zip(*[PcIdxList[iDayMap] for iDayMap in idxDayMap])]
                for i, idx in enumerate(np.argsort(cCv)):
                    if isPcList[idx]:
                        col = [0, 0, 0]
                    else:
                        col = [.7, .7, .7]
                    a[0].plot([s[idx] + i for s in S_concat], color=col)
                    a[1].plot([s[idx] + i for s in S_concat], color=col)

            zScoreRankedRec = [[], []]
            for iRec, (iDayTemp,iDayMap) in enumerate(zip(idxDayTemp, idxDayMap)):

                # f_, a_ = plt.subplots(1, 4)
                # for i, idx in enumerate(np.argsort(cCv)):
                #     a_[iRec].plot([s+i/10 for s in sConvTauList[iDayTemp][idx]])

                # get pair numbers and tau corr
                tauVecTmp = tauVecSingleList[iDayTemp]
                idxTauPairsTmp = tauPairsSingleList[iDayTemp]

                # to remove nearby cells
                isSeparated = isSeparatedList[iDayTemp]

                # #create filter for Nan's
                isAnyNan = [np.isnan(t[0]) for t in tauVecTmp]

                # filter out cells that are rarely active
                meanRate = [np.mean([s > 0 for s in c]) for c in sConvTauList[iDayTemp]]
                rate = [meanRate[p[0]] > thresh and meanRate[p[1]] > thresh for p in idxTauPairsTmp]

                # create filter for spatial separation and significance and non nan correlation
                # filter = [True if (isSep and t[1]<0.05 and not np.isnan(t[0])) else False for isSep, t in zip(isSeparated, tauVecTmp)]
                tauFilter = [sep and r and not iNan for sep, iNan, r in zip(isSeparated, isAnyNan, rate)]
                meanRate = rate = isAnyNan = isSeparated = []

                # filter tau values and pair numbers
                tauVec = [t[0] for t,f in zip(tauVecTmp, tauFilter) if f]
                idxTauPairs = [idx for idx,f in zip(idxTauPairsTmp, tauFilter) if f]

                keptCells = np.unique([x for pair in idxTauPairs for x in pair])
                svmCoef = [cCv[iCell] for iCell in keptCells]
                nKeptCells = len(keptCells)
                # old method for above
                # tau_vec_tmp = [scstat.kendalltau(c1, c2)[0] if isSep else np.nan
                #                for isSep, (c1, c2) in zip(isSeparated, itertools.combinations(sTau, 2))]
                #
                # tau_vec = [t for t in tau_vec_tmp if not np.isnan(t)]
                #
                # idxTauPairs = [(x1, x2) for [x1, x2], t in
                #                zip(itertools.combinations(range(len(sTau)), 2), tau_vec_tmp)
                #                if not np.isnan(t)]

                if doPerc:
                    # to cut pairs in segments
                    idxSort = np.argsort(tauVec)[::-1]
                    nSegPairs = 20
                    nPairsPerSeg = len(idxSort) // nSegPairs

                    #cut ordered tau pairs in nSegPairs segments
                    sortedSegPairList = [[idxTauPairs[n] for n in idxSort[i * nPairsPerSeg:(i + 1) * nPairsPerSeg]]
                                      for i in range(nSegPairs)]

                    # since the nb of pairs is not always a multiple of nSegPairs, add last ones to last segment
                    sortedSegPairList[nSegPairs-1].extend([idxTauPairs[n] for n in idxSort[nSegPairs * nPairsPerSeg:]])

                    # evaluate how many times each cell is in each segment
                    # segPairsCount = [[np.sum([np.any([p[0] == i, p[1] == i]) for p in pairsSeg]) for i in range(nCells)]
                    #                  for pairsSeg in sortedSegPairList]
                    segPairsCountTmp = [[np.sum([np.any([p[0] == i, p[1] == i]) for p in pairsSeg]) for i in keptCells]
                                     for pairsSeg in sortedSegPairList]

                    #normalize for the number of pairs a cell is in
                    segPairsCount = [x/np.sum(segPairsCountTmp,0) for x in segPairsCountTmp]

                else:
                    nBins = len(tauBins)+1
                    binList = [*range(nBins)]

                    # sort pairs in different segments
                    equalSegIdx = np.digitize(tauVec, tauBins) # a vector containing the bin number for each pair

                    # for each segment, list the pairs it contains
                    equalSegPairList = [[idxTauPairs[iPair] for iPair in np.where(equalSegIdx == iBin)[0]]
                                        for iBin in binList]

                    aBin.plot(tauPlot, [len(x) for x in equalSegPairList], '.')

                    # evaluate how many times each cell is in each segment
                    segPairsCountTmp = [[np.sum([np.any([p[0] == i, p[1] == i]) for p in pairsSeg]) for i in keptCells]
                                     for pairsSeg in equalSegPairList]

                    #normalize for the number of pairs a cell is in
                    segPairsCount = [x/np.sum(segPairsCountTmp,0) for x in segPairsCountTmp]

                    if 'RCT' in sessListTemp[iDayTemp]: col = 'r'
                    elif 'CYL' in sessListTemp[iDayTemp]: col='b'

                    zScore = [(s - np.mean(s)) / np.std(s) for s in segPairsCount] # zscore for each tau, bin across all cells the particpation in that tau bin
                    zScoreRanked = [np.transpose(zScore)[idx] for idx in np.argsort([cCv[k] for k in keptCells])]
                    aRankedZscore[i_day][iRec].imshow(np.transpose(zScoreRanked), aspect='auto', vmin=-1, vmax=1)
                    aRankedZscore[i_day][iRec].set_title(sessListTemp[iDayTemp])

                    SegCoefIdx = np.digitize([c / np.std(cCv) for c in cCv], [-1,-0.5, 0, 0.5, 1]) #find bins for cCv (all cells)
                    svmSeg = [[] for _ in range(12)] # pre-allocate (each [] for a bin created line above)
                    for iCell, iSeg in enumerate([SegCoefIdx[k] for k in keptCells]):
                        svmSeg[iSeg].append(np.transpose(zScore)[iCell])
                    if 'RCT' in sessListTemp[iDayTemp]:
                        zScoreRankedRec[0].append([np.mean(s, 0) if np.size(s) > 1 else [np.nan for _ in binList] for s in svmSeg])
                    elif 'CYL' in sessListTemp[iDayTemp]:
                        zScoreRankedRec[1].append([np.mean(s, 0) if np.size(s) > 1 else [np.nan for _ in binList] for s in svmSeg])
                    aRankedZscore2[i_day][iRec].imshow([np.mean(s, 0) if np.size(s) > 1 else [np.nan for _ in binList] for s in svmSeg],
                                                       aspect='auto', vmin=-0.5, vmax=0.5)
                    aRankedZscore2[i_day][iRec].set_title(sessListTemp[iDayTemp])


                    for ii, iTauSeg in enumerate(segPairsCount):
                        zscore = (iTauSeg - np.mean(iTauSeg)) / np.std(iTauSeg)
                        for i in range(np.min([10,len(zscore)])):
                            idxNeg = np.argsort(zscore)[::-1][i]
                            if zscore[idxNeg] > 0:
                                if 'RCT' in sessListTemp[iDayTemp]:
                                    aTauSvmScat[iAni][0].scatter(cCv[keptCells[idxNeg]],tauPlot[ii], s=1, color=col, alpha=0.5)
                                elif 'CYL' in sessListTemp[iDayTemp]:
                                    aTauSvmScat[iAni][1].scatter(cCv[keptCells[idxNeg]], tauPlot[ii], s=1, color=col,
                                                                 alpha=0.5)

                    if i_day > 6:
                        zscoreNeg = (segPairsCount[1] - np.mean(segPairsCount[1])) / np.std(segPairsCount[1])
                        zscorePos = (segPairsCount[-1] - np.mean(segPairsCount[-1])) / np.std(segPairsCount[-1])
                        for i in range(10):
                            idxPos = np.argsort(zscorePos)[::-1][i]
                            idxNeg = np.argsort(zscoreNeg)[::-1][i]
                            # a[1].plot([s[keptCells[idx]] + i for s in S_concat])
                            a[0].plot([s[keptCells[idxNeg]] + np.where([x == keptCells[idxNeg] for x in np.argsort(cCv)])[0][0]
                                       for s in S_concat], color=col, alpha=0.8)
                            a[1].plot([s[keptCells[idxPos]] + np.where([x == keptCells[idxPos] for x in np.argsort(cCv)])[0][0]
                                       for s in S_concat], color=col, alpha=0.8)



                    aTauWeight[iAni][0].scatter(svmCoef, segPairsCount[0], s=1, alpha=0.3, color=col)
                    aTauWeight[iAni][1].scatter(svmCoef, segPairsCount[1], s=1, alpha=0.3, color=col)
                    aTauWeight[iAni][2].scatter(svmCoef, segPairsCount[-2], s=1, alpha=0.3, color=col)
                    aTauWeight[iAni][3].scatter(svmCoef, segPairsCount[-1], s=1, alpha=0.3, color=col)

                if 'RCT' in sessListTemp[iDayTemp]:
                    for iIdx, idx in enumerate([0, 1, 2, nBins-3, nBins-2, nBins-1]):
                        tauWeight_toPlot[iIdx][0].extend(svmCoef)
                        tauWeight_toPlot[iIdx][1].extend(segPairsCount[idx])

                    corrTauWeightDay[0].append([np.corrcoef(pairCount, svmCoef)[1, 0] for pairCount in segPairsCount])
                    aTauDay.flat[i_day].plot(segPairsCount[0], svmCoef, 'r.')
                    segWeightCount = [[np.mean([pairCount[x]
                                                for i_x, x in enumerate(np.argsort(svmCoef)[::-1])
                                                if i_x > nKeptCells - (nKeptCells // 5) * i_perc])
                                       for i_perc in [1, 5]]
                                      for pairCount in segPairsCount]

                else:
                    aTauWeight[iAni][0].scatter(svmCoef, segPairsCount[0], s=1, alpha=0.3, color='k')
                    aTauWeight[iAni][1].scatter(svmCoef, segPairsCount[1], s=1, alpha=0.3, color='k')
                    aTauWeight[iAni][2].scatter(svmCoef, segPairsCount[-2], s=1, alpha=0.3, color='k')
                    aTauWeight[iAni][3].scatter(svmCoef, segPairsCount[-1], s=1, alpha=0.3, color='k')

                    corrTauWeightDay[1].append([np.corrcoef(pairCount, svmCoef)[1, 0] for pairCount in segPairsCount])
                    aTauDay.flat[i_day].plot(segPairsCount[0], svmCoef, 'b.')
                    segWeightCount = [[np.mean([pairCount[x]
                                                for i_x, x in enumerate(np.argsort(svmCoef))
                                                if i_x > nKeptCells - (nKeptCells // 5) * i_perc])
                                       for i_perc in [1, 5]]
                                      for pairCount in segPairsCount]
                SegWeightRatioDay.append([s[0]/s[1] for s in segWeightCount])

            zScoreRankedDay.append([np.nanmean(z, 0) for z in zScoreRankedRec])
            SegWeightRatio.append(np.nanmean(SegWeightRatioDay, 0))
            corrTauWeight[0].append(np.nanmean(corrTauWeightDay[0], 0))
            corrTauWeight[1].append(np.nanmean(corrTauWeightDay[1], 0))

        else:
            SegWeightRatio.append(np.nan)
            for i in range(2):
                corrTauWeight[i].append(np.nan)

    print(corrTauWeight)
    print(SegWeightRatio)

    aZScoreRanked[iAni][0].imshow(np.transpose(np.nanmean(zScoreRankedDay, 0)[0]))
    aZScoreRanked[iAni][1].imshow(np.transpose(np.nanmean(zScoreRankedDay, 0)[1]))

    zScoreRankedAni.append([np.transpose(np.nanmean(zScoreRankedDay, 0)[0]),
                            np.transpose(np.nanmean(zScoreRankedDay, 0)[1])])

    for iIdx, val in enumerate(tauWeight_toPlot):
        H, xedges, yedges = np.histogram2d(val[0], val[1], bins=15)
        # aTauWeight2d[iAni][iIdx].hist2d(tauWeight_toPlot[iIdx][0], tauWeight_toPlot[iIdx][1], 20)
        aTauWeight2d[iAni][iIdx].pcolormesh(xedges, yedges, np.transpose([h / np.sum(h) for h in H]))
    # fTau, aTau = plt.subplots(1,2)
    # aTau[0].plot(SegWeightRatio, '.')
    # aTau[1].plot(corrTauWeight[0], 'r.')
    # aTau[1].plot(corrTauWeight[1], 'b.')

    corrTauWeightAll.append(corrTauWeight)
    SegWeightRatioAll.append(SegWeightRatio)

Title_perf = ['current decode current _best',
              'current decode current _cv']

x_day = np.array([1, 2, 3, 8, 9, 10, 15, 16, 17])

f,a = plt.subplots(1,2)
a[0].imshow(np.mean([z[0] for z in zScoreRankedAni], 0))
a[1].imshow(np.mean([z[1] for z in zScoreRankedAni], 0))

if saveOutput:
    print('saving...')
    with open(save_path, 'wb') as f:
        pickle.dump([animalList, corrTauWeightAll, SegWeightRatioAll, tauBins], f, pickle.HIGHEST_PROTOCOL)
    print('saved')
else:
    print('not saving')

print('done done done')

plt.show(block=True)