# import matplotlib as mpl
# # mpl.use('TkAgg')
# from mpl_toolkits.mplot3d import Axes3D

import matplotlib.pyplot as plt
import numpy as np
import pickle
from scipy.ndimage import gaussian_filter1d
from sklearn import decomposition, manifold
import pandas as pd

path_all = '/Volumes/GoogleDrive/My Drive/Fenton_Lab/CellPress_CircleSquare/MNS_Analysis/Data/CircleSquare/'
pathFigure = '/Volumes/GoogleDrive/My Drive/Fenton_Lab/CellPress_CircleSquare/MNS_Analysis/Figures/CircleSquare' \
             '/savedFigs/IsoMap_Simon/ '

TimeBinSize_AVC = 1000
TimeBinSize_TAU = 1000

threshMult = 2
thresh = 0.01 * threshMult * (TimeBinSize_TAU / 1000)

saveFig = False

sigma = 1  # Kernel width

n_neighbors = 5
target_dim = 3
doStd = True
doPTI = False

saveAni = ['M29']

animalList = ['M29']
# ['M20', 'M29','M34','M19', 'M39']
nAni = len(animalList)

indic = [-1, 0, 1]
envList = ['HMC1', 'HMC2', 'CYL1', 'CYL2', 'RCT1',
           'RCT2']  # environment HMC = Home Cage (-1), CYL = Cylinder (1), RCT = Rectangle (0)
indicColors = [[0, 0, 0], [0.5, 0.5, 0.5, 0.5], [0, 0, 1], [0, 0, 0.5], [1, 0, 0], [0.5, 0, 0]]

tauBins = [-0.05, -0.04, -0.03, -0.025, -0.02, -0.0175, -0.015, -0.0125, -0.01, -0.005,
           0.005, 0.01, 0.02, 0.03, 0.04, 0.06,
           0.1]  # in this code only looking at most negative and most positive, i.e, <-0.05 and >+0.1
tauPlot = [tauBins[0] - (tauBins[1] - tauBins[0]) / 2]
tauPlot.extend([b + d / 2 for b, d in zip(tauBins, np.diff(tauBins))])
tauPlot.extend([tauBins[-1] - (tauBins[-2] - tauBins[-1]) / 2])

meanIsoAll = []
meanPcaAll = []
meanAvAll = []

ratioDiscardList = [20, 10]  # [20,10,4,2] # discarding 1/ratioDiscardList %
discardPop = 'neg'  # kc: Keep Cells, pos: Positive Cells Pairs, neg: Negative Cell Pairs
discard = True  # if True: discard the 1/ratioDiscard % most (anti-)correlated if False: only keep the 1/ratioDiscard
# % most (anti-)correlated

if doPTI:
    loadStr = 'Zlinear'
else:
    loadStr = 'TempCorr'

target_dim = 3  # dimension to which we want to reduce dimensionality
fAvProj, aAvProj = plt.subplots(1, 2)

for iAni, animal in enumerate(animalList):

    print(animal)

    with open(path_all + 'analysisFiles_Simon/' + loadStr + '_List_' +
              str(int(TimeBinSize_AVC)) + 'msAV_' + str(int(TimeBinSize_TAU)) + 'msKT'
                                                                                '_CircleSquare_' + animal + '.file',
              'rb') as f:
        outputList = pickle.load(f)

    # from cnmfe_temporal_analysis
    # listsToSave = [[day_real_list, day_count_list, sess_list, S_conv_avc_list, S_conv_tau_list, isSeparatedList,
    #                 tauPairsSingleList, tauVecSingleList, trainingList],
    #                [trackList, [], aMaxLocList],
    #                [splitList, iDayList]]

    # from cnmfe_zscore_analysis listsToSave = [[day_real_list, day_count_list, sess_list, S_conv_avc_list,
    # [zConvTauList, S_conv_tau_list, sConvTauPcList], isSeparatedList, tauPairsSingleList, [zTauVecSingleList,
    # tauVecSingleList], trainingList, xBinnedPcList, yBinnedPcList], [trackList, [zBinnedRateMapsList,
    # zRateMapsList], aMaxLocList], [atnSplitList, iDayList]]
    if doPTI:
        dayCountList = outputList[0][1]
        sessList = outputList[0][2]
        sTauList = outputList[0][4][0]
        isSeparatedList = outputList[0][5]
        tauPairsSingleList = outputList[0][6]
        tauVecSingleList = outputList[0][7][0]
        trainingList = outputList[0][8]
    else:
        dayCountList = outputList[0][1]
        sessList = outputList[0][2]
        sTauList = outputList[0][4]
        isSeparatedList = outputList[0][5]
        tauPairsSingleList = outputList[0][6]
        tauVecSingleList = outputList[0][7]
        trainingList = outputList[0][8]

    outputList = []

    figSvm = plt.figure(figsize=(8, 8))
    if target_dim < 4:
        figScatter = plt.figure(figsize=(8, 8))

    meanIsoAni = []
    meanPcaAni = []
    meanAvAni = []

    for iDay, day in enumerate(range(1, 10)):
        print([iDay, day])
        sessIdx = np.where([d == day for d in dayCountList])[0]

        if len(sessIdx) < 5:
            meanIsoAni.append([])
            meanPcaAni.append([])
            meanAvAni.append([])
            continue

        sConcat = []
        sessIndic = []
        sessIndicFull = []
        meanIsoDay = [[] for _ in ratioDiscardList]
        meanPcaDay = [[] for _ in ratioDiscardList]
        meanAvDay = [[] for _ in ratioDiscardList]

        idxPosList = [[] for _ in ratioDiscardList]
        idxNegList = [[] for _ in ratioDiscardList]
        keptCellsList = []

        for iSessIdx in sessIdx:

            if 'HMC' in sessList[iSessIdx][:3]:
                sessIndic.append((-1) * np.ones_like(sTauList[iSessIdx][0]))
                sessIndicFull.append([sessList[iSessIdx] for s in sTauList[iSessIdx][0]])
                sConcat.append(gaussian_filter1d(sTauList[iSessIdx], sigma, 1))
                # posDay.append([Xlist[iSessIdx], Ylist[iSessIdx]])

            elif 'CYL' in sessList[iSessIdx][:3]:
                sessIndic.append(np.ones_like(sTauList[iSessIdx][0]))
                sessIndicFull.append([sessList[iSessIdx] for s in sTauList[iSessIdx][0]])
                sConcat.append(gaussian_filter1d(sTauList[iSessIdx], sigma, 1))
                # posDay.append([Xlist[iSessIdx], Ylist[iSessIdx]])

            elif 'RCT' in sessList[iSessIdx][:3]:
                sessIndic.append(np.zeros_like(sTauList[iSessIdx][0]))
                sessIndicFull.append([sessList[iSessIdx] for s in sTauList[iSessIdx][0]])
                sConcat.append(gaussian_filter1d(sTauList[iSessIdx], sigma, 1))
                # posDay.append([Xlist[iSessIdx], Ylist[iSessIdx]])

            # FIND PAIRS OVER-REPRESENTED IN CERTAIN TAU CORR
            # get pair numbers and tau corr
            tauVecTmp = tauVecSingleList[iSessIdx]
            idxTauPairsTmp = tauPairsSingleList[iSessIdx]

            isSeparated = isSeparatedList[iSessIdx]  # to remove nearby cells
            isAnyNan = [np.isnan(t[0]) for t in tauVecTmp]  # create filter for Nan's

            # filter out cells that are rarely active
            meanRate = [np.mean([s > 0 for s in c]) for c in sTauList[iSessIdx]]
            rate = [meanRate[p[0]] > thresh and meanRate[p[1]] > thresh for p in idxTauPairsTmp]

            # create filter for spatial separation and significance and non nan correlation
            tauFilter = [sep and r and not iNan for sep, iNan, r in zip(isSeparated, isAnyNan, rate)]
            meanRate = rate = isAnyNan = isSeparated = []  # delete variable from memory

            # filter tau values and pair numbers
            tauVec = [t[0] for t, f in zip(tauVecTmp, tauFilter) if f]
            idxTauPairs = [idx for idx, f in zip(idxTauPairsTmp, tauFilter) if f]

            keptCells = np.unique([x for pair in idxTauPairs for x in pair])
            nKeptCells = len(keptCells)

            nBins = len(tauBins) + 1
            binList = [*range(nBins)]
            # sort pairs in different segments
            equalSegIdx = np.digitize(tauVec, tauBins)  # a vector containing the bin number for each pair

            # for each segment, list the cell pairs it contains
            equalSegPairList = [[idxTauPairs[iPair] for iPair in np.where(equalSegIdx == iBin)[0]]
                                for iBin in binList]

            # evaluate how many times each cell is in each segment
            segPairsCountTmp = [[np.sum([np.any([p[0] == i, p[1] == i]) for p in pairsSeg]) for i in keptCells]
                                for pairsSeg in equalSegPairList]

            # normalize for the number of pairs a cell is in
            segPairsCount = [x / np.sum(segPairsCountTmp, 0) for x in segPairsCountTmp]

            keptCellsList.append(keptCells)

            for iR, ratio in enumerate(ratioDiscardList):
                nDiscard = len(keptCells) // ratio  # number of cells discarded

                if 'HMC' not in sessList[iSessIdx][:3]:
                    idxPosList[iR].append([keptCells[idx] for idx in np.argsort(segPairsCount[-1])[::-1][:nDiscard]
                                           if segPairsCount[-1][idx] > 0])
                    idxNegList[iR].append([keptCells[idx] for idx in np.argsort(segPairsCount[0])[::-1][:nDiscard]
                                           if segPairsCount[0][idx] > 0])

        # sConcatFilt = [[s for iCell, s in enumerate(sRec) if iCell not in np.unique(idxNegList)]
        #                for sRec in sConcat]

        keptCellsAll = np.unique([x_ for x in keptCellsList for x_ in x])

        for iR, ratio in enumerate(ratioDiscardList):
            if discardPop == 'pos':
                discardCells = np.unique([x_ for x in idxPosList[iR] for x_ in x])
                strCode = 'percPos_noHMC'

            elif discardPop == 'neg':
                discardCells = np.unique([x_ for x in idxNegList[iR] for x_ in x])
                strCode = 'percNeg_noHMC'

            elif discardPop == 'kc':
                discardCells = []
                strCode = ''

            elif discardPop == 'randPos':
                randLen = len(np.unique([x_ for x in idxNegList[iR] for x_ in x]))
                discardCells = np.random.choice(keptCellsAll, randLen, replace=False)
                strCode = ''

            if discard:
                sConcatFilt = [[s for iCell, s in enumerate(sRec)
                                if (iCell in keptCellsAll) and (iCell not in discardCells)]
                               for sRec in sConcat]

                print(len([iCell for iCell, s in enumerate(sConcat[0])
                           if (iCell in keptCellsAll) and (iCell not in discardCells)]))
            else:
                sConcatFilt = [[s for iCell, s in enumerate(sRec)
                                if (iCell in keptCellsAll) and (iCell in discardCells)]
                               for sRec in sConcat]

                print(len([iCell for iCell, s in enumerate(sConcat[0])
                           if (iCell in keptCellsAll) and (iCell in discardCells)]))

            if len([iCell for iCell, s in enumerate(sConcat[0]) if
                    (iCell in keptCellsAll) and (iCell in discardCells)]) < 5:
                continue

            if doStd:
                sAll = np.transpose(np.concatenate(sConcatFilt, 1)) / np.std(
                    np.transpose(np.concatenate(sConcatFilt, 1)), 0)
            else:
                sAll = np.transpose(np.concatenate(sConcatFilt, 1))

            sessIndicAll = np.concatenate(sessIndic, 0)
            sessIndicFullAll = np.concatenate(sessIndicFull, 0)

            if doPTI:
                # to handle negative numbers
                data_to_use = np.sign(sAll.copy()) * np.sqrt(sAll.copy() * np.sign(sAll.copy()))
                # data_to_use = sAll.copy()
            else:
                data_to_use = np.sqrt(sAll.copy())

            # do Isomap transformation
            iso_instance = manifold.Isomap(n_neighbors=n_neighbors, n_components=target_dim, path_method='D')
            dataIso = iso_instance.fit_transform(data_to_use)

            # do PCA transformation
            pca1 = decomposition.PCA(target_dim)
            dataPca = pca1.fit_transform(data_to_use)

            ax = figSvm.add_subplot(3, 3, iDay + 1)
            for iIndic, envIndic in enumerate(envList):
                projDataIso = dataIso[[envIndic in s for s in sessIndicFullAll], :]
                meanIsoDay[iR].append(np.mean(projDataIso, 0))

                projDataPca = dataPca[[envIndic in s for s in sessIndicFullAll], :]
                meanPcaDay[iR].append(np.mean(projDataPca, 0))

                projDataAV = data_to_use[[envIndic in s for s in sessIndicFullAll], :]
                meanAvDay[iR].append(np.mean(projDataAV, 0))

        meanIsoAni.append(meanIsoDay)
        meanPcaAni.append(meanPcaDay)
        meanAvAni.append(meanAvDay)

    meanIsoAll.append(meanIsoAni)
    meanPcaAll.append(meanPcaAni)
    meanAvAll.append(meanAvAni)

    if animal in saveAni and target_dim == 2 and saveFig and False:
        figScatter.savefig(pathFigure + 'IsomapProj_' + str(target_dim) + 'D_' + animal + '.eps', format='eps')

    print(strCode)

pairList = [[(2, 3), (4, 5)], [(2, 4), (2, 5), (3, 4), (3, 5)]]
colorList = ['k', 'g', 'm']

# exporting data
df = pd.DataFrame()
metricStr = ['isoMap', 'pca', 'full']

f, a = plt.subplots(2, len(ratioDiscardList))

for iR, ratio in enumerate(ratioDiscardList):

    if discard:
        strTitle = 'discardCells-' + str(int(100 / ratio)) + strCode
    else:
        strTitle = 'keptCells-' + str(int(100 / ratio)) + strCode

    for iMetric, metric in enumerate([meanIsoAll, meanPcaAll, meanAvAll]):
        z = [[[np.mean([np.sqrt(np.sum((m[iR][x1] - m[iR][x2]) ** 2)) for x1, x2 in pairs]) for pairs in pairList]
              if (len(m) > 0 and len(m[iR]) > 0) else [np.nan, np.nan] for m in mAni] for mAni in metric]
        zRatio = [[zDay[0] / zDay[1] for zDay in zAni] for zAni in z]
        print(z)

        for zAni in zRatio:
            a[0][iR].plot([*range(9)], zAni, '.', color=colorList[iMetric])
        a[0][iR].errorbar([*range(9)], np.nanmean(zRatio, 0), np.nanstd(zRatio, 0) / np.sqrt(nAni),
                          color=colorList[iMetric])

        dataBoxPlot = [[z for zAni in zRatio for z in zAni[0 + iWk:3 + iWk] if np.isfinite(z)] for iWk in [0, 3, 6]]
        tagData = [[[animalList[iAni], iWk / 3, iDay] for iAni, zAni in enumerate(zRatio) for iDay, z in
                    enumerate(zAni[0 + iWk:3 + iWk]) if np.isfinite(z)] for iWk in [0, 3, 6]]
        a[1][iR].boxplot(dataBoxPlot, positions=[0 + iMetric, 6 + iMetric, 12 + iMetric])

        if iMetric == 0:
            # export data
            if discard:
                typeStr = '_discard_' + str(int(100 / ratio)) + discardPop
            else:
                typeStr = '_keep_' + str(int(100 / ratio)) + discardPop
            dataStr = [metricStr[iMetric] + typeStr]  # create index string
            df = df.append(pd.DataFrame(data=np.concatenate(dataBoxPlot), columns=dataStr).transpose())
            # dataStr = [s + '_' + metricStr[iMetric] + typeStr for s in ['ani', 'week', 'day']]  # create index string
            # df = df.append(pd.DataFrame(data=np.concatenate(tagData), columns=dataStr).transpose())

    a[1][iR].set_xlim([-1, 15])
    a[1][iR].set_ylim([-0.2, 2])
    a[0][iR].set_title(strTitle)
    print(strCode)

if saveFig:
    if discard:
        f.savefig(pathFigure + 'DimRed_discardPercPop_' + discardPop + '.eps', format='eps')
        df.to_csv(pathFigure + 'dataDimRed_discard_' + discardPop + '.csv')
    else:
        f.savefig(pathFigure + 'DimRed_keepPercPop_' + discardPop + '.eps', format='eps')
        df.to_csv(pathFigure + 'dataDimRed_keep_' + discardPop + '.csv')

    plt.savefig(pathFigure + 'QuantDimRed_dimRed_' + str(target_dim) + 'D' + strCode + '.eps', format='eps')
    if False:  # target_dim == 2:
        fAvProj.savefig(pathFigure + 'QuantDimRed_CS_2DAverageProj_meanEnv_wk3.eps', format='eps')

plt.show(block=True)
