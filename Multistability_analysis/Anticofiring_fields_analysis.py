import numpy as np
import pickle
from scipy.ndimage import gaussian_filter1d
from scipy.io import savemat

path_all = '/Volumes/GoogleDrive/My Drive/Fenton_Lab/CellPress_CircleSquare/MNS_Analysis/Data/CircleSquare/'

TimeBinSize_AVC = 1000
TimeBinSize_TAU = 1000

threshMult = 2
thresh = 0.01 * threshMult * (TimeBinSize_TAU / 1000)

sigma = 1  # Kernel width

doStd = True
doPTI = False
expt = 'CircleSquare'
saveOutput = False

saveAni = []
animalList = ['M29','M34','M35','M39']
# ['M20', 'M29','M34','M19', 'M39']
nAni = len(animalList)

indic = [-1, 0, 1]
envList = ['HMC1', 'HMC2', 'CYL1', 'CYL2', 'RCT1',
           'RCT2']  # environment HMC = Home Cage (-1), CYL = Cylinder (1), RCT = Rectangle (0)
indicColors = [[0, 0, 0], [0.5, 0.5, 0.5, 0.5], [0, 0, 1], [0, 0, 0.5], [1, 0, 0], [0.5, 0, 0]]

tauBins = [-0.05, 0.1]  # in this code only looking at most negative and most positive, i.e, <-0.05 and >+0.1
tauPlot = [tauBins[0] - (tauBins[1] - tauBins[0]) / 2]
tauPlot.extend([b + d / 2 for b, d in zip(tauBins, np.diff(tauBins))])
tauPlot.extend([tauBins[-1] - (tauBins[-2] - tauBins[-1]) / 2])

ratioDiscardList = [4]  # [20,10,4,2] # discarding 1/ratioDiscardList %
discardPop = 'kc'  # kc: Keep Cells, pos: Positive Cells Pairs, neg: Negative Cell Pairs
discard = True  # if True: discard the 1/ratioDiscard % most (anti-)correlated if False: only keep the 1/ratioDiscard
# % most (anti-)correlated

if doPTI:
    loadStr = 'Zlinear'
else:
    loadStr = 'TempCorr'

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
    overlap_I = []
    overlap_II = []
    overlap_hmc = []
    overlap_between = []

    for iDay, day in enumerate(range(1, 10)):
        print([iDay, day])
        # need to add a condition to skip if day is not present
        sessIdx = np.where([d == day for d in dayCountList])[0]

        if len(sessIdx) < 1:
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

                # if 'HMC' not in sessList[iSessIdx][:3]:
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

            IsoToSave = [data_to_use, sAll, sessIndicAll, sessIndicFullAll, tauVec, idxTauPairs, tauVecTmp,
                         idxTauPairsTmp, keptCellsAll, discardCells, sConcat, sessIndic, sessIndicFull,idxPosList, idxNegList]


            intersection_I = np.intersect1d(idxNegList[0][2], idxNegList[0][4])
            intersection_II = np.intersect1d(idxNegList[0][1], idxNegList[0][3])
            intersection_hmc = np.intersect1d(idxNegList[0][0], idxNegList[0][5])
            intersection_between = np.intersect1d(intersection_I, intersection_II)

            rep_I = np.size(intersection_I) / np.maximum(np.size(idxNegList[0][2]), 1) * 100
            rep_II = np.size(intersection_II) / np.maximum(np.size(idxNegList[0][3]), 1) * 100
            rep_hmc = np.size(intersection_hmc) / np.maximum(np.size(idxNegList[0][0]), 1) * 100
            rep_between = np.size(intersection_between) / np.maximum(np.maximum(np.size(idxNegList[0][2]), np.size(idxNegList[0][3])),1) * 100

            print('intersection_I', intersection_I)
            print('intersection_II', intersection_II)
            print('intersection_hmc', intersection_hmc)

            print('intersection_I', rep_I)
            print('intersection_II', rep_II)
            print('intersection_hmc', rep_hmc)
            print('intersection_between', rep_between)

            overlap_I.append(rep_I)
            overlap_II.append(rep_II)
            overlap_hmc.append(rep_hmc)
            overlap_between.append(rep_between)

    overlap_tosave = [overlap_I, overlap_II, overlap_hmc, overlap_between]

    mdic_overlap = {"Overlap":overlap_tosave}

    pathSave_overlap = path_all + 'analysisFiles_Simon/' + 'Overlap_anticofiring' + '_' + animal + str(int(100 / ratio)) + '.mat'
    savemat(pathSave_overlap, mdic_overlap)

print('done done done')
