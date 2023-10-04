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

# update registration protocol

pathAll = '/Users/eliott/Dropbox/dDocuments/NYU/Fenton/MNS_ANALYSIS/Data/CircleSquare/analysisFiles/'
pathFigure = '/Users/eliott/Dropbox/dDocuments/NYU/Fenton/MNS_ANALYSIS/Figures/CircleSquare/savedFigs/'

BinNb = 12  # number of bins for place maps
nsep = 18
sampRate = 10  # Hz

AnimalList = ['M19']
nAni = len(AnimalList)

dayPlot = [1, 2, 3, 8, 9, 10, 15, 16, 17]

mapC_all = []
mapC_PCs = []

PercPcRec_all = []
PercPcDay_all = []

aLinCorr = plt.subplots(2, 2)[1]

cmap = plt.get_cmap('viridis')
cmap.set_bad('w', 1.)

aActivity = plt.subplots(2, 2)[1]
aHistActivity = plt.subplots(1, nAni)[1]

percStableDaySepAnimals = []
percTotStableDaySepAnimals = []
matchedDaysepAnimals = []
percDayStableDaySepAnimals = []
percDayTotStableDaySepAnimals = []
activityZScoreAnimals = []

for iAnimal, animal in enumerate(AnimalList):

    print(animal)

    loadpath = pathAll + 'MapCorr_' + str(int(BinNb)) + '_Bins_' + animal + '.file'

    with open(loadpath, 'rb') as f:
        output = pickle.load(f)
    #
    # [maps_corr, maps_corr_pc, day_count_pair, day_real_pair, sess_pair, day_sep, same_diff,
    #  LinMapsCorr, LinMapsCorrPC] = output[0]

    # [dayRealList, dayCountList, sessList, rateMapsList,
    #  PcIdxList, LinMapsList, occList, cohList, pInfoList] = output[1]

    dayRealList = output[1][0]
    dayCountList = output[1][1]
    sessList = output[1][2]
    rateMapsList = output[1][3]
    PcIdxList = output[1][4]
    occList = output[1][6]
    # cohList = output[1][7]
    # pInfoList = output[1][8]

    sTracked = output[2][2]

    output = []

    with open(pathAll + animal + '_assignments_cs.file', 'rb') as f:
        assignments = pickle.load(f)

    dayAni = np.unique(dayRealList)

    aActivity.flat[0].plot(dayPlot, [np.sum(np.isfinite(assignments[:, i])) for i in range(9)], '.-')
    aActivity.flat[1].plot(dayAni, [np.mean([np.sum(p) for p, d in zip(PcIdxList, dayRealList) if d == iDay]) for iDay in dayAni], '.-')
    aActivity.flat[2].plot(dayAni, [np.mean([np.mean(p) for p, d in zip(PcIdxList, dayRealList) if d == iDay]) for iDay in dayAni], '.-')
    aActivity.flat[3].plot(dayAni, [np.mean([sampRate * np.mean(s) for s, d in zip(sTracked, dayRealList) if d == iDay]) for iDay in dayAni], '.-')

    aActivity.flat[0].set_title('total # cells isolated')
    aActivity.flat[1].set_title('mean # PCs/rec')
    aActivity.flat[2].set_title('mean % cells PCs/rec')
    aActivity.flat[3].set_title('mean spike rate')


    plt.figure()
    activityMean = np.array([np.array(sampRate * np.mean(s, 1)) for s in sTracked])
    activityZScore = [(actMean - np.mean(actMean)) / np.std(actMean) for actMean in activityMean]
    activityZScoreDay = [np.concatenate(activityZScore[i*4:4+i*4]) for i in range(9)]
    activityZScoreAnimals.append(activityZScoreDay)

    for i_plot, actMean in enumerate(activityZScoreDay):
        aHistActivity[iAnimal].hist(actMean, np.linspace(-1.5, 4, 101),
                                    histtype='step', cumulative=True, density=True,
                                    color=[i_plot / 9, i_plot / 9, 1 - i_plot / 9])

    PcSameDayStability = []
    PcDiffDayStability = []
    PcOverlap = []

    # aPcCriteria = plt.subplots(9, 4)[1]
    for iDay, dayCount in enumerate(range(9)):

        idxDay = np.where([d == dayCount + 1 for d in dayCountList])[0]

        pairSame = [[idx_pair for idx_pair in itertools.combinations(idxDay, 2) if
                     (sessList[idx_pair[0]][:3] == sessList[idx_pair[1]][:3]) == iBool]
                    for iBool in [True, False]]

        PcSameDayStability_tmp = []
        for idxPair in pairSame[0]:
            PcSameDayStability_tmp.append(np.mean([PcIdxList[idxPair[0]][idx] for idx in np.where(PcIdxList[idxPair[1]])[0]]))
            PcSameDayStability_tmp.append(np.mean([PcIdxList[idxPair[1]][idx] for idx in np.where(PcIdxList[idxPair[0]])[0]]))

        PcDiffDayStability_tmp = []
        for idxPair in pairSame[1]:
            PcDiffDayStability_tmp.append(np.mean([PcIdxList[idxPair[0]][idx] for idx in np.where(PcIdxList[idxPair[1]])[0]]))
            PcDiffDayStability_tmp.append(np.mean([PcIdxList[idxPair[1]][idx] for idx in np.where(PcIdxList[idxPair[0]])[0]]))

        # to evaluate overlap
        isPcEnv = [np.all([PcIdxList[idx] for idx in iIdxPair], 0) for iIdxPair in pairSame[0]]
        # PcOverlap.append(np.mean([np.sum(np.all(isPcEnv, 0)) / x for x in np.sum(isPcEnv, 1)]))
        # np.mean([PcIdxList[pairSame[0][1][0]][idx] for idx in np.where(isPcEnv[0])[0]])

        PcDiffDayStability.append(PcDiffDayStability_tmp)
        PcSameDayStability.append(PcSameDayStability_tmp)

        for iRec, idxRec in enumerate(idxDay):
            rateMaps = rateMapsList[idxRec]
            occ = occList[idxRec]
            # pInfoTest = pInfoList[idxRec]
            # cohTest = cohList[idxRec]
            # PC_idx_tmp = [i for i, v in enumerate(PC_test) if v]
            # aPcCriteria[iDay][iRec].plot(p_info, coh, '.')
            # aPcCriteria[iDay][iRec].plot([0, 0.03], [0.5, 0.5])
            # aPcCriteria[iDay][iRec].plot([0.0005, 0.0005], [0, 1])

        if iDay>6 and np.all(np.sum(isPcEnv,1)>0):

            idxPc = np.concatenate([np.random.choice(np.where(isPcEnv[ii])[0], np.min([10, len(np.where(isPcEnv[ii])[0])]),
                                                    replace=False) for ii in [0, 1]])
            aPcDay = plt.subplots(4, len(idxPc))[1]
            for iRec, idxRec in enumerate(idxDay):
                for iCell, idxCell in enumerate(idxPc):

                    p_info = place_info_content(occList[idxRec], rateMapsList[idxRec][idxCell])

                    aPcDay[iRec][iCell].imshow(nan_gaussian_filter(rateMapsList[idxRec][idxCell], 0.5), cmap=cmap)
                    aPcDay[iRec][iCell].axis('off')
                    if PcIdxList[idxRec][idxCell]:
                        aPcDay[iRec][iCell].set_title('PC ')
                        # aPcDay[iRec][iCell].set_title('PC ' + '{:.4f}'.format(p_info))
                    # else:
                        # aPcDay[iRec][iCell].set_title('{:.4f}'.format(p_info))
    percDayStable = []
    percDayTot = []
    DayStableDaySep = []
    DayStableDayPair = []

    for nComb, dayCountComb in enumerate(itertools.combinations(range(9), 2)):

        DayStableDaySep.append(dayCountComb[1] - dayCountComb[0])
        DayStableDayPair.append(dayCountComb)

        # find recording number for each day of the pair of days (dayCountComb)
        idxDayComb = [np.where([d == dayCount + 1 for d in dayCountList])[0] for dayCount in dayCountComb]

        # for each day find the pair of rec in the same env within a day(?)
        pairSameComb = [[idx_pair for idx_pair in itertools.combinations(idxDay, 2) for env in ['RCT', 'CYL']
                          if sessList[idx_pair[0]][:3] == sessList[idx_pair[1]][:3]
                         and sessList[idx_pair[0]][:3] == env]
                        for idxDay in idxDayComb]

        isPcEnv = [[np.all([PcIdxList[idx] for idx in iIdxPair], 0) for iIdxPair in pairSame]
                   for pairSame in pairSameComb]

        both_dual = np.all(~np.isnan(assignments[:, dayCountComb]), 1)
        # n_matched.append(sum(both_dual))
        CellList = [[n_cell for n_cell in assignments[both_dual, day_id].astype(int)]
                    for day_id in dayCountComb]

        PC_assign = [[[isPcEnv[iDay][env][n_cell]
                      for n_cell in CellList[iDay]]
                     for env in range(2)]
                     for iDay in range(2)]

        percDayStable.append(np.mean([np.mean([PC_assign[0][env][idx] for idx in np.where(PC_assign[1][env])[0]]) for env in range(2)]))
        percDayTot.append(np.mean([np.sum([PC_assign[0][env][idx] for idx in np.where(PC_assign[1][env])[0]])/np.sum(isPcEnv[1][env]) for env in range(2)]))

        # matchFind = np.where(assignments[:, iDayCount] == 1)[0][0]

    percDayStableDaySep = [[np.nanmean([np.nanmean(p) for idx, p in enumerate(percDayStable)
                                         if DayStableDaySep[idx] == i_day
                                        if np.min(DayStableDayPair[idx]) >= day_start])
                            for i_day in range(9)]
                          for day_start in [0, 3, 6]]

    percTotDayStableDaySep = [[np.nanmean([np.nanmean(p) for idx, p in enumerate(percDayTot)
                                        if DayStableDaySep[idx] == i_day
                                        if np.min(DayStableDayPair[idx]) >= day_start])
                            for i_day in range(9)]
                          for day_start in [0, 3, 6]]

    percDayStableDaySepAnimals.append(percDayStableDaySep)
    percDayTotStableDaySepAnimals.append(percTotDayStableDaySep)

    # a = plt.subplots(1, 2)[1]
    # for i in range(3):
    #     a[0].plot(percDayStableDaySep[i])
    #     a[1].plot(percTotDayStableDaySep[i])
    #
    # plt.figure()
    # plt.plot(np.nanmean(PcSameDayStability, 1))
    # plt.plot(np.nanmean(PcDiffDayStability, 1))
    # plt.plot(PcOverlap)

    dayRealPair = []
    daySep = []
    sameDiff = []
    dayCountPair = []
    perc_stable = []
    perc_tot = []
    n_matched =[]

    for n_comb, i_comb in enumerate(itertools.combinations(range(len(dayCountList)), 2)):

        dayReal = tuple([dayRealList[i] for i in i_comb])
        sessId = tuple([sessList[i] for i in i_comb])
        dayCount = tuple([dayCountList[i]-1 for i in i_comb])

        daySep.append(np.abs(dayReal[1] - dayReal[0]))

        if dayReal[1] - dayReal[0] < 0:
            print('day 2 < day 1 !')

        sameDiff.append(sessId[0][:3] == sessId[1][:3])
        dayCountPair.append(dayCount)

        both_dual = np.all(~np.isnan(assignments[:, dayCount]), 1)
        n_matched.append(sum(both_dual))
        CellList = [[n_cell for n_cell in assignments[both_dual, day_id].astype(int)]
                    for day_id in dayCount]

        PC_assign = [[PcIdxList[rec_id][n_cell]
                      for n_cell in CellList[i_rec]] for i_rec, rec_id in enumerate(i_comb)]

        perc_stable.append(np.mean([PC_assign[0][idx] for idx in np.where(PC_assign[1])[0]]))
        perc_tot.append(np.sum([PC_assign[0][idx] for idx in np.where(PC_assign[1])[0]])/np.sum(PcIdxList[i_comb[1]]))

    # for mean in day sep but with different start time
    percStableDaySep = [[[np.nanmean([np.nanmean(p) for idx, p in enumerate(perc_stable)
                                        if sameDiff[idx] == ibool if daySep[idx] == i_day
                                        if np.min(dayCountPair[idx]) >= day_start])
                            for i_day in range(nsep)]
                           for ibool in [True, False]]
                          for day_start in [0, 3, 6]]

    percTotStableDaySep = [[[np.nanmean([np.nanmean(p) for idx, p in enumerate(perc_tot)
                                        if sameDiff[idx] == ibool if daySep[idx] == i_day
                                        if np.min(dayCountPair[idx]) >= day_start])
                            for i_day in range(nsep)]
                           for ibool in [True, False]]
                          for day_start in [0, 3, 6]]

    matchedDaySep = [[[np.nanmean([np.nanmean(p) for idx, p in enumerate(n_matched)
                                        if sameDiff[idx] == ibool if daySep[idx] == i_day
                                        if np.min(dayCountPair[idx]) >= day_start])
                            for i_day in range(nsep)]
                           for ibool in [True, False]]
                          for day_start in [0, 3, 6]]

    # a = plt.subplots(1, 3)[1]
    # for i in range(3):
    #     a[0].plot(percStableDaySep[i][0])
    # for i in range(3):
    #     a[0].plot(percStableDaySep[i][1])
    #
    # for i in range(3):
    #     a[1].plot(percTotStableDaySep[i][0])
    # for i in range(3):
    #     a[1].plot(percTotStableDaySep[i][1])
    #
    # for i in range(3):
    #     a[2].plot(matchedDaySep[i][0])

    percStableDaySepAnimals.append(percStableDaySep)
    percTotStableDaySepAnimals.append(percTotStableDaySep)
    matchedDaysepAnimals.append(matchedDaySep)

plt.figure()
iDiff = 0
iDayMin = 0
wkDayIdx = [[0], [1], [2], [5, 6, 7, 8, 9], [12, 13, 14, 15, 16]]
barCenterSep = [[0, 8, 12, 20, 24], [1, 9, 13, 21, 25]]

for iDiff in range(2):
    plt.bar(barCenterSep[iDiff],
            np.mean([[np.mean([m[iDayMin][iDiff][iIdx] for iIdx in idxList]) for idxList in wkDayIdx]
                     for m in percStableDaySepAnimals], 0))

    plt.errorbar(barCenterSep[iDiff],
                 np.mean([[np.mean([m[iDayMin][iDiff][iIdx] for iIdx in idxList]) for idxList in wkDayIdx]
                          for m in percStableDaySepAnimals], 0),
                 np.std([[np.mean([m[iDayMin][iDiff][iIdx] for iIdx in idxList]) for idxList in wkDayIdx]
                         for m in percStableDaySepAnimals], 0) / np.sqrt(3),
                 ls='', color='k')
    plt.title('%PCs remaining PCs')

a = plt.subplots(1, 2)[1]
for i in range(3):
    a[0].errorbar([*range(nsep)], np.mean(percStableDaySepAnimals, 0)[i][0], np.std(percStableDaySepAnimals, 0)[i][0]/np.sqrt(3))
for i in range(3):
    a[0].errorbar([*range(nsep)], np.mean(percStableDaySepAnimals, 0)[i][1], np.std(percStableDaySepAnimals, 0)[i][1]/np.sqrt(3))

for i in range(3):
    a[1].errorbar([*range(nsep)], np.mean(percTotStableDaySepAnimals, 0)[i][0], np.std(percTotStableDaySepAnimals, 0)[i][0]/np.sqrt(3))
for i in range(3):
    a[1].errorbar([*range(nsep)], np.mean(percTotStableDaySepAnimals, 0)[i][1], np.mean(percTotStableDaySepAnimals, 0)[i][1]/np.sqrt(3))

a = plt.subplots(1, 2)[1]
for i in range(3):
    a[0].errorbar([0, 1, 2, 7, 8, 9, 14, 15, 16], np.mean(percDayStableDaySepAnimals, 0)[i],
                  np.std(percDayStableDaySepAnimals, 0)[i]/np.sqrt(3))

    a[1].errorbar([0, 1, 2, 7, 8, 9, 14, 15, 16], np.mean(percDayTotStableDaySepAnimals, 0)[i],
                  np.std(percDayTotStableDaySepAnimals, 0)[i] / np.sqrt(3))

a = plt.subplots(1, 2)[1]
activityZScoreConcat = [np.concatenate([r[i] for r in activityZScoreAnimals]) for i in range(len(activityZScoreAnimals[0]))]
for i_plot, actMean in enumerate(activityZScoreConcat):
    a[0].hist(actMean, np.linspace(-1.5, 4, 51),
             histtype='step', cumulative=False, density=True,
             color=[i_plot / 9, i_plot / 9, 1 - i_plot / 9])

    a[1].hist(actMean, np.linspace(-1.5, 4, 301),
             histtype='step', cumulative=True, density=True,
             color=[i_plot / 9, i_plot / 9, 1 - i_plot / 9])

# for i in range(3):
#     a[2].plot(np.mean(matchedDaysepAnimals, 0)[i][0])


plt.show(block=True)