import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import pickle
import itertools
import time
import scipy.stats as scstat
from sklearn import linear_model
import pandas as pd
import random

path_all = '/Users/eliott/Dropbox/dDocuments/NYU/Fenton/MNS_ANALYSIS/Data/CircleSquare/'
pathFigure = '/Users/eliott/Dropbox/dDocuments/NYU/Fenton/MNS_ANALYSIS/Figures/CircleSquare/savedFigs/z/'
TimeBinSize_AVC = 10 * 1000
TimeBinSize_TAU = 1000
mSize = 1
BinNb = 12
SaveFig = True

animalList = ['M19','M20', 'M29','M35','M34', 'M39']
nAni = len(animalList)

plotEveryRecCorr = False
day_thresh = -1
thresh = 0.02*(TimeBinSize_TAU/1000)

zTauAniAll = []
rTauAniAll = []
mapCoAniAll = []
zTauPcAniAll = []
rTauPcAniAll = []
mapCoPcAniAll = []
zTFitList = [[], []]
zTCorrList = [[], []]
zMFitList = [[], []]
zMCorrList = [[], []]
zCorr = []
zStabilityCorr = []

fAllAni, aAllAni = plt.subplots(nrows=2, ncols=2, figsize=(9, 9)) #zMap scatter plot
fZTAllAni, aZTAllAni = plt.subplots(nrows=2, ncols=2, figsize=(9, 9)) #rTau scatter plot
fZM,aZM = plt.subplots(nAni, 5)
fZT,aZT = plt.subplots(nAni, 5)
fZStabilityCorr, aZStabilityCorr = plt.subplots(1,1)
fTauLogDist, aTauLogDist = plt.subplots(2,nAni)

for iAni, animal in enumerate(animalList):

    print(animal)

    zMCorr = []
    zMCorrPc = []
    zTCorr = []
    zTCorrPc = []
    zAllDay = [[],[]]

    analStr = 'zLin'
    # analStr = 'RandzLin'
    # analStr = 'z'
    # analStr = 'Randz'

    with open(path_all + 'analysisFiles/' + 'Zlinear_List_' +
              str(int(TimeBinSize_AVC)) + 'msAV_' + str(int(TimeBinSize_TAU)) + 'msKT' +
              '_CircleSquare_' + animal +'.file', 'rb') as f:
        outputList = pickle.load(f)

    # outputList = [[day_real_list, day_count_list, sess_list, S_conv_avc_list, S_conv_tau_list, isSeparatedList,
    #                 tauPairsSingleList, tauVecSingleList],
    #                [Sector_list, Shock_list]]
    dayCountListTau = outputList[0][1]
    sessListTau = outputList[0][2]
    sConv = outputList[0][4]
    isSeparatedList = outputList[0][5]
    tauPairsSingleList = outputList[0][6]
    tau = outputList[0][7]
    outputList = []

    sConvTauList = sConv[1] # for rate calculation
    zTauVecList = tau[0]  # [0] for z
    regTauVecList= tau[1]  # [1] for regular tau


    loadpath = path_all + 'analysisFiles/' + 'MapCorr_' + str(int(BinNb)) + '_Bins_' + animal + '.file'
    with open(loadpath, 'rb') as f:
        output = pickle.load(f)

    # to_save = [[maps_corr, [maps_corr_pc_any, maps_corr_pc_all, maps_corr_pc_not], day_count_pair, day_real_pair,
    #             sess_pair, day_sep, same_diff,
    #             LinMapsCorr, LinMapsCorrPC, cellCount, trainingPair],
    #            [day_real_list, day_count_list, sess_list, rate_maps_list, PC_idx_list, LinMapsList, OccList,
    #             cohList, pInfoList, trainingList],
    #            [X_tracked_list, Y_tracked_list, s_tracked_list],
    #            [nRepeat]]
    dayCountListMap = output[1][1]
    sessListMap = output[1][2]
    rateMapList = output[1][3]
    PcIdxList = output[1][4]
    output = []

    zTauAni = []
    rTauAni = []
    mapCoAni = []
    zTauPcAni = []
    rTauPcAni = []
    mapCoPcAni = []
    zCorrAni = [[], [], [], []]

    fZScatter, aZScatter = plt.subplots(2, 6)
    fZScatter2, aZScatter2 = plt.subplots(2, 6)
    for iDay, day in enumerate(range(4,10)):
        zDayList = []
        zDayList2 = [[], []]
        zFilterList = []
        pcList = []
        for session in ['RCT1', 'RCT2', 'CYL1', 'CYL2']:

            mapRecIdxList = np.where([d == day and session in iSess for d, iSess in zip(dayCountListMap, sessListMap)])[0]
            tauRecIdxList = np.where([d == day and session in iSess for d, iSess in zip(dayCountListTau, sessListTau)])[0]

            if len(mapRecIdxList) == 0 or len(tauRecIdxList) == 0:
                continue

            for tauRecIdx,mapRecIdx in zip(tauRecIdxList,mapRecIdxList):

                zTauVec = zTauVecList[tauRecIdx]
                zDayList.append([z[0] for z in zTauVec])
                rTauVec = regTauVecList[tauRecIdx]
                tauPairs = tauPairsSingleList[tauRecIdx]
                meanRate = [np.mean([s > 0 for s in c]) for c in sConvTauList[tauRecIdx]]
                isSep = isSeparatedList[tauRecIdx]

                PcIdx = PcIdxList[mapRecIdx]
                rateMaps = rateMapList[mapRecIdx]

                #create rate filters for tau and z
                rate = [meanRate[p[0]]>thresh and meanRate[p[1]]>thresh for p in tauPairs]
                # filter for z
                zNan = [np.isnan(t[0]) for t in zTauVec]
                zAntiCorr = [t[0]<2 for t in rTauVec]
                zTauFilter = [i and r and a and not n for i,n,r,a in zip(isSep, zNan, rate, zAntiCorr)]
                zFilterList.append(zTauFilter)
                #filter for tau
                rNan = [np.isnan(t[0]) for t in rTauVec]
                # rAntiCorr = [t[0]<-0.01 for t in rTauVec]
                rTauFilter = [i and r and not n for i,n,r in zip(isSep, rNan, rate)]

                # filter tau and z
                zTauVecFiltered = [t[0] for t, fr, fz in zip(zTauVec, rTauFilter, zTauFilter) if fz] #(fr and fz)]
                rTauVecFiltered = [t[0] for t, fr, fz in zip(rTauVec, rTauFilter, zTauFilter) if fz] #(fr and fz)]
                tauPairsFiltered = [p for p,fr,fz in zip(tauPairs, rTauFilter, zTauFilter) if fz] #(fr and fz)]

                # split according to PC or not PC
                zTauVecPc = [t for t, p in zip(zTauVecFiltered, tauPairsFiltered) if (PcIdx[p[0]] or PcIdx[p[1]])]
                rTauVecPc = [t for t, p in zip(rTauVecFiltered, tauPairsFiltered) if (PcIdx[p[0]] or PcIdx[p[1]])]
                tauPairsPc = [p for p in tauPairsFiltered if (PcIdx[p[0]] or PcIdx[p[1]])]

                zTauVecNonPc = [t for t, p in zip(zTauVecFiltered, tauPairsFiltered) if (not PcIdx[p[0]] and not PcIdx[p[1]])]
                rTauVecNonPc = [t for t, p in zip(rTauVecFiltered, tauPairsFiltered) if (not PcIdx[p[0]] and not PcIdx[p[1]])]
                tauPairsNonPc = [p for p in tauPairsFiltered if (not PcIdx[p[0]] and not PcIdx[p[1]])]

                pcList.append(PcIdx)

                # compute rate maps Correlation
                mask = ~np.isnan(rateMaps[0].flatten())
                mapCo = [np.corrcoef(rateMaps[p[0]].flatten()[mask], rateMaps[p[1]].flatten()[mask])[0, 1]
                         for p in tauPairsFiltered]
                mapCoPc = [np.corrcoef(rateMaps[p[0]].flatten()[mask], rateMaps[p[1]].flatten()[mask])[0, 1]
                         for p in tauPairsPc]
                mapCoNonPc = [np.corrcoef(rateMaps[p[0]].flatten()[mask], rateMaps[p[1]].flatten()[mask])[0, 1]
                         for p in tauPairsNonPc]

                # extend to have list for all recs for this animal
                # tauAni.extend(tauVecFiltered)
                zTauAni.extend(zTauVecNonPc)
                zTauPcAni.extend(zTauVecPc)
                rTauAni.extend(rTauVecNonPc)
                rTauPcAni.extend(rTauVecPc)
                # mapCoAni.extend(mapCo)
                mapCoAni.extend(mapCoNonPc)
                mapCoPcAni.extend(mapCoPc)

                aZT[iAni][0].scatter(zTauVecNonPc, rTauVecNonPc, s=1, alpha=0.3)
                regr = linear_model.LinearRegression()
                regr.fit(np.reshape(zTauVecNonPc, [-1, 1]), rTauVecNonPc)
                x1 = np.min(zTauVecNonPc)
                x2 = np.max(zTauVecNonPc)
                aZT[iAni][0].plot([x1, x2], [x1 * regr.coef_ + regr.intercept_, x2 * regr.coef_ + regr.intercept_])
                zTCorr.append(np.corrcoef(zTauVecNonPc, rTauVecNonPc)[0, 1])

                aZT[iAni][1].scatter(zTauVecPc, rTauVecPc, s=1, alpha=0.3)
                regr = linear_model.LinearRegression()
                regr.fit(np.reshape(zTauVecPc, [-1, 1]), rTauVecPc)
                x1 = np.min(zTauVecPc)
                x2 = np.max(zTauVecPc)
                aZT[iAni][1].plot([x1, x2], [x1 * regr.coef_ + regr.intercept_, x2 * regr.coef_ + regr.intercept_])
                zTCorrPc.append(np.corrcoef(zTauVecPc, rTauVecPc)[0, 1])

                # plot with regression line and compute correlation
                # non PC
                aZM[iAni][0].scatter(zTauVecNonPc, mapCoNonPc, s=1, alpha=0.3)
                if len(zTauVecNonPc)>5 and len(mapCoNonPc)>5:
                    tauVal = [t for t, x in zip(zTauVecNonPc,mapCoNonPc) if np.isfinite(x)]
                    mapVal = [x for x in mapCoNonPc if np.isfinite(x)]
                    regr = linear_model.LinearRegression()
                    regr.fit(np.reshape(tauVal, [-1, 1]),mapVal)
                    x1 = np.min(tauVal)
                    x2 = np.max(tauVal)
                    aZM[iAni][0].plot([x1, x2], [x1* regr.coef_ + regr.intercept_, x2* regr.coef_ + regr.intercept_])
                    zMCorr.append(np.corrcoef(zTauVecNonPc,mapCoNonPc)[0,1])

                # PC
                aZM[iAni][1].scatter(zTauVecPc, mapCoPc, s=1, alpha=0.3)
                if len(zTauVecPc)>5 and len(mapCoPc)>5:
                    tauVal = [t for t, x in zip(zTauVecPc,mapCoPc) if np.isfinite(x)]
                    regr = linear_model.LinearRegression()
                    regr.fit(np.reshape(tauVal, [-1, 1]),[x for x in mapCoPc if np.isfinite(x)])
                    x1 = np.min(tauVal)
                    x2 = np.max(tauVal)
                    aZM[iAni][1].plot([x1, x2], [x1* regr.coef_ + regr.intercept_, x2* regr.coef_ + regr.intercept_])
                    zMCorrPc.append(np.corrcoef(zTauVecPc, mapCoPc)[0,1])

        #scatter plot of z over session within a day to evaluate stability
        if len(zDayList)>3:

            zFilter = [np.all(z) for z in zip(*zFilterList)]
            # print(np.mean(zFilter))
            for sameDiff, idxPairList in enumerate([[[0,1],[2,3]],[[0,2],[0,3],[1,2],[1,3]]]):
                for idxPair in idxPairList:

                    altFilter = [f1 and f2 for f1,f2 in zip(zFilterList[idxPair[0]], zFilterList[idxPair[1]])]
                    pcs = [p1 or p2 for p1, p2 in zip(pcList[idxPair[0]], pcList[idxPair[1]])]
                    nPcs = [(not p1) and (not p2) for p1, p2 in zip(pcList[idxPair[0]], pcList[idxPair[1]])]

                    tauPairsAltFilt = [p for p,fz in zip(tauPairs, altFilter) if fz]

                    z1 = [z1 for z1, z2, f in zip(zDayList[idxPair[0]], zDayList[idxPair[1]], altFilter) if f]
                    z2 = [z2 for z1, z2, f in zip(zDayList[idxPair[0]], zDayList[idxPair[1]], altFilter) if f]

                    z1Pc = [t for t, p in zip(z1, tauPairsAltFilt) if (pcs[p[0]] or pcs[p[1]])]
                    z2Pc = [t for t, p in zip(z2, tauPairsAltFilt) if (pcs[p[0]] or pcs[p[1]])]
                    z1NPc = [t for t, p in zip(z1, tauPairsAltFilt) if (nPcs[p[0]] and nPcs[p[1]])]
                    z2NPc = [t for t, p in zip(z2, tauPairsAltFilt) if (nPcs[p[0]] and nPcs[p[1]])]

                    if sameDiff == 0:
                        zAllDay[0].append([z1NPc, z2NPc, z1Pc, z2Pc])
                    else:
                        zAllDay[1].append([z1NPc, z2NPc, z1Pc, z2Pc])

                    if len(z1Pc) > 2:
                        regr = linear_model.LinearRegression()
                        regr.fit(np.reshape(z1Pc, [-1, 1]), z2Pc)
                        aZScatter2[sameDiff][iDay].scatter(z1Pc, z2Pc, s=1, color='r', alpha=0.5)
                        x1 = np.min(z1Pc)
                        x2 = np.max(z1Pc)
                        aZScatter2[sameDiff][iDay].plot([x1, x2], [x1* regr.coef_ + regr.intercept_, x2* regr.coef_ + regr.intercept_], 'r')
                        if sameDiff == 0:
                            zCorrAni[1].append(np.corrcoef(z1Pc, z2Pc)[0,1])
                        else:
                            zCorrAni[3].append(np.corrcoef(z1Pc, z2Pc)[0, 1])

                    if len(z1NPc)>2:
                        regr = linear_model.LinearRegression()
                        regr.fit(np.reshape(z1NPc, [-1, 1]), z2NPc)
                        aZScatter2[sameDiff][iDay].scatter(z1NPc, z2NPc, s=1, color='k', alpha=0.5)
                        x1 = np.min(z1NPc)
                        x2 = np.max(z1NPc)
                        aZScatter2[sameDiff][iDay].plot([x1, x2], [x1* regr.coef_ + regr.intercept_, x2* regr.coef_ + regr.intercept_], 'k')
                        if sameDiff == 0:
                            zCorrAni[0].append(np.corrcoef(z1NPc, z2NPc)[0, 1])
                        else:
                            zCorrAni[2].append(np.corrcoef(z1NPc, z2NPc)[0, 1])

                    # print(np.mean(altFilter))
                    regr = linear_model.LinearRegression()
                    regr.fit(np.reshape(z1, [-1, 1]), z2)
                    aZScatter[sameDiff][iDay].scatter(z1, z2, s=1)
                    x1 = np.min(z1)
                    x2 = np.max(z1)
                    aZScatter[sameDiff][iDay].plot([x1, x2], [x1* regr.coef_ + regr.intercept_, x2* regr.coef_ + regr.intercept_])
                    aZScatter[sameDiff][iDay].set_ylim([-1, 1])
                    aZScatter[sameDiff][iDay].set_xlim([-1, 1])

    # empty list to save memory
    tau = dayCountListTau = sConvTauList = tauPairsSingleList = isSeparatedList = sessListTau = sConv = []

    zCorr.append(zCorrAni)

    c = [[np.corrcoef([z for zRec in zAllDay[isDiff] for z in zRec[0]], [z for zRec in zAllDay[isDiff] for z in zRec[1]])[0][1],
          np.corrcoef([z for zRec in zAllDay[isDiff] for z in zRec[2]], [z for zRec in zAllDay[isDiff] for z in zRec[3]])[0][1]]
         for isDiff in [0,1]]
    aZStabilityCorr.plot([i for _ in c for i in _], '.')
    zStabilityCorr.append([i for _ in c for i in _])

    # plot average of z/Map correlations
    aZM[iAni][2].bar([0, 1], [np.mean(zMCorr), np.mean(zMCorrPc)])
    aZM[iAni][2].errorbar([0, 1], [np.mean(zMCorr), np.mean(zMCorrPc)], [np.std(zMCorr), np.std(zMCorrPc)], ls = '', ecolor = 'k')

    # plot average of z/Tau correlations
    aZT[iAni][2].bar([0, 1], [np.mean(zTCorr), np.mean(zTCorrPc)])
    aZT[iAni][2].errorbar([0, 1], [np.mean(zTCorr), np.mean(zTCorrPc)], [np.std(zTCorr), np.std(zTCorrPc)], ls = '', ecolor = 'k')

    # extend to have list for animals
    zTauAniAll.extend(zTauAni)
    rTauAniAll.extend(rTauAni)
    mapCoAniAll.extend(mapCoAni)
    zTauPcAniAll.extend(zTauPcAni)
    rTauPcAniAll.extend(rTauPcAni)
    mapCoPcAniAll.extend(mapCoPcAni)

    # compute and plot regression of animal (all recs) - Non PC - z/Map
    tauVal = [t for t, x in zip(zTauAni, mapCoAni) if np.isfinite(x)]
    mapVal = [x for x in mapCoAni if np.isfinite(x)]
    regr = linear_model.LinearRegression()
    regr.fit(np.reshape(tauVal, [-1, 1]), mapVal)
    x1 = np.min(tauVal)
    x2 = np.max(tauVal)
    for ax in [[aZM[iAni][3], aZM[iAni][4]], aAllAni[0]]:
        ax[0].scatter(zTauAni, mapCoAni, s=1, color='k', alpha=0.2)
        ax[0].plot([x1, x2], [x1* regr.coef_ + regr.intercept_, x2* regr.coef_ + regr.intercept_], 'm--')

    # then compute regression and corr for each intervals - Non PC - z/Map
    fitListTmp = []
    corrListTmp = []
    for tauInterval in [[-2,2], [-2, 0], [0, 2]]:
        tSeg = [t for t in tauVal if tauInterval[0] < t <= tauInterval[1]]
        mSeg = [m for t, m in zip(tauVal, mapVal) if tauInterval[0] < t <= tauInterval[1]]
        fFit = np.polyfit(tSeg, mSeg, 1)
        fitListTmp.append(fFit)
        corrListTmp.append(np.corrcoef(tSeg, mSeg)[0,1])
    zMFitList[0].append(fitListTmp)
    zMCorrList[0].append(corrListTmp)

    # compute and plot regression of animal (all recs) - PCs - z/Map
    tauVal = [t for t, x in zip(zTauPcAni, mapCoPcAni) if np.isfinite(x)]
    mapVal = [x for x in mapCoPcAni if np.isfinite(x)]
    regr = linear_model.LinearRegression()
    regr.fit(np.reshape(tauVal, [-1, 1]), mapVal)
    x1 = np.min(tauVal)
    x2 = np.max(tauVal)
    for ax in [[aZM[iAni][3], aZM[iAni][4]], aAllAni[0]]:
        ax[1].scatter(zTauPcAni, mapCoPcAni, s=1, color='k', alpha=0.2)
        ax[1].plot([x1, x2], [x1* regr.coef_ + regr.intercept_, x2* regr.coef_ + regr.intercept_], 'm--')

    # then compute regression and corr for each intervals - PCs - z/Map
    fitListTmp = []
    corrListTmp = []
    for tauInterval in [[-2,2], [-2, 0], [0, 2]]:
        tSeg = [t for t in tauVal if tauInterval[0] < t <= tauInterval[1]]
        mSeg = [m for t, m in zip(tauVal, mapVal) if tauInterval[0] < t <= tauInterval[1]]
        f = np.polyfit(tSeg, mSeg, 1)
        fitListTmp.append(f)
        corrListTmp.append(np.corrcoef(tSeg, mSeg)[0, 1])
    zMFitList[1].append(fitListTmp)
    zMCorrList[1].append(corrListTmp)

    # compute and plot regression of animal (all recs) - Non PC - z/Tau
    regr = linear_model.LinearRegression()
    regr.fit(np.reshape(zTauAni, [-1, 1]), rTauAni)
    x1 = np.min(zTauAni)
    x2 = np.max(zTauAni)
    for ax in [[aZT[iAni][3], aZT[iAni][4]], aZTAllAni[0]]:
        ax[0].scatter(zTauAni, rTauAni, s=1, color='k', alpha=0.2)
        ax[0].plot([x1, x2], [x1 * regr.coef_ + regr.intercept_, x2 * regr.coef_ + regr.intercept_], 'm--')

    # then compute regression and corr for each intervals - Non PC - z/Map
    fitListTmp = []
    corrListTmp = []
    for tauInterval in [[-2,2], [-2, 0], [0, 2]]:
        tSeg = [t for t in zTauAni if tauInterval[0] < t <= tauInterval[1]]
        mSeg = [m for t, m in zip(zTauAni, rTauAni) if tauInterval[0] < t <= tauInterval[1]]
        f = np.polyfit(tSeg, mSeg, 1)
        fitListTmp.append(f)
        corrListTmp.append(np.corrcoef(tSeg, mSeg)[0,1])
    zTFitList[0].append(fitListTmp)
    zTCorrList[0].append(corrListTmp)

    # compute and plot regression of animal (all recs) - PCs - z/Tau
    regr = linear_model.LinearRegression()
    regr.fit(np.reshape(zTauPcAni, [-1, 1]), rTauPcAni)
    x1 = np.min(zTauPcAni)
    x2 = np.max(zTauPcAni)
    for ax in [[aZT[iAni][3], aZT[iAni][4]], aZTAllAni[0]]:
        ax[1].scatter(zTauPcAni, rTauPcAni, s=1, color='k', alpha=0.2)
        ax[1].plot([x1, x2], [x1 * regr.coef_ + regr.intercept_, x2 * regr.coef_ + regr.intercept_], 'm--')

    # then compute regression and corr for each intervals - Non PC - z/Map
    fitListTmp = []
    corrListTmp = []
    for tauInterval in [[-2,2], [-2, 0], [0, 2]]:
        tSeg = [t for t in zTauPcAni if tauInterval[0] < t <= tauInterval[1]]
        mSeg = [m for t, m in zip(zTauPcAni, rTauPcAni) if tauInterval[0] < t <= tauInterval[1]]
        f = np.polyfit(tSeg, mSeg, 1)
        fitListTmp.append(f)
        corrListTmp.append(np.corrcoef(tSeg, mSeg)[0,1])
    zTFitList[1].append(fitListTmp)
    zTCorrList[1].append(corrListTmp)

    # plot z log distribution
    aTauLogDist[0][iAni].hist(zTauAni, 100, log=True)
    aTauLogDist[1][iAni].hist(zTauPcAni, 100, log=True)
    aTauLogDist[0][iAni].set_title(animal+'_nonPCs')
    aTauLogDist[1][iAni].set_title(animal+'_PCs')
    fTauLogDist.suptitle('tau distribution')

# plot correlation of z btw recordings
aZCorr = plt.subplots(1,nAni)[1]
for iAni in range(nAni):
    aZCorr[iAni].boxplot([[x for x in xList if np.isfinite(x)] for xList in zCorr[iAni]])

# plot 2d hist
aZTAllAni[1][0].hist2d(zTauAniAll,rTauAniAll,100)
aZTAllAni[1][1].hist2d(zTauPcAniAll,rTauPcAniAll,100)

# compute and plot linear regression - non PC - zTau/rTau
regr = linear_model.LinearRegression()
zTauVal = [zt for zt, rt in zip(zTauAniAll, rTauAniAll) if np.isfinite(zt) and np.isfinite(rt)]
rTauVal = [rt for zt, rt in zip(zTauAniAll, rTauAniAll) if np.isfinite(zt) and np.isfinite(rt)]
regr.fit(np.reshape(zTauVal, [-1, 1]),rTauVal)
x1 = np.min(zTauVal)
x2 = np.max(zTauVal)
aZTAllAni[0][0].plot([x1, x2], [x1* regr.coef_ + regr.intercept_, x2* regr.coef_ + regr.intercept_], 'r')

# compute and plot linear regression - PC - zTau/rTau
regr = linear_model.LinearRegression()
zTauVal = [zt for zt, rt in zip(zTauPcAniAll, rTauPcAniAll) if np.isfinite(zt) and np.isfinite(rt)]
rTauVal = [rt for zt, rt in zip(zTauPcAniAll, rTauPcAniAll) if np.isfinite(zt) and np.isfinite(rt)]
regr.fit(np.reshape(zTauVal, [-1, 1]),rTauVal)
x1 = np.min(zTauVal)
x2 = np.max(zTauVal)
aZTAllAni[0][1].plot([x1, x2], [x1* regr.coef_ + regr.intercept_, x2* regr.coef_ + regr.intercept_], 'r')

# plot reg at different tau segments - Non PC - z/Tau
rr = np.linspace(np.min(zTauAniAll),np.max(zTauAniAll),100)
for tauInterval in [[-2,2], [-2, 0], [0, 2]]:
    tSeg = [t for t in zTauAniAll if tauInterval[0]<t<=tauInterval[1]]
    mSeg = [m for t,m in zip(zTauAniAll,rTauAniAll) if tauInterval[0]<t<=tauInterval[1]]
    f = np.polyfit(tSeg,mSeg,1)
    aZTAllAni[1][0].plot(rr, [r*f[0]+f[1] for r in rr], 'r-.')

# plot reg at different tau segments - Non PC - z/Tau
rr = np.linspace(np.min(zTauPcAniAll),np.max(zTauPcAniAll),100)
for tauInterval in [[-2,2], [-2, 0], [0, 2]]:
    tSeg = [t for t in zTauPcAniAll if tauInterval[0]<t<=tauInterval[1]]
    mSeg = [m for t,m in zip(zTauPcAniAll,rTauPcAniAll) if tauInterval[0]<t<=tauInterval[1]]
    f = np.polyfit(tSeg,mSeg,1)
    aZTAllAni[1][1].plot(rr, [r*f[0]+f[1] for r in rr], 'r-.')

# compute and plot linear regression - non PC - z/Map
regr = linear_model.LinearRegression()
tauVal = [t for t, x in zip(zTauAniAll, mapCoAniAll) if np.isfinite(x)]
mapVal = [x for x in mapCoAniAll if np.isfinite(x)]
regr.fit(np.reshape(tauVal, [-1, 1]),mapVal)
x1 = np.min(tauVal)
x2 = np.max(tauVal)
aAllAni[0][0].plot([x1, x2], [x1* regr.coef_ + regr.intercept_, x2* regr.coef_ + regr.intercept_], 'r')

# same at different tau segments
rr = np.linspace(np.min(tauVal),np.max(tauVal),100)
for tauInterval in [[-2,2], [-2, 0], [0, 2]]:
    tSeg = [t for t in tauVal if tauInterval[0]<t<=tauInterval[1]]
    mSeg = [m for t,m in zip(tauVal,mapVal) if tauInterval[0]<t<=tauInterval[1]]
    f = np.polyfit(tSeg,mSeg,1)
    aAllAni[1][0].plot(rr, [r*f[0]+f[1] for r in rr], 'r-.')

# plot 2d hist
aAllAni[1][0].hist2d(tauVal,mapVal,100)

# compute and plot linear regression - PCs
regr = linear_model.LinearRegression()
tauVal = [t for t, x in zip(zTauPcAniAll, mapCoPcAniAll) if np.isfinite(x)]
mapVal = [x for x in mapCoPcAniAll if np.isfinite(x)]
regr.fit(np.reshape(tauVal, [-1, 1]), mapVal)
x1 = np.min(tauVal)
x2 = np.max(tauVal)
aAllAni[0][1].plot([x1, x2], [x1* regr.coef_ + regr.intercept_, x2* regr.coef_ + regr.intercept_], 'r')

# same at different tau segments
rr = np.linspace(np.min(tauVal),np.max(tauVal),100)
for tauInterval in [[-2,2], [-2, 0], [0, 2]]:
    tSeg = [t for t in tauVal if tauInterval[0]<t<=tauInterval[1]]
    mSeg = [m for t,m in zip(tauVal,mapVal) if tauInterval[0]<t<=tauInterval[1]]
    f = np.polyfit(tSeg,mSeg,1)
    aAllAni[1][1].plot(rr, [r*f[0]+f[1] for r in rr], 'r-.')

# plot 2d hist
aAllAni[1][1].hist2d(tauVal,mapVal,100)

fMapSeg, aMapSeg = plt.subplots(2,2)
titles = [['fit z/Map nonPC','corr z/Map nonPC'],['fit z/Map PCs','corr z/Map PCs']]
for iIsPc in range(2):
    for iAni in range(nAni):
        aMapSeg[0][iIsPc].plot([f[0] for f in zMFitList[iIsPc][iAni]], '.')
        aMapSeg[1][iIsPc].plot([f for f in zMCorrList[iIsPc][iAni]], '.')
    aMapSeg[0][iIsPc].set_title(titles[iIsPc][0])
    aMapSeg[1][iIsPc].set_title(titles[iIsPc][1])

df=pd.DataFrame([[zAni[0] for zAni in zSameDiff] for zSameDiff in zMCorrList] , columns=animalList, index=['nPC', 'PC'])
df.to_csv(pathFigure + analStr + '_z-MapCorr.csv')


fTauSeg, aTauSeg = plt.subplots(2,2)
titles = [['fit z/Tau nonPC','corr z/Tau nonPC'],['fit z/Tau PCs','corr z/Tau PCs']]
for iIsPc in range(2):
    for iAni in range(nAni):
        aTauSeg[0][iIsPc].plot([f[0] for f in zTFitList[iIsPc][iAni]], '.')
        aTauSeg[1][iIsPc].plot([f for f in zTCorrList[iIsPc][iAni]], '.')
    aTauSeg[0][iIsPc].set_title(titles[iIsPc][0])
    aTauSeg[1][iIsPc].set_title(titles[iIsPc][1])

df=pd.DataFrame([[zAni[0] for zAni in zSameDiff] for zSameDiff in zTCorrList], columns=animalList, index=['nPC', 'PC'])
df.to_csv(pathFigure + analStr + '_z-TauCorr.csv')

aZM[0][0].set_title('non PCs')
aZM[0][1].set_title('Place Cells')
aZM[0][3].set_title('non PCs')
aZM[0][4].set_title('Place Cells')
aZM[0][2].set_title('avg+/-std corrCoef')
fZM.suptitle('mapCorr vs z for each animal')

aZT[0][0].set_title('non PCs')
aZT[0][1].set_title('Place Cells')
aZT[0][3].set_title('non PCs')
aZT[0][4].set_title('Place Cells')
aZT[0][2].set_title('avg+/-std corrCoef')
fZT.suptitle('tay vs z for each animal')

aAllAni[0][0].set_title('zT vs map - nPCs')
aAllAni[0][1].set_title('zT vs map - PCs')

aZTAllAni[0][0].set_title('zT vs rT - nPCs')
aZTAllAni[0][1].set_title('zT vs rT - PCs')

for ax in aAllAni[0]:
    ax.set_xlim([-0.55,1.05])
    ax.set_ylim([-0.3, 1.05])

for ax in aZTAllAni[0]:
    ax.set_xlim([-0.55, 1.05])
    ax.set_ylim([-0.3, 1.05])

for ax in aMapSeg.flat:
    ax.set_ylim([0, 0.6])
#
# for ax in a:
#     ax.set_xlim([-0.2, 1])
#     ax.set_ylim([-0.2, 1])
# aTauSeg[0][0].set_ylim([0, 1.5])
# aTauSeg[0][1].set_ylim([0, 1.5])
aTauSeg[1][0].set_ylim([-1, 1])
aTauSeg[1][1].set_ylim([-1, 1])

aZStabilityCorr.errorbar([0,1,2,3], np.mean(zStabilityCorr,0),np.std(zStabilityCorr,0), color='k')
aZStabilityCorr.set_ylim([-0.2, 0.5])
df = pd.DataFrame(zStabilityCorr, index=animalList, columns=['nPC-Same', 'PC-Same','nPC-Diff', 'PC-Diff'])
df.to_csv(pathFigure + analStr + '_PCoStability.csv')

if SaveFig:
    fMapSeg.savefig(pathFigure + analStr + 'MapCorr_allTau.eps', format='eps')
    fTauSeg.savefig(pathFigure+ analStr + 'TauCorr_allTau.eps', format='eps')
    fZStabilityCorr.savefig(pathFigure + analStr + 'StabilityCorr_allTau.eps', format='eps')
    fAllAni.savefig(pathFigure + analStr + 'MapCorr_scatter.eps', format='eps')
    fZTAllAni.savefig(pathFigure + analStr + 'TauCorr_scatter.eps', format='eps')
    fAllAni.savefig(pathFigure + analStr + 'MapCorr_scatter.png', format='png')
    fZTAllAni.savefig(pathFigure + analStr + 'TauCorr_scatter.png', format='png')
    fTauLogDist.savefig(pathFigure + 'zTauDistribution.eps', format='eps')

print('done done done')
plt.show(block=True)