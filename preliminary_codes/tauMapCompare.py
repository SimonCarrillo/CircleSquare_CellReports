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

path_all = '/Volumes/GoogleDrive/My Drive/Fenton_Lab/CellPress_CircleSquare/MNS_Analysis/Data/CircleSquare/'
pathFigure = '/Volumes/GoogleDrive/My Drive/Fenton_Lab/CellPress_CircleSquare/MNS_Analysis/Figures/CircleSquare/savedFigs_Simon/'
TimeBinSize_AVC = 1 * 1000
TimeBinSize_TAU = 1 * 1000
mSize = 1
BinNb = 10

animalList = ['M34']
nAni = len(animalList)

SaveFig = True
day_thresh = -1
thresh = 0.04*(TimeBinSize_TAU/1000)

tauAniAll = []
mapCoAniAll = []
tauPcAniAll = []
mapCoPcAniAll = []
fitList = [[], []]
corrList = [[], []]

fAllAni, aAllAni = plt.subplots(nrows=2, ncols=2, figsize=(9, 9))
fTauSeg, aTauSeg = plt.subplots(2,2)
fTauLogDist, aTauLogDist = plt.subplots(2,nAni)
fEachAni,aEachAni = plt.subplots(nAni, 5)

for iAni, animal in enumerate(animalList):

    print(animal)

    corr = []
    corrPc = []

    with open(path_all + 'analysisFiles_Simon/' + 'TempCorr_List_' +
              str(int(TimeBinSize_AVC)) + 'msAV_' + str(int(TimeBinSize_TAU)) + 'msKT'
              '_CircleSquare_' + animal +'.file', 'rb') as f:
        outputList = pickle.load(f)

    # outputList = [[day_real_list, day_count_list, sess_list, S_conv_avc_list, S_conv_tau_list, isSeparatedList,
    #                 tauPairsSingleList, tauVecSingleList],
    #                [Sector_list, Shock_list]]
    dayCountListTau = outputList[0][1]
    sessListTau = outputList[0][2]
    sConvTauList = outputList[0][4]
    isSeparatedList = outputList[0][5]
    tauPairsSingleList = outputList[0][6]
    tauVecSingleList = outputList[0][7]
    outputList = []

    loadpath = path_all + 'analysisFiles_Simon/' + 'MapCorr_' + str(int(BinNb)) + '_Bins_' + animal + '.file'
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

    tauAni = []
    mapCoAni = []
    tauPcAni = []
    mapCoPcAni = []

    for day in range(4,10):

        for session in ['RCT1', 'RCT2', 'CYL1', 'CYL2']:

            mapRecIdxList = np.where([d == day and session in iSess for d, iSess in zip(dayCountListMap, sessListMap)])[0]
            tauRecIdxList = np.where([d == day and session in iSess for d, iSess in zip(dayCountListTau, sessListTau)])[0]

            if len(mapRecIdxList) == 0 or len(tauRecIdxList) == 0:
                continue

            for tauRecIdx,mapRecIdx in zip(tauRecIdxList,mapRecIdxList):

                tauVec = tauVecSingleList[tauRecIdx]
                tauPairs = tauPairsSingleList[tauRecIdx]
                meanRate = [np.mean([s > 0 for s in c]) for c in sConvTauList[tauRecIdx]]
                isSep = isSeparatedList[tauRecIdx]

                PcIdx = PcIdxList[mapRecIdx]
                rateMaps = rateMapList[mapRecIdx]

                #create filters for tau
                tNan = [np.isnan(t[0]) for t in tauVec]
                rate = [meanRate[p[0]]>thresh and meanRate[p[1]]>thresh for p in tauPairs]
                antiCorr = [t[0]<2 for t in tauVec]
                tauFilter = [i and r and a and not n for i,n,r,a in zip(isSep, tNan, rate, antiCorr)]

                # filter tau
                tauVecFiltered = [t[0] for t,f in zip(tauVec, tauFilter) if f]
                tauPairsFiltered = [p for p,f in zip(tauPairs, tauFilter) if f]

                tauVecPc = [t for t, p in zip(tauVecFiltered, tauPairsFiltered) if (PcIdx[p[0]] or PcIdx[p[1]])]
                tauPairsPc = [p for p in tauPairsFiltered if (PcIdx[p[0]] or PcIdx[p[1]])]

                tauVecNonPc = [t for t, p in zip(tauVecFiltered, tauPairsFiltered) if (not PcIdx[p[0]] and not PcIdx[p[1]])]
                tauPairsNonPc = [p for p in tauPairsFiltered if (not PcIdx[p[0]] and not PcIdx[p[1]])]

                mask = ~np.isnan(rateMaps[0].flatten())
                mapCo = [np.corrcoef(rateMaps[p[0]].flatten()[mask], rateMaps[p[1]].flatten()[mask])[0, 1]
                         for p in tauPairsFiltered]
                mapCoPc = [np.corrcoef(rateMaps[p[0]].flatten()[mask], rateMaps[p[1]].flatten()[mask])[0, 1]
                         for p in tauPairsPc]
                mapCoNonPc = [np.corrcoef(rateMaps[p[0]].flatten()[mask], rateMaps[p[1]].flatten()[mask])[0, 1]
                         for p in tauPairsNonPc]

                # tauAni.extend(tauVecFiltered)
                tauAni.extend(tauVecNonPc)
                tauPcAni.extend(tauVecPc)
                # mapCoAni.extend(mapCo)
                mapCoAni.extend(mapCoNonPc)
                mapCoPcAni.extend(mapCoPc)


                aEachAni[iAni][0].scatter(tauVecNonPc, mapCoNonPc, s=1, alpha=0.3)
                if len(tauVecNonPc)>5 and len(mapCoNonPc)>5:
                    tauVal = [t for t, x in zip(tauVecNonPc,mapCoNonPc) if np.isfinite(x)]
                    mapVal = [x for x in mapCoNonPc if np.isfinite(x)]
                    regr = linear_model.LinearRegression()
                    regr.fit(np.reshape(tauVal, [-1, 1]),mapVal)
                    x1 = np.min(tauVal)
                    x2 = np.max(tauVal)
                    aEachAni[iAni][0].plot([x1, x2], [x1* regr.coef_ + regr.intercept_, x2* regr.coef_ + regr.intercept_], 'm--')
                    corr.append(np.corrcoef(tauVecNonPc,mapCoNonPc)[0,1])

                aEachAni[iAni][1].scatter(tauVecPc, mapCoPc, s=1, alpha=0.3)
                if len(tauVecPc)>5 and len(mapCoPc)>5:
                    tauVal = [t for t, x in zip(tauVecPc,mapCoPc) if np.isfinite(x)]
                    regr = linear_model.LinearRegression()
                    regr.fit(np.reshape(tauVal, [-1, 1]),[x for x in mapCoPc if np.isfinite(x)])
                    x1 = np.min(tauVal)
                    x2 = np.max(tauVal)
                    aEachAni[iAni][1].plot([x1, x2], [x1* regr.coef_ + regr.intercept_, x2* regr.coef_ + regr.intercept_], 'm--')
                    corrPc.append(np.corrcoef(tauVecPc, mapCoPc)[0,1])


    tauVecSingleList = dayCountListTau = sConvTauList = tauPairsSingleList = isSeparatedList = sessListTau = []

    aEachAni[iAni][2].bar([0, 1], [np.mean(corr), np.mean(corrPc)])
    aEachAni[iAni][2].errorbar([0, 1], [np.mean(corr), np.mean(corrPc)], [np.std(corr), np.std(corrPc)], ls = '', ecolor = 'k')

    tauAniAll.extend(tauAni)
    mapCoAniAll.extend(mapCoAni)
    tauPcAniAll.extend(tauPcAni)
    mapCoPcAniAll.extend(mapCoPcAni)

    tauVal = [t for t, x in zip(tauAni, mapCoAni) if np.isfinite(x)]
    mapVal = [x for x in mapCoAni if np.isfinite(x)]
    regr = linear_model.LinearRegression()
    regr.fit(np.reshape(tauVal, [-1, 1]), mapVal)
    x1 = np.min(tauVal)
    x2 = np.max(tauVal)

    fitListTmp = []
    corrListTmp = []
    for tauInterval in [[-2, 2], [-1, -0.05], [-0.05, 0], [0, 0.1], [0.1, 0.2], [0.2, 1]]:
        tSeg = [t for t in tauVal if tauInterval[0] < t <= tauInterval[1]]
        mSeg = [m for t, m in zip(tauVal, mapVal) if tauInterval[0] < t <= tauInterval[1]]
        fFit = np.polyfit(tSeg, mSeg, 1)
        fitListTmp.append(fFit)
        corrListTmp.append(np.corrcoef(tSeg, mSeg)[0,1])
    fitList[0].append(fitListTmp)
    corrList[0].append(corrListTmp)

    for ax in [[aEachAni[iAni][3], aEachAni[iAni][4]], aAllAni[0]]:
        ax[0].scatter(tauAni, mapCoAni, s=1, color='k', alpha=0.2)
        ax[0].plot([x1, x2], [x1* regr.coef_ + regr.intercept_, x2* regr.coef_ + regr.intercept_], 'm--')

    tauVal = [t for t, x in zip(tauPcAni, mapCoPcAni) if np.isfinite(x)]
    mapVal = [x for x in mapCoPcAni if np.isfinite(x)]
    regr = linear_model.LinearRegression()
    regr.fit(np.reshape(tauVal, [-1, 1]), mapVal)
    x1 = np.min(tauVal)
    x2 = np.max(tauVal)

    fitListTmp = []
    corrListTmp = []
    for tauInterval in [[-2, 2], [-1, -0.05], [-0.05, 0], [0, 0.1], [0.1, 0.2], [0.2, 1]]:
        tSeg = [t for t in tauVal if tauInterval[0] < t <= tauInterval[1]]
        mSeg = [m for t, m in zip(tauVal, mapVal) if tauInterval[0] < t <= tauInterval[1]]
        fFit = np.polyfit(tSeg, mSeg, 1)
        fitListTmp.append(fFit)
        corrListTmp.append(np.corrcoef(tSeg, mSeg)[0, 1])
    fitList[1].append(fitListTmp)
    corrList[1].append(corrListTmp)

    for ax in [[aEachAni[iAni][3], aEachAni[iAni][4]], aAllAni[0]]:
        ax[1].scatter(tauPcAni, mapCoPcAni, s=1, color='k', alpha=0.2)
        ax[1].plot([x1, x2], [x1* regr.coef_ + regr.intercept_, x2* regr.coef_ + regr.intercept_], 'm--')

    aTauLogDist[0][iAni].hist(tauAni, 100, log=True)
    aTauLogDist[1][iAni].hist(tauPcAni, 100, log=True)
    aTauLogDist[0][iAni].set_title(animal+'_nonPCs')
    aTauLogDist[1][iAni].set_title(animal+'_PCs')
    fTauLogDist.suptitle('tau distribution')

regr = linear_model.LinearRegression()
tauVal = [t for t, x in zip(tauAniAll, mapCoAniAll) if np.isfinite(x)]
mapVal = [x for x in mapCoAniAll if np.isfinite(x)]
regr.fit(np.reshape(tauVal, [-1, 1]),mapVal)
x1 = np.min(tauVal)
x2 = np.max(tauVal)
aAllAni[0][0].plot([x1, x2], [x1* regr.coef_ + regr.intercept_, x2* regr.coef_ + regr.intercept_], 'r')

# rr = np.linspace(np.min(tauVal),np.max(tauVal),100)
# f = np.polyfit(tauVal,mapVal,1)
# aAllAni[0][0].plot(rr, [r*f[0]+f[1] for r in rr], '--')
# f = np.polyfit(tauVal,mapVal,2)
# aAllAni[0][0].plot(rr, [(r**2)*f[0]+r*f[1]+f[2] for r in rr], '--')
# f = np.polyfit(tauVal,mapVal,3)
# aAllAni[0][0].plot(rr, [(r**3)*f[0]+(r**2)*f[1]+r*f[2]+f[3] for r in rr], '--')

# for tauInterval in [[-2, 2], [-1,-0.05], [-0.05,0], [0,0.1], [0.1,1]]:
#     tSeg = [t for t in tauVal if tauInterval[0]<t<=tauInterval[1]]
#     mSeg = [m for t,m in zip(tauVal,mapVal) if tauInterval[0]<t<=tauInterval[1]]
#     f = np.polyfit(tSeg,mSeg,1)
#     print(len(tSeg))
#     print(f)
#     aAllAni[1][0].plot(rr, [r*f[0]+f[1] for r in rr], 'r-.')

aAllAni[1][0].hist2d(tauVal,mapVal,100)

regr = linear_model.LinearRegression()
tauVal = [t for t, x in zip(tauPcAniAll, mapCoPcAniAll) if np.isfinite(x)]
mapVal = [x for x in mapCoPcAniAll if np.isfinite(x)]
regr.fit(np.reshape(tauVal, [-1, 1]), mapVal)
x1 = np.min(tauVal)
x2 = np.max(tauVal)
aAllAni[0][1].plot([x1, x2], [x1* regr.coef_ + regr.intercept_, x2* regr.coef_ + regr.intercept_], 'r')

aAllAni[0][0].set_ylim([-0.3, 1.05])
aAllAni[0][1].set_ylim([-0.3, 1.05])
aAllAni[0][0].set_xlim([-0.3, 1.05])
aAllAni[0][1].set_xlim([-0.3, 1.05])

# rr = np.linspace(np.min(tauVal),np.max(tauVal),100)
# f = np.polyfit(tauVal,mapVal,1)
# aAllAni[0][1].plot(rr, [r*f[0]+f[1] for r in rr])
# f = np.polyfit(tauVal,mapVal,2)
# aAllAni[0][1].plot(rr, [(r**2)*f[0]+r*f[1]+f[2] for r in rr])
# f = np.polyfit(tauVal,mapVal,3)
# aAllAni[0][1].plot(rr, [(r**3)*f[0]+(r**2)*f[1]+r*f[2]+f[3] for r in rr])
#
# for tauInterval in [[-2, 2], [-1,-0.05], [-0.05,0], [0,0.1], [0.1,0.2], [0.2,1]]:
#     tSeg = [t for t in tauVal if tauInterval[0]<t<=tauInterval[1]]
#     mSeg = [m for t,m in zip(tauVal,mapVal) if tauInterval[0]<t<=tauInterval[1]]
#     f = np.polyfit(tSeg,mSeg,1)
#     print(len(tSeg))
#     print(f)
#     aAllAni[1][1].plot(rr, [r*f[0]+f[1] for r in rr], 'r-.')

aAllAni[1][1].hist2d(tauVal,mapVal,100)

titleList = [['fit tau/map nPC', 'fit tau/map PCs'],['corr tau/map nPC', 'corr tau/map PCs']]
for iIsPc in range(2):
    for iAni in range(nAni):
        aTauSeg[0][iIsPc].plot([f[0] for f in fitList[iIsPc][iAni]], '.')
        aTauSeg[1][iIsPc].plot([f for f in corrList[iIsPc][iAni]], '.')
        aTauSeg[0][iIsPc].set_title(titleList[0][iIsPc])
        aTauSeg[1][iIsPc].set_title(titleList[1][iIsPc])

df=pd.DataFrame([[cAni[0] for cAni in cIsPC] for cIsPC in corrList] , columns=animalList, index=['nPC', 'PC'])
df.to_csv(pathFigure + 'tau-MapCorr.csv')

aEachAni[0][0].set_title('non PCs')
aEachAni[0][1].set_title('Place Cells')
aEachAni[0][3].set_title('non PCs (all Recs)')
aEachAni[0][4].set_title('Place Cells (all Recs)')
aEachAni[0][2].set_title('avg+/-std corrCoef')
fEachAni.suptitle('mapCorr vs Tau for each animal')

aAllAni[0][0].set_title('non PCs')
aAllAni[0][1].set_title('Place Cells')
fAllAni.suptitle('mapCorr vs Tau for all animals')

# for ax in aAllAni.flat:
#     ax.set_xlim([-0.2, 1])
#     ax.set_ylim([-0.2, 1])
#
# for ax in a:
#     ax.set_xlim([-0.2, 1])
#     ax.set_ylim([-0.2, 1])
aTauSeg[0][0].set_ylim([0, 1.5])
aTauSeg[0][1].set_ylim([0, 1.5])
aTauSeg[1][0].set_ylim([0, .6])
aTauSeg[1][1].set_ylim([0, .6])

if SaveFig:
    fTauSeg.savefig(pathFigure+'TauMapPerTauSegment.eps',format='eps')
    fAllAni.savefig(pathFigure+'TauMapScatter.png')
    fTauLogDist.savefig(pathFigure+'TauDistribution.eps', format='eps')

print('done done done')
plt.show(block=True)