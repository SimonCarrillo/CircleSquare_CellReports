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

TimeBinSize_AVC = 10 * 1000
TimeBinSize_TAU = 1000

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

animalList = ['M19', 'M29', 'M35', 'M20', 'M34', 'M39']

strTimeScale = str(int(TimeBinSize_AVC)) + 'msAV_' + str(int(TimeBinSize_TAU)) + 'msKT'

nAni = len(animalList)

doPTI = True

if doPTI:
    pathFigure = '/Users/eliott/Dropbox/dDocuments/NYU/Fenton/MNS_ANALYSIS/Figures/' + expt + '/savedFigs/z/' # for z
else:
    pathFigure = '/Users/eliott/Dropbox/dDocuments/NYU/Fenton/MNS_ANALYSIS/Figures/' + expt + '/savedFigs/'

print(doPTI)
print(pathFigure)

netCoMean = []
netCoNorm = []
netCoZ = []
fMean,aMean = plt.subplots(2,nAni)

for iAni, animal in enumerate(animalList):

    print(animal)

    if doPTI:
        with open(pathAll + 'analysisFiles/' + 'Zlinear_List_' +
                  str(int(TimeBinSize_AVC)) + 'msAV_' + str(int(TimeBinSize_TAU)) + 'msKT'
                  '_CircleSquare_' + animal +'.file', 'rb') as f:
            outputComb = pickle.load(f)
    else:
        with open(pathAll + 'analysisFiles/' + 'TempCorr_List_' +
                  str(int(TimeBinSize_AVC)) + 'msAV_' + str(int(TimeBinSize_TAU)) + 'msKT'
                  '_CircleSquare_' + animal +'.file', 'rb') as f:
            outputComb = pickle.load(f)

    if doPTI:
        sConvTauList = outputComb[0][4][1]
        tauVecSingleList = outputComb[0][7][0]
        # [dayRealList, _ , _, sConvAvcList, [zConvTauList, sConvtauList], isSeparatedList,
        #   tauPairsSingleList, [zTauVecSingleList, tauVecSingleList], trainingList] = outputComb[0]
    else:
        sConvTauList = outputComb[0][4]
        tauVecSingleList = outputComb[0][7]
        # [dayRealList, _ , _, sConvAvcList, sConvtauList, isSeparatedList,
        #  tauPairsSingleList, tauVecSingleList, trainingList] = outputComb[0]
    dayCountList = outputComb[0][1]
    sessList = outputComb[0][2]
    isSeparatedList = outputComb[0][5]
    tauPairsSingleList = outputComb[0][6]

    outputComb = []

    daySave = []

    f, a = plt.subplots(9, 4)
    fDist, aDist = plt.subplots(3, 3)
    fRand, aRand = plt.subplots(3, 3)

    netCoAniMean = [[] for _ in range(9)]
    netCoAniNorm = [[] for _ in range(9)]

    for iDay, day in enumerate(range(1,10)):

        if day not in dayCountList: continue

        daySave.append(day)

        # find index of non-HMC recording from the 'day' defined above
        recIdx = np.where([d == day and 'HMC' not in s for d,s in zip(dayCountList, sessList)])[0]

        for iRec, rIdx in enumerate(recIdx):

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
            # # # filter out cells that are rarely active
            # meanRate = [np.mean([s > 0 for s in c]) for c in sConvTau]
            # rate = [meanRate[p[0]] > thresh and meanRate[p[1]] > thresh for p in tauPairs]

            # create filter for spatial separation and significance and non nan correlation
            tauFilter = [sep and not iNan for sep, iNan in zip(isSepDay, isAnyNan)]
            meanRate = rate = isAnyNan = isSeparated = []

            # filter tau values and pair numbers
            tauVec = [t for t, f in zip(tauVec, tauFilter) if f]
            idxTauPairs = [idx for idx, f in zip(tauPairs, tauFilter) if f]

            netCo = np.sum([[s1 * s2 * t for s1, s2 in zip(sConvTau[idP[0]], sConvTau[idP[1]])]
                            for idP, t in zip(idxTauPairs, tauVec)], 0)

            netCoRand = []
            for _ in range(100):
                tRand = tauVec.copy()
                np.random.shuffle(tRand)
                netCoRand.extend(np.sum([[s1 * s2 * t for s1, s2 in zip(sConvTau[idP[0]], sConvTau[idP[1]])]
                                     for idP, t in zip(idxTauPairs, tRand)], 0))

            netCoAniMean[iDay].append(np.mean(netCo))
            netCoAniNorm[iDay].append(np.mean(netCo) / np.std(netCoRand))

            a[iDay][iRec].plot(netCo)
            aDist.flat[iDay].hist(netCo, 50, density=True, histtype='step')
            aRand.flat[iDay].hist(netCoRand, 50, density=True, histtype='step')
            aDist.flat[iDay].set_title('PTI')
            # print(iDay)
            # print(netCoAniMean)
            # print(netCoAniNorm)

    # netCoAll.append(netCoAni)
    aMean[0][iAni].errorbar([*range(9)], [np.mean(n) for n in netCoAniMean], [np.std(n) for n in netCoAniMean])
    aMean[1][iAni].errorbar([*range(9)], [np.mean(n) for n in netCoAniNorm], [np.std(n) for n in netCoAniNorm])
    netCoMean.append(netCoAniMean)

    m = np.mean([n_ for nD in netCoAniMean for n_ in nD])
    s = np.std([n_ for nD in netCoAniMean for n_ in nD])
    netCoZ.append([[(n_-m)/s for n_ in nD] for nD in netCoAniMean])

    netCoNorm.append(netCoAniNorm)

print(netCoMean)
netCoMeanBarPlot = [[n_ for n in netCoMean for nDay in n[iWeek:iWeek+3] for n_ in nDay] for iWeek in [0,3,6]]
netCoNormBarPlot = [[n_ for n in netCoNorm for nDay in n[iWeek:iWeek+3] for n_ in nDay] for iWeek in [0,3,6]]
netCoZBarPlot = [[n_ for n in netCoZ for nDay in n[iWeek:iWeek+3] for n_ in nDay] for iWeek in [0,3,6]]
tagNetCo = [[[animalList[iAni],iWeek/3, iDay] for iAni, n in enumerate(netCoZ) for iDay, nDay in enumerate(n[iWeek:iWeek+3]) for n_ in nDay] for iWeek in [0,3,6]]


f,a = plt.subplots(1,3)
a[0].boxplot(netCoMeanBarPlot)
a[1].boxplot(netCoZBarPlot)
a[2].boxplot(netCoNormBarPlot)

a[0].set_title('mean')
a[1].set_title('mean - zScored')
a[2].set_title('Normed mean')

if doPTI:
    f.savefig(pathFigure + 'PTI_netCo.eps', format='eps')
else:
    f.savefig(pathFigure+'netCo.eps', format='eps')

# exporting data
df = pd.DataFrame()
# add data to the panda
dataStr = ['mean'] #create index string
df = df.append(pd.DataFrame(data=np.concatenate(netCoMeanBarPlot),columns=dataStr).transpose())

dataStr = ['mean-zScored'] #create index string
df = df.append(pd.DataFrame(data=np.concatenate(netCoZBarPlot),columns=dataStr).transpose())

dataStr = ['normedMean'] #create index string
df = df.append(pd.DataFrame(data=np.concatenate(netCoNormBarPlot),columns=dataStr).transpose())

dataStr = [s + '_tag' for s in ['ani', 'week', 'day']] #create index string
df = df.append(pd.DataFrame(data=np.concatenate(tagNetCo),columns=dataStr).transpose())

if doPTI:
    df.to_csv(pathFigure + 'PTI_dataNetworkCoherence.csv')
else:
    df.to_csv(pathFigure + 'dataNetworkCoherence.csv')

plt.show(block=True)