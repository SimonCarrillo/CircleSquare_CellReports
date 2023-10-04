import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as scstat
import pickle
from sklearn import linear_model
import pandas as pd

expt = 'CircleSquare'

BinNb = 12  # for place maps
TimeBinSize = 1 * 1000
TimeBinSize_TAU = 1000
x_day = np.array([1, 2, 3, 8, 9, 10, 15, 16, 17])

threshMult = 2
doPerc = False

pathAll = '/Users/eliott/Dropbox/dDocuments/NYU/Fenton/MNS_ANALYSIS/Data/' + expt + '/'
pathFigure = '/Users/eliott/Dropbox/dDocuments/NYU/Fenton/MNS_ANALYSIS/Figures/' + expt + '/savedFigs/'

# TAU/WEIGHTS
if doPerc:
    pathLoad = pathAll + 'analysisFiles/' + 'CS_Decoder_SVM_' + str(
        int(TimeBinSize)) + 'msBins_' + str(
        int(TimeBinSize_TAU)) + 'msKT_perc_thresh' + str(int(threshMult))+ '.file'
else:
    pathLoad = pathAll + 'analysisFiles/' + 'CS_Decoder_SVM_' + str(
        int(TimeBinSize)) + 'msBins_' + str(
        int(TimeBinSize_TAU)) + 'msKT_seg_thresh' + str(int(threshMult))+ '.file'

print(pathLoad)

with open(pathLoad, 'rb') as f:
    output = pickle.load(f)
animalList, corrTauWeightAll, SegWeightRatioAll, tauBins = output

perfSum = []

# TAU/WEIGHTS
corrTauWeighPick = []
corrTauWeighPickAll = []
RatioSegWeightPick = []
RatioSegWeightPickAll = []
RatioSegWeightPickWeek = []

nAni = len(animalList)

fSeg, aSeg = plt.subplots(nAni, 2)

# TAU/WEIGHTS
nSeg = len(SegWeightRatioAll[0][7])
if doPerc:
    xPlot = [*range(nSeg)]
    idPick = [[0, 1], [9, 10], [18, 19]]
    pickList = ['<10%', '45-55%', '>90%']
else:
    xPlot = [tauBins[0] - (tauBins[1] - tauBins[0]) / 2]
    xPlot.extend([b + d / 2 for b, d in zip(tauBins, np.diff(tauBins))])
    xPlot.extend([tauBins[-1] - (tauBins[-2] - tauBins[-1]) / 2])

    if TimeBinSize_TAU == 300:
        idPick = [np.where([x < -0.025 for x in xPlot])[0][-1:],
                  np.where([- 0.012 < x < -0.005 for x in xPlot])[0][-1:],
                  np.where([x == 0 for x in xPlot])[0],
                  np.where([0.1 < x < 0.2 for x in xPlot])[0],
                  np.where([0.2 < x for x in xPlot])[0]]
        pickList = ['<-0.025', '-0.01', '0', '0.1-0.2', '>0.2']

    elif TimeBinSize_TAU == 1000:
        idPick = [np.where([x < -0.05 for x in xPlot])[0],
                  np.where([- 0.025 < x < -0.015 for x in xPlot])[0],
                  np.where([x == 0 for x in xPlot])[0],
                  np.where([0.1 < x < 0.2 for x in xPlot])[0],
                  np.where([0.2 < x for x in xPlot])[0]]
        pickList = ['<-0.05', '-0.02', '0', '0.1-0.2', '>0.2']
        print(xPlot)
        print(idPick)

for iAni, animal in enumerate(animalList):
    for iDay in range(9):
        if np.size(corrTauWeightAll[iAni][0][iDay])>1:
            aSeg[iAni][0].plot(xPlot, SegWeightRatioAll[iAni][iDay], color=[1 - iDay / 8, iDay / 8, iDay / 8])
            aSeg[iAni][1].plot(xPlot, corrTauWeightAll[iAni][0][iDay], 'r')
            aSeg[iAni][1].plot(xPlot, corrTauWeightAll[iAni][1][iDay], 'k')

    aSeg[iAni][0].plot(xPlot, [1 for x in range(nSeg)], 'y--')
    aSeg[iAni][1].plot(xPlot, [0 for x in range(nSeg)], 'y--')

    # remove nans
    cEnv1 = [-c if (np.size(c) > 1 and np.sum(np.isfinite(c))>1) else [np.nan for _ in range(9)] for c in corrTauWeightAll[iAni][0]]
    cEnv2 = [c if (np.size(c) > 1 and np.sum(np.isfinite(c))>1) else [np.nan for _ in range(9)] for c in corrTauWeightAll[iAni][1]]
    ratio = [np.array(c) if (np.size(c) > 1 and np.sum(np.isfinite(c)) > 1) else [] for c in SegWeightRatioAll[iAni]]

    #take average of first 2, last 2, middle 2 and then average over days
    cPick = [np.nanmean([[np.nanmean([cc[i] for i in iId]) for iId in idPick] for cc in cEnv1], 0),
             np.nanmean([[np.nanmean([cc[i] for i in iId]) for iId in idPick] for cc in cEnv2], 0)]
    cPickAll = [[[np.nanmean([cc[i] for i in iId]) for iId in idPick] for cc in cEnv1],
             [[np.nanmean([cc[i] for i in iId]) for iId in idPick] for cc in cEnv2]]

    rPickAll = [[np.nanmean([cc[i] for i in iId]) for iId in idPick] if len(cc) > 0 else [np.nan for x in range(nSeg)] for cc in ratio]
    rPick = np.nanmean([[np.nanmean([cc[i] for i in iId]) for iId in idPick] for cc in ratio if len(cc) > 0],0)
    rPickWeek = [np.nanmean([[np.nanmean([cc[i] for i in iId]) for iId in idPick] for cc in ratioWeek if len(cc)>0],0)
                 for ratioWeek in [[ratio[x] for x in days] for days in [[0,1,2],[3,4,5],[6,7,8]]]]

    # save for plotting and average calculation
    corrTauWeighPick.append(cPick)
    corrTauWeighPickAll.append(cPickAll)
    RatioSegWeightPick.append(rPick)
    RatioSegWeightPickAll.append(rPickAll)
    RatioSegWeightPickWeek.append(rPickWeek)


fTauAverage = plt.figure(figsize=[12.8,9.6])
ax1 = fTauAverage.add_subplot(413)
ax2 = fTauAverage.add_subplot(414)
ax3 = fTauAverage.add_subplot(241)
ax4 = fTauAverage.add_subplot(242)
axBox1 = fTauAverage.add_subplot(243)
axBox2 = fTauAverage.add_subplot(244)

nPick = len(RatioSegWeightPick[0])

ratioCorrBoxPlot = [[rDay[iSeg] for rAni in RatioSegWeightPickAll for rDay in rAni if np.isfinite(rDay[iSeg])] for iSeg in range(5)]
tagRatioCorr = [[[animalList[iAni], iDay, pickList[iSeg]] for iAni, rAni in enumerate(RatioSegWeightPickAll) for iDay,rDay in enumerate(rAni) if np.isfinite(rDay[iSeg])] for iSeg in range(5)]
axBox1.boxplot(ratioCorrBoxPlot)

envList = ['RCT', 'CYL'] #got from SvmDecodeCs_tauSegments.py (line 349, 362)
corrWeightBoxPlot = [[cDay[iSeg] for cAni in corrTauWeighPickAll for iEnv in range(2) for cDay in cAni[iEnv] if np.isfinite(cDay[iSeg])] for iSeg in range(nPick)]
tagCorrWeight = [[[animalList[iAni], envList[iEnv], iDay, pickList[iSeg]] for iAni, cAni in enumerate(corrTauWeighPickAll) for iEnv in range(2) for iDay, cDay in enumerate(cAni[iEnv]) if np.isfinite(cDay[iSeg])] for iSeg in range(nPick)]
axBox2.boxplot(corrWeightBoxPlot, positions=[2*x-0.3 for x in range(nPick)])
# axBox2.boxplot([[cDay[iSeg] for cAni in corrTauWeighPickAll for cDay in cAni[0]] for iSeg in range(nPick)],
#                positions=[2*x-0.3 for x in range(nPick)])
# axBox2.boxplot([[cDay[iSeg] for cAni in corrTauWeighPickAll for cDay in cAni[1]] for iSeg in range(nPick)],
#                positions=[2*x+0.3 for x in range(nPick)])
axBox2.set_xlim([-1,9])

for iDay in range(9):
    if len([s[0] for s in SegWeightRatioAll if np.size(s[0]) > 1])>0:
        ax1.plot(xPlot, np.nanmean([s[iDay] for s in SegWeightRatioAll if np.size(s[iDay]) > 1], 0),
                            '.-', color=[1 - iDay / 8, iDay / 8, iDay / 8])

for SegWeightRatioAni in SegWeightRatioAll:
    ax2.plot(xPlot, np.nanmean([s for s in SegWeightRatioAni if np.size(s) > 1], 0), '.-')

ax1.set_title('average per day')
ax2.set_title('average per ani')

ax3.bar([*range(nPick)], np.mean(RatioSegWeightPick,0))
ax3.errorbar([*range(nPick)], np.mean(RatioSegWeightPick,0), np.std(RatioSegWeightPick,0)/np.sqrt(nAni), color='k', ls = '')

ax4.bar([x-0.15 for x in range(nPick)], np.mean(corrTauWeighPick,0)[0], 0.2)
ax4.bar([x+0.15 for x in range(nPick)], np.mean(corrTauWeighPick,0)[1], 0.2)
ax4.errorbar([x-0.15 for x in range(nPick)], np.mean(corrTauWeighPick,0)[0],
                     np.std(corrTauWeighPick,0)[0]/np.sqrt(6), color='k', ls = '')
ax4.errorbar([x+0.15 for x in range(nPick)], np.mean(corrTauWeighPick,0)[1],
                     np.std(corrTauWeighPick,0)[1]/np.sqrt(6), color='k', ls = '')

# plot
for rPickAni in RatioSegWeightPickWeek:
    ax3.plot([x - 0.15 for x in range(nPick)], rPickAni[0], 'k.')
    ax3.plot([x for x in range(nPick)], rPickAni[1], 'k.')
    ax3.plot([x + 0.15 for x in range(nPick)], rPickAni[2], 'k.')

for rPickAni in RatioSegWeightPick:
    # ax3.plot(rPickAni, 'k.')
    ax3.plot(rPickAni, '--', color=[0.7, 0.7, 0.7])
    ax3.set_xticks([])
    ax3.set_xticks([*range(nPick)])
    ax3.set_xticklabels(pickList)


for cPickAni in corrTauWeighPick:
    ax4.plot([x-0.15 for x in range(nPick)], cPickAni[0], 'r.')
    ax4.plot([x+0.15 for x in range(nPick)], cPickAni[1], 'k.')
    ax4.set_xticks([])
    ax4.set_xticks([*range(nPick)])
    ax4.set_xticklabels(pickList)

#exporting data
df = pd.DataFrame()

dataStr = ['ratioCorr']  # create index string
df = df.append(pd.DataFrame(data=np.concatenate(ratioCorrBoxPlot), columns=dataStr).transpose())
dataStr = [s + '_tag_ratioCorr' + s for s in ['ani', 'day', 'pick']]  # create index string
df = df.append(pd.DataFrame(data=np.concatenate(tagRatioCorr), columns=dataStr).transpose())

dataStr = ['CorrWeight']  # create index string
df = df.append(pd.DataFrame(data=np.concatenate(corrWeightBoxPlot), columns=dataStr).transpose())
dataStr = [s + '_tag_CorrWeight' + s for s in ['ani', 'env', 'day', 'pick']]  # create index string
df = df.append(pd.DataFrame(data=np.concatenate(tagCorrWeight), columns=dataStr).transpose())

df.to_csv(pathFigure + 'dataSvmTauSegments.csv')

fTauAverage.savefig(pathFigure + 'TauSegmentAverage.eps', format='eps')
fSeg.savefig(pathFigure + 'TauSegmentAllAni.eps', format='eps')

print('done done done')
plt.show(block=True)