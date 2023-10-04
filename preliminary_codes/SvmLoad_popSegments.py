import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as scstat
import pickle
from sklearn import linear_model
import pandas as pd

expt = 'CircleSquare'

TimeBinSize = 1 * 1000
x_day = np.array([1, 2, 3, 8, 9, 10, 15, 16, 17])

dayPick = 8
doPerc = False
nSeg = 5

pathAll = '/Volumes/GoogleDrive/My Drive/Fenton_Lab/CellPress_CircleSquare/MNS_ANALYSIS/Data/' + expt + '/'
pathFigures = '/Volumes/GoogleDrive/My Drive/Fenton_Lab/CellPress_CircleSquare/MNS_ANALYSIS/Figures/' + expt + '/savedFigs_Simon/'

path_load = pathAll + 'analysisFiles/' + 'CS_Decoder_SVM_' + str(
    int(TimeBinSize)) + 'msBins_' + 'popSegmentsTESTany_'+ str(int(nSeg)) +'seg' '.file'

with open(path_load, 'rb') as f:
    output = pickle.load(f)
performance, weights, perfTopCells, ratio, coefStd, animalList = output
#weights [animals, noCV/CV, day]

ratioPercPc = ratio[0]
ratioPercHighRate = ratio[1]
ratioPercLowRate = ratio[2]

perfSum = []
perfSumAllVals = []

animalList = animalList

nAni = len(animalList)

aWeightCpm = plt.subplots(nAni+1, 9)[1]
aDist = plt.subplots(nAni, 3)[1]
aDotDensity = plt.subplots(1, 1)[1]
fSeg, aSeg = plt.subplots(nAni, 2)

weightCorr = [[[] for iDay in range(9)] for iAni in animalList]
weightsAll = [[[],[]] for iDay in range(9)]

print(TimeBinSize)

# to export data
df = pd.DataFrame()

for iAni, animal in enumerate(animalList):

    perfTopCellsReshaped = [[pp[iSeg] if len(pp) > iSeg else np.nan for pp in perfTopCells[iAni]] for iSeg in range(nSeg)]
    for iLine, perf in enumerate(perfTopCellsReshaped):
        if iLine < 5:
            lineColor = [0.5, iLine / 5, 1 - iLine / 5]
        else:
            lineColor = [0.5 + (iLine - 5) / 10, 1, 0]
        aDist[iAni][0].plot(x_day, perf, '.-', color=lineColor)
        aDist[iAni][0].set_ylim([0.45, 1])
        aDotDensity.plot(x_day+iLine/15, perf, '.', color=lineColor)

    perfSumAni = [np.reshape(perf, [-1, 3])/np.nanmean(perf[0:3]) for perf in perfTopCellsReshaped]
    # perfSumAni = [np.array([perf[0:3] - np.nanmean(perf[0:3]),
    #                         perf[3:6] - np.nanmean(perf[0:3]),
    #                         perf[6:9] - np.nanmean(perf[3:6])]) for perf in perfTopCellsReshaped]
    for iLine, perf in enumerate(perfSumAni):
        if iLine < 5:
            lineColor = [0.5, iLine / 5, 1 - iLine / 5]
        else:
            lineColor = [0.5 + (iLine - 5) / 10, 1, 0]
        aDist[iAni][1].errorbar([1, 2, 3], np.nanmean(perf, 1), np.nanstd(perf, 1), color=lineColor)
        # a2[iAni][1].set_ylim([0, 1.4])
        aDist[iAni][1].set_xlim([0.8, 3.2])

    for i_plot, iCoef in enumerate([np.abs(w) for w in weights[iAni][1]]):
        if len(iCoef)>1:
            aDist[iAni][2].hist(iCoef, 100, histtype='step', cumulative=True,
                        density=True, color=[1 - i_plot / 8, i_plot / 8, i_plot / 8])


    perfSum.append(np.nanmean(perfSumAni, 2))
    perfSumAllVals.append(perfSumAni)

    with open(pathAll + 'analysisFiles/' + 'AssignmentsByPairs_' + animal + '.file', 'rb') as f:
        [assignments, [assignmentsByPairs, idPairs]] = pickle.load(f)

    for i_day in range(9):
        dayPair = (i_day, dayPick)

        if np.diff(dayPair)[0] == 0:
            findSelf = np.where([dayPair[0] == iId[0] for iId in idPairs])[0]
            if len(findSelf)>0:
                nMatched = len(assignmentsByPairs[findSelf[0]][0])+len(assignmentsByPairs[findSelf[0]][2])
            else:
                findSelf = np.where([dayPair[0] == iId[1] for iId in idPairs])[0]
                nMatched = len(assignmentsByPairs[findSelf[0]][1]) + len(assignmentsByPairs[findSelf[0]][3])
            CellList = [[*range(nMatched)] for x in range(2)]
            f1Score = 1
        else:
            # will raise error if pair is not found
            if dayPair[0] < dayPair[1]:
                findAssign = np.where([dayPair == iPair for iPair in idPairs])[0][0]
                assign = assignmentsByPairs[findAssign]
                CellList = [assign[0], assign[1]]
            else:
                findAssign = np.where([dayPair[::-1] == iPair for iPair in idPairs])[0][0]
                assign = assignmentsByPairs[findAssign]
                CellList = [assign[1], assign[0]]
            nMatched = len(assign[0])
            f1Score = assign[4]['f1_score']

        # check for misalignment
        if f1Score < 0.3: continue

        if i_day == dayPick: continue

        if np.all([len(weights[iAni][1][dayPair[i]])>1 for i in [0, 1]]):
            W = [[weights[iAni][1][dayPair[i]][n] for n in CellList[i]] for i in [0, 1]]
            # W = [[weights[iAni][0][day_pair[i]][n] for n in CellList[i]] for i in [0, 1]]

            aWeightCpm[iAni][i_day].plot(W[0], W[1], '.')
            aWeightCpm[nAni][i_day].scatter(W[0], W[1], alpha=0.3, s=1)
            aWeightCpm[iAni][i_day].plot([np.min(W), np.max(W)], [np.min(W), np.max(W)], 'y--')

            weightsAll[i_day][0].extend(W[0])
            weightsAll[i_day][1].extend(W[1])
            weightCorr[iAni][i_day].append(np.corrcoef((W[0], W[1]))[0,1])

            regr = linear_model.LinearRegression()
            regr.fit(np.reshape(W[0], [-1, 1]), W[1])
            pred_plot = regr.predict([[np.min(W[0])], [np.max(W[0])]])

            c, p = scstat.pearsonr(W[0], W[1])
            aWeightCpm[iAni][i_day].plot([[np.min(W[0])], [np.max(W[0])]], pred_plot, 'r')
            aWeightCpm[iAni][i_day].text(np.min(W), 0.1*np.min(W)+ 0.9*np.max(W),
                                  'c=' + '{:.2f}'.format(c) + ' p=' + '{:.3f}'.format(p))

            aWeightCpm[iAni][i_day].set_aspect('equal', 'box')
        else:
            aWeightCpm[iAni][i_day].set_aspect('equal', 'box')

fPerfSeg, aPerfSeg = plt.subplots(2, 3)
zipped = zip(np.nanmean(perfSum,0), np.nanstd(perfSum,0), np.sum(~np.isnan(perfSum),0))
for iLine, (m, s, n) in enumerate(zipped):
    if iLine<5:
        lineColor = [0.5, iLine / 5, 1 - iLine / 5]
    else:
        lineColor = [0.5+(iLine-5) / 10, 1, 0]
    aPerfSeg[0][0].errorbar([1,2,3], m, s/np.sqrt(n[0]), color=lineColor)

for iWeek in range(3):
    aPerfSeg[0][1].errorbar([*range(nSeg)], np.nanmean(perfSum,0)[:,iWeek], np.nanstd(perfSum,0)[:,iWeek]/np.sqrt(nAni))

aPerfSeg[0][2].plot([np.mean([pAni[iDay] for pAni in perfTopCells if len(pAni[iDay])>1],0) for iDay in range(9)])
aPerfSeg[0][2].errorbar([*range(9)], np.nanmean(performance, 0)[1], np.nanstd(performance, 0)[1],
                            color='k', fmt='-.')

for iSeg in range(nSeg):
    perfSegBoxPlot = [[pSeg[iSeg] for pAni in perfTopCells for pSeg in pAni[0 + iWk:3 + iWk] if np.isfinite(pSeg[0])]
                      for iWk in [0, 3, 6]]
    tagPerfSegBoxPlot = [[[animalList[iAni], iWk, iDay] for iAni, pAni in enumerate(perfTopCells) for iDay, pSeg in enumerate(pAni[0 + iWk:3 + iWk]) if np.isfinite(pSeg[0])]
                      for iWk in [0, 3, 6]] # to record animal, day and week when exporting data
    aPerfSeg[1][2].boxplot(perfSegBoxPlot, positions=[2 + 0.8 * iSeg, 9 + 0.8 * iSeg, 16 + 0.8 * iSeg])
    aPerfSeg[1][1].boxplot([[pSeg[iSeg] for pAni in perfTopCells for pSeg in pAni[0 + iWk:3 + iWk] if np.isfinite(pSeg[0])]
                            for iWk in [0, 3, 6]], positions=[2 + 4 * iSeg, 3 + 4 * iSeg, 4 + 4 * iSeg])
    aPerfSeg[1][0].boxplot([[p for pAni in perfSumAllVals for p in pAni[iSeg][iWk] if np.isfinite(p)] for iWk in range(3)],
                           positions = [2 + 0.8 * iSeg, 9 + 0.8 * iSeg, 16 + 0.8 * iSeg])

    dataStr = ['perfperSeg_seg' + str(iSeg) + '_week' + str(iWk) for iWk in range(3)]  # create index string
    df = df.append(pd.DataFrame(data=perfSegBoxPlot, index=dataStr))
    dataStr = ['TagPerfperSeg_seg' + str(iSeg) + '_week' + str(iWk) for iWk in range(3)]  # create index string
    df = df.append(pd.DataFrame(data=tagPerfSegBoxPlot, index=dataStr))

aPerfSeg[1][2].set_xlim([0,24])
aPerfSeg[1][1].set_xlim([0,42])
aPerfSeg[1][0].set_xlim([0,20])

# PLOT: DECODING PERF, COEF STD AND %PC
fSummary, aSummary = plt.subplots(3, 4)
Title_perf = ['current decode current _best',
              'current decode current _cv']

# plot decoding perf as line
for i_plot in range(2):
    aSummary.flat[i_plot].plot(x_day, np.ones_like(x_day), 'y--')
    aSummary.flat[i_plot].plot(x_day, 0.5*np.ones_like(x_day), 'y--')
    aSummary.flat[i_plot].errorbar(x_day, np.nanmean(performance, 0)[i_plot], np.nanstd(performance, 0)[i_plot],
                            color='k', fmt='-.')

    for i_ani_perf in performance:
        aSummary.flat[i_plot].plot(x_day, i_ani_perf[i_plot], '.')

    aSummary.flat[i_plot].set_ylim([0, 1.2])
    aSummary.flat[i_plot].set_title(Title_perf[i_plot])
    aSummary.flat[i_plot].set_ylim([0.4, 1.1])

# plot decoding perf as Bars
for i_plot in range(2):
    aSummary.flat[i_plot + 4].plot(x_day, np.ones_like(x_day), 'y--')
    aSummary.flat[i_plot + 4].plot(x_day, 0.5 * np.ones_like(x_day), 'y--')

    for i_ani_perf in performance:
        aSummary.flat[i_plot+4].plot(x_day, i_ani_perf[i_plot], 'k.')

    # average each week for each animal
    perfWeekAverg = [np.nanmean(np.reshape(p[i_plot], [3, 3]), 1) for p in performance]

    # make bar and errorbars
    aSummary.flat[i_plot+4].bar([2, 9, 16], np.mean(perfWeekAverg, 0))
    aSummary.flat[i_plot+4].errorbar([2, 9, 16], np.mean(perfWeekAverg, 0),
                              np.std(perfWeekAverg, 0) / np.sqrt(nAni),
                              ls='', color='k')

    aSummary.flat[i_plot+4].set_ylim([0, 1.2])
    aSummary.flat[i_plot+4].set_title(Title_perf[i_plot])
    aSummary.flat[i_plot+4].set_ylim([0.4, 1.1])

# plot dcoding performance
perfBoxPlot = [[p for pAni in performance for p in pAni[1][0 + iWk:3 + iWk] if np.isfinite(p)] for iWk in [0, 3, 6]]
tagPerfBoxPlot = [[[animalList[iAni], iWk, iDay] for iAni, pAni in enumerate(performance)
                   for iDay, p in enumerate(pAni[1][0 + iWk:3 + iWk]) if np.isfinite(p)] for iWk in [0, 3, 6]] # to record animal, day and week when exporting data
aSummary.flat[7].boxplot(perfBoxPlot, positions = [2, 9, 16], notch=False)
aSummary.flat[7].set_title(Title_perf[1])
aSummary.flat[7].plot(x_day, np.ones_like(x_day), 'y--')
aSummary.flat[7].plot(x_day, 0.5 * np.ones_like(x_day), 'y--')

#export data
dataStr = ['decodePerf_week'+ str(iWk) for iWk in range(3)] #create index string
df = df.append(pd.DataFrame(data=perfBoxPlot,index=dataStr))
dataStr = ['decodePerfTag_week'+ str(iWk) for iWk in range(3)] #create index string
df = df.append(pd.DataFrame(data=tagPerfBoxPlot,index=dataStr))

aSummary.flat[2].plot(x_day, np.ones_like(x_day), 'y--')  # plot dashed line at 1
# plot line average
# aSummary.flat[2].errorbar(x_day, np.nanmean(ratioPercPc, 0), np.nanstd(ratioPercPc, 0)/np.sqrt(len(animalList)))

# plot a point for each day for each animal
for r in ratioPercPc:
    aSummary.flat[2].plot(x_day, r, 'k.')

# average each week for each animal
ratioPercPcWeekAverg = [np.nanmean(np.reshape(r, [3,3]),1) for r in ratioPercPc]

# make bar and errorbars
aSummary.flat[2].bar([2, 9, 16], np.mean(ratioPercPcWeekAverg, 0))
aSummary.flat[2].errorbar([2, 9, 16], np.mean(ratioPercPcWeekAverg, 0), np.std(ratioPercPcWeekAverg, 0)/np.sqrt(nAni),
             ls='', color='k')

ratioPercPcBoxPlot = [[r for rAni in ratioPercPc for r in rAni[0+iWk:3+iWk] if np.isfinite(r)] for iWk in [0,3,6]]
ratioPercHighRateBoxPlot = [[r for rAni in ratioPercHighRate for r in rAni[0+iWk:3+iWk] if np.isfinite(r)] for iWk in [0,3,6]]
ratioPercLowRateBoxPlot = [[r for rAni in ratioPercLowRate for r in rAni[0+iWk:3+iWk] if np.isfinite(r)] for iWk in [0,3,6]]
tagBoxPlot = [[[animalList[iAni], iWk, iDay] for iAni, rAni in enumerate(ratioPercPc)
               for iDay, r in enumerate(rAni[0+iWk:3+iWk]) if np.isfinite(r)] for iWk in [0,3,6]] # to record animal, day and week when exporting data

aSummary.flat[6].boxplot(ratioPercPcBoxPlot, positions=[2, 9, 16], notch = False)
aSummary.flat[8].boxplot(ratioPercHighRateBoxPlot, positions=[2, 9, 16], notch = False)
aSummary.flat[9].boxplot(ratioPercLowRateBoxPlot, positions=[2, 9, 16], notch = False)

aSummary.flat[2].set_title('ratio %PC in top 20 vs overall')
aSummary.flat[6].set_title('ratio %PC in top 20 vs overall')
aSummary.flat[8].set_title('ratio %HR in top 20 vs overall')
aSummary.flat[9].set_title('ratio %LR in top 20 vs overall')

#export data
dataStr = ['ratioPercPc_week'+ str(iWk) for iWk in range(3)] #create index string
df = df.append(pd.DataFrame(data=ratioPercPcBoxPlot,index=dataStr))

dataStr = ['ratioHighRate_week'+ str(iWk) for iWk in range(3)] #create index string
df = df.append(pd.DataFrame(data=ratioPercHighRateBoxPlot,index=dataStr))

dataStr = ['ratioLowRate_week'+ str(iWk) for iWk in range(3)] #create index string
df = df.append(pd.DataFrame(data=ratioPercLowRateBoxPlot,index=dataStr))

dataStr = ['tagRatio_week'+ str(iWk) for iWk in range(3)] #create index string
df = df.append(pd.DataFrame(data=tagBoxPlot,index=dataStr))


coefStdNormed = [[i_s/np.nanmean(s[0:3]) for i_s in s] for s in coefStd]
aSummary.flat[3].plot(x_day, np.ones_like(x_day), 'y--')
aSummary.flat[3].errorbar(x_day, np.nanmean(coefStdNormed, 0), np.nanstd(coefStdNormed, 0)/np.sqrt(nAni))
for c in coefStdNormed:
    aSummary.flat[3].plot(x_day, c, 'k.')
aSummary.flat[3].set_title('coef dist std - normed to 1st week')

plt.figure()
weightCorrBoxPlot = [[w for wAni in weightCorr for wDay in wAni[0+iWk:3+iWk] for w in wDay] for iWk in [0,3,6]]
tagWCorr = [[[animalList[iAni], iWk, iDay] for iAni, wAni in enumerate(weightCorr) for iDay, wDay in enumerate(wAni[0+iWk:3+iWk]) for w in wDay] for iWk in [0,3,6]]
plt.boxplot(weightCorrBoxPlot)
plt.savefig(pathFigures+'WeightCorr_' + str(int(TimeBinSize)) + 'msBins_dayPick'+ str(int(dayPick)) +'.eps', format='eps')

#export data
dataStr = ['WeightCorr_' + str(int(TimeBinSize)) + 'msBins_dayPick'+ str(int(dayPick)) + '_week'+ str(iWk) for iWk in range(3)] #create index string
df = df.append(pd.DataFrame(data=weightCorrBoxPlot,index=dataStr))
dataStr = ['tagWeightCorr_' + str(int(TimeBinSize)) + 'msBins_dayPick'+ str(int(dayPick)) + '_week'+ str(iWk) for iWk in range(3)] #create index string
df = df.append(pd.DataFrame(data=tagWCorr,index=dataStr))

# f,a = plt.subplots(2,3)
# for iweek in range(3):
#
#     w = [[w for wday in weightsall[0+iweek*3:3+iweek*3] for w in wday[i]] for i in range(2)]
#
#     a[1][iweek].hist2d(w[0], w[1], 50, normed=true, vmax=3)
#     a[1][iweek].set_ylim([-1, 1])
#     a[1][iweek].set_xlim([-1, 1])
#
#     a[0][iweek].scatter(w[0], w[1], s=1, alpha=0.3)
#     a[0][iweek].plot([np.min(w), np.max(w)], [np.min(w), np.max(w)], 'y--')
#     regr = linear_model.linearregression()
#     regr.fit(np.reshape(w[0], [-1, 1]), w[1])
#     pred_plot = regr.predict([[np.min(w[0])], [np.max(w[0])]])
#     c, p = scstat.pearsonr(w[0], w[1])
#     a[0][iweek].plot([[np.min(w[0])], [np.max(w[0])]], pred_plot, 'r')
#     a[0][iweek].text(np.min(w), 0.1 * np.min(w) + 0.9 * np.max(w),
#                                  'c=' + '{:.3f}'.format(c) + ' p=' + '{:.3f}'.format(p))
#
#     a[0][iweek].set_ylim([-2, 2])
#     a[0][iweek].set_xlim([-2, 2])

# f.savefig(pathFigures+'Weights_' + str(int(TimeBinSize)) + 'msBins_dayPick'+ str(int(dayPick)) +'.eps', format='eps')
#
# fPerfSeg.savefig(pathFigures+'popSegmentsSummary_' + str(int(TimeBinSize)) + 'msBins_'+ str(int(nSeg)) + 'seg.eps', format='eps')
# fSummary.savefig(pathFigures+'decodeSummary_' + str(int(TimeBinSize)) + 'msBins.eps', format='eps')

df.to_csv(pathFigures + 'dataDecoding_popSegments.csv')

print('done done done')
plt.show(block=True)