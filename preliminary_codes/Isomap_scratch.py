import matplotlib as mpl
mpl.use('TkAgg')
from mpl_toolkits.mplot3d import Axes3D

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pickle
import itertools
import time
from scipy.ndimage import gaussian_filter1d
from sklearn import decomposition, manifold, neighbors
from scipy.spatial.distance import pdist
from ripser import ripser as tda
import matplotlib.gridspec as gridspec
from Pos_file_processing import *
import pandas as pd

path_all = '/Volumes/GoogleDrive/My Drive/Fenton_Lab/CellPress_CircleSquare/MNS_Analysis/Data/CircleSquare/'
TimeBinSize_AVC = 1000 # not in use
TimeBinSize_TAU = 1000

# to remvoe cells with low activity
threshMult = 2
thresh = 0.01* threshMult*(TimeBinSize_TAU/1000)

saveFig = True

sigma = 1  # Kernel width (in s) to smooth data

n_neighbors = 5
target_dim = 3
h2MaxDt = 700

preThrsh = False
Bettithrsh = True # treshold nb of data point to h2MaxDt
doBetti = True
doStd = True
doDiscardCells = False  # keep false otherwise use other script Isomap_discardPercPop.py
doPTI = False

if doPTI:
    pathFigure = '/Volumes/GoogleDrive/My Drive/Fenton_Lab/CellPress_CircleSquare/MNS_Analysis/Figures/CircleSquare/savedFigs/z/IsoMap_Simon/'
else:
    pathFigure = '/Volumes/GoogleDrive/My Drive/Fenton_Lab/CellPress_CircleSquare/MNS_Analysis/Figures/CircleSquare/savedFigs/IsoMap_Simon/'

animalList = ['M19']
nAni = len(animalList)

envList = ['HMC1','HMC2', 'CYL1','CYL2', 'RCT1','RCT2']
indicColors  = [[0,0,0],[0.5,0.5,0.5,0.5], [0,0,1],[0,0,0.5], [1,0,0],[0.5,0,0]]

tauBins = [-0.05, -0.04, -0.03, -0.025, -0.02, -0.0175, -0.015, -0.0125, -0.01, -0.005,
           0.005, 0.01, 0.02, 0.03, 0.04, 0.06, 0.1]
tauPlot = [tauBins[0] - (tauBins[1] - tauBins[0]) / 2]
tauPlot.extend([b + d / 2 for b, d in zip(tauBins, np.diff(tauBins))])
tauPlot.extend([tauBins[-1] - (tauBins[-2] - tauBins[-1]) / 2])

meanIsoAll = []
meanPcaAll = []
meanAvAll  = []

if doPTI:
    loadStr = 'Zlinear'
else:
    loadStr = 'TempCorr'

if doDiscardCells:
    ratioDiscard = 2
    discard = 'kc'
    print(discard)

if target_dim == 2:
    fAvProj, aAvProj = plt.subplots(1, 2)
    fAvProjOrig, aAvProjOrig = plt.subplots(1, 2)

for iAni, animal in enumerate(animalList):

    print(animal)

    with open(path_all + 'analysisFiles_Simon/' + loadStr + '_List_' +
              str(int(TimeBinSize_AVC)) + 'msAV_' + str(int(TimeBinSize_TAU)) + 'msKT'
              '_CircleSquare_' + animal +'.file', 'rb') as f:
        outputList = pickle.load(f)

    # from cnmfe_temporal_analysis
    # listsToSave = [[day_real_list, day_count_list, sess_list, S_conv_avc_list, S_conv_tau_list, isSeparatedList,
    #                 tauPairsSingleList, tauVecSingleList, trainingList],
    #                [trackList, [], aMaxLocList],
    #                [splitList, iDayList]]

    # from cnmfe_zscore_analysis
    # listsToSave = [[day_real_list, day_count_list, sess_list, S_conv_avc_list, [zConvTauList, S_conv_tau_list, sConvTauPcList], isSeparatedList,
    #                 tauPairsSingleList, [zTauVecSingleList,tauVecSingleList], trainingList, xBinnedPcList, yBinnedPcList],
    #                [trackList, [zBinnedRateMapsList, zRateMapsList], aMaxLocList],
    #                [atnSplitList, iDayList]]
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
    if target_dim<4:
        figScatter = plt.figure(figsize=(8, 8))

    if target_dim == 10 and doBetti:
        figBetti = plt.figure(figsize=(10, 8))

    meanIsoAni = []
    meanPcaAni = []
    meanAvAni  = []

    for iDay, day in enumerate(range(1,10)):
        print([iDay,day])
        sessIdx = np.where([d==day for d in dayCountList])[0]

        if len(sessIdx)<5:
            meanIsoAni.append([])
            meanPcaAni.append([])
            meanAvAni.append([])
            continue

        sConcat = []
        sessIndic = []
        sessIndicFull = []
        meanIsoDay = []
        meanPcaDay = []
        meanAvDay = []

        idxPosList = []
        idxNegList = []
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

            if doDiscardCells:
                # FIND PAIRS OVER-REPRESENTED IN CERTAIN TAU CORR
                # get pair numbers and tau corr
                tauVecTmp = tauVecSingleList[iSessIdx]
                idxTauPairsTmp = tauPairsSingleList[iSessIdx]

                isSeparated = isSeparatedList[iSessIdx] # to remove nearby cells
                isAnyNan = [np.isnan(t[0]) for t in tauVecTmp] #create filter for Nan's

                # filter out cells that are rarely active
                meanRate = [np.mean([s > 0 for s in c]) for c in sTauList[iSessIdx]]
                rate = [meanRate[p[0]] > thresh and meanRate[p[1]] > thresh for p in idxTauPairsTmp]

                # create filter for spatial separation and significance and non nan correlation
                tauFilter = [sep and r and not iNan for sep, iNan, r in zip(isSeparated, isAnyNan, rate)]
                meanRate = rate = isAnyNan = isSeparated = [] # delete variable from memory

                # filter tau values and pair numbers
                tauVec = [t[0] for t, f in zip(tauVecTmp, tauFilter) if f]
                idxTauPairs = [idx for idx, f in zip(idxTauPairsTmp, tauFilter) if f]

                keptCells = np.unique([x for pair in idxTauPairs for x in pair])
                nKeptCells = len(keptCells)

                nBins = len(tauBins) + 1
                binList = [*range(nBins)]
                # sort pairs in different segments
                equalSegIdx = np.digitize(tauVec, tauBins)  # a vector containing the bin number for each pair

                # for each segment, list the pairs it contains
                equalSegPairList = [[idxTauPairs[iPair] for iPair in np.where(equalSegIdx == iBin)[0]]
                                    for iBin in binList]

                # evaluate how many times each cell is in each segment
                segPairsCountTmp = [[np.sum([np.any([p[0] == i, p[1] == i]) for p in pairsSeg]) for i in keptCells]
                                    for pairsSeg in equalSegPairList]

                # normalize for the number of pairs a cell is in
                segPairsCount = [x / np.sum(segPairsCountTmp, 0) for x in segPairsCountTmp]

                keptCellsList.append(keptCells)

                nDiscard = len(keptCells)//ratioDiscard

                if 'HMC' not in sessList[iSessIdx][:3]:
                    idxPosList.append([keptCells[idx] for idx in np.argsort(segPairsCount[-1])[::-1][:nDiscard]
                                       if segPairsCount[-1][idx]>0])
                    idxNegList.append([keptCells[idx] for idx in np.argsort(segPairsCount[0])[::-1][:nDiscard]
                                       if segPairsCount[0][idx]>0])

        # sConcatFilt = [[s for iCell, s in enumerate(sRec) if iCell not in np.unique(idxNegList)]
        #                for sRec in sConcat]

        if doDiscardCells:
            keptCellsAll = np.unique([x_ for x in keptCellsList for x_ in x])
            if discard == 'pos':
                discardCells = np.unique([x_ for x in idxPosList for x_ in x])
            elif discard == 'neg':
                discardCells = np.unique([x_ for x in idxNegList for x_ in x])
            elif discard == 'kc':
                discardCells = []

            sConcatFilt = [[s for iCell, s in enumerate(sRec)
                            if (iCell in keptCellsAll) and (iCell not in discardCells)]
                           for sRec in sConcat]

            if discard == 'kc':
                strCode = 'keptCells'
            elif discard == 'pos':
                strCode = 'keptCells-' + str(int(100 / ratioDiscard)) + 'percPos_noHMC'
            elif discard == 'neg':
                strCode = 'keptCells-' + str(int(100 / ratioDiscard)) + 'percNeg_noHMC'

            print(len([iCell for iCell, s in enumerate(sConcat[0])
                       if (iCell in keptCellsAll) and (iCell not in discardCells)]))

        else:
            sConcatFilt = [[s for iCell, s in enumerate(sRec)]
                           for sRec in sConcat]

        # sConcatFilt = [[s for iCell, s in enumerate(sRec)
        #                 if iCell in keptCellsAll]
        #                for sRec in sConcat]


        if doStd:
            sAll = np.transpose(np.concatenate(sConcatFilt, 1))/np.std(np.transpose(np.concatenate(sConcatFilt, 1)),0)
        else:
            sAll = np.transpose(np.concatenate(sConcatFilt, 1))

        # remove cells with NaNs (there was only one in one day (M29))
        if np.sum(np.any(np.isnan(sAll), 1)) / len(sAll) > 0:
            print(np.sum(np.any(np.isnan(sAll), 1)) / len(sAll))
            sAll = np.array([[sCell for sCell, filt in zip(sTime, (np.any(np.isnan(sAll), 0))) if not filt] for sTime in sAll])
            print(np.sum(np.any(np.isnan(sAll), 1)) / len(sAll))

        sessIndicAll = np.concatenate(sessIndic, 0)
        sessIndicFullAll = np.concatenate(sessIndicFull, 0)

        #thresh to remove points close to zero
        # sAll = [s for s in sAll if np.mean(s)>0.1]
        # sessIndicAll = [sess for s, sess in zip(sAll, sessIndicAll) if np.mean(s)>0.1]
        # sessIndicFullAll = [sess for s, sess in zip(sAll, sessIndicFullAll) if np.mean(s)>0.1]

        if preThrsh:
            print('neighbor thresholding 1')
            # a) find number of neighbors of each point within radius of 1st percentile of all
            # pairwise dist.
            dist = pdist(sAll, 'euclidean')
            rad = np.percentile(dist, 1)
            neigh = neighbors.NearestNeighbors()
            neigh.fit(sAll)
            num_nbrs = [*map(len, neigh.radius_neighbors(X=sAll, radius=rad, return_distance=False))]

            # b) threshold out points with low density
            thrsh_prcnt = 5
            threshold = np.percentile(num_nbrs, thrsh_prcnt)
            while np.mean([num_nbrs > threshold]) < 0.8:
                threshold -= 1
                print(threshold)
            thrshRates = sAll[num_nbrs > threshold]
            thrshIndic = sessIndicAll[num_nbrs > threshold]
            sAll = thrshRates
            sessIndicAll = thrshIndic


        if doPTI:
            data_to_use = np.sign(sAll.copy())*np.sqrt(sAll.copy()*np.sign(sAll.copy())) # to handle negative numbers
            # data_to_use = sAll.copy()
        else:
            data_to_use = np.sqrt(sAll.copy())

        # do Isomap transformation
        iso_instance = manifold.Isomap(n_neighbors, target_dim)
        dataIso = iso_instance.fit_transform(data_to_use)

        # do PCA transformation
        pca1 = decomposition.PCA(target_dim)
        dataPca = pca1.fit_transform(data_to_use)

        if target_dim == 10 and doBetti:
            if Bettithrsh:
                print('neighbor thresholding 2')
                # a) find number of neighbors of each point within radius of 1st percentile of all
                # pairwise dist.
                dist = pdist(dataIso, 'euclidean')
                rad = np.percentile(dist, 1)
                neigh = neighbors.NearestNeighbors()
                neigh.fit(dataIso)
                num_nbrs = [*map(len, neigh.radius_neighbors(X=dataIso, radius=rad, return_distance=False))]

                # b) threshold out points with low density
                thrsh_prcnt = 20
                threshold = np.percentile(num_nbrs, thrsh_prcnt)
                thrsh_rates = dataIso[num_nbrs > threshold]
                dataBetti = thrsh_rates
            else:
                dataBetti = dataIso

            print('computing Betti numbers...')
            results = {'h0': [], 'h1': [], 'h2': []}
            #Betti
            # H0 & H1
            H1_rates = dataBetti
            print('h0-1')
            barcodes = tda(H1_rates, maxdim=1, coeff=2)['dgms']
            results['h0'] = barcodes[0]
            results['h1'] = barcodes[1]
            if len(dataBetti) > h2MaxDt:
                print('shortening data for Betti nb 2')
                idx = np.random.choice(np.arange(len(dataBetti)), h2MaxDt, replace=False)
                H2_rates = dataBetti[idx]
            else:
                H2_rates = dataBetti
            print('h2')
            barcodes = tda(H2_rates, maxdim=2, coeff=2)['dgms']
            results['h2'] = barcodes[2]
            print('done')

            plot_barcode = True
            if plot_barcode:
                print('plotting...')
                col_list = ['r', 'g', 'm', 'c']
                h0, h1, h2 = results['h0'], results['h1'], results['h2']
                # replace the infinity bar (-1) in H0 by a really large number
                h0[~np.isfinite(h0)] = 100
                # Plot the longest barcodes only
                plot_prcnt = [99, 98, 90]  # order is h0, h1, h2
                to_plot = []
                for curr_h, cutoff in zip([h0, h1, h2], plot_prcnt):
                    bar_lens = curr_h[:, 1] - curr_h[:, 0]
                    if len(curr_h)>0:
                        plot_h = curr_h[bar_lens > np.percentile(bar_lens, cutoff)]
                        to_plot.append(plot_h)

                # gs = gridspec.GridSpec(3, 4)
                for curr_betti, curr_bar in enumerate(to_plot):
                    # ax = fig.add_subplot(gs[curr_betti, :])
                    ax = figBetti.add_subplot(3, 9, curr_betti*9+ iDay +1)
                    for i, interval in enumerate(reversed(curr_bar)):
                        ax.plot([interval[0], interval[1]], [i, i], '.-', color=col_list[curr_betti],
                                lw=1.5)
                    # ax.set_xlim([0, xlim])
                    # ax.set_xticks([0, xlim])
                    ax.set_ylim([-1, len(curr_bar)])
                    ax.set_xlim([-1, 100])
                    # ax.set_yticks([])
                print('plotted')

        # 2D plot
        if target_dim == 2:
            meanEnv = [[], [], [], [], [], []]
            ax = figScatter.add_subplot(3, 3, iDay+1)

            for iDimRed, dimRed in enumerate([dataIso, dataPca, data_to_use]):

                for iIndic, envIndic in enumerate(envList): #['HMC1','HMC2', 'CYL1','CYL2', 'RCT1','RCT2']
                    proj_data = dimRed[[envIndic in s for s in sessIndicFullAll], :]
                    meanEnv[iIndic] = np.mean(proj_data, 0)

                    # scatter plot for Isomap with average vector
                    if iDimRed == 0:
                        ax.plot([0, np.mean(proj_data, 0)[0]], [0, np.mean(proj_data, 0)[1]], '*-', color=indicColors[iIndic],
                                lw=2)
                        ax.scatter(proj_data[:, 0], proj_data[:, 1],
                                   s=1, alpha=0.2, color=indicColors[iIndic])

                # plotting 2D projection
                if iDimRed < 1 and iDay>-1:
                    # rotate vectors to align them
                    hmc = np.mean(meanEnv[:2], 0)
                    cyl = np.mean(meanEnv[2:4],0)

                    #create rotation array to put hmc to the right
                    h = np.sqrt(hmc[0]**2+hmc[1]**2) #to make the rotation vector len=1
                    rot = np.array([[hmc[0], hmc[1]], [-hmc[1], hmc[0]]])/h

                    # rotate cyl and extract sign
                    cylRot = np.sum(cyl * rot, 1)
                    doFlip = [1, np.sign(cylRot[1])]

                    for iEnv in range(3):
                        # to plot each individual session as vector
                        for iVec in meanEnv[0+2*iEnv:2+2*iEnv]:
                            newProj = np.sum(iVec * rot, 1)*doFlip
                            aAvProj[iDimRed].plot([0, newProj[0]], [0, newProj[1]],
                                                  color=indicColors[2 * iEnv])

                        # to plot average of each pair of 2 sessions as vector
                        # newProj = np.sum(np.mean(meanEnv[0+2*iEnv:2+2*iEnv], 0) * rot, 1)*doFlip
                        # aAvProj[iDimRed].plot([0, newProj[0]], [0, newProj[1]], color=indicColors[2*iEnv])




        if target_dim == 3:
            ax = figScatter.add_subplot(3, 3, iDay+1, projection='3d')
            for iIndic, envIndic in enumerate(envList):
                proj_data = dataIso[[envIndic in s for s in sessIndicFullAll], :]
                ax.scatter(proj_data[:, 0], proj_data[:, 1], proj_data[:, 2],
                           s=3, alpha=0.2, edgecolor='face', color=indicColors[iIndic])
                ax.plot([0,np.mean(proj_data,0)[0]], [0,np.mean(proj_data,0)[1]], [0,np.mean(proj_data,0)[2]],
                        '*-', color=indicColors[iIndic], lw = 5)
                # ax.set_title('1 STD')
                # ax.plot(proj_data[:, 0], proj_data[:, 1], proj_data[:, 2], '.--', color=cols[iIndic])

        # # separate RCT and CYL
        dataSvm = [d for d, s in zip(dataIso, sessIndicAll) if s >= 0]
        sessSvm = [s for s in sessIndicAll if s >= 0]
        perfIso, coef, i = SVM_decode(np.array(dataSvm), np.array(sessSvm), do_cv = False, cv_param=[1e-10, 10000, 0])
        # perf, coef = SVM_decode(np.array(data_to_use), np.array(sessSvm), cv_param=[1e-10, 10000, 0])
        # print([perfIso, perf])

        ax = figSvm.add_subplot(3, 3, iDay + 1)
        for iIndic, envIndic in enumerate(envList):
            projDataIso = dataIso[[envIndic in s for s in sessIndicFullAll], :]
            meanIsoDay.append(np.mean(projDataIso, 0))

            projDataPca = dataPca[[envIndic in s for s in sessIndicFullAll], :]
            meanPcaDay.append(np.mean(projDataPca, 0))

            projDataAV = data_to_use[[envIndic in s for s in sessIndicFullAll], :]
            meanAvDay.append(np.mean(projDataAV, 0))

            ax.plot(np.sum(projDataIso * coef, 1), iIndic*np.ones_like(np.sum(projDataIso * coef, 1)), '.',
                     ms = 1, color=indicColors[iIndic])

        meanIsoAni.append(meanIsoDay)
        meanPcaAni.append(meanPcaDay)
        meanAvAni.append(meanAvDay)

    meanIsoAll.append(meanIsoAni)
    meanPcaAll.append(meanPcaAni)
    meanAvAll.append(meanAvAni)

    if target_dim<4 and saveFig:
        figScatter.savefig(pathFigure + 'IsomapProj_' +str(target_dim)+ 'D_' + animal +'.eps', format='eps')

    if not doDiscardCells: strCode = 'no discard'

    if doBetti and target_dim==10 and saveFig: figBetti.savefig(pathFigure + 'IsomapBetti_' + animal +'.eps', format='eps')

    print(strCode)

f,a  = plt.subplots(1,2)
pairList = [[(2,3),(4,5)], [(2,4),(2,5),(3,4),(3,5)]]
colorList = ['k', 'g', 'm']

metricStr = ['isoMap', 'pca', 'full']

#exporting data
df = pd.DataFrame()

for iMetric, metric in enumerate([meanIsoAll, meanPcaAll, meanAvAll]):
    z = [[[np.mean([np.sqrt(np.sum((m[x1]-m[x2])**2)) for x1,x2 in pairs]) for pairs in pairList]
          if len(m)>0 else [np.nan, np.nan] for m in mAni] for mAni in metric]
    zRatio = [[zDay[0]/zDay[1] for zDay in zAni] for zAni in z]
    print(z)

    for zAni in zRatio:
        a[0].plot([*range(9)], zAni, '.', color=colorList[iMetric])
    a[0].errorbar([*range(9)], np.nanmean(zRatio,0), np.nanstd(zRatio,0)/np.sqrt(nAni), color=colorList[iMetric])

    dataBoxPlot = [[z for zAni in zRatio for z in zAni[0 + iWk:3 + iWk] if np.isfinite(z)] for iWk in [0, 3, 6]]
    tagData = [[[animalList[iAni], iWk/3, iDay] for iAni, zAni in enumerate(zRatio) for iDay, z in enumerate(zAni[0 + iWk:3 + iWk]) if np.isfinite(z)] for iWk in [0, 3, 6]]
    a[1].boxplot(dataBoxPlot, positions=[0+iMetric, 6+iMetric, 12+iMetric])

    #export data
    # dataStr = [metricStr[iMetric]]  # create index string
    # df = df.append(pd.DataFrame(data=np.concatenate(dataBoxPlot), columns=dataStr).transpose())
    # dataStr = ['tag_' + s + '_' + metricStr[iMetric] for s in ['ani', 'week', 'day']]  # create index string
    # df = df.append(pd.DataFrame(data=np.concatenate(tagData), columns=dataStr).transpose())

a[1].set_xlim([-1, 15])
f.suptitle(strCode)

# df.to_csv(pathFigure + 'dataDimRed_kc.csv')

print(strCode)

if saveFig:
    plt.savefig(pathFigure + 'QuantDimRed_dimRed_'+str(target_dim)+'D'+ strCode + '.eps', format='eps')
    if False: #target_dim == 2:
        fAvProj.savefig(pathFigure + 'QuantDimRed_CS_2DAverageProj_meanEnv_wk3.eps', format='eps')
    print('saved')
else:
    print('not saved')

plt.show(block=True)