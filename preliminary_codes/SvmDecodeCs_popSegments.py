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
#import caiman.base.rois
import scipy
from scipy.io import loadmat
import pickle
from sklearn import svm
import time

#parameters
expt = 'CircleSquare'
frame = 'Arena'
saveOutput = True

BinNb = 12  # for place maps
TimeBinSize = 1 * 1000
min_cell_nb = 15
nPercReshape = 5
x_day = np.array([1, 2, 3, 8, 9, 10, 15, 16, 17])

print(TimeBinSize)

# SVM parameters
cvRep = 100
maxIter = 2000000
verbose = 0
cvParams = [1e-10, maxIter, verbose]

pathAll = '/Volumes/GoogleDrive/My Drive/Fenton_Lab/CellPress_CircleSquare/MNS_ANALYSIS/Data/' + expt + '/'
ExcelPath = pathAll + 'FilenameMatch.xlsx'

savePath = pathAll + 'analysisFiles_Simon/' + 'CS_Decoder_SVM_' + str(
    int(TimeBinSize)) + 'msBins_' + 'popSegments_'+str(int(nPercReshape)) +'seg.file'

animalList = ['M19', 'M20', 'M29', 'M34', 'M35', 'M39']

performance = []
weights = []
performanceTopCells = []
pcRatioAll = []
highRateRatioAll = []
lowRateRatioAll = []
StdAll = []

aPc_all = plt.subplots(2, 1)[1]

for animal in animalList:

    aBin = plt.subplots(1, 1)[1]

    performanceTopCellsAni = []

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

    # pre-allocate
    day_count_load_order = []
    pcRatio = []
    highRateRatio = []
    lowRateRatio = []
    coefAll = []
    performance_tmp = [[] for i in range(2)]
    weights_tmp = [[] for i in range(2)]
    SegWeightRatio = []

    # create figures
    f_svm, a_svm = plt.subplots(9, 7)
    fTauDay, aTauDay = plt.subplots(3, 3)

    for i_day, day in enumerate(Days):

        # load data
        xlOut, fileNames, physio = loadPhysiology(day, animal, xlAll, pathAll, expt, frame)
        VideoName, DayCount, Sessions, doAnalysis, Training, atn = xlOut
        MnsTSFileName, TrackerFileName = fileNames
        S_all, A_tmp, rec_idx = physio

        # in case not all day in Days are loaded
        day_count_load_order.append(DayCount)

        sess_indic = []
        s_all = []

        # extract data for the day
        for i_vid, vid in enumerate(VideoName):

            if doAnalysis[i_vid] > 0:

                # extract physiology from longer file
                s_tmp = np.transpose(S_all[rec_idx[i_vid][0]:rec_idx[i_vid][1]] / np.std(S_all, 0))
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

            start = time.time()

            #concatenate data of the day
            sess_concat_tmp = np.concatenate(sess_indic, 0)
            S_concat_tmp = np.transpose(np.concatenate(s_all, 1))

            S_concat = S_concat_tmp[np.any(~np.isnan(S_concat_tmp), 1)]
            sess_concat = sess_concat_tmp[np.any(~np.isnan(S_concat_tmp), 1)]

            # self training nth recording with cross validation
            pCv, cCv,  i__ = SVM_decode(S_concat, sess_concat, do_cv=True, cv_rep=cvRep, cv_param=[1e-10, 100000, 1])
            performance_tmp[1].append(pCv)
            weights_tmp[1].append(cCv)

            # EVALUATE ROLE OF EACH SEGMENT OF CELLS
            coefRank = np.argsort(np.abs(cCv))[::-1]
            nReshape = len(coefRank) // nPercReshape
            coefRankSegmented = np.reshape(coefRank[:(len(coefRank)//nReshape)*nReshape], [-1, nReshape])

            performanceTopCells_tmp = []
            for iSeg in range(np.min([len(coefRankSegmented), nPercReshape])):
                idx = coefRankSegmented[iSeg]
                pSeg, cSeg, i__ = SVM_decode(S_concat[:, idx], sess_concat, do_cv=True, cv_rep=cvRep, cv_param=[1e-10, 100000, 0])
                performanceTopCells_tmp.append(pSeg)
            performanceTopCellsAni.append(performanceTopCells_tmp)

            # finding the index of place cells
            nCells = len(cCv)
            idx_day = np.where([d == day for d in dayRealListMap])[0]
            pair_same = [idx_pair for idx_pair in itertools.combinations(idx_day, 2) if
                         sessListMap[idx_pair[0]][:3] == sessListMap[idx_pair[1]][:3]]
            isPcEnv = [np.any([PcIdxList[idx] for idx in iIdxPair], 0) for iIdxPair in pair_same]

            rate = np.sum(S_concat,0)
            highRateThresh = np.percentile(np.sum(S_concat, 0), 80)
            lowRateThresh = np.percentile(np.sum(S_concat, 0), 20)
            isHighRate = [r > highRateThresh for r in rate]
            isLowRate = [r < lowRateThresh for r in rate]

            # compute %cell that are PCs (in any env) in the top 20% of coef or on average
            pcRatio.append([np.mean([np.any(isPcEnv, 0)[x]
                                     for i_x, x in enumerate(np.argsort(np.abs(cCv)))
                                     if i_x > nCells - (nCells // 5) * i_perc])
                            for i_perc in [1, 5]])

            highRateRatio.append([np.mean([isHighRate[x]
                                     for i_x, x in enumerate(np.argsort(np.abs(cCv)))
                                     if i_x > nCells - (nCells // 5) * i_perc])
                            for i_perc in [1, 5]])

            lowRateRatio.append([np.mean([isLowRate[x]
                                     for i_x, x in enumerate(np.argsort(np.abs(cCv)))
                                     if i_x > nCells - (nCells // 5) * i_perc])
                            for i_perc in [1, 5]])

            # save coefs
            coefAll.append(cCv)

            # self training nth recording
            clf = svm.LinearSVC(tol=1e-10, max_iter=100000, verbose=1)
            clf.fit(S_concat, sess_concat)
            a_svm[i_day][2].plot([np.sum(s * clf.coef_[0]) for s in S_concat], '.-')
            a_svm[i_day][2].plot(-clf.intercept_ * np.ones_like(sess_concat))
            a_svm[i_day][2].set_title(clf.score(S_concat, sess_concat))
            performance_tmp[0].append(clf.score(S_concat, sess_concat))
            weights_tmp[0].append(clf.coef_[0])
            a_svm[i_day][6].plot(clf.coef_[0])

            # for plotting example only
            idx_train = np.concatenate([np.random.choice(np.where(sess_concat == i_bool)[0],
                                 int(np.sum(sess_concat == i_bool)*0.7), replace=False) for i_bool in [0, 1]], 0)
            idx_test = np.array([x for x in range(len(sess_concat)) if x not in idx_train])
            clf = svm.LinearSVC(tol=1e-10, max_iter=100000, verbose=1)
            clf.fit(S_concat[idx_train], sess_concat[idx_train])
            a_svm[i_day][3].plot([np.sum(s * clf.coef_[0]) for s in S_concat])
            a_svm[i_day][3].plot(-clf.intercept_ * np.ones_like(sess_concat))
            a_svm[i_day][3].set_title(np.round(pCv,2))

        else:
            pcRatio.append([np.nan, np.nan])
            lowRateRatio.append([np.nan, np.nan])
            highRateRatio.append([np.nan, np.nan])
            coefAll.append([np.nan])
            performanceTopCellsAni.append([np.nan])
            for i in range(2):
                weights_tmp[i].append([np.nan])
                performance_tmp[i].append(np.nan)

    # fTau, aTau = plt.subplots(1,2)
    # aTau[0].plot(SegWeightRatio, '.')
    # aTau[1].plot(corrTauWeight[0], 'r.')
    # aTau[1].plot(corrTauWeight[1], 'b.')

    a_svm[0][6].set_title("cv'd coefs")

    aPc = plt.subplots(1, 4)[1]
    aPc[0].plot(x_day, pcRatio, '.-', color='k')
    aPc[0].plot(x_day, highRateRatio, '.-', color='r')
    aPc[0].plot(x_day, lowRateRatio, '.-', color='b')
    aPc[0].set_title('% PC in top 20 vs overall')
    aPc[1].plot(x_day, [c[0] / c[1] for c in pcRatio], '.-', color='k')
    aPc[1].plot(x_day, [c[0] / c[1] for c in highRateRatio], '.-', color='r')
    aPc[1].plot(x_day, [c[0] / c[1] for c in lowRateRatio], '.-', color='b')
    aPc[1].set_title('ratio')

    coefAllAbs = [np.abs(c) for c in coefAll]

    for i_plot, iCoef in enumerate(coefAllAbs):
        if len(iCoef)>1:
            aPc[2].hist(iCoef, 100, histtype='step', cumulative=True,
                        density=True, color=[1 - i_plot / 10, i_plot / 10, i_plot / 10])
    aPc[2].set_title('coef dist')

    aPc[3].plot([np.std(iCoef) for iCoef in coefAllAbs], '.-')
    aPc[3].set_title('dist std')

    aPc_all[0].plot(x_day, [c[0] / c[1] for c in pcRatio], '.-', color='k')
    aPc_all[0].plot(x_day, [c[0] / c[1] for c in highRateRatio], '.-', color='r')
    aPc_all[0].plot(x_day, [c[0] / c[1] for c in lowRateRatio], '.-', color='b')
    aPc_all[1].plot([np.std(iCoef) for iCoef in coefAllAbs])

    pcRatioAll.append([c[0] / c[1] for c in pcRatio])
    highRateRatioAll.append([c[0] / c[1] for c in highRateRatio])
    lowRateRatioAll.append([c[0] / c[1] for c in lowRateRatio])
    StdAll.append([np.std(iCoef) for iCoef in coefAllAbs])

    performance.append(performance_tmp)
    weights.append(weights_tmp)
    performanceTopCells.append(performanceTopCellsAni)

    plt.figure()
    perfTopCellsReshapedAni = [[pp[iSeg] if len(pp) > iSeg else np.nan for pp in performanceTopCellsAni]
                            for iSeg in range(nPercReshape)]

    for i_plot, p in enumerate(perfTopCellsReshapedAni):
        plt.plot(p, '.-.', color=[1 - i_plot / nPercReshape, i_plot / nPercReshape, i_plot / nPercReshape])

    plt.plot(performance_tmp[0], 'k-.')
    plt.suptitle('perf per' + str(nPercReshape) + '% segment of cells')

    # print('tmp saving...')
    # with open(save_path, 'wb') as f:
    #     pickle.dump([performance, weights, performanceTopCells, CoefDevAll, StdAll, animalList], f, pickle.HIGHEST_PROTOCOL)
    # print('tmp saved')

Title_perf = ['current decode current _best',
              'current decode current _cv']

performance_av = np.mean(performance, 0)

x_day = np.array([1, 2, 3, 8, 9, 10, 15, 16, 17])

a = plt.subplots(2, 2)[1]

for i_plot, i_idx in enumerate([0, 1]):
    a.flat[i_plot].errorbar(x_day, np.mean(performance, 0)[i_idx], np.std(performance, 0)[i_idx]/2,
                            color='k', fmt='-.')

    for i_ani_perf in performance:
        a.flat[i_plot].plot(x_day, i_ani_perf[i_idx], '.')

    a.flat[i_plot].set_ylim([0.4, 1])
    a.flat[i_plot].set_title(Title_perf[i_idx])

if saveOutput:
    print('saving...')
    with open(savePath, 'wb') as f:
        pickle.dump([performance, weights, performanceTopCells, [pcRatioAll,highRateRatioAll, lowRateRatioAll], StdAll, animalList], f, pickle.HIGHEST_PROTOCOL)
    print('saved')
else:
    print('not saving')

plt.figure()
plt.errorbar(x_day, np.nanmean(pcRatioAll, 0), np.nanstd(pcRatioAll, 0)/np.sqrt(len(animalList)))
for c in pcRatioAll:
    plt.plot(x_day, c, 'k.')
plt.suptitle('ratio %PC in top 20 vs overall')

plt.figure()
plt.errorbar(x_day, np.nanmean(highRateRatioAll, 0), np.nanstd(highRateRatioAll, 0)/np.sqrt(len(animalList)))
for c in highRateRatioAll:
    plt.plot(x_day, c, 'k.')
plt.suptitle('ratio HR in top 20 vs overall')

plt.figure()
plt.errorbar(x_day, np.nanmean(lowRateRatioAll, 0), np.nanstd(lowRateRatioAll, 0)/np.sqrt(len(animalList)))
for c in lowRateRatioAll:
    plt.plot(x_day, c, 'k.')
plt.suptitle('ratio LR in top 20 vs overall')

StdAllNormed = [[i_s/np.nanmean(s[0:3]) for i_s in s] for s in StdAll]

plt.figure()
plt.errorbar(x_day, np.nanmean(StdAllNormed, 0), np.nanstd(StdAllNormed, 0)/np.sqrt(len(animalList)))
for c in StdAllNormed:
    plt.plot(x_day, c, 'k.')
plt.suptitle('coef dist std - normed to 1st week')

aPc_all[0].set_title('ratio %PC in top 20 vs overall')
aPc_all[1].set_title('coef dist std')

print('done done done')

plt.show(block=True)