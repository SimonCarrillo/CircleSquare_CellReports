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
from sklearn import svm

expt = 'CircleSquare'
frame = 'Arena'

BinNb = 12  # for place maps
TimeBinSize = 1 * 1000
min_cell_nb = 15
LoadAssignments = True
SaveAssignments = False

pathAll = '/Users/eliott/Dropbox/dDocuments/NYU/Fenton/MNS_ANALYSIS/Data/' + expt + '/'
ExcelPath = pathAll + 'FilenameMatch.xlsx'

btwPerfMat  = []
selfPerfMat = []
fMatBtw, aMatBtw = plt.subplots(2,3)
fMatSelf, aMatSelf = plt.subplots(2,3)

for iAni, animal in enumerate(['M19', 'M20', 'M29', 'M34','M35', 'M39']):

    print(animal)

    # process excel file
    xlread = pd.ExcelFile(ExcelPath)
    xlAll = pd.read_excel(xlread, animal)
    xlAll = xlAll[pd.notna(xlAll.videoname)]

    Days = np.unique(list(xlAll.day))
    n_sess = len(xlAll)

    nsep = np.max(Days) - np.min(Days)
    nday = len(Days)

    dayCountList = []
    coefDev = []
    coefAll = []

    sConcat = []
    sessConcat = []

    performance_tmp = [[] for i in range(4)]

    with open(pathAll + 'analysisFiles/' + 'AssignmentsByPairs_' + animal + '.file', 'rb') as f:
        [assignments, [assignmentsByPairs, idPairs]] = pickle.load(f)

    for i_day, day in enumerate(Days):

        xlOut, fileNames, physio = loadPhysiology(day, animal, xlAll, pathAll, expt, frame)
        VideoName, dayCount, Sessions, doAnalysis, training, atnSplit = xlOut
        MnsTSFileName, TrackerFileName = fileNames
        S_all = physio[0]
        rec_idx = physio[2]

        dayCountList.append(dayCount)

        sess_indic = []
        s_all = []

        for i_vid, (vid, sess) in enumerate(zip(VideoName, Sessions)):

            if doAnalysis[i_vid] > 0:
                # extract physiology from longer file AND normalize
                s_tmp = np.transpose(S_all[rec_idx[i_vid][0]:rec_idx[i_vid][1]]/np.std(S_all, 0))

                # load video time file and adjust for delay to start MNS
                MnsFrameCount, MnsTS = read_timestamps(MnsTSFileName[i_vid])

                # correction for lost frames during tiff conversion
                if sess[:3] == 'WIR' or sess[:3] == 'PTC' or len(s_tmp[0])>3333:
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
                    print('')
                    # sess_indic.append(np.zeros_like(S_binned_avc[0]))
                elif 'CYL' in Sessions[i_vid][:3]:
                    sess_indic.append(np.ones_like(sBinned[0]))
                    s_all.append(sBinned)
                elif 'RCT' in Sessions[i_vid][:3]:
                    sess_indic.append(np.zeros_like(sBinned[0]))
                    s_all.append(sBinned)

        if len(sess_indic) > 0:
            sess_concat_tmp = np.concatenate(sess_indic, 0)
            S_concat_tmp = np.transpose(np.concatenate(s_all, 1))
            sConcat.append(S_concat_tmp[np.any(~np.isnan(S_concat_tmp), 1)])
            sessConcat.append(sess_concat_tmp[np.any(~np.isnan(S_concat_tmp), 1)])
        else:
            sConcat.append([])
            sessConcat.append([])

    # # load PCs to get Pc_idx
    # with open(pathAll + 'analysisFiles/MapCorr_' + str(int(BinNb)) + '_Bins_' + animal + '.file', 'rb') as f:
    #     output = pickle.load(f)
    # PcDay_real_list = output[1][0]
    # PcSess_list = output[1][2]
    # PcRate_maps_list = output[1][3]
    # PcIdx_list = output[1][4]

    # load registration
    with open(pathAll + 'analysisFiles/' + 'AssignmentsByPairs_' + animal + '.file', 'rb') as f:
        [assignments, [assignmentsByPairs, idPairs]] = pickle.load(f)

    # [day_real_list, day_count_list, sess_list, rate_maps_list,
    #  PC_idx_list, LinMapsList, OccList, cohList, pInfoList] = output[1]

    selfPerf = []
    btwPerf = []
    dayPairList = []

    for nComb, iComb in enumerate(itertools.permutations(range(len(dayCountList)), 2)):

        # correct for 0-based vs 1-based indexing discrepancy
        dayPair = tuple([dayCountList[i]-1 for i in iComb])

        # check for non-extracted day
        if np.any([len(sConcat[i])<1 for i in iComb]): continue

        sConcatPair = [np.transpose(sConcat[i]) for i in iComb]
        sessConcatPair = [sessConcat[i] for i in iComb]

        if np.diff(dayPair)[0] == 0:
            nMatched = len(sConcatPair[0])
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

        sConcatMatched = [[iS[n] for n in iAssign] for iS, iAssign in zip(sConcatPair, CellList)]

        # idxTrain =
        # idxTest  =

        # train on day iComb[0] no CV
        clf = svm.LinearSVC(tol=1e-10, max_iter=10000, verbose=1)
        clf.fit(np.transpose(sConcatMatched[0]), sessConcatPair[0])

        # test on day iComb[1]
        btwPerf.append(clf.score(np.transpose(sConcatMatched[1]), sessConcatPair[1]))

        # test on self w/ CV
        selfPerfRand = []
        for iT in range(10):
            idxTrain = np.concatenate([
                np.random.choice(np.where(sessConcatPair[0] == i_bool)[0],
                                 int(np.sum(sessConcatPair[0] == i_bool) * 0.7), replace=False)
                for i_bool in [0, 1]], 0)
            idxTest = np.array([x for x in range(len(sessConcatPair[0])) if x not in idxTrain])

            clf = svm.LinearSVC(tol=1e-10, max_iter=10000, verbose=1)
            clf.fit(np.transpose(sConcatMatched[0])[idxTrain], sessConcatPair[0][idxTrain])
            selfPerfRand.append(clf.score(np.transpose(sConcatMatched[0])[idxTest], sessConcatPair[0][idxTest]))

        selfPerf.append(np.mean(selfPerfRand))

        dayPairList.append(dayPair)

    btwPerfMatAni = []
    selfPerfMatAni = []
    for dTrain in range(9):
        btwPerfMatTmp = []
        selfPerfMatTmp = []
        for dTest in range(9):

            if (dTrain, dTest) not in dayPairList:
                btwPerfMatTmp.append(np.nan)
                selfPerfMatTmp.append(np.nan)
                continue

            findIdx = np.where([idx == (dTrain, dTest) for idx in dayPairList])[0][0]
            btwPerfMatTmp.append(btwPerf[findIdx])
            selfPerfMatTmp.append(selfPerf[findIdx])

        btwPerfMatAni.append(btwPerfMatTmp)
        selfPerfMatAni.append(selfPerfMatTmp)

    zPlot = aMatBtw.flat[iAni].imshow(btwPerfMatAni, vmax=1, vmin=0.5)
    plt.colorbar(zPlot, ax=aMatBtw.flat[iAni])
    aMatBtw.flat[iAni].set_ylabel('day train')
    aMatBtw.flat[iAni].set_xlabel('day test')
    aMatBtw.flat[iAni].set_title(animal)

    zPlot = aMatSelf.flat[iAni].imshow(selfPerfMatAni, vmax=1, vmin=0.5)
    plt.colorbar(zPlot, ax=aMatSelf.flat[iAni])
    aMatSelf.flat[iAni].set_ylabel('day train')
    aMatSelf.flat[iAni].set_xlabel('day test')
    aMatSelf.flat[iAni].set_title(animal)


    btwPerfMat.append(btwPerfMatAni)
    selfPerfMat.append(selfPerfMatAni)

f,a = plt.subplots(1,2)
zPlot0 = a[0].imshow(np.nanmean(btwPerfMat,0), vmax=0.75, vmin=0.5)
zPlot1 = a[1].imshow(np.nanmean(selfPerfMat,0), vmax=0.75, vmin=0.5)

a[0].set_title('Decoding perf')
a[1].set_title("self decoding - CV'd control")

plt.colorbar(zPlot0, ax=a[0])
plt.colorbar(zPlot1, ax=a[1])
a[0].set_ylabel('day train')
a[0].set_xlabel('day test')
a[1].set_ylabel('day train+test')
a[1].set_xlabel('day matched')

plt.show(block=True)