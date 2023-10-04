import matplotlib as mpl
mpl.use('TkAgg')

import itertools
import numpy as np
import scipy.stats as scstat
import pandas as pd
from Pos_file_processing import *
import pickle
import time
from scipy.io import savemat

# analyze recordings then save
# temporal analysis for both AV (activity vectors) and Tau (correlation Kendall-Tau)

expt = 'CircleSquare'  # here you can change to ATN (neophobia) or PlaceAvoidance (APA)

animal = 'M39' #['M39','M35','M34','M29','M19','M20']

TimeBinSize_AVC = 1 * 1000
TimeBinSize_TAU = 1 * 1000
min_cell_nb = 15
saveOutput = True
saveOutput2 = False
CombComputation = False

pathAll = '/Volumes/GoogleDrive/My Drive/Fenton_Lab/CellPress_CircleSquare/MNS_Analysis/Data/' + expt + '/'
ExcelPath = pathAll + 'FilenameMatch.xlsx'

mapC_sep_all = []
cv_sep_all = []
tau_sep_all = []
mapC_acr_all = []
cv_acr_all = []
tau_acr_all = []

print(animal)
print(expt)

# process excel file
xlread = pd.ExcelFile(ExcelPath)
xlAll = pd.read_excel(xlread, animal)
xlAll = xlAll[pd.notna(xlAll.videoname)]

Days = np.unique(list(xlAll.day))
print(xlAll.session)
print(Days)

day_real_list = []
iDayList = []
day_count_list = []
sess_list = []
rate_maps_list = []
S_conv_avc_list = []
S_conv_tau_list = []
PC_idx_list = []
templates = []
day_count_load_order = []
trackList = []
Alist = []
aMaxLocList = []
tauVecSingleList = []
pvalueVecSingleList = []
tauPairsSingleList = []
isSeparatedList = []
trainingList = []
atnSplitList = []
aDistList = []

for iDay, day in enumerate(Days):

    print(iDay)

    xlOut, fileNames, physio = loadPhysiology(day, animal, xlAll, pathAll, expt, 'Room')
    VideoName, DayCount, Sessions, doAnalysis, Training, atnSplit = xlOut
    MnsTSFileName, TrackerFileName = fileNames
    S_all, A_tmp, rec_idx = physio

    # in case not all day in Days are loaded
    day_count_load_order.append(DayCount)

    for i_vid, vid in enumerate(VideoName):

        if doAnalysis[i_vid] > 0:
            print(vid)
            # extract physiology from longer file
            s_tmp = np.transpose(S_all[rec_idx[i_vid][0]:rec_idx[i_vid][1]])

            # load video time file and adjust for delay to start MNS
            MnsFrameCount, MnsTS = read_timestamps(MnsTSFileName[i_vid])

            # correction for lost frames during tiff conversion
            if len(s_tmp[0])<=3333:
                MnsTS_corrected_for_tiff = MnsTS[1:1 + 3 * ((len(MnsTS) - 1) // 3)]
            else:
                if len(s_tmp[0])<=6666:
                    MnsTS_corrected_for_tiff = MnsTS[1:10000] + MnsTS[10001:10001 + 3 * (
                            (len(MnsTS) - 10001) // 3)]

                elif len(s_tmp[0])<=9999:
                    MnsTS_corrected_for_tiff = MnsTS[1:10000] + MnsTS[10001:20000] + MnsTS[20001:20001 + 3 * (
                                (len(MnsTS) - 20001) // 3)]

                elif len(s_tmp[0])<=13332:
                    MnsTS_corrected_for_tiff = MnsTS[1:10000] + MnsTS[10001:20000] + MnsTS[20001:30000] + MnsTS[30001:30001 + 3 * (
                                (len(MnsTS) - 30001) // 3)]

                elif len(s_tmp[0])<=16665:
                    MnsTS_corrected_for_tiff = MnsTS[1:10000] + MnsTS[10001:20000] + MnsTS[20001:30000] + MnsTS[30001:40000] + MnsTS[40001:40001 + 3 * (
                                (len(MnsTS) - 40001) // 3)]

                elif len(s_tmp[0])<= 19998:
                    MnsTS_corrected_for_tiff = MnsTS[1:10000] + MnsTS[10001:20000] + MnsTS[20001:30000] + MnsTS[30001:40000] + MnsTS[40001:50000] + MnsTS[50001:50001 + 3 * (
                                (len(MnsTS) - 50001) // 3)]

                elif len(s_tmp[0])<= 23331:
                    MnsTS_corrected_for_tiff = MnsTS[1:10000] + MnsTS[10001:20000] + MnsTS[20001:30000] + MnsTS[30001:40000] + MnsTS[40001:50000] + MnsTS[50001:60000] + MnsTS[60001:60001 + 3 * (
                                (len(MnsTS) - 60001) // 3)]

                elif len(s_tmp[0])<= 26664:
                    MnsTS_corrected_for_tiff = MnsTS[1:10000] + MnsTS[10001:20000] + MnsTS[20001:30000] + MnsTS[30001:40000] + MnsTS[40001:50000] + MnsTS[50001:60000] + MnsTS[60001:70000] + MnsTS[70001:70001 + 3 * (
                                (len(MnsTS) - 70001) // 3)]
                else:
                    raise Exception('file length not supported')

            # correction for delay between scope and camera
            MnsTS_delay_adjusted = [ts + 100 for ts in MnsTS_corrected_for_tiff]

            # correction of Mns timestamps for re sampling
            MnsTS_resampld = np.min(np.reshape(MnsTS_delay_adjusted, [-1, 3]), 1)  # used for binning later

            if len(MnsTS_resampld) == len(s_tmp[0]):
                print('Matched length')
            else:
                print([len(MnsTS_resampld), len(s_tmp[0])])

            # load tracker files and exclude bad frames and subsample
            TrackFrameCount, TrackTS, Xpos, Ypos, Sector, Shock = Pos_process(TrackerFileName[i_vid])
            Xpos_corr = [x if x + y > 0 else np.nan for x, y in zip(Xpos, Ypos)]
            Ypos_corr = [y if x + y > 0 else np.nan for x, y in zip(Xpos, Ypos)]

            # match position to physiology, right?
            matchMns = np.digitize(TrackTS, MnsTS_resampld)
            test_time_diff = [(u - MnsTS_resampld[v - 1]) < 130 for u, v in zip(TrackTS, matchMns)]
            Xpos_resampld = [np.nanmean([Xpos_corr[idx] for idx in np.where(matchMns == i)[0]]) # no usage of test time diff
                             for i in range(1, 1+len(MnsTS_resampld))]
            Ypos_resampld = [np.nanmean([Ypos_corr[idx] for idx in np.where(matchMns == i)[0]])
                             for i in range(1, 1+len(MnsTS_resampld))]
            Sector_resampld = [np.nanmean([Sector[idx] for idx in np.where(matchMns == i)[0]])
                               for i in range(1, 1+len(MnsTS_resampld))]
            Shock_resampld = [np.nanmean([Shock[idx] for idx in np.where(matchMns == i)[0]])
                              for i in range(1, 1+len(MnsTS_resampld))]

            # pre processing for rate maps computation
            pos_mns = [[x, y] for x, y in zip(Xpos_resampld, Ypos_resampld)]
            trck_and_speed = [~np.isnan(x[0]) for x in pos_mns]
            # pos_mns = [[Xpos[i - 1], Ypos[i - 1], v[i - 2]] for i in match_Mns]
            # trck_and_speed = [(x[2] > speed_tresh and ~np.isnan(x[0])) for x in pos_mns]

            # this is where temporal and PCs analyses diverge, here the nan's are kept
            X_tracked = [x[0] if check else np.nan for x, check in zip(pos_mns, trck_and_speed)]
            Y_tracked = [x[1] if check else np.nan for x, check in zip(pos_mns, trck_and_speed)]
            Sector_tracked = [x if check else np.nan for x, check in zip(Sector_resampld, trck_and_speed)]
            Shock_tracked = [x if check else np.nan for x, check in zip(Shock_resampld, trck_and_speed)]
            # do we need X_tracked anymore? would be useful for speed thresh?
            # s_tracked = [[spike if check else np.nan for spike, check in zip(s, trck_and_speed)] for s in s_tmp]

            # binning for calcium trace (w/o taking tracker into account, see older code for that)
            # AVC
            maxTimeBin = np.ceil(np.max(MnsTS_resampld)) + 1 + TimeBinSize_AVC
            minTimeBin = np.floor(np.min(MnsTS_resampld)) - 1
            TimeBins = np.arange(minTimeBin, maxTimeBin, TimeBinSize_AVC)
            MnsMatchTimeBins = np.digitize(MnsTS_resampld, TimeBins)

            sBinned_avc = [
                [np.mean([s[ii] for ii in np.where(MnsMatchTimeBins == i + 1)[0]]) for i in range(len(TimeBins) - 1)
                 if np.size(np.where(MnsMatchTimeBins == i + 1)) > 0]  # making sure there is no empty bins
                for s in s_tmp]

            if np.sum([np.sum([np.isnan(c) for c in Ct]) for Ct in sBinned_avc])>0:
                raise ValueError('AV matrix contains nans')

            XBinned = [np.nanmean([X_tracked[ii] for ii in np.where(MnsMatchTimeBins == i + 1)[0]])
                       for i in range(len(TimeBins) - 1)
                       if np.size(np.where(MnsMatchTimeBins == i + 1)) > 0]  # making sure there is no empty bins
            YBinned = [np.nanmean([Y_tracked[ii] for ii in np.where(MnsMatchTimeBins == i + 1)[0]])
                       for i in range(len(TimeBins) - 1)
                       if np.size(np.where(MnsMatchTimeBins == i + 1)) > 0]  # making sure there is no empty bins

            SectorBinned = [np.nanmean([Sector_tracked[ii] for ii in np.where(MnsMatchTimeBins == i + 1)[0]])
                            for i in range(len(TimeBins) - 1)
                            if np.size(np.where(MnsMatchTimeBins == i + 1)) > 0]  # making sure there is no empty bins
            ShockBinned = [np.nanmean([Shock_tracked[ii] for ii in np.where(MnsMatchTimeBins == i + 1)[0]])
                           for i in range(len(TimeBins) - 1)
                           if np.size(np.where(MnsMatchTimeBins == i + 1)) > 0]  # making sure there is no empty bins

            # TAU
            maxTimeBin = np.ceil(np.max(MnsTS_resampld)) + 1 + TimeBinSize_TAU
            TimeBins = np.arange(minTimeBin, maxTimeBin, TimeBinSize_TAU)
            MnsMatchTimeBins = np.digitize(MnsTS_resampld, TimeBins)
            sBinned_tau = [
                [np.mean([s[ii] for ii in np.where(MnsMatchTimeBins == i + 1)[0]]) for i in range(len(TimeBins) - 1)
                 if np.size(np.where(MnsMatchTimeBins == i + 1))>0] # making sure there is no empty bins
                for s in s_tmp]

            if np.sum([np.sum([np.isnan(c) for c in Ct]) for Ct in sBinned_avc]) > 0:
                raise ValueError('AV matrix _for tau calc_ contains nans')

            #if REC2 ATN: concatenate with previous and rerun analysis and replace in all lists
            if 'REC1' in Sessions[i_vid]:
                XBinnedPrev = XBinned
                YBinnedPrev = YBinned
                SectorBinnedPrev = SectorBinned
                ShockBinnedPrev = ShockBinned
                sBinned_tauPrev = sBinned_tau
                sBinned_avcPrev = sBinned_avc
                continue

            if 'REC2' in Sessions[i_vid]:
                XBinned = XBinnedPrev + XBinned
                YBinned = YBinnedPrev + YBinned
                SectorBinned = SectorBinnedPrev + SectorBinned
                ShockBinned = ShockBinnedPrev + ShockBinned
                sBinned_tau = np.concatenate([sBinned_tauPrev, sBinned_tau], 1)
                sBinned_avc = np.concatenate([sBinned_avcPrev, sBinned_avc], 1)

            # compute tau for this single recording only (as opposed to calculated on a pair of recordings)
            tauVecSingle = [scstat.kendalltau(c1, c2) for c1, c2 in itertools.combinations(sBinned_tau, 2)]

            # keep pairs id
            tauPairsSingle = [(x1, x2) for x1, x2 in itertools.combinations(range(len(sBinned_tau)), 2)]

            # to remove nearby cells
            aMaxLoc = [np.unravel_index(np.argmax(iA), [752, 480]) for iA in A_tmp] # find location of max
            aDist = [np.linalg.norm([x1 - x2 for x1,x2 in zip(u, v)]) for u,v in itertools.combinations(aMaxLoc, 2)] # compute euclidian distance
            isSeparated = [u>30 for u in aDist] # see if separated

            day_real_list.append(day)
            iDayList.append(iDay)
            day_count_list.append(DayCount)
            trainingList.append(Training)
            sess_list.append(Sessions[i_vid])
            S_conv_avc_list.append(sBinned_avc)
            S_conv_tau_list.append(sBinned_tau)
            trackList.append([XBinned,YBinned,SectorBinned,ShockBinned])
            # WHY remove nan's like that???? should have been fixed by the line making sure there are no empty bins
            # error will be raised in case there are nan's in sBinned_avc or sBinned_tau
            # S_conv_avc_list.append([[s_i for s_i in s if ~np.isnan(s_i)] for s in sBinned_avc])
            # # S_conv_tau_list.append([[s_i for s_i in s if ~np.isnan(s_i)] for s in sBinned_tau])
            # Sector_list.append([s_i for s_i in SectorBinned if ~np.isnan(s_i)])
            # Shock_list.append([s_i for s_i in ShockBinned if ~np.isnan(s_i)])
            Alist.append(A_tmp)
            aMaxLocList.append(aMaxLoc)
            tauVecSingleList.append(tauVecSingle)
            tauPairsSingleList.append(tauPairsSingle)
            isSeparatedList.append(isSeparated)
            atnSplitList.append(atnSplit)
            aDistList.append(aDist)

disttosave = [tauVecSingleList, aDistList]

with open(pathAll + 'analysisFiles_Simon/' +
          'Dist_' + '_' + animal +
          '.file', 'wb') as f:
    pickle.dump(disttosave, f, pickle.HIGHEST_PROTOCOL)

mdic_dist = {"Dist": disttosave}

pathSave_dist = pathAll + 'analysisFiles_Simon/' + 'Dist'+ '_' + animal +'.mat'
savemat(pathSave_dist, mdic_dist)

print('saved!')