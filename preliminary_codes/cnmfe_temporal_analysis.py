import matplotlib as mpl
mpl.use('TkAgg')

import itertools
import numpy as np
import scipy.stats as scstat
import pandas as pd
from Pos_file_processing import *
import pickle
import time

# analyze recordings then save
# temporal analysis for both AV (activity vectors) and Tau (correlation Kendall-Tau)

expt = 'PlaceAvoidance'  # here you can change to ATN (neophobia) or PlaceAvoidance (APA)

animal = 'M20' #['M39','M35','M34','M29','M19','M20']

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
tauPairsSingleList = []
isSeparatedList = []
trainingList = []
atnSplitList = []

for iDay, day in enumerate(Days):

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

if saveOutput:
    print('saving...')
    listsToSave = [[day_real_list, day_count_list, sess_list, S_conv_avc_list, S_conv_tau_list, isSeparatedList,
                    tauPairsSingleList, tauVecSingleList, trainingList],
                   [trackList, [], aMaxLocList],
                   [atnSplitList, iDayList]]

    with open(pathAll + 'analysisFiles_Simon/' +
              'TempCorr_List_' + str(int(TimeBinSize_AVC)) +
              'msAV_' + str(int(TimeBinSize_TAU)) + 'msKT' +
              '_' + expt + '_' + animal +
              '.file', 'wb') as f:
        pickle.dump(listsToSave, f, pickle.HIGHEST_PROTOCOL)
    print('saved!')

else:
    print('NOT SAVING')


if CombComputation:

    with open(pathAll + 'analysisFiles_Simon/' + 'AssignmentsByPairs_' + animal + '.file', 'rb') as f:
        [assignments, [assignmentsByPairs, idPairs]] = pickle.load(f)


    day_sep = []
    same_diff = []
    is_HMC = []
    day_count_pair = []
    day_real_pair = []
    sess_pair = []
    trainingPair = []

    cv = []
    tauVecAllComb = []
    isSeparatedAllComb = []
    ratePairAllComb = []

    trackCompletion = 0
    trackIncrement = 0.1
    totLoop = np.sum([1 for iComb in itertools.combinations(range(len(day_count_list)), 2)])

    for nComb, iComb in enumerate(itertools.combinations(range(len(day_count_list)), 2)):

        if nComb/totLoop>trackCompletion:
            trackCompletion += trackIncrement
            print(nComb/totLoop)

        day_real = tuple([day_real_list[i] for i in iComb])
        day_count = tuple([int(day_count_list[i] - 1) for i in iComb])
        sess_id = tuple([sess_list[i] for i in iComb])
        iDayPair = tuple([iDayList[i] for i in iComb])  # for finding cell pair assignments

        # not needed anymore because assignments are done separately
        # recording number in case not all recording are loaded, otherwise identical to day_count
        # day_count_loaded = tuple([np.where([d-1 == i_d for d in day_count_load_order])[0][0] for i_d in day_count])

        is_HMC_tmp = (sess_id[0][:3] == 'HMC') | (sess_id[1][:3] == 'HMC')

        # former method
        # both_dual = np.all(~np.isnan(assignments[:, day_count_loaded]), 1)
        # #both_dual = np.all(~np.isnan(assignments[:, day_count]), 1)
        # n_matched = sum(both_dual)
        # CellList = [[n_cell for n_cell in assignments[both_dual, day_id].astype(int)]
        #             for day_id in day_count_loaded]

        if np.diff(iDayPair)[0]==0:
            nMatched = len(S_conv_tau_list[iComb[0]])
            CellList = [[*range(nMatched)] for x in range(2)]
            f1Score = 1
        else:
            findAssign = np.where([iDayPair == iPair for iPair in idPairs])[0][0] #will raise error if pair is not found
            assign = assignmentsByPairs[findAssign]
            nMatched = len(assign[0])
            CellList = [assign[0], assign[1]]
            f1Score = assign[4]['f1_score']

        if f1Score<0.3:
            print(iComb)
            print(day_real)
            print(sess_id)
            print('SKIPPING - LOW f1 - '+ str(f1Score))
            # pairsSkipped.append(day_real)

            print(nMatched)

        if (not is_HMC_tmp) and (nMatched > min_cell_nb) and (f1Score>0.3):  # & (day_idx[0]>2)& (day_idx[1]>2):
            is_HMC.append(is_HMC_tmp)
            day_count_pair.append(day_count)
            day_real_pair.append(day_real)
            sess_pair.append(sess_id)
            day_sep.append(np.abs(day_real[1] - day_real[0]))
            same_diff.append(sess_id[0][:3] == sess_id[1][:3])
            trainingPair.append(tuple([trainingList[i] for i in iComb]))

            # mask = ~np.any([np.isnan(rate_maps_list[idx][0].flatten()) for idx in iComb], 0)
            # print('keeping {} pixels out of{} '.format(sum(mask), BinSize ** 2))
            #
            # Rmaps = [[rate_maps_list[rec_id][n_cell]
            #           for n_cell in CellList[i_rec]] for i_rec, rec_id in enumerate(iComb)]
            #
            # PC_assign = [[PC_idx_list[rec_id][n_cell]
            #               for n_cell in CellList[i_rec]] for i_rec, rec_id in enumerate(iComb)]
            #
            # PC_any = np.any(PC_assign, 0)
            # Rmaps_pc = [[m for m, v in zip(maps, PC_any) if v]
            #             for maps in Rmaps]
            #
            # maps_corr.append([np.corrcoef([m.flatten()[mask] for m in M])[0, 1] for M in zip(*Rmaps)])

            # Correlations
            # extract traces for matched cells
            C_cv = [np.array([S_conv_avc_list[rec_id][n_cell] for n_cell in CellList[i_rec]])
                     for i_rec, rec_id in enumerate(iComb)]

            C_tau = [np.array(
                [S_conv_tau_list[rec_id][n_cell] for n_cell in CellList[i_rec]])
                for i_rec, rec_id in enumerate(iComb)]

            # C_tau_pc = [[c for c, v in zip(c_tau, PC_any) if v]
            #                 for c_tau in C_tau]

            # compute AV correlation
            cv_tmp = np.corrcoef(C_cv[0].transpose(), C_cv[1].transpose())
            # cv_cross_tmp = cv_tmp[:C_cv[0].shape[1], -C_cv[1].shape[1]:].flatten()
            # cv_same_tmp = cv_tmp[:C_cv[0].shape[1], :C_cv[0].shape[1]].flatten()

            # to remove pairs closeby
            # aPairComb = [[Alist[rec][iCell] for iCell in CellList[iRec]] for iRec, rec in enumerate(iComb)]
            # aMaxLocComb = [[np.unravel_index(np.argmax(iA), [752, 480]) for iA in aRec] for aRec in aPairComb] # find location of max
            aMaxLocComb = [[aMaxLocList[rec][iCell] for iCell in CellList[iRec]] for iRec, rec in enumerate(iComb)]
            aDistComb = [[np.linalg.norm([x1 - x2 for x1,x2 in zip(u, v)]) for u,v in itertools.combinations(iAMax, 2)] for iAMax in aMaxLocComb] # compute euclidian distance
            isSeparatedComb = [u>30 and v>30 for u,v in zip(*aDistComb)] # see if separated

            # compute tau for all pairs
            # start = time.time()
            tauVecComb = np.stack([[scstat.kendalltau(c1, c2)
                                   for c1,c2 in itertools.combinations(C, 2)]
                                  for C in C_tau])

            ratePair = [[[np.mean([s > 0 for s in c1]), np.mean([s > 0 for s in c2])]
                         for c1,c2 in itertools.combinations(C, 2)] for C in C_tau]

            # tauVecComb2 = np.stack([[scstat.kendalltau(c1[(c1+c2)>0], c2[(c1+c2)>0])
            #                        for c1,c2 in itertools.combinations(C, 2)]
            #                       for C in C_tau])

            # end = time.time(); print(end - start)

            # cv.append(cv_tmp)
            cv.append([cv_tmp, [C_cv[0].shape[1], C_cv[1].shape[1]]])
            # cv[1].append([c for c in cv_same_tmp if ~np.isnan(c)])
            # cv[2].append([c for c in cv_cross_tmp if ~np.isnan(c)])
            tauVecAllComb.append(tauVecComb)
            isSeparatedAllComb.append(isSeparatedComb)
            ratePairAllComb.append(ratePair)
            # tau_vec_sep.append(tau_vec_tmp_separated)

if saveOutput2:
    print('saving...')

    cvToSave = [cv, day_count_pair, day_real_pair, sess_pair, day_sep, same_diff, trainingPair]

    tauToSave = [tauVecAllComb, day_count_pair, day_real_pair, sess_pair, day_sep, same_diff,
                 isSeparatedAllComb, trainingPair, ratePairAllComb]

    combToSave = [cv, tauVecAllComb, day_count_pair, day_real_pair, sess_pair, day_sep, same_diff,
                  isSeparatedAllComb, trainingPair, ratePairAllComb]

    with open(pathAll + 'analysisFiles_Simon/' +
              'TempCorr_KT_Comb_' + str(int(TimeBinSize_TAU)) + 'msKT_' +
               expt + '_' + animal +
              '.file', 'wb') as f:
        pickle.dump(tauToSave, f, pickle.HIGHEST_PROTOCOL)

    with open(pathAll + 'analysisFiles_Simon/' +
              'TempCorr_AV_Comb_' + str(int(TimeBinSize_AVC)) + 'msAV_' +
               expt + '_' + animal +
              '.file', 'wb') as f:
        pickle.dump(cvToSave, f, pickle.HIGHEST_PROTOCOL)

    with open(pathAll + 'analysisFiles_Simon/' +
              'TempCorr_Comb_' + str(int(TimeBinSize_AVC)) +
              'msAV_' + str(int(TimeBinSize_TAU)) + 'msKT' +
              '_' + expt + '_' + animal +
              '.file', 'wb') as f:
        pickle.dump(combToSave, f, pickle.HIGHEST_PROTOCOL)

    print('saved!')

else:
    print('NOT SAVING')

print('done done done')

