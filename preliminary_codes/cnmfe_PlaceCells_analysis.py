import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import itertools
import numpy as np
import pandas as pd
from Pos_file_processing import *
import pickle

# FILE DESCRIPTION:
# Input: excel file,
# Compute regular and linearized place maps and determine place cells (using shuffled timestamps)
# compute map correlations, regular and linearized, for PCs and all Cells
# Save above mentioned metrics as well as matched position and calcium activity vectors

# to do :
# compute speed, discuss with Andre

expt = 'CircleSquare'
frame = 'Arena' # 'Room' or 'Arena'

animal = ['M35','M20', 'M29','M34','M19', 'M39']

sampRate = 10  # Hz

BinNb = 6  # number of bins for place maps
BinNB_Lin = 30
nRepeat = 1000  # for randomizing for calculating coherence and information content

pathAll = '/Volumes/GoogleDrive/My Drive/Fenton_Lab/CellPress_CircleSquare/MNS_Analysis/Data/' + expt + '/'

if expt == 'PlaceAvoidance':
    save_path = pathAll + 'analysisFiles_Simon/' + 'MapCorr_' + str(int(BinNb)) + '_Bins_' + frame + '_' + animal + '.file'
else:
    save_path = pathAll + 'analysisFiles_Simon/' + 'MapCorr_' + str(int(BinNb)) + '_Bins_' + animal + '.file'

# process excel file
ExcelPath = pathAll + 'FilenameMatch.xlsx'
xlread = pd.ExcelFile(ExcelPath)
xlAll = pd.read_excel(xlread, animal)
xlAll = xlAll[pd.notna(xlAll.videoname)]

Days = [int(x) for x in np.unique(list(xlAll.day))]
nsep = np.max(Days) - np.min(Days)

print(animal)
print(Days)
print(nRepeat)

# create variables
day_real_list = []  # day number in real time e.g. for Mon-Tues-Wed on 2 sequential weeks: [1,2,3,8,9,10]
day_count_list = [] # sequential day number e.g. for Mon-Tues-Wed on 2 sequential weeks: [1,2,3,4,5,6]
sess_list = []
rate_maps_list = []
LinMapsList = []
PC_idx_list = []
day_count_load_order = []
X_tracked_list = []
Y_tracked_list = []
s_tracked_list = []
OccList = []
cohList = []
pInfoList = []
trainingList = []

for day in Days:

    print(day)

    # extract data
    xlOut, fileNames, physio = loadPhysiology(day, animal, xlAll, pathAll, expt, frame)
    VideoName, DayCount, Sessions, doAnalysis, Training, atn = xlOut
    MnsTSFileName, TrackerFileName = fileNames
    S_all = physio[0]
    rec_idx = physio[2]

    # in case not all day in Days are loaded
    day_count_load_order.append(DayCount)

    for i_vid, vid in enumerate(VideoName):

        if doAnalysis[i_vid]>0 and 'HMC' not in Sessions[i_vid]:
            # extract physiology from longer file
            s_tmp = np.transpose(S_all[rec_idx[i_vid][0]:rec_idx[i_vid][1]])

            # load video time file and adjust for delay to start MNS and conversion
            MnsFrameCount, MnsTS = read_timestamps(MnsTSFileName[i_vid])

            # correction for lost frames during tiff conversion and division by 3
            if len(s_tmp[0])<=3333:
                MnsTS_corrected_for_tiff = MnsTS[1:1 + 3 * ((len(MnsTS) - 1) // 3)]

            elif len(s_tmp[0])<=6666:
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
            MnsTS_resampld = np.min(np.reshape(MnsTS_delay_adjusted, [-1, 3]), 1)

            if len(MnsTS_resampld) == len(s_tmp[0]):
                print('Matched length')
            else:
                print([len(MnsTS_resampld), len(s_tmp[0])])

            # load tracker files and exclude bad frames and subsample
            TrackFrameCount, TrackTS, Xpos, Ypos, Sector, Shock = Pos_process(TrackerFileName[i_vid])
            Xpos_corr = [x if x + y > 0 else np.nan for x, y in zip(Xpos, Ypos)]
            Ypos_corr = [y if x + y > 0 else np.nan for x, y in zip(Xpos, Ypos)]

            # skip this recording if no good frame
            if np.all([(np.isnan(x) and np.isnan(y)) for x,y in zip(Xpos_corr,Ypos_corr)]):
                print('RECORDING CORRUPTED')
                print(day)
                print(Sessions[i_vid])
                print(TrackerFileName[i_vid])
                continue

            # match position to physiology
            matchMns = np.digitize(TrackTS, MnsTS_resampld)
            test_time_diff = [(u - MnsTS_resampld[v - 1]) < 130 for u, v in zip(TrackTS, matchMns)]
            Xpos_resampld = [np.nanmean([Xpos_corr[idx] for idx in np.where(matchMns == i)[0] if test_time_diff[idx]])
                             for i in range(1, len(MnsTS_resampld))]
            Ypos_resampld = [np.nanmean([Ypos_corr[idx] for idx in np.where(matchMns == i)[0] if test_time_diff[idx]])
                             for i in range(1, len(MnsTS_resampld))]

            # compute animal speed

            # pre processing for rate maps computation
            pos_mns = [[x, y] for x, y in zip(Xpos_resampld, Ypos_resampld)]
            trck_and_speed = [~np.isnan(x[0]) for x in pos_mns]
            # pos_mns = [[Xpos[i - 1], Ypos[i - 1], v[i - 2]] for i in match_Mns]
            # trck_and_speed = [(x[2] > speed_tresh and ~np.isnan(x[0])) for x in pos_mns]

            # this is where temporal and PCs analyses diverge, here the nan's are NOT kept
            X_tracked = [x[0] for x, check in zip(pos_mns, trck_and_speed) if check]
            Y_tracked = [x[1] for x, check in zip(pos_mns, trck_and_speed) if check]
            s_tracked = [[spike for spike, check in zip(s, trck_and_speed) if check] for s in s_tmp]

            #if REC2 ATN: concatenate with previous and rerun analysis and replace in all lists
            if 'REC1' in Sessions[i_vid]:
                xTrackedPrev = X_tracked
                yTrackedPrev = Y_tracked
                sTrackedPrev = s_tracked
                continue

            if 'REC2' in Sessions[i_vid]:
                X_tracked = xTrackedPrev + X_tracked
                Y_tracked = yTrackedPrev + Y_tracked
                s_tracked = np.concatenate([sTrackedPrev, s_tracked], 1)

            # COMPUTING THE RATE MAPS
            occ = np.histogram2d(X_tracked, Y_tracked, bins=BinNb)[0]/sampRate
            rate_maps_temp = [np.histogram2d(X_tracked, Y_tracked, bins=BinNb, weights=s)[0]
                              / occ for s in s_tracked]

            # linearizing with polar coordinates
            Ycenter = Y_tracked - (np.nanmax(Y_tracked) + np.nanmin(Y_tracked)) / 2
            Xcenter = X_tracked - (np.nanmax(X_tracked) + np.nanmin(X_tracked)) / 2
            alpha = np.arctan2(Ycenter, Xcenter)

            L_occ = np.histogram(alpha, bins=BinNB_Lin)[0]
            LinMapsTemp = [np.histogram(alpha, bins=BinNB_Lin, weights=s)[0] / L_occ
                 for s in s_tracked]

            # analyze rate maps
            p_info = [place_info_content(occ, m) for m in rate_maps_temp]
            coh = [neighbor_sum(m) for m in rate_maps_temp]

            # Randomizing
            sTrackedCopy = s_tracked.copy()
            sTrackedTimePairs = list(zip(*sTrackedCopy))

            pInfoShuffled = []
            cohShuffled = []
            for iRand in range(nRepeat):
                np.random.shuffle(sTrackedTimePairs)
                sShuffled = list(zip(*sTrackedTimePairs))

                RateMapsShuffled = [np.histogram2d(X_tracked, Y_tracked, bins=BinNb, weights=s)[0]
                                    / occ for s in sShuffled]
                pInfoShuffled.append([place_info_content(occ, m) for m in RateMapsShuffled])

                cohShuffled.append([neighbor_sum(m) for m in RateMapsShuffled])

            pInfoTest = []
            for p, z in zip(p_info, zip(*pInfoShuffled)):
                Out = [p, (p - np.mean(z)) / np.std(z)]
                pInfoTest.append([Out, Out[0]>1 and Out[1]>1.96])

            cohTest = []
            for p, z in zip(coh, zip(*cohShuffled)):
                Out = [p, (p - np.mean(z)) / np.std(z)]
                testTmp = Out[0] > 0.5 and Out[1] > 1.96
                cohTest.append([Out, testTmp])

            # determine whether cells are place cells
            PC_test = [c[0][1] > 1.96 and p[0][1] > 1.96 for c, p in zip(cohTest, pInfoTest)]

            day_real_list.append(day)
            day_count_list.append(DayCount)
            trainingList.append(Training)
            sess_list.append(Sessions[i_vid])
            rate_maps_list.append(rate_maps_temp)
            LinMapsList.append(LinMapsTemp)
            # S_conv_avc_list.append([[s_i for s_i in s if ~np.isnan(s_i)] for s in sBinned_avc])
            # S_conv_tau_list.append([[s_i for s_i in s if ~np.isnan(s_i)] for s in sBinned_tau])
            PC_idx_list.append(PC_test)
            X_tracked_list.append(X_tracked)
            Y_tracked_list.append(Y_tracked) # 10Hz
            s_tracked_list.append(s_tracked)
            cohList.append(cohTest)
            pInfoList.append(pInfoTest)
            OccList.append(occ)


with open(pathAll + 'analysisFiles_Simon/' + 'AssignmentsByPairs_' + animal + '.file', 'rb') as f:
    [assignments, [assignmentsByPairs, idPairs]] = pickle.load(f)


day_sep = []
same_diff = []
is_HMC = []
day_count_pair = []
day_real_pair = []
sess_pair = []
cellCount = []
trainingPair = []

maps_corr = []
maps_corr_pc_all = [] # a place cell in both environment
maps_corr_pc_any = [] # a place cell in either environment
maps_corr_pc_not = [] # a place cell in neither environment
LinMapsCorr = []
LinMapsCorrPC = []

# for displaying analysis progress on screen
trackCompletion = 0
trackIncrement = 0.1
tot_loop = np.sum([1 for i_comb in itertools.combinations(range(len(day_count_list)), 2)])

pairsSkipped = []

for n_comb, i_comb in enumerate(itertools.combinations(range(len(day_count_list)), 2)):

    # for displaying analysis progress on screen
    if n_comb/tot_loop>trackCompletion:
        trackCompletion += trackIncrement
        print(n_comb/tot_loop)

    day_real = tuple([day_real_list[i] for i in i_comb])
    day_count = tuple([int(day_count_list[i] - 1) for i in i_comb])
    sess_id = tuple([sess_list[i] for i in i_comb])

    # recording number in case not all recording are loaded, otherwise identical to day_count
    day_count_loaded = tuple([np.where([d-1 == i_d for d in day_count_load_order])[0][0] for i_d in day_count])

    is_HMC_tmp = (sess_id[0][:3] == 'HMC') | (sess_id[1][:3] == 'HMC')

    # former matching algorithm
    # both_dual = np.all(~np.isnan(assignments[:, day_count_loaded]), 1)
    # CellList2 = [[n_cell for n_cell in assignments[both_dual, day_id].astype(int)]
    #             for day_id in day_count_loaded]

    # matching cells across days
    if np.diff(day_count_loaded)[0]==0:
        nMatched = len(s_tracked_list[i_comb[0]])
        CellList = [[*range(nMatched)] for x in range(2)]
        f1Score = 1
    else:
        findAssign = np.where([day_count_loaded == iPair for iPair in idPairs])[0][0] #will raise error if pair is not found
        assign = assignmentsByPairs[findAssign]
        nMatched = len(assign[0])
        CellList = [assign[0], assign[1]]
        f1Score = assign[4]['f1_score']

    if f1Score<0.3:
        print(i_comb)
        print(day_real)
        print(sess_id)
        print('SKIPPING - LOW f1 - '+ str(f1Score))
        pairsSkipped.append(day_real)

    if (not is_HMC_tmp) and (f1Score>0.3):  # & (day_idx[0]>2)& (day_idx[1]>2):
        is_HMC.append(is_HMC_tmp)
        day_count_pair.append(day_count)
        day_real_pair.append(day_real)
        sess_pair.append(sess_id)
        day_sep.append(np.abs(day_real[1] - day_real[0]))
        same_diff.append(sess_id[0][:3] == sess_id[1][:3])
        trainingPair.append(tuple([trainingList[i] for i in i_comb]))


        mask = ~np.any([np.isnan(rate_maps_list[idx][0].flatten()) for idx in i_comb], 0)
        # maskLin = ~np.any([np.isnan(LinMapsList[idx][0]) for idx in i_comb], 0)
        # print('keeping {} pixels out of {} '.format(sum(mask), BinNb ** 2))

        Rmaps = [[rate_maps_list[rec_id][n_cell]
                  for n_cell in CellList[i_rec]] for i_rec, rec_id in enumerate(i_comb)]

        LinRmaps = [[LinMapsList[rec_id][n_cell]
                    for n_cell in CellList[i_rec]] for i_rec, rec_id in enumerate(i_comb)]

        PC_assign = [[PC_idx_list[rec_id][n_cell]
                      for n_cell in CellList[i_rec]] for i_rec, rec_id in enumerate(i_comb)]

        PC_any = np.any(PC_assign, 0)
        PC_all = np.all(PC_assign, 0)
        PC_not = [not v for v in PC_any]

        cellCount.append([nMatched, np.sum(PC_any), np.sum(PC_all)])

        Rmaps_pc_any = [[m for m, v in zip(maps, PC_any) if v]
                        for maps in Rmaps]
        Rmaps_pc_all = [[m for m, v in zip(maps, PC_all) if v]
                        for maps in Rmaps]
        Rmaps_pc_not = [[m for m, v in zip(maps, PC_not) if v]
                        for maps in Rmaps]
        LinRmapsPc = [[lm for lm, v in zip(maps, PC_all) if v]
                      for maps in LinRmaps]

        # computing correlations
        maps_corr.append([np.corrcoef([m.flatten()[mask] for m in M])[0, 1] for M in zip(*Rmaps)])
        maps_corr_pc_all.append([np.corrcoef([m.flatten()[mask] for m in M])[0, 1] for M in zip(*Rmaps_pc_all)])
        maps_corr_pc_any.append([np.corrcoef([m.flatten()[mask] for m in M])[0, 1] for M in zip(*Rmaps_pc_any)])
        maps_corr_pc_not.append([np.corrcoef([m.flatten()[mask] for m in M])[0, 1] for M in zip(*Rmaps_pc_not)])

        # compute ensemble correlation
        LinRmapsNorm = [np.transpose(np.transpose(l) / np.max(l, 1)) for l in LinRmaps]
        # exclude cell with no activity
        maskLin = ~np.any(np.any(np.isnan(LinRmapsNorm), 2), 0)
        # save
        LinMapsCorr.append(np.corrcoef([LM[maskLin].flatten() for LM in LinRmapsNorm])[0, 1])

        # same for PCs
        if np.sum(PC_all) > 1:
            LinRmapsPcNorm = [np.transpose(np.transpose(l) / np.max(l, 1)) for l in LinRmapsPc]
            maskLinPc = ~np.any(np.any(np.isnan(LinRmapsPcNorm), 2), 0)
            LinMapsCorrPC.append(np.corrcoef([LM[maskLinPc].flatten() for LM in LinRmapsPcNorm])[0, 1])
        else:
            LinMapsCorrPC.append(np.nan)


#save:
# [result of correlation]
# [analysis of each recording]
# [X_tracked_list, Y_tracked_list, s_tracked_list]
# [nRepeat]

to_save = [ [maps_corr, [maps_corr_pc_any, maps_corr_pc_all, maps_corr_pc_not], day_count_pair, day_real_pair, sess_pair, day_sep, same_diff,
            LinMapsCorr, LinMapsCorrPC, cellCount, trainingPair],
           [day_real_list, day_count_list, sess_list, rate_maps_list, PC_idx_list, LinMapsList, OccList,
            cohList, pInfoList, trainingList],
           [X_tracked_list, Y_tracked_list, s_tracked_list],
           [nRepeat]]

with open(save_path, 'wb') as f:
    pickle.dump(to_save, f, pickle.HIGHEST_PROTOCOL)

# SAME VS DIFF
mapC_env_all = [[[m for idx, v in enumerate(same_diff) for m in maps_corr[idx] if v == ibool if day_sep[idx] == i_day]
            for i_day in range(nsep)]
            for ibool in [True, False]]

mapC_env_PC = [[[m for idx, v in enumerate(same_diff) for m in maps_corr_pc_any[idx] if v == ibool if day_sep[idx] == i_day]
            for i_day in range(nsep)]
            for ibool in [True, False]]

# f, a = plt.subplots(1, 2)
#
# a[0].hist([m for m in mapC_env_all[0][0] if not np.isnan(m)], 100, histtype='step', cumulative=False, density=True)
# a[0].hist([m for m in mapC_env_all[1][0] if not np.isnan(m)], 100, histtype='step', cumulative=False, density=True)
#
# a[1].hist([m for m in mapC_env_PC[0][0] if not np.isnan(m)], 100, histtype='step', cumulative=False, density=True)
# a[1].hist([m for m in mapC_env_PC[1][0] if not np.isnan(m)], 100, histtype='step', cumulative=False, density=True)
#
# f, a = plt.subplots(1, 2)
#
# a[0].hist([m for m in mapC_env_all[0][7] if not np.isnan(m)], 100, histtype='step', cumulative=False, density=True)
# a[0].hist([m for m in mapC_env_all[1][7] if not np.isnan(m)], 100, histtype='step', cumulative=False, density=True)
#
# a[1].hist([m for m in mapC_env_PC[0][7] if not np.isnan(m)], 100, histtype='step', cumulative=False, density=True)
# a[1].hist([m for m in mapC_env_PC[1][7] if not np.isnan(m)], 100, histtype='step', cumulative=False, density=True)

print('done')

print(str(len(pairsSkipped)) + ' rec pairs skipped')
print(str(len(set(pairsSkipped))) + ' day pairs skipped')
print(set(pairsSkipped))

plt.show(block=True)