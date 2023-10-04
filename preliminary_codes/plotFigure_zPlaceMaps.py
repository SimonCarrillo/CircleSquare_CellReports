import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Pos_file_processing import *
import pickle
import itertools

pathAll = '/Volumes/GoogleDrive/My Drive/Fenton_Lab/CellPress_CircleSquare/MNS_Analysis/Data/CircleSquare/'
pathFigure = '/Volumes/GoogleDrive/My Drive/Fenton_Lab/CellPress_CircleSquare/MNS_Analysis/Figures/CircleSquare/savedFigs_Simon/z/'
BinNb = 10  # number of bins for place maps
smoothFilter = 1/3

AnimalList = ['M20']#, 'M29','M20','M35','M19']
nAni = len(AnimalList)

randomize = False
linearZ = True

strRand = 'Rand' if randomize else ''
strLin = 'linear' if linearZ else ''
strSaveName = strRand + 'Z' + strLin
TimeBinSize_AVC = 10 * 1000
TimeBinSize_TAU = 1 * 1000

for iAni, animal in enumerate(AnimalList):

    print(animal)

    loadpath = pathAll + 'analysisFiles_Simon/' + 'MapCorr_' + str(int(BinNb)) + '_Bins_' + animal + '.file'

    with open(loadpath, 'rb') as f:
        output = pickle.load(f)

    # load analyzed data
    # [maps_corr, maps_corr_pc, day_count_pair, day_real_pair, sess_pair, day_sep, same_diff,
    #  LinMapsCorr, LinMapsCorrPC, cellCount] = output[0]
    # [maps_corr_pc_any, maps_corr_pc_all, maps_corr_pc_not] = maps_corr_pc
    # [day_real_list, dayCountList, sessList, rateMapList, PC_idx_list, LinMapsList, OccList,
    #  cohList, pInfoList] = output[1]

    dayCountList = output[1][1]
    sessList = output[1][2]
    rateMapList = output[1][3]
    PcIdxList = output[1][4]
    cohList = output[1][7]
    pInfoList = output[1][8]
    output = None
    # pInfo:  [[p, (p - np.mean(shuffle)) / np.std(shuffle)], [>1 and >1.96]

    zLoadpath = pathAll + 'analysisFiles_Simon/' + strSaveName + '_List_' + str(int(TimeBinSize_AVC)) + 'msAV_' + str(int(TimeBinSize_TAU)) + 'msKT_CircleSquare_' + animal + '.file'
    with open(zLoadpath, 'rb') as f:
        output = pickle.load(f)

    zBinnedRateMapsList, zRateMapsList = output[1][1]
    zDayCountList = output[0][1]
    zSessList = output[0][2]
    output = None

    day = 7
    print(day)

    idxDay = np.where([(d == day) and ('HMC' not in sess) for d, sess in zip(dayCountList, sessList)])[0]
    sess = [sessList[x] for x in idxDay]
    maps = [rateMapList[x] for x in idxDay]
    PcIdx = [PcIdxList[x] for x in idxDay]
    coh = [cohList[x] for x in idxDay]
    pInfo = [pInfoList[x] for x in idxDay]

    # get rec idx for same environments
    pairSame = [[idx for idx, s in enumerate(sess) if strEnv in s] for strEnv in ['RCT','CYL']]

    # check if PC in both env
    isPcEnv = np.any([np.all([PcIdx[idx] for idx in iIdxPair], 0) for iIdxPair in pairSame], 0)
    # isPcEnv = np.any([np.any([[not pc for pc in PcIdx[idx]] for idx in iIdxPair], 0) for iIdxPair in pairSame], 0)
    isPcIdx = np.where(isPcEnv)[0]

    # same for the zScored data
    zIdxDay = np.where([(d == day) and ('HMC' not in sess) for d, sess in zip(zDayCountList, zSessList)])[0]
    zBinnedMaps = [zBinnedRateMapsList[x] for x in zIdxDay]
    zMaps = [zRateMapsList[x] for x in zIdxDay]
    zSess = [zSessList[x] for x in zIdxDay]
    zPairSame = [[idx for idx, s in enumerate(zSess) if strEnv in s] for strEnv in ['RCT', 'CYL']]

    cmap = mpl.cm.viridis
    cmap.set_bad('w', 1.)

    for iRec in range(len(idxDay)):

        f,a = plt.subplots(12,15, figsize=(16, 14))

        for iCol in range(15):

            for iRow in range(4):

                if iCol + 15*iRow + 30>len(isPcIdx): continue
                idx = isPcIdx[iCol + 15*iRow + 20]  #idx = isPcIdx[iCol + 15*iRow + 45]

                m = maps[iRec][idx]
                p = pInfo[iRec][idx]
                c = coh[iRec][idx]
                pc = PcIdx[iRec][idx]

                a[0+3*iRow][iCol].imshow(nan_gaussian_filter(m, smoothFilter), cmap=cmap, vmin=-0.3, vmax=0.3)
                a[0+3*iRow][iCol].axis('off')
                cTest = 'x' if c[0][1] > 1.96 else 'o'
                pTest = 'x' if p[0][1] > 1.96 else 'o'
                pc = 'P' if pc else 'N'
                a[0+3*iRow][iCol].set_title('%i%s %.2f' % (idx,pc,c[0][0]) + cTest + '_' + '%.2f' % (p[0][0]) + pTest,
                             fontsize=6)

                a[1+3*iRow][iCol].imshow(zBinnedMaps[iRec][idx], vmin=-0.3, vmax=0.3)
                a[2+3*iRow][iCol].imshow(zMaps[iRec][idx], vmin=-0.3, vmax=0.3)

        strFig = animal + 'PC_day' + str(day) + '_' + sess[iRec]

        f.suptitle(strFig)
        f.savefig(pathFigure + strFig + '_2.eps', format='eps', dpi=300)




print('done done done')

plt.show(block=True)