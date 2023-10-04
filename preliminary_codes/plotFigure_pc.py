import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Pos_file_processing import *
import pickle
import itertools

path_all = '/Volumes/GoogleDrive/My Drive/Fenton_Lab/CellPress_CircleSquare/MNS_Analysis/Data/CircleSquare/'
pathFigure = '/Volumes/GoogleDrive/My Drive/Fenton_Lab/CellPress_CircleSquare/MNS_Analysis/Figures/CircleSquare/savedFigs_Simon/'
BinNb = 10  # number of bins for place maps
smoothFilter = 1/3

AnimalList = ['M19']#, 'M29','M20','M35','M19']
nAni = len(AnimalList)

for iAni, animal in enumerate(AnimalList):

    print(animal)

    loadpath = path_all + 'analysisFiles_Simon/' + 'MapCorr_' + str(int(BinNb)) + '_Bins_' + animal + '.file'

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

    # sTracked = output[2][2]
    output = None

    # pInfo:  [[p, (p - np.mean(shuffle)) / np.std(shuffle)], [>1 and >1.96]

    day = 7
    print(day)

    idxDay = np.where([d == day for d in dayCountList])[0]
    sess = [sessList[x] for x in idxDay]
    maps = [rateMapList[x] for x in idxDay]
    PcIdx = [PcIdxList[x] for x in idxDay]
    coh = [cohList[x] for x in idxDay]
    pInfo = [pInfoList[x] for x in idxDay]

    # get rec idx for same environments
    pairSame = [[idx for idx, s in enumerate(sess) if strEnv in s] for strEnv in ['RCT','CYL']]

    # check if PC in both env
    # isPcEnv = np.any([np.all([PcIdx[idx] for idx in iIdxPair], 0) for iIdxPair in pairSame], 0)
    isPcEnv = np.any([np.any([[not pc for pc in PcIdx[idx]] for idx in iIdxPair], 0) for iIdxPair in pairSame], 0)
    isPcIdx = np.where(isPcEnv)[0]

    cmap = mpl.cm.viridis
    cmap.set_bad('w', 1.)

    for iRec in range(len(idxDay)):

        f,a = plt.subplots(8,15, figsize=(16, 14))

        for ax, idx in zip(a.flat, isPcIdx):

            m = maps[iRec][idx]
            p = pInfo[iRec][idx]
            c = coh[iRec][idx]
            pc = PcIdx[iRec][idx]

            ax.imshow(nan_gaussian_filter(m, smoothFilter), cmap=cmap)
            ax.axis('off')
            cTest = 'x' if c[0][1] > 1.96 else 'o'
            pTest = 'x' if p[0][1] > 1.96 else 'o'
            pc = 'P' if pc else 'N'
            ax.set_title('%i%s %.2f' % (idx,pc,c[0][0]) + cTest + '_' + '%.2f' % (p[0][0]) + pTest,
                         fontsize=6)

        strFig = animal + 'NPC_day' + str(day) + '_' + sess[iRec]

        f.suptitle(strFig)

        f.savefig(pathFigure + strFig + '.eps', format='eps', dpi=300)

        mapCorr = [[x for x in range(len(maps[0]))]]
        sessPairList = ['CellID']

    for iPair in itertools.combinations(range(len(idxDay)), 2):
        sessPair = '%s_%s' % (sess[iPair[0]], sess[iPair[1]])
        print(sessPair)
        mapPair = [maps[iRec] for iRec in iPair]
        mask = ~np.any([np.isnan(m[0].flatten()) for m in maps], 0)

        mapCorr.append([np.corrcoef([m.flatten()[mask] for m in M])[0, 1] for M in zip(*mapPair)])
        sessPairList.append(sessPair)

    df = pd.DataFrame(np.transpose(mapCorr), columns=sessPairList)
    df.to_csv(pathFigure + 'data' + strFig + '.csv')




print('done done done')

plt.show(block=True)