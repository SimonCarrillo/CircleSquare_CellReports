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
TimeBinSize_AVC = 1 * 1000
TimeBinSize_TAU = 1 * 1000

saveFig = True
doPTI = True

if doPTI:
    pathFigure = '/Volumes/GoogleDrive/My Drive/Fenton_Lab/CellPress_CircleSquare/MNS_Analysis/Figures/CircleSquare/savedFigs_Simon/z/'
else:
    pathFigure = '/Volumes/GoogleDrive/My Drive/Fenton_Lab/CellPress_CircleSquare/MNS_Analysis/Figures/CircleSquare/savedFigs_Simon/'

animalList = ['M19']
nAni = len(animalList)


if doPTI:
    loadStr = 'Zlinear'
else:
    loadStr = 'TempCorr'

f,a = plt.subplots(2,3)

for iAni, animal in enumerate(animalList):

    print(animal)

    # load data
    with open(path_all + 'analysisFiles_Simon/' + loadStr + '_List_' +
              str(int(TimeBinSize_AVC)) + 'msAV_' + str(int(TimeBinSize_TAU)) + 'msKT'
              '_CircleSquare_' + animal +'.file', 'rb') as fFile:
        outputList = pickle.load(fFile)

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

    # compute average of first half for each cell for each recording
    meanFirstHalf = [[np.nanmean(sCell[:150]) for sCell in sRec] for sRec, sess in zip(sTauList, sessList)
                      if 'HMC' not in sess]
    # compute average of second half for each cell for each recording
    meanSecondHalf = [[np.nanmean(sCell[150:]) for sCell in sRec] for sRec, sess in zip(sTauList, sessList)
                      if 'HMC' not in sess]

    # compute difference between halves for each cell for each recording
    diffHalf = [first - second for iRec in range(len(meanFirstHalf))
                for first, second in zip(meanFirstHalf[iRec], meanSecondHalf[iRec])]

    # randomize each cell time activity
    randS = [[np.random.permutation(sCell) for sCell in sRec] for sRec, sess in zip(sTauList, sessList)
             if 'HMC' not in sess]

    # same as above for ramdomized data
    meanFirstHalfRand = [[np.nanmean(sCell[:150]) for sCell in sRec] for sRec in randS]
    meanSecondHalfRand = [[np.nanmean(sCell[150:]) for sCell in sRec] for sRec in randS]
    diffHalfRand = [first - second for iRec in range(len(meanFirstHalf))
                for first, second in zip(meanFirstHalfRand[iRec], meanSecondHalfRand[iRec])]

    # plot histogram for each animal
    a.flat[iAni].hist([x for x in diffHalf if not np.isnan(x)], np.linspace(-0.02, 0.02, 100))
    a.flat[iAni].hist([x for x in diffHalfRand if not np.isnan(x)], np.linspace(-0.02,0.02,100), histtype='step')
    a.flat[iAni].set_title(animal)

# save plot
if saveFig:
    f.savefig(pathFigure + 'evaluateRateHalves.eps', format='eps')

plt.show(block=True)