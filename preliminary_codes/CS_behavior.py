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
from sklearn import linear_model

path_all = '/Users/eliott/Dropbox/dDocuments/NYU/Fenton/MNS_ANALYSIS/Data/CircleSquare/analysisFiles/'
pathFigure = '/Users/eliott/Dropbox/dDocuments/NYU/Fenton/MNS_ANALYSIS/Figures/CircleSquare/behavior/'
BinNb = 12  # number of bins for place maps
TimeBinSize_AVC = 5 * 1000
TimeBinSize_TAU = 300
nsep = 18
sampRate = 10  # Hz

xAxisDays = [1, 2, 3, 8, 9, 10, 15, 16, 17]

AnimalList = ['M19', 'M20', 'M29', 'M34', 'M35', 'M39']

nAni = len(AnimalList)

occAll = [[[] for rec in range (4)] for day in range(9)]

f1Sec, a1Sec = plt.subplots(2,3)

for iAni, animal in enumerate(AnimalList):

    print(animal)

    loadpath = path_all + 'MapCorr_' + str(int(BinNb)) + '_Bins_' + animal + '.file'

    with open(loadpath, 'rb') as f:
        output = pickle.load(f)

    # load analyzed data
    OccList = output[1][6]
    dayCountList = output[1][1]
    sessList = output[1][2]
    xTrackedList = output[2][0]
    yTrackedList = output[2][1]

    Output = None

    fTrace, axTrace = plt.subplots(9, 4, figsize=(9, 10))
    fOcc, axOcc = plt.subplots(9, 4, figsize=(9, 10))

    # plot behavior
    dayTrack = -1
    for iIdx, occ in enumerate(OccList):
        iDay = dayCountList[iIdx]-1

        isCyl = 'CYL' in sessList[iIdx] # is it a CYL recording
        isSecond = '2' in sessList[iIdx] # is it a CYL2 or RCT2 recording

        recTrack = isSecond + 2*isCyl

        axOcc[iDay][recTrack].imshow(occ, vmin=0, vmax=3)
        axTrace[iDay][recTrack].plot(xTrackedList[iIdx], yTrackedList[iIdx])

        axOcc[iDay][recTrack].set_title(sessList[iIdx])
        axTrace[iDay][recTrack].set_title(sessList[iIdx])

        occAll[iDay][recTrack].append(occ)

    for ax in axTrace.flat:
        ax.axis('off')

    for ax in axOcc.flat:
        ax.axis('off')

    fOcc.suptitle(animal)
    fTrace.suptitle(animal)

    fOcc.subplots_adjust(hspace=0.5)
    fTrace.subplots_adjust(hspace=0.5)

    fOcc.savefig(pathFigure + animal + 'traces.eps', format='eps')
    fTrace.savefig(pathFigure + animal + 'dwellMaps.eps', format='eps')

    dist1Sec = []
    for xtracked,yTracked in zip(xTrackedList,yTrackedList):
        # bin animal position (and export bin identity!)
        [occTmp, x_edge, y_edge, loc] = scstat.binned_statistic_2d(xtracked, yTracked, [], bins=BinNb,
                                                                   statistic='count', expand_binnumbers=True)

        # compute distance between animal position a second later (10 timebins)
        dist1Sec.extend([np.linalg.norm([loc[0][i + 0] - loc[0][i + 10], loc[1][i + 0] - loc[1][i + 10]])
                         for i in range(len(loc[0]) - 10)])

    a1Sec.flat[iAni].hist(dist1Sec,np.linspace(-0.5,18.5,20), density=True)
    a1Sec.flat[iAni].set_title(animal)

    print('la')

fMeanOcc, aMeanOcc = plt.subplots(9, 4, figsize=(9, 10))
for iDay in range(9):
    for iRec in range(4):
        aMeanOcc[iDay][iRec].imshow(np.mean(occAll[iDay][iRec], 0), vmin=0, vmax=3)
        aMeanOcc[iDay][iRec].axis('off')

aMeanOcc[0][0].set_title('RCT1')
aMeanOcc[0][1].set_title('RCT2')
aMeanOcc[0][2].set_title('CYL1')
aMeanOcc[0][3].set_title('CYL2')
fMeanOcc.savefig(pathFigure + 'averageDwellMaps.eps', format='eps')
f1Sec.savefig(pathFigure + 'distance1Sec.eps', format='eps')

plt.show(block=True)
