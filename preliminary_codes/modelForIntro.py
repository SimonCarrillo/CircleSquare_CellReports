import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib as mpl
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

path_all = '/Users/eliott/Dropbox/dDocuments/NYU/Fenton/MNS_ANALYSIS/Data/CircleSquare/'
pathFigure = '/Users/eliott/Dropbox/dDocuments/NYU/Fenton/MNS_ANALYSIS/Figures/CircleSquare/savedFigs/'
BinNb = 12  # number of bins for place maps
sampRate = 10  # Hz
cvRep = 100
nMinutes = 240
nSteps = sampRate * 60 * nMinutes #sample rate * 60 seconds * minutes
minSpeed = 0.1
maxSpeed = 1 #10 samples per seconds, max speed ~ 10 per s = 10/10 per sample ~ 0.3 per sample
speedChange = [-0.1,0,0.1]
dayPick = 5

iDayCountToMatch = 0 # to match PCs across days

xAxisDays = [1, 2, 3, 8, 9, 10, 15, 16, 17]

animalList = ['M19', 'M35', 'M29','M20','M35','M39']
nAni = len(animalList)

envList = ['RCT1', 'RCT2', 'CYL1', 'CYL2']

cmap = mpl.cm.viridis
cmap.set_bad('w', 1.)
smoothFilter = 1/3

bigMapsList = [[], [], [], []]

print(str(nMinutes)+ ' minutes')
print('day pick = ' + str(dayPick))

# # past results
# f,a = plt.subplots(1,1)
# A = [0.991865104,0.980224479,0.971194792 ,0.924270833, 0.864492708]
# a.bar([0,1,2,3,4],A)
# f.savefig(pathFigure + 'modelResults.eps')

#load maps from all animals
for iAni, animal in enumerate(animalList):

    loadpath = path_all + 'analysisFiles/' + 'MapCorr_' + str(int(BinNb)) + '_Bins_' + animal + '.file'

    with open(loadpath, 'rb') as f:
        output = pickle.load(f)

    # load analyzed data
    # [maps_corr, maps_corr_pc, day_count_pair, day_real_pair, sess_pair, day_sep, same_diff,
    #  LinMapsCorr, LinMapsCorrPC, cellCount] = output[0]
    # [maps_corr_pc_any, maps_corr_pc_all, maps_corr_pc_not] = maps_corr_pc
    [day_real_list, day_count_list, sess_list, rate_maps_list, PC_idx_list, LinMapsList, OccList,
     cohList, pInfoList] = output[1]
    # sTracked = output[2][2]

    recList = np.where([d==dayPick and 'HMC' not in s for d,s in zip(day_count_list, sess_list)])[0]

    print(animal)
    print(len(recList)==4)

    PcList = np.any([PC_idx_list[idx] for idx in recList], 0)

    if len(recList)==4:
        for recIdx in recList:

            #filter rate maps
            rateMapsTmp = [rM for rM, pc in zip(rate_maps_list[recIdx], PcList) if pc]

            #norm
            # rateMaps = [rM/np.nanmax(rM) if np.nanmax(rM)>0 else rM for rM in rateMapsTmp]
            #no norm
            rateMaps = rateMapsTmp

            # make larger map
            bigMaps = [np.concatenate([np.concatenate([rateMaps[0+iMap*4],rateMaps[1+iMap*4]],0),
                                     np.concatenate([rateMaps[2+iMap*4],rateMaps[3+iMap*4]],0)],1)
                      for iMap in range(len(rateMaps)//4)]

            if envList[0] in sess_list[recIdx]:
                bigMapsList[0].extend(bigMaps)
            elif envList[1] in sess_list[recIdx]:
                bigMapsList[1].extend(bigMaps)
            elif envList[2] in sess_list[recIdx]:
                bigMapsList[2].extend(bigMaps)
            elif envList[3] in sess_list[recIdx]:
                bigMapsList[3].extend(bigMaps)

print(len(bigMapsList[0]))

f,aMask = plt.subplots(1,4)

limSet = [[6,6], [6,12], [12,12], [12,24], [24,12], [24,24]]

for iLim, lims in enumerate(limSet):

    print(lims)

    fOcc, aOcc = plt.subplots(2,4)
    fTraces, aTraces = plt.subplots(2, 4)

    sAll = []
    sessIndic = []

    for iEnv in range(4):

        mask = np.any(np.isnan(bigMapsList[iEnv]),0)

        if iLim==0: aMask[iEnv].imshow(mask)

        xLim = lims[0]
        yLim = lims[1]
        speed = minSpeed
        pos = [np.random.randint(xLim), np.random.randint(yLim)]

        [_, _, _, loc] = scstat.binned_statistic_2d([pos[0]], [pos[1]], [],
                                                    bins=[x for x in range(24)], statistic='count',
                                                    expand_binnumbers=True)
        while mask[loc[0][0] - 1][loc[1][0] - 1]:
            pos = [np.random.randint(xLim), np.random.randint(yLim)]
            [_, _, _, loc] = scstat.binned_statistic_2d([pos[0]], [pos[1]], [],
                                                        bins=[x for x in range(24)], statistic='count',
                                                        expand_binnumbers=True)

        posList = [pos]

        for s in range(nSteps):
            dX = np.random.choice([-1,0,1],2)
            speed += np.random.choice(speedChange,1)[0]
            if speed>maxSpeed: speed = maxSpeed
            if speed<minSpeed: speed = minSpeed
            newPos = pos + dX * speed
            if newPos[0]>xLim: newPos[0] = xLim
            if newPos[1]>yLim: newPos[1] = yLim
            if newPos[0]<0: newPos[0] = 0
            if newPos[1]<0: newPos[1] = 0

            [_, _, _, loc] = scstat.binned_statistic_2d([newPos[0]], [newPos[1]], [],
                                                        bins=[x for x in range(24)], statistic='count',
                                                        expand_binnumbers=True)

            if mask[loc[0][0]-1][loc[1][0]-1]: newPos = pos
            posList.append(newPos)
            pos = newPos

        aOcc[0][iEnv].plot([x[0] for x in posList], [x[1] for x in posList])
        aOcc[0][iEnv].set_ylim([0, 24])
        aOcc[0][iEnv].set_xlim([0, 24])
        aOcc[0][iEnv].set_aspect('equal', adjustable='box')

        [occTmp, x_edge, y_edge, loc] = scstat.binned_statistic_2d([x[0] for x in posList], [x[1] for x in posList], [],
                                                                   bins=[x for x in range(24)], statistic = 'count',
                                                                   expand_binnumbers = True)

        aOcc[1][iEnv].imshow(occTmp, origin='lower')
        aOcc[1][iEnv].axis('off')


        if iLim==0:
            f,a = plt.subplots(6,14, figsize=[17,10])
            for ax,m in zip(a.flat,bigMapsList[iEnv]):
                ax.imshow(nan_gaussian_filter(m, smoothFilter), cmap=cmap)
                ax.axis('off')
            f.savefig(pathFigure + 'bigMap_day' + str(dayPick) + '_' + envList[iEnv] + '.eps', fmt='eps', dpi=350)

        muRegular = [[r[loc[0, i] - 1][loc[1, i] - 1] / sampRate
                       for i in range(len(loc[0]))] for r in bigMapsList[iEnv]]

        muBinned = [[np.mean(s[0+iBin*sampRate:sampRate+iBin*sampRate]) for iBin in range(nSteps//sampRate)]
            for s in muRegular]


        for iCell, mu in enumerate(muRegular):
            aTraces[0][iEnv].plot([m+iCell/10 for m in mu])
        for iCell, mu in enumerate(muBinned):
            aTraces[1][iEnv].plot([m+iCell/10 for m in mu])

        if iLim == len(limSet)-1:
            fOcc.savefig(pathFigure + 'bigMapOcc_day' + str(dayPick) + '_' + str(nMinutes) + 'mn.eps',
                         fmt='eps', dpi=250)
            fTraces.savefig(pathFigure + 'bigMapActivity_day' + str(dayPick) + '_' + str(nMinutes) + 'mn.eps',
                         fmt='eps')

        if iEnv<2:
            sessIndic.append(np.ones_like(muBinned[0]))
        else:
            sessIndic.append(np.zeros_like(muBinned[0]))

        sAll.append(muBinned)

    sessConcat = np.concatenate(sessIndic, 0)
    muConcat = np.transpose(np.concatenate(sAll, 1))

    pCv, cCv, i__ = SVM_decode(muConcat, sessConcat, do_cv=True, cv_rep=cvRep, cv_param=[1e-10, 100000, 0])

    f,a = plt.subplots(1,2)
    a[0].imshow(np.transpose(muConcat), aspect='auto', vmax= 0.05)
    a[1].imshow([np.transpose(muConcat)[x] for x in np.argsort(cCv)[::-1]], aspect='auto', vmax= 0.05)


    print(pCv)

print()
print('done done done')
# load maps from all animal
#
# extract place cells
#
# combine cells
#
# define random walk
#
# model of activity
#
# find peak

plt.show(block=True)