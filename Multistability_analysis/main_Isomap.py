import os

import matplotlib as mpl

mpl.use('TkAgg')
from mpl_toolkits.mplot3d import Axes3D

import matplotlib.pyplot as plt
import numpy as np
import pickle
from sklearn import manifold
import pandas as pd
from scipy.io import savemat

TimeBinSize_AVC = 1000
TimeBinSize_TAU = 1000

n_neighbors = 5
target_dim = 10
doStd = True
doPTI = True
expt = 'CircleSquare'
saveOutput = True

if expt == 'CircleSquare':
    path_all = '/Volumes/GoogleDrive/My Drive/Fenton_Lab/CellPress_CircleSquare/MNS_Analysis/Data/CircleSquare/'
elif expt == 'PlaceAvoidance':
    path_all = '/Volumes/GoogleDrive/My Drive/Fenton_Lab/CellPress_CircleSquare/MNS_Analysis/Data/PlaceAvoidance/'

saveAni = ['M19']

animalList = ['M19']
nAni = len(animalList)

indic = [-1, 0, 1]
envList = ['HMC1', 'HMC2', 'CYL1', 'CYL2', 'RCT1',
           'RCT2']  # environment HMC = Home Cage (-1), CYL = Cylinder (1), RCT = Rectangle (0)
indicColors = [[0, 0, 0], [0.5, 0.5, 0.5, 0.5], [0, 0, 1], [0, 0, 0.5], [1, 0, 0], [0.5, 0, 0]]

ratioDiscardList = [4]  # [20,10,4,2] # discarding 1/ratioDiscardList %
discardPop = 'kc'  # kc: Keep Cells, pos: Positive Cells Pairs, neg: Negative Cell Pairs
discard = True  # if True: discard the 1/ratioDiscard % most (anti-)correlated if False: only keep the 1/ratioDiscard
# % most (anti-)correlated

if doPTI:
    loadStr = 'Zlinear'
else:
    loadStr = 'TempCorr'

for iAni, animal in enumerate(animalList):

    IsoToSave = []

    for iDay, day in enumerate(range(1, 10)):

        for iR, ratio in enumerate(ratioDiscardList):

            if os.path.exists(path_all + 'analysisFiles_Simon/' +
                              'CellDiscarded_' + discardPop + 'Ratio' + str(int(ratio)) + loadStr + str(
                int(TimeBinSize_TAU)) + expt + '_' + animal + '_' + 'Day' + str(int(iDay)) + '.file'):

                with open(path_all + 'analysisFiles_Simon/' + 'CellDiscarded_' + discardPop + 'Ratio' + str(
                        int(ratio)) + loadStr + str(int(TimeBinSize_TAU)) + expt + '_' + animal + '_' + 'Day' + str(
                    int(iDay)) + '.file', 'rb') as f:
                    [data_to_use, sAll, sessIndicAll, sessIndicFullAll, tauVec, idxTauPairs, tauVecTmp,
                     idxTauPairsTmp, keptCellsList, discardCells, sConcat, sessIndic, sessIndicFull] = pickle.load(f)

                print(animal, 'Day', iDay, 'DiscardRatio', ratio, discardPop)

                # do Isomap transformation
                iso_instance = manifold.Isomap(n_neighbors=n_neighbors, n_components=target_dim, path_method='D')
                dataIso = iso_instance.fit_transform(data_to_use)

                # ## figures
                # figScatter3d = plt.figure(figsize=(8, 8))
                # ax = Axes3D(figScatter3d)
                #
                # for iIndic, envIndic in enumerate(envList):
                #     proj_data = dataIso[[envIndic in s for s in sessIndicFullAll], :]
                #     ax.scatter(proj_data[:, 0], proj_data[:, 1], proj_data[:, 2],
                #                s=3, alpha=0.2, edgecolor='face', color=indicColors[iIndic])
                #     ax.plot([0, np.mean(proj_data, 0)[0]], [0, np.mean(proj_data, 0)[1]], [0, np.mean(proj_data, 0)[2]],
                #             '*-', color=indicColors[iIndic], lw=5)
                #     ax.set_xlabel('Component 1')
                #     ax.set_ylabel('Component 2')
                #     ax.set_zlabel('Component 3')

                IsoToSaveTmp = [data_to_use, dataIso, sAll, sessIndicAll, sessIndicFullAll, tauVec, idxTauPairs,
                                tauVecTmp, idxTauPairsTmp, envList, keptCellsList, discardCells]
            else:
                IsoToSaveTmp = []
                continue

        IsoToSave.append(IsoToSaveTmp)

    if saveOutput:
        print('saving...')
    
        with open(path_all + 'analysisFiles_Simon/' +
                  'ISO_' + str(int(target_dim)) + '_' + 'CellDiscarded_' + discardPop + 'Ratio' + str(
            int(ratio)) + loadStr + str(
            int(TimeBinSize_TAU)) + expt + '_' + animal + '_' +
                  '.file', 'wb') as f:
            pickle.dump(IsoToSave, f, pickle.HIGHEST_PROTOCOL)
    
        pathSave_ISO = path_all + 'analysisFiles_Simon/' + 'ISO_' + str(
            int(target_dim)) + '_' + 'CellDiscarded_' + discardPop + 'Ratio' + str(
            int(ratio)) + loadStr + str(int(TimeBinSize_TAU)) + expt + '_' + animal + '_' + '.mat'
    
        mdic_ISO = {"ISO": IsoToSave}
        savemat(pathSave_ISO, mdic_ISO)
    
        print('saved!')

print('done done done')
