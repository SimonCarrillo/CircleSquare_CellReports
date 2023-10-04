import matplotlib as mpl
mpl.use('TkAgg')

import pickle
from scipy.io import savemat

expt = 'CircleSquare'

path_all = '/Volumes/GoogleDrive/My Drive/Fenton_Lab/CellPress_CircleSquare/MNS_Analysis/Data/CircleSquare/'
pathFigure = '/Volumes/GoogleDrive/My Drive/Fenton_Lab/CellPress_CircleSquare/MNS_Analysis/Figures/CircleSquare/savedFigs/IsoMap_Simon/'

TimeBinSize_AVC = 1000
TimeBinSize_TAU = 1000

threshMult = 2
thresh = 0.01* threshMult*(TimeBinSize_TAU/1000)

saveOutput = True

doStd = True
doPTI = True

saveAni = ['M19']

animalList = saveAni  #['M19']#,'M20', 'M29','M34','M19', 'M39']
nAni = len(animalList)

indic = [-1, 0, 1]
envList = ['HMC1','HMC2', 'CYL1','CYL2', 'RCT1','RCT2'] # environment HMC = Home Cage (-1), CYL = Cylinder (1), RCT = Rectangle (0)

if doPTI:
    loadStr = 'Zlinear'
else:
    loadStr = 'TempCorr'

for iAni, animal in enumerate(animalList):

    print(animal)

    with open(path_all + 'analysisFiles_Simon/' + loadStr + '_List_' +
              str(int(TimeBinSize_AVC)) + 'msAV_' + str(int(TimeBinSize_TAU)) + 'msKT'
              '_CircleSquare_' + animal +'.file', 'rb') as f:
        outputList = pickle.load(f)

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

if saveOutput:
    print('saving...')

    pathSave = path_all + 'analysisFiles_Simon/' + loadStr + '_List_' + str(int(TimeBinSize_AVC)) + 'msAV_' + str(int(TimeBinSize_TAU)) + 'msKT' +'_CircleSquare_' + animal + '.mat'

    mdic = {'dayCountList': dayCountList,'sessList':sessList, 'sTauList':sTauList,'isSeparatedList':isSeparatedList,'tauPairsSingleList':tauPairsSingleList,'tauVecSingleList':tauVecSingleList,'trainingList':trainingList}
    savemat(pathSave, mdic)

    print('saved!')

else:
    print('NOT SAVING')

print('done done done')