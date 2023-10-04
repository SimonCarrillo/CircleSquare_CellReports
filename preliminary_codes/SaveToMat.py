import numpy as np
import pickle
from scipy.io import savemat

# analyze recordings then save
# temporal analysis for both AV and Tau

expt = 'CircleSquare'

BinNb = 10  # for place maps
TimeBinSize_AVC = 1 * 1000
TimeBinSize_TAU = 1000
min_cell_nb = 15
LoadAssignments = True
SaveAssignments = False

pathAll = '/Volumes/GoogleDrive/My Drive/Fenton_Lab/CellPress_CircleSquare/MNS_Analysis/Data/' + expt + '/'

for animal in ['M19','M20','M29']:

    print(animal)

    pathSave = pathAll + 'analysisFiles_Simon/' + str(int(TimeBinSize_TAU)) + 'ms' + '_' + expt + '_' + animal + '.mat'

    # with open(path_all + animal + '_assignments_cs.file', 'rb') as f:
    #     assignments = pickle.load(f)

    with open(pathAll + 'analysisFiles_Simon/' +
              'TempCorr_List_' + str(int(TimeBinSize_AVC)) +
              'msAV_' + str(int(TimeBinSize_TAU)) + 'msKT' + '_' + expt + '_' + animal +
              '.file', 'rb') as f:
        outputList = pickle.load(f)

    dayRealList = outputList[0][1]
    sessList = outputList[0][2]
    # sAvcList = outputList[0][3]
    sTauList = outputList[0][4]

    output = []

    Days = np.unique(dayRealList)

    sAll = []
    sessIndicAll = []
    posAll = []
    dayAll = []

    for i_day, day in enumerate(Days):

        print(day)
        sessIdx = np.where([d == day for d in dayRealList])[0]

        sessIndic = []
        sDay = []
        posDay = []

        for iSessIdx in sessIdx:

            if 'HMC' in sessList[iSessIdx][:3]:
                sessIndic.append((-1)*np.ones_like(sTauList[iSessIdx][0]))
                sDay.append(sTauList[iSessIdx])
                # posDay.append([Xlist[iSessIdx], Ylist[iSessIdx]])

            elif 'CYL' in sessList[iSessIdx][:3]:
                sessIndic.append(np.ones_like(sTauList[iSessIdx][0]))
                sDay.append(sTauList[iSessIdx])
                # posDay.append([Xlist[iSessIdx], Ylist[iSessIdx]])

            elif 'RCT' in sessList[iSessIdx][:3]:
                sessIndic.append(np.zeros_like(sTauList[iSessIdx][0]))
                sDay.append(sTauList[iSessIdx])
                # posDay.append([Xlist[iSessIdx], Ylist[iSessIdx]])

        sAll.append(np.transpose(np.concatenate(sDay, 1)))
        sessIndicAll.append(np.concatenate(sessIndic, 0))
        dayAll.append(day)
        # posAll.append(np.transpose(np.concatenate(posDay, 1)))

    if animal == 'M34':
        sAll[0] = sAll[0][:-1]
        sessIndicAll[0] = sessIndicAll[0][:-1]

    mdic = {"spike": sAll, "pos": posAll, "sess": sessIndicAll, 'days':dayAll}
    savemat(pathSave, mdic)
