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

# to do :
# check ordering of files
# compute speed
# check PCness with new method (new cnmfe + subsampling)
# do random corr for place maps and other metrics

expt = 'ATN'

animal = 'M39'

winsizeList = [10, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250]
# winsizeList=[100,110]
flowItList = [9]
flowLevelsList = [5]

path_all = '/Users/eliott/Dropbox/dDocuments/NYU/Fenton/MNS_ANALYSIS/Data/' + expt + '/'
ExcelPath = path_all + 'FilenameMatch.xlsx'

# process excel file
xlread = pd.ExcelFile(ExcelPath)
xl_all = pd.read_excel(xlread, animal)
xl_all = xl_all[pd.notna(xl_all.videoname)]

print(animal)
Days = np.unique(list(xl_all.day))

A_animal = []
templates = []

print('la')

for day in Days: #[Days[d] for d in [0,1,2, 3,4,5, 6,7,8]]: # 0,1,2, 3,4,5, 6,7,8
    xl = xl_all[xl_all.day == day]
    date = np.unique(xl.date)[0]
    print(date)
    MatPath = path_all + 'mat/' + animal + '_' + date + '_out_1b_ss3' + '.mat'
    template_path = path_all + 'templates/' + animal + '_' + date + '_template.mat'

    # load physiology
    f = h5py.File(MatPath, 'r')
    A_tmp = np.array(f['output']['A'])
    A_animal.append(scipy.sparse.csc_matrix(A_tmp.transpose()))


    # load template
    x = loadmat(template_path)
    templates.append(x['template2'])
    f = []
    x = []

accuracy = [[] for x in flowItList]
overallMatch = [[] for x in flowItList]
numberMatched = [[] for x in flowItList]
fScore = [[] for x in flowItList]

for iFlowIt, flowIt in enumerate(flowItList):

    for iFlowLevels, flowLevels in enumerate(flowLevelsList):

        accuracyWin = []
        numberMatchedWin = []
        overallMatchWin = []
        fScoreWin = []

        for winsize in winsizeList:
            print(flowIt)
            print(flowLevels)
            print(winsize)
            accuracy_tmp = []
            numberMatched_tmp = []

            # assignments = caiman.base.rois.register_multisession(A_animal, (480, 752), align_flag=False, templates=templates, winsize=winsize)[1]

            assign = []
            for iPair in itertools.combinations(range(len(A_animal)), 2):
                print(iPair)
                assign.append(caiman.base.rois.register_ROIs(A_animal[iPair[0]], A_animal[iPair[1]], (480, 752),
                                                     template1=templates[iPair[0]], template2=templates[iPair[1]],
                                                     align_flag=True, use_opt_flow=True,
                                                     max_dist=15, winsize=winsize,
                                                     flow_levels=flowLevels, flow_iterations = flowIt))

            # for iDay in range(len(A_animal)-1):
            #     print(iDay)
            #     assign.append(caiman.base.rois.register_ROIs(A_animal[iDay], A_animal[iDay+1], (480, 752),
            #                                          template1=templates[iDay], template2=templates[iDay-1],
            #                                          align_flag=True, use_opt_flow=True,
            #                                          max_dist=15, winsize=winsize,
            #                                          flow_levels=flowLevels, flow_iterations = flowIt))

            accuracyWin.append([len(iAssign[0])/(len(iAssign[0])+len(iAssign[2])) for iAssign in assign])
            numberMatchedWin.append([len(iAssign[0]) for iAssign in assign])
            # overallMatchWin.append(np.sum(np.all(np.isfinite(assignments),1))/len(assignments))
            fScoreWin.append([a[4]['f1_score'] for a in assign])

        accuracy[iFlowIt].append(accuracyWin)
        numberMatched[iFlowIt].append(numberMatchedWin)
        overallMatch[iFlowIt].append(overallMatchWin)
        fScore[iFlowIt].append(fScoreWin)


metricName = ['accuracy', 'numberMatched', 'fScore']

for iMetric, metric in enumerate([accuracy, numberMatched, fScore]):
    f, a = plt.subplots(1,len(flowItList))
    if len(flowItList)<2:
        for iFlowLevels, flowLevels in enumerate(flowLevelsList):
            if len(flowLevelsList)<2:
                a.plot(winsizeList, metric[iFlowIt][iFlowLevels])
                a.errorbar(winsizeList, np.mean(metric[iFlowIt][iFlowLevels],1),
                           np.std(metric[iFlowIt][iFlowLevels],1), color='k')
            else:
                a.errorbar(winsizeList, np.mean(metric[iFlowIt][iFlowLevels], 1),
                           np.std(metric[iFlowIt][iFlowLevels], 1))
        a.set_title(' flowIt =' + str(flowIt))
        a.legend(['flowLvl ' + str(x) for x in flowLevelsList])
        # a.set_ylim([0, 1])
    else:
        for iFlowIt, flowIt in enumerate(flowItList):
            for iFlowLevels, flowLevels in enumerate(flowLevelsList):
                if len(flowLevelsList) < 2:
                    a[iFlowIt].plot(winsizeList, metric[iFlowIt][iFlowLevels])
                    a[iFlowIt].errorbar(winsizeList, np.mean(metric[iFlowIt][iFlowLevels], 1),
                                        np.std(metric[iFlowIt][iFlowLevels], 1), color='k')
                else:
                    a[iFlowIt].errorbar(winsizeList, np.mean(metric[iFlowIt][iFlowLevels], 1),
                                        np.std(metric[iFlowIt][iFlowLevels], 1))
            a[iFlowIt].set_title(' flowIt =' + str(flowIt))
            a[iFlowIt].legend(['flowLvl ' + str(x) for x in flowLevelsList])
            # a[iFlowIt].set_ylim(yLims)
    f.suptitle(animal + ' ' + metricName[iMetric])

# a = plt.subplots(1,3)[1]
# a[0].plot(overallMatch)
# a[1].plot(winsizeList, np.mean(accuracy,1))
# a[2].plot(winsizeList, np.mean(numberMatched,1))
#
# print(accuracy)
# print(numberMatched)

print('done done done')

plt.show(block=True)
