import matplotlib as mpl
# mpl.use('TkAgg')
import matplotlib.pyplot as plt

import itertools
import h5py
import numpy as np
import pandas as pd
import caiman.base.rois
import scipy
from scipy.io import loadmat
import pickle

expt = 'CircleSquare'

animal = 'M20'

winsize = 75
flowIt = 7
flowLvl = 5

path_all = '/Volumes/GoogleDrive/My Drive/Fenton_Lab/CellPress_CircleSquare/MNS_Analysis/Data/' + expt + '/'
ExcelPath = path_all + 'FilenameMatch.xlsx'

pathSave = path_all + 'analysisFiles_Simon/' + 'AssignmentsByPairs_' + animal + '.file'

# process excel file
xlread = pd.ExcelFile(ExcelPath)
xl_all = pd.read_excel(xlread, animal)
xl_all = xl_all[pd.notna(xl_all.videoname)]

print(animal)
Days = np.unique(list(xl_all.day))

A_animal = []
templates = []

for day in Days: #[Days[u] for u in [0,1,3,4,5,6,7,8]]:
    xl = xl_all[xl_all.day == day]
    date = np.unique(xl.date)[0]
    print(date)
    MatPath = path_all + 'mat/' + animal + '_' + date + '_out_1b_ss3' + '.mat'
    templatePath = path_all + 'templates/' + animal + '_' + date + '_template.mat'

    # load physiology
    f = h5py.File(MatPath, 'r')
    A_tmp = np.array(f['output']['A'])
    A_animal.append(scipy.sparse.csc_matrix(A_tmp.transpose()))


    # load template
    x = loadmat(templatePath)
    templates.append(x['template2'])
    f = []
    x = []

assignments = caiman.base.rois.register_multisession(A_animal, (480, 752),
                                                     align_flag=True, templates=templates,
                                                     winsize=winsize, max_dist = 15)[1]

print('Registration 1 done')

nDay = len(Days)

assignmentsByPairs = []
idPairs = []
for iPair in itertools.combinations(range(nDay),2):
    print(iPair)
    assignmentsByPairs.append(caiman.base.rois.register_ROIs(A_animal[iPair[0]], A_animal[iPair[1]], (480, 752),
                                                         template1 = templates[iPair[0]],
                                                         template2 = templates[iPair[1]],
                                                         align_flag=True, use_opt_flow = True,
                                                         max_dist = 15, winsize=winsize, plot_results=True,
                                                             flow_levels = flowLvl, flow_iterations = flowIt))
    idPairs.append(iPair)
    # idPairs.append([Days[iPair[0]], Days[iPair[1]]])

print('Registration 2 done')

with open(pathSave, 'wb') as f:
    pickle.dump([assignments,[assignmentsByPairs,idPairs]], f, pickle.HIGHEST_PROTOCOL)

print(pathSave)

print('done done done')

plt.show(block=True)