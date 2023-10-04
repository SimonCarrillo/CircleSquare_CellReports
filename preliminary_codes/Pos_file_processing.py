# Tools module by Eliott

def neighbor_sum(rate_map):
    import numpy as np
    import scipy.ndimage

    mask = np.isnan(rate_map)
    rate_map[mask] = 0

    conv_rate_map = scipy.ndimage.convolve(rate_map, np.ones([3, 3]), mode='constant')

    conv_rate_map[mask] = np.nan
    rate_map[mask] = np.nan
    out_tmp = np.corrcoef(rate_map[~mask].flat, conv_rate_map[~mask].flat)[0, 1]

    # rate_map = (rate_map) / np.std(rate_map[~np.isnan(rate_map)])
    # out_tmp = []
    # for index, val in np.ndenumerate(rate_map):
    #     if not np.isnan(val):
    #         map_tmp = rate_map.copy()
    #         map_tmp[index] = 0
    #         map_tmp[np.isnan(rate_map)] = 0
    #         sum_neigh = scipy.ndimage.convolve(map_tmp, np.ones([3, 3]), mode='constant')[index]
    #
    #         map_norm = np.ones(rate_map.shape)
    #         map_norm[index] = 0
    #         map_norm[np.isnan(rate_map)] = 0
    #         norm = scipy.ndimage.convolve(map_norm, np.ones([3, 3]), mode='constant')[index]
    #         if norm > 0:
    #             out_tmp.append(abs(val - sum_neigh/norm))

    # print(out_tmp)
    return np.mean(out_tmp)


def place_info_content(occ, rate_map):
    import numpy as np
    T = np.sum(occ)  # total time spent
    R = (np.sum(occ[~np.isnan(rate_map)] * rate_map[~np.isnan(rate_map)])) / T  # mean rate
    info_map = (occ / T) * (rate_map / R) * np.log2(rate_map / R)

    return np.sum(info_map[~np.isnan(info_map)])


def nan_gaussian_filter(map, sig):
    import scipy.ndimage
    import numpy as np

    V = map.copy()
    V[np.isnan(map)] = 0
    V_filtered = scipy.ndimage.gaussian_filter(V, sig)

    W = 0 * map.copy() + 1
    W[np.isnan(map)] = 0
    W_filtered = scipy.ndimage.gaussian_filter(W, sig)

    mask = W_filtered < 0.4

    ratio = V_filtered / W_filtered
    out = scipy.ndimage.gaussian_filter(ratio, sig)
    out[mask] = np.nan
    return out


def read_timestamps(filename):
    file = open(filename, 'r')
    columns = (line.split(sep='\t') for line in file)
    next(columns)  # skip first line
    z = list(zip(*columns))
    frameNum = list(map(int, z[1]))
    sysClock = list(map(int, z[2]))
    file.close()
    return frameNum, sysClock


def Pos_process(filename):
    file = open(filename, 'r')
    uu = (line for line in file)
    arg = (str.split(sep='\t') for str in uu if not (str[0] == '%' or str[0] == '\t'))
    z = list(zip(*arg))

    NumFrame = list(map(int, z[0]))
    TimeFrame = list(map(int, z[1]))
    Xpos = list(map(int, z[2]))
    Ypos = list(map(int, z[3]))
    Sector = list(map(int, z[4]))
    Shock = list(map(int, z[5]))
    file.close()

    return NumFrame, TimeFrame, Xpos, Ypos, Sector, Shock


def center_gravity(rate_map):
    import numpy as np
    y = 0
    x = 0
    w = 0

    for index, val in np.ndenumerate(rate_map):
        if not (np.isnan(val)):
            x += index[0] * val
            y += index[1] * val
            w += val

    return x / w, y / w


def center_high_rate(occ, rate_map):
    import numpy as np
    y = 0
    x = 0
    w = 0

    T = np.sum(occ)  # total time spent
    R = (np.sum(occ[~np.isnan(rate_map)] * rate_map[~np.isnan(rate_map)])) / T  # mean rate

    thresh = R + 4 * np.nanstd(rate_map.flat)
    for index, val in np.ndenumerate(rate_map):
        if not (np.isnan(val)) and val > thresh:
            x += index[0] * val
            y += index[1] * val
            w += val

    if w == 0:
        return np.nan, np.nan
    else:
        return x / w, y / w


def center_high_rate_contiguous(rm):
    import numpy as np
    import scipy.ndimage

    rate_map = nan_gaussian_filter(rm.copy(), 0.5)

    thresh = 4 * np.nanstd(rate_map.flat)
    rate_map[rate_map > thresh] = 1
    rate_map[rate_map <= thresh] = 0

    mask = np.isnan(rate_map)
    rate_map[mask] = 0
    conv_rate_map = scipy.ndimage.convolve(rate_map, np.ones([3, 3]), mode='constant')
    conv_rate_map[mask] = np.nan
    rate_map[mask] = np.nan

    map_norm = np.ones(rate_map.shape)
    map_norm[mask] = 0
    norm = scipy.ndimage.convolve(map_norm, np.ones([3, 3]), mode='constant')
    norm[mask] = np.inf

    conv_rate_map == norm - 1

    y = 0
    x = 0
    w = 0

    for index, val in np.ndenumerate(zip(conv_rate_map, rate_map)):
        if val:
            x += index[0] * val
            y += index[1] * val
            w += val

    return x / w, y / w


def conv_nan(x, w):
    import numpy as np
    half_w = int(np.floor(w / 2))

    a = []
    for idx in range(len(x)):
        min_b = max(0, idx - half_w)
        max_b = min(len(x), idx + half_w + 1)
        a.append(np.nanmean(x[min_b:max_b]))

    return a


def fisherZ(r):
    import numpy as np
    z = 0.5 * (np.log(1 + r) - np.log(1 - r))

    return z


def loadPhysiology(day, animal, xlAll, pathAll, expt, frame):
    import numpy as np
    import h5py
    # set up path
    xl = xlAll[xlAll.day == day]
    date = np.unique(xl.date)[0]
    print(date)
    TrackerPath = pathAll + 'tracker_files/' + date + '/'
    MnsTSPath = pathAll + 'timestamps/' + date + '/'
    MatPath = pathAll + 'mat/' + animal + '_' + date + '_out_1b_ss3' + '.mat'

    # extract recording info
    VideoName = list(xl.videoname)
    Sessions = list(xl.session)
    DayCount = list(xl.day_idx)[0]
    TrackerName = list(xl.trackername)
    doAnalysis = list(xl.analyze)
    Training = list(xl.training)[0]
    if expt == 'PlaceAvoidance' and frame == 'Arena':
        TrackerName = list(xl.arena)

    if expt == 'ATN':
        AtnSplit = [list(xl.start), list(xl.soundStart),
                    list(xl.foodStart), list(xl.foodEnd),
                    list(xl.soundEnd), list(xl.end)]
    else:
        AtnSplit = []

    MnsTSFileName = [MnsTSPath + s + '_timestamp.dat' for s in VideoName]
    TrackerFileName = [TrackerPath + str(s) for s in TrackerName]

    # check ordering ?!

    # evaluate length of each recording to later extract from longer file
    rec_idx = []
    rec_len = 0
    for fn in MnsTSFileName:
        VideoFrameNum, sysClock = read_timestamps(fn)
        nExtraTiff = (len(VideoFrameNum) - 1) // 10000
        # because of segmentation : frame 2 to 10001 (first and last frame always lost) then remaining modulo 3
        rec_idx.append([rec_len, rec_len + ((len(VideoFrameNum) - 1 - 10000 * nExtraTiff) // 3) + 3333 * nExtraTiff])
        rec_len += ((len(VideoFrameNum) - 1 - 10000 * nExtraTiff) // 3) + 3333 * nExtraTiff

    # load physiology
    f = h5py.File(MatPath, 'r')
    S_all = np.array(f['output']['S'])
    # C_all = np.array(f['output']['C'])
    A_tmp = np.array(f['output']['A'])

    xlOut = [VideoName, DayCount, Sessions, doAnalysis, Training, AtnSplit]
    fileNames = [MnsTSFileName, TrackerFileName]
    physio = [S_all, A_tmp, rec_idx]

    return [xlOut, fileNames, physio]


def SVM_decode(s, sess, do_cv, cv_rep, cv_param):
    import numpy as np
    from sklearn import svm
    from sklearn import metrics

    if not cv_rep == 100:
        print('cv rep = ' + str(cv_rep))

    if do_cv:
        # train & test w/ CV
        perf_tmp = []
        coef_tand = []
        acc_all = []
        roc_tmp = []

        for iT in range(cv_rep):

            n_train = int(np.min([np.sum(sess == i_bool) * (2 / 3) for i_bool in [0, 1]]))
            # n_train = int(np.min([np.sum(sess == i_bool) * (2/3) for i_bool in [-1, 1]]))

            idx_train = np.concatenate([
                np.random.choice(np.where(sess == i_bool)[0], n_train, replace=False)
                for i_bool in [0, 1]], 0)
            # for i_bool in [-1, 1]], 0)

            idx_test = np.array([x for x in range(len(sess)) if x not in idx_train])

            clf = svm.LinearSVC(tol=cv_param[0], max_iter=cv_param[1], verbose=cv_param[2])
            #clf = svm.SVC(probability=False, kernel="rbf", C=2.8, gamma=.0073, verbose=10)
            clf.fit(s[idx_train], sess[idx_train])
            perf_tmp.append(clf.score(s[idx_test], sess[idx_test]))
            coef_tand.append(clf.coef_[0])
            acc_tmp = np.zeros_like(sess) * np.nan
            acc_tmp[idx_test] = (clf.predict(s[idx_test])-sess[idx_test])**2
            #acc_tmp[idx_test] = (clf.predict(s[idx_test]) == sess[idx_test])

            try:
                roc_tmp.append(metrics.roc_auc_score(clf.predict(s[idx_test]), sess[idx_test]))
            except ValueError:
                roc_tmp.append(np.nan)
            pass

            acc_all.append(acc_tmp)

        perf = np.nanmean(perf_tmp)
        coef = np.nanmean(coef_tand, 0)
        acc = np.nanmean(acc_all, 0)
        intercept = []
        roc_all = np.nanmean(roc_tmp)

    else:
        # train & test
        clf = svm.LinearSVC(tol=cv_param[0], max_iter=cv_param[1], verbose=cv_param[2])
        clf.fit(s, sess)
        perf = clf.score(s, sess)
        coef = clf.coef_[0]
        intercept = clf.intercept_
        acc = []
        try:
            roc_all = metrics.roc_auc_score(s, sess)
        except ValueError:
            roc_all = np.nan

    return [perf, coef, roc_all, [intercept, acc]]


def runsTest(series, s0):
    import numpy as np
    runs = 0
    for i in range(len(series)):
        if (series[i] > s0 > series[i - 1]) or (series[i] < s0 < series[i - 1]):
            runs += 1

    n1 = np.sum([val > s0 for val in series])
    n2 = np.sum([val < s0 for val in series])

    runs_exp = ((2 * n1 * n2) / (n1 + n2)) + 1
    stan_dev = np.sqrt((2 * n1 * n2 * (2 * n1 * n2 - n1 - n2)) / (((n1 + n2) ** 2) * (n1 + n2 - 1)))

    z = (runs - runs_exp) / stan_dev

    return z
