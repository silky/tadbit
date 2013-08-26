"""
06 Aug 2013


"""
from pytadbit.imp.imp_modelling import generate_3d_models
import numpy as np

from scipy.optimize import anneal

global COUNT
COUNT = 0


def grid_search(zscores=None, upfreq_range=(0, 1, 0.1), lowfreq_range=(-1, 0, 0.1),
                scale_range=(0.005, 0.005, 0.001), maxdist_range=(400, 1500, 100),
                resolution=None, values=None, n_models=500,
                cutoff=300, n_keep=100, n_cpus=1, close_bins=1, verbose=True):
    count = 0
    if type(maxdist_range) == tuple:
        maxdist_step = maxdist_range[2]
        maxdist_arange = range(maxdist_range[0], maxdist_range[1] + maxdist_step,
                               maxdist_step)
    else:
        maxdist_arange = maxdist_range
    if type(lowfreq_range) == tuple:
        lowfreq_step = lowfreq_range[2]
        lowfreq_arange = np.arange(lowfreq_range[0],
                                   lowfreq_range[1] + lowfreq_step / 2,
                                   lowfreq_step)
    else:
        lowfreq_arange = lowfreq_range
    if type(upfreq_range) == tuple:
        upfreq_step = upfreq_range[2]
        upfreq_arange = np.arange(upfreq_range[0],
                                  upfreq_range[1] + upfreq_step / 2,
                                  upfreq_step)
    else:
        upfreq_arange = upfreq_range
    if type(scale_range) == tuple:
        scale_step = scale_range[2]
        scale_arange = np.arange(scale_range[0],
                                  scale_range[1] + scale_step / 2,
                                  scale_step)
    else:
        scale_arange = scale_range
    results = np.empty((len(scale_arange), len(maxdist_arange),
                        len(upfreq_arange), len(lowfreq_arange)))
    for w, scale in enumerate(scale_arange):
        for x, maxdist in enumerate(maxdist_arange):
            for y, upfreq in enumerate(upfreq_arange):
                for z, lowfreq in enumerate(lowfreq_arange):
                    tmp = {'kforce'   : 5,
                           'lowrdist' : 100,
                           'maxdist'  : maxdist,
                           'upfreq'   : upfreq,
                           'lowfreq'  : lowfreq,
                           'scale'    : scale}
                    tdm = generate_3d_models(zscores, resolution, n_models,
                                             n_keep, config=tmp, n_cpus=n_cpus,
                                             values=values,
                                             close_bins=close_bins)
                    count += 1
                    if verbose:
                        print '%5s  ' % (count), upfreq, lowfreq, maxdist, scale,
                    try:
                        result = tdm.correlate_with_real_data(cutoff=cutoff)[0]
                        if verbose:
                            print result
                        results[w, x, y, z] = result
                    except:
                        print 'ERROR'
    return results, scale_arange, maxdist_arange, upfreq_arange, lowfreq_arange


def to_optimize(params, zscores, resolution, values, n_models, n_keep,
                n_cpus=1):
    upfreq, lowfreq, maxdist = params
    tmp = {'kforce'   : 5,
           'lowrdist' : 100,
           'maxdist'  : maxdist,
           'upfreq'   : upfreq,
           'lowfreq'  : lowfreq}
    tdm = generate_3d_models(zscores, resolution, n_models, n_keep,
                             config=tmp, n_cpus=n_cpus, values=values)
    global COUNT
    COUNT += 1
    print '%5s  ' % (COUNT), params,
    try:
        result = tdm.correlate_with_real_data(cutoff=200)[0]
        print result
        return 1. - result
    except:
        print 'ERROR'
        return 1.0


def optimize(zscores, resolution, values):
    zscvals = sorted(reduce(lambda x, y: x+y,
                            [zscores[z].values() for z in zscores]))
    # lower bound must be higher than percentil 10 of zscores
    lzsc = zscvals[int(len(zscvals)*0.1)]
    # upper bound must be lower than percentil 90 of zscores
    uzsc = zscvals[int(len(zscvals)*0.9)]
    #
    print [(0.,uzsc),(lzsc,0.),(400, 1500)]
    #
    print anneal(to_optimize, (uzsc/2, lzsc/2, 700),
                 args=(zscores, resolution,
                       values, 500, 100, 8),
                 lower=(0, lzsc, 400), upper=(uzsc, 0, 2000), full_output=True)
    # print fmin_tnc(to_optimize, (uzsc/2, lzsc/2, 700), args=(zscores, resolution,
    #                                                    values, 500, 100, 8),
    #                bounds=((0.,uzsc),(lzsc,0.),(400, 2000)),
    #                approx_grad=True, epsilon=.01)
    # print fmin_l_bfgs_b(to_optimize, (uzsc/2, lzsc/2, 700), args=(zscores, resolution,
    #                                                               values, 8, 4, 8),
    #                     bounds=((0.,uzsc),(lzsc,0.),(100, 2000)),
    #                     approx_grad=True, epsilon=[.05,.05,10])
    # print fmin_slsqp(to_optimize, (uzsc/2, lzsc/2, 700), args=(zscores, resolution,
    #                                                            values, 8, 4, 8),
    #                  bounds=[(0.,uzsc),(lzsc,0.),(100, 2000)], epsilon=0.01)


def main():

    from pytadbit import Chromosome
    crm = '2R'
    crmbit = Chromosome('2R')

    for xnam in ['TR2', 'TR1', 'BR']:
        crmbit.add_experiment(xnam, resolution=10000, 
                              xp_handler='/home/fransua/db/hi-c/corces_dmel/10Kb/{0}/{0}_{1}_10Kb.txt'.format(crm, xnam))

    exp = crmbit.experiments['TR1'] + crmbit.experiments['TR2'] + crmbit.experiments['BR']

    start = 190
    end   = 295

    matrix = exp.get_hic_matrix()
    end += 1

    new_matrix = [[] for _ in range(end-start)]
    for i in xrange(start, end):
        for j in xrange(start, end):
            new_matrix[i - start].append(matrix[i][j])

    tmp = Chromosome('tmp')
    tmp.add_experiment('exp1', xp_handler=[new_matrix],
                       resolution=exp.resolution)

    exp2 = tmp.experiments[0]
    exp2.normalize_hic(method='visibility')
    exp2.get_hic_zscores(remove_zeros=True)
    values = [[float('nan') for _ in xrange(exp2.size)]
              for _ in xrange(exp2.size)]
    for i in xrange(exp2.size):
        # zeros are rows or columns having a zero in the diagonal
        if i in exp2._zeros:
            continue
        for j in xrange(i + 1, exp2.size):
            if j in exp2._zeros:
                continue
            if (not exp2.hic_data[0][i * exp2.size + j] 
                or not exp2.hic_data[0][i * exp2.size + j]):
                continue
            try:
                values[i][j] = (exp2.hic_data[0][i * exp2.size + j] /
                                exp2.wght[0][i * exp2.size + j])
                values[j][i] = (exp2.hic_data[0][i * exp2.size + j] /
                                exp2.wght[0][i * exp2.size + j])
            except ZeroDivisionError:
                values[i][j] = 0.0
                values[j][i] = 0.0


    optimize(exp2._zscores, exp.resolution, values)

    results = grid_search(upfreq_range=(0,1), lowfreq_range=(-1,0), freq_step=0.1,
                          zscores=exp2._zscores, resolution=exp2.resolution,
                          values=values,maxdist_range=(400, 1400),
                          n_cpus=8, n_models=10, n_keep=2)



    from matplotlib import pyplot as plt
    result, max_dist_arange, upfreq_arange, lowfreq_arange = results
    x = [i for i in max_dist_arange for j in upfreq_arange for k in lowfreq_arange]
    y = [j for i in max_dist_arange for j in upfreq_arange for k in lowfreq_arange]
    z = [k for i in max_dist_arange for j in upfreq_arange for k in lowfreq_arange]
    col = [result[i,j,k] for i in range(len(max_dist_arange)) for j in range(len(upfreq_arange)) for k in range(len(lowfreq_arange))]
    
    fig = plt.figure()
    from mpl_toolkits.mplot3d import Axes3D
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel('maxdist')
    ax.set_ylabel('upfreq')
    ax.set_zlabel('lowfreq')
    lol = ax.scatter(x, y, z, c=col, s=100, alpha=0.9)
    cbar = fig.colorbar(lol)
    cbar.ax.set_ylabel('Correlation value')
    plt.show()
