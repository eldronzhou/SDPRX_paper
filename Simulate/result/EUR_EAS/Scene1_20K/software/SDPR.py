#!/usr/bin/env python
import numpy as np
from scipy import stats
import cPickle as pickle
import gzip
import gibbs, ld
import argparse
import sys, util
import pandas as pd


def SDPR_gibbs(beta_margin, N, ld_boundaries, ref_ld_mat, mcmc_samples, 
    burn, max_cluster, save_mcmc, VS=True):
    M = max_cluster
    trace = {'alpha':[], 'num_cluster':[], 'beta':np.zeros(shape=(mcmc_samples, len(beta_margin))),
        'suffstats':[], 'h2':[]}

    # initialize
    state = gibbs.initial_state(data=beta_margin, N=N, num_clusters=M, a0k=.1, b0k=.1)
    #state['suffstats'] = gibbs.update_suffstats(state)
    #state['cluster_var'] = gibbs.sample_sigma2(state, VS=True)
    state['cluster_var'] = {0:0, 1:.165/(len(beta_margin)*5.6e-4)} 
    state['pi'] = {0:1-5.6e-4, 1:5.6e-4}

    for i in range(mcmc_samples):
        # update everything
        gibbs.gibbs(state, ld_boundaries=ld_boundaries, ref_ld_mat=ref_ld_mat, VS=VS)

        if (i > burn):
            trace['h2'].append(state['h2']*state['eta']**2)

        # record the result
        trace['beta'][i,] = state['beta']*state['eta']
        # trace['pi'][i,:] = np.array(state['pi'].values())
        # trace['cluster_var'][i,:] = np.array(state['cluster_var'].values())
        # trace['alpha'].append(state['alpha'])
        # trace['num_cluster'].append( np.sum(np.array(state['pi'].values()) > .0001) )
        # trace['suffstats'].append(state['suffstats'])

        util.progressBar(value=i+1, endvalue=mcmc_samples)

    # calculate posterior average
    poster_mean = np.mean(trace['beta'][burn:mcmc_samples], axis=0)
    
    print 'h2: ' + str(np.median(trace['h2']))

    if save_mcmc is not None:
        f = gzip.open(args.save_mcmc, 'wb')
        pickle.dump(trace, f, protocol=2)
        f.close()

    return poster_mean

def pipeline(args):
    
    # sanity check

    if args.bfile is not None and args.load_ld is not None:
        raise ValueError('Both --bfile and --load_ld flags were set. \
            Please use only one of them.')

    if args.bfile is None and args.load_ld is None:
        raise ValueError('Both --bfile and --load_ld flags were not set. \
            Please use one of them.')

    N = args.N

    print('Load summary statistics from {}'.format(args.ss))
    ss = pd.read_table(args.ss)
    # beta_margin = np.array(np.sign(ss['BETA']) * abs(stats.norm.ppf(ss['P'] / 2.0)) / np.sqrt(N))
    beta_margin = np.array(ss['T_STAT']) / np.sqrt(N)
    assert np.all(np.isreal(beta_margin)), 'Something wrong with summary stats.'

    # to be correct flip 
    # flip = pd.read_table("/ysm-gpfs/pi/zhao/gz222/UKB_simulate/flip.txt", header=None)
    # flip_idx = np.array(flip[0])-1
    # beta_margin[flip_idx] = -beta_margin[flip_idx]

    if args.load_ld is not None:
        print('Load pre-computed reference LD from {}'.format(args.load_ld))
        f = gzip.open(args.load_ld, 'r')
        ld_dict = pickle.load(f)
        f.close()
        ld_boundaries = ld_dict[0]
        ref_ld_mat = ld_dict[1]
    else:
        print('Calculating reference LD. May take ~ 2 hours ...')
        ld_boundaries = ld.parse_ld_boundaries(min_block_wd=100, ss=ss, block_path=args.block)
        ref_ld_mat = Parallel(n_jobs=args.threads)(delayed(ld.calc_ref_ld)(i, ref_path=args.bfile, 
                ld_boundaries=ld_boundaries) for i in range(len(ld_boundaries))) 
        if args.save_ld is not None:
            print('Save computed reference LD to {}'.format(args.save_ld))
            f = gzip.open(args.save_ld, 'wb')
            pickle.dump([ld_boundaries, ref_ld_mat], f, protocol=2)
            f.close()

    print('Start MCMC ...')
    res = SDPR_gibbs(beta_margin=beta_margin, N=N, ld_boundaries=ld_boundaries, ref_ld_mat=ref_ld_mat, 
                 mcmc_samples=args.mcmc_samples, burn=args.burn, max_cluster=args.M, save_mcmc=args.save_mcmc, VS=args.VS)

    print('Done!\nWrite output to {}'.format(args.out+'.txt'))
    # res[flip_idx] = -res[flip_idx]
    ss['post_beta'] = res / np.sqrt(2*ss['A1_FREQ']*(1-ss['A1_FREQ']))
    ss.to_csv(args.out+'.txt', columns=['ID', 'A1', 'BETA' ,'post_beta'], sep="\t", index=False)


parser = argparse.ArgumentParser(prog='SDPR',
                                formatter_class=argparse.RawDescriptionHelpFormatter,
                                description="Version 0.0.1 Test Only")

parser.add_argument('--ss', type=str, required=True,
                        help='Path to cleaned summary statistics. e.g. /home/tutor/myss.txt')

parser.add_argument('--block', type=str, default=None,
                        help='Path to LD block. Prefix for .bed file output by ldetect. e.g. /home/ldetect/bed/EUR-chr')

parser.add_argument('--bfile', type=str, default=None,
                        help='Path to reference LD file. Prefix for plink .bed/.bim/.fam.')

parser.add_argument('--N', type=int, default=None, required=True,
                        help='Number of individuals in summary statistic sile.')

parser.add_argument('--M', type=int, default=20,
                        help='Maximum number of normal components in Truncated Dirichlet Process.')

parser.add_argument('--VS', type=bool, default=True, 
                        help='Whether to perform variable selection.')

parser.add_argument('--threads', type=int, default=1, 
                        help='Number of Threads used.')

parser.add_argument('--seed', type=int, 
                        help='Specify the seed for numpy random number generation.')

parser.add_argument('--mcmc_samples', type=int, default=1500,
                        help='Specify the total number of iterations in MCMC.')

parser.add_argument('--burn', type=int, default=200,
                        help='Specify the total number of iterations to be discarded before \
                        Markov Chain approached the stationary distribution.')

parser.add_argument('--save_ld', type=str, default=None,
                        help='Prefix of the location to save calculated LD Reference file \
                        in pickled and gzipped format.')

parser.add_argument('--load_ld', type=str, default=None,
                        help='Prefix of the location to load calculated LD Reference file \
                        in pickled and gzipped format.')

parser.add_argument('--save_mcmc', type=str, default=None,
                        help='Prefix of the location to save intermediate output of MCMC \
                        in pickled and gzipped format.')

parser.add_argument('--out', type=str, required=True,
                        help='Prefix of the location for the output tab deliminated .txt file.')

def main():
    if sys.version_info[0] != 2:
        print('ERROR: SDPR currently does not support Python 3')
        sys.exit(1)
    pipeline(parser.parse_args())

if __name__ == '__main__':
    main()
