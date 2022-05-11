#!/usr/bin/env python
import numpy as np
from scipy import stats
import cPickle as pickle
import gzip
import gibbs, ld
import argparse
import sys, util
import pandas as pd


def SDPRX_gibbs(beta_margin1, beta_margin2, N1, N2, rho, ld_boundaries, ref_ld_mat1, ref_ld_mat2, mcmc_samples, 
    burn, max_cluster, save_mcmc, n_threads, VS=True):
    M = [1, 1000, 1000, 1000]
    trace = {'alpha':[], 'num_cluster':[], 'beta1':np.zeros(shape=(mcmc_samples, len(beta_margin1))),
        'beta2':np.zeros(shape=(mcmc_samples, len(beta_margin2))),
        'suffstats':[], 'h2_1':[], 'h2_2':[]}

    # initialize
    state = gibbs.initial_state(data1=beta_margin1, data2=beta_margin2, ld_boundaries=ld_boundaries, M=M, N1=N1, N2=N2, a0k=.5, b0k=.5)
    state['suffstats'] = gibbs.update_suffstats(state)
    state['cluster_var'] = gibbs.sample_sigma2(state, rho=rho, VS=True)

    for i in range(mcmc_samples):
        # update everything
        gibbs.gibbs_stick_break(state, rho=rho, ld_boundaries=ld_boundaries, ref_ld_mat1=ref_ld_mat1, 
            ref_ld_mat2=ref_ld_mat2, n_threads=n_threads, VS=VS)

        if (i > burn):
            trace['h2_1'].append(state['h2_1']*state['eta']**2)
            trace['h2_2'].append(state['h2_2']*state['eta']**2)

	print 'h2_1: ' + str(state['h2_1']*state['eta']**2) + 'h2_2: ' + str(state['h2_2']*state['eta']**2)

        # record the result
        trace['beta1'][i,] = state['beta1']*state['eta']
        trace['beta2'][i,] = state['beta2']*state['eta']
        # trace['pi'][i,:] = np.array(state['pi'].values())
        # trace['cluster_var'][i,:] = np.array(state['cluster_var'].values())
        # trace['alpha'].append(state['alpha'])
        # trace['num_cluster'].append( np.sum(np.array(state['pi'].values()) > .0001) )
        # trace['suffstats'].append(state['suffstats'])

        util.progressBar(value=i+1, endvalue=mcmc_samples)

    # calculate posterior average
    poster_mean1 = np.mean(trace['beta1'][burn:mcmc_samples], axis=0)
    poster_mean2 = np.mean(trace['beta2'][burn:mcmc_samples], axis=0)
    
    print 'h2_1: ' + str(np.median(trace['h2_1'])) + ' h2_2: ' + str(np.median(trace['h2_2']))

    print state['pi_pop']

    if save_mcmc is not None:
        f = gzip.open(args.save_mcmc, 'wb')
        pickle.dump(trace, f, protocol=2)
        f.close()

    return poster_mean1, poster_mean2

def pipeline(args):
    
    # sanity check

    if args.bfile is not None and args.load_ld is not None:
        raise ValueError('Both --bfile and --load_ld flags were set. \
            Please use only one of them.')

    if args.bfile is None and args.load_ld is None:
        raise ValueError('Both --bfile and --load_ld flags were not set. \
            Please use one of them.')

    N1 = args.N1; N2 = args.N2

    print('Load summary statistics from {}'.format(args.ss1))
    ss1 = pd.read_table(args.ss1)
    print('Load summary statistics from {}'.format(args.ss2))
    ss2 = pd.read_table(args.ss2)
    # beta_margin = np.array(np.sign(ss['BETA']) * abs(stats.norm.ppf(ss['P'] / 2.0)) / np.sqrt(N))
    beta_margin1 = np.array(ss1['T_STAT']) / np.sqrt(N1)
    beta_margin2 = np.array(ss2['T_STAT']) / np.sqrt(N2)
    assert np.all(np.isreal(beta_margin1)), 'Something wrong with summary stats 1.'
    assert np.all(np.isreal(beta_margin2)), 'Something wrong with summary stats 2.'

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
        ref_ld_mat1 = ld_dict[1]
        ref_ld_mat2 = ld_dict[2]
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
    res1, res2 = SDPRX_gibbs(beta_margin1=beta_margin1, beta_margin2=beta_margin2, N1=N1, N2=N2, rho=args.rho, ld_boundaries=ld_boundaries, ref_ld_mat1=ref_ld_mat1, 
                 ref_ld_mat2=ref_ld_mat2, mcmc_samples=args.mcmc_samples, 
                 burn=args.burn, max_cluster=args.M, save_mcmc=args.save_mcmc, n_threads=args.threads, VS=args.VS)

    print('Done!\nWrite output to {}'.format(args.out+'.txt'))
    # res[flip_idx] = -res[flip_idx]
    ss1['post_beta'] = res1 / np.sqrt(2*ss1['A1_FREQ']*(1-ss1['A1_FREQ']))
    ss1.to_csv(args.out+'_1.txt', columns=['ID', 'A1', 'BETA' ,'post_beta'], sep="\t", index=False)
    ss2['post_beta'] = res2 / np.sqrt(2*ss2['A1_FREQ']*(1-ss2['A1_FREQ']))
    ss2.to_csv(args.out+'_2.txt', columns=['ID', 'A1', 'BETA' ,'post_beta'], sep="\t", index=False)


parser = argparse.ArgumentParser(prog='SDPR',
                                formatter_class=argparse.RawDescriptionHelpFormatter,
                                description="Version 0.0.1 Test Only")

parser.add_argument('--ss1', type=str, required=True,
                        help='Path to cleaned summary statistics 1. e.g. /home/tutor/myss.txt')

parser.add_argument('--ss2', type=str, required=True,
                        help='Path to cleaned summary statistics 2. e.g. /home/tutor/myss.txt')

parser.add_argument('--block', type=str, default=None,
                        help='Path to LD block. Prefix for .bed file output by ldetect. e.g. /home/ldetect/bed/EUR-chr')

parser.add_argument('--bfile', type=str, default=None,
                        help='Path to reference LD file. Prefix for plink .bed/.bim/.fam.')

parser.add_argument('--N1', type=int, default=None, required=True,
                        help='Number of individuals in summary statistic sile 1.')

parser.add_argument('--N2', type=int, default=None, required=True,
                        help='Number of individuals in summary statistic sile 2.')

parser.add_argument('--rho', type=float, default=0, required=True,
                        help='Transethnic genetic correlation.')

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
