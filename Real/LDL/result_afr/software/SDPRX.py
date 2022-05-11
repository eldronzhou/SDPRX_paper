#!/usr/bin/env python
import numpy as np
from scipy import stats, linalg
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

    state['a'] = 0.1; state['c'] = 1
    state['A1'] = [np.linalg.solve(ref_ld_mat1[j]+state['a']*np.identity(ref_ld_mat1[j].shape[0]), ref_ld_mat1[j]) for j in range(len(ld_boundaries))]
    state['B1'] = [np.dot(ref_ld_mat1[j], state['A1'][j]) for j in range(len(ld_boundaries))]
    state['A2'] = [np.linalg.solve(ref_ld_mat2[j]+state['a']*np.identity(ref_ld_mat2[j].shape[0]), ref_ld_mat2[j]) for j in range(len(ld_boundaries))]
    state['B2'] = [np.dot(ref_ld_mat2[j], state['A2'][j]) for j in range(len(ld_boundaries))]

    for i in range(mcmc_samples):
        # update everything
        gibbs.gibbs_stick_break(state, rho=rho, ld_boundaries=ld_boundaries, ref_ld_mat1=ref_ld_mat1, 
            ref_ld_mat2=ref_ld_mat2, n_threads=n_threads, VS=VS)

        if (i > burn):
            trace['h2_1'].append(state['h2_1']*state['eta']**2)
            trace['h2_2'].append(state['h2_2']*state['eta']**2)

	#if (i % 100 == 0):
	 #   print 'h2_1: ' + str(state['h2_1']*state['eta']**2) + 'h2_2: ' + str(state['h2_2']*state['eta']**2) + ' max_beta1: ' + str(np.max(state['beta1']*state['eta'])) + ' max_beta2: ' + str(np.max(state['beta2']*state['eta']))

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
    
    print 'h2_1: ' + str(np.median(trace['h2_1'])) + ' h2_2: ' + str(np.median(trace['h2_2'])) + ' max_beta1: ' + str(np.max(poster_mean1)) + ' max_beta2: ' + str(np.max(poster_mean2))

    print state['pi_pop']

    if save_mcmc is not None:
        f = gzip.open(args.save_mcmc, 'wb')
        pickle.dump(trace, f, protocol=2)
        f.close()

    return poster_mean1, poster_mean2

def pipeline(args):
    
    # sanity check

    N1 = args.N1; N2 = args.N2

    print('Load summary statistics from {}'.format(args.ss1))
    ss1 = pd.read_table(args.ss1, sep=" ")
    print('Load summary statistics from {}'.format(args.ss2))
    ss2 = pd.read_table(args.ss2)

    valid = pd.read_table(args.valid, header=None)[1]

    common = list(set(ss1.SNP) & set(ss2.SNP) & set(valid))
    ss1 = ss1[ss1.SNP.isin(common)]
    ss2 = ss2[ss2.SNP.isin(common)]
    ref_ld_mat1 = []; ref_ld_mat2 = []; ld_boundaries = []
    A1 = []; SNP = []; beta_margin1 = []; beta_margin2 = []
    sz = []
    left = 0

    f = gzip.open(args.load_ld + '/chr_' + str(args.chr) +'.gz', 'r')
    ld_dict = pickle.load(f)
    f.close()
    
    snps = ld_dict[0]; a1 = ld_dict[1]; a2 = ld_dict[2]
    ref_boundary = ld_dict[3]; ref1 = ld_dict[4]; ref2 = ld_dict[5]

    ref = pd.DataFrame({'SNP':snps, 'A1':a1, 'A2':a2})
    tmp_ss1 = pd.merge(ref, ss1, on="SNP", how="left")
    tmp_ss2 = pd.merge(ref, ss2, on="SNP", how="left")

    for i in range(len(ref_boundary)):
	tmp_blk_ss1 = tmp_ss1.iloc[ref_boundary[i][0]:ref_boundary[i][1]]
	tmp_blk_ss2 = tmp_ss2.iloc[ref_boundary[i][0]:ref_boundary[i][1]]
	tmp_beta1 = tmp_ss1.iloc[ref_boundary[i][0]:ref_boundary[i][1]]['BETA'] / np.sqrt(tmp_ss1.iloc[ref_boundary[i][0]:ref_boundary[i][1]]['N'])
	tmp_beta2 = tmp_ss2.iloc[ref_boundary[i][0]:ref_boundary[i][1]]['BETA'] / np.sqrt(tmp_ss2.iloc[ref_boundary[i][0]:ref_boundary[i][1]]['N'] )
	idx1_ss1 = np.logical_and(tmp_blk_ss1['A1_x'] == tmp_blk_ss1['A1_y'], tmp_blk_ss1['A2_x'] == tmp_blk_ss1['A2_y'])
	idx2_ss1 = np.logical_and(tmp_blk_ss1['A1_x'] == tmp_blk_ss1['A2_y'], 
		tmp_blk_ss1['A2_x'] == tmp_blk_ss1['A1_y'])
	tmp_beta1[idx2_ss1] = -tmp_beta1[idx2_ss1] 
	idx1_ss2 = np.logical_and(tmp_blk_ss2['A1_x'] == tmp_blk_ss2['A1_y'], 
		tmp_blk_ss2['A2_x'] == tmp_blk_ss2['A2_y'])
	idx2_ss2 = np.logical_and(tmp_blk_ss2['A1_x'] == tmp_blk_ss2['A2_y'], 
		tmp_blk_ss2['A2_x'] == tmp_blk_ss2['A1_y'])
	tmp_beta2[idx2_ss2] = -tmp_beta2[idx2_ss2] 
	idx1 = np.logical_or(idx1_ss1, idx2_ss1)
	idx2 = np.logical_or(idx1_ss2, idx2_ss2)
	idx = np.logical_and(idx1, idx2)
	if np.sum(idx) == 0:
	    continue
	ref_ld_mat1.append(ref1[i][idx,:][:,idx])
	ref_ld_mat2.append(ref2[i][idx,:][:,idx])
	ld_boundaries.append([left, left+np.sum(idx)])
	sz.extend(list(tmp_blk_ss1[idx].N))
	SNP.extend(list(tmp_blk_ss1[idx].SNP)); A1.extend(list(tmp_blk_ss1[idx].A1_x))
	beta_margin1.extend(list(tmp_beta1[idx])); beta_margin2.extend(list(tmp_beta2[idx]))
	left += np.sum(idx)

    print str(args.chr) + ':' + str(len(beta_margin1))
    
    sz = np.array(sz)
    ss2 = pd.read_csv('/ysm-gpfs/pi/zhao-data/gz222/SDPR_revise/LDL/summ_stats/Mc_LDL.txt', '\t')
    snp2 = np.array(ss2['rsid'])
    
    for i in range(len(ld_boundaries)):
	start = ld_boundaries[i][0]
	end = ld_boundaries[i][1]
	N = sz[start:end]
	tmp = np.minimum(np.array([N]*len(N)), np.array([N]*len(N)).T)*1.0 / (1.1*np.array([N]*len(N)) * np.array([N]*len(N)).T)
	np.fill_diagonal(tmp, 1.0/N)
	ref_ld_mat1[i] = ref_ld_mat1[i]*tmp
	SNP_sub = SNP[start:end]; sz_sub = sz[start:end]
	idx1 = np.where(np.logical_and(np.isin(SNP_sub, snp2, invert=True), sz_sub<100e3))[0]
	idx2 = np.where(np.logical_and(np.isin(SNP_sub, snp2), sz_sub<150e3))[0]
	for m in range(len(idx1)):
	    for n in range(len(idx2)):
		ref_ld_mat1[i][idx1[m], idx2[n]] = 0
		ref_ld_mat1[i][idx2[n], idx1[m]] = 0						
		
    for i in range(len(ld_boundaries)):				
	eig = linalg.eigvalsh(ref_ld_mat1[i])
	if (min(eig) < 0):
	    ref_ld_mat1[i] -= 1.1*min(eig)*np.identity(ref_ld_mat1[i].shape[0])

    print('Start MCMC ...')
    res1, res2 = SDPRX_gibbs(beta_margin1=np.array(beta_margin1)/args.c1, beta_margin2=np.array(beta_margin2)/args.c2, N1=N1, N2=N2, rho=args.rho, ld_boundaries=ld_boundaries, ref_ld_mat1=ref_ld_mat1, 
                 ref_ld_mat2=ref_ld_mat2, mcmc_samples=args.mcmc_samples, 
                 burn=args.burn, max_cluster=args.M, save_mcmc=args.save_mcmc, n_threads=args.threads, VS=args.VS)

    print('Done!\nWrite output to {}'.format(args.out+'.txt'))
    
    out1 = pd.DataFrame({'SNP':SNP, 'A1':A1, 'post_beta':res1})
    out1.to_csv(args.out+'_1.txt', columns=['SNP', 'A1', 'post_beta'], sep="\t", index=False)
    out2 = pd.DataFrame({'SNP':SNP, 'A1':A1, 'post_beta':res2})
    out2.to_csv(args.out+'_2.txt', columns=['SNP', 'A1', 'post_beta'], sep="\t", index=False)


parser = argparse.ArgumentParser(prog='SDPR',
                                formatter_class=argparse.RawDescriptionHelpFormatter,
                                description="Version 0.0.1 Test Only")

parser.add_argument('--ss1', type=str, required=True,
                        help='Path to cleaned summary statistics 1. e.g. /home/tutor/myss.txt')

parser.add_argument('--ss2', type=str, required=True,
                        help='Path to cleaned summary statistics 2. e.g. /home/tutor/myss.txt')

parser.add_argument('--block', type=str, default=None,
                        help='Path to LD block. Prefix for .bed file output by ldetect. e.g. /home/ldetect/bed/EUR-chr')

parser.add_argument('--valid', type=str, default=None,
                        help='Path to valid .bim file..')

parser.add_argument('--N1', type=int, default=None, required=True,
                        help='Number of individuals in summary statistic sile 1.')

parser.add_argument('--N2', type=int, default=None, required=True,
                        help='Number of individuals in summary statistic sile 2.')

parser.add_argument('--chr', type=int, default=None, required=True,
	                        help='Chromosome.')

parser.add_argument('--c1', type=float, default=1.0,
	                                help='C1.')

parser.add_argument('--c2', type=float, default=1.0,
	                                help='C2.')

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
