import parse_genet
import argparse
import numpy as np, pandas as pd
from scipy import stats, linalg
from stick_break import *

parser = argparse.ArgumentParser(prog='SDPR',
                                formatter_class=argparse.RawDescriptionHelpFormatter,
                                description="Version 0.0.1 Test Only")

parser.add_argument('--ss', type=str, required=True,
                        help='Path to cleaned summary statistics. e.g. /home/tutor/myss.txt')

parser.add_argument('--N', type=int, required=True,
                        help='Sample size')

parser.add_argument('--out', type=str, required=True,
                        help='Prefix of the location for the output tab deliminated.')

parser.add_argument('--chrom', type=int, required=True,
	                        help='which chromosome.')

parser.add_argument('--M', type=int, default=1000, required=True, help="Maxium of truncated components")

parser.add_argument('--a', type=float, default=0.1, required=False, help="Initial value of a")

parser.add_argument('--threads', type=int, default=1, required=False, help="Threads")

parser.add_argument('--c', type=float, default=1, required=False, help="Initial value of c")

if __name__ == '__main__':

    args = parser.parse_args()

    final_h2 = []; SNP = []; A1 = []; final_postermean = []; final_pi = []
    for chrom in range(args.chrom, args.chrom+1):
        beta_margin = []; ref_ld_mat = [];
        ref_dict = parse_genet.parse_ref('/ysm-gpfs/pi/zhao/gz222/UKB_real/ref/SDPR_shrink/ldblk_1kg_eur/snpinfo_1kg_hm3', int(chrom)) 
        vld_dict = parse_genet.parse_bim('/ysm-gpfs/pi/zhao/gz222/UKB_real/genotype/Ukb_imp_v2_hm3', int(chrom)) 
        sst_dict = parse_genet.parse_sumstats(ref_dict, vld_dict, args.ss,
                                              args.N) 
        ld_blk, blk_size = parse_genet.parse_ldblk('/ysm-gpfs/pi/zhao/gz222/UKB_real/ref/SDPR_shrink/ldblk_1kg_eur/', sst_dict, int(chrom)) 
        SNP.extend(sst_dict['SNP']); A1.extend(sst_dict['A1']); beta_margin.extend(sst_dict['BETA'])
        ref_ld_mat.extend(ld_blk)
        ref_ld_mat = [ref_ld_mat[i] for i in range(len(ref_ld_mat)) if np.shape(ref_ld_mat[i])[0] != 0]
        ld_boundaries = []; add_on = 0
        for i in range(len(ref_ld_mat)):
            left = 0+add_on; right = np.shape(ref_ld_mat[i])[0]+add_on; add_on = right;
            if np.shape(ref_ld_mat[i])[0] > 0:
                ld_boundaries.append([left, right])
        beta_margin = np.array(beta_margin)
        
        mcmc_samples = 1000; M = args.M
        trace = {'pi': np.zeros(shape=(mcmc_samples, M)), 'cluster_var': np.zeros(shape=(mcmc_samples, M)),
                'alpha':[], 'num_cluster':[], 'beta':np.zeros(shape=(mcmc_samples, len(beta_margin))),
                'suffstats':[], 'h2':[], 'sigma_e2':[]}
	accept = []

        state = initial_state(data=beta_margin/args.c, num_clusters=M, ld_boundaries=ld_boundaries,
            partition=np.array([0]*len(beta_margin)), N=args.N, a0k=.5, b0k=.5)
        state['suffstats'] = [update_suffstats(state, j) for j in range(len(np.unique(state['partition'])))]
        state['cluster_var'] = [sample_sigma2(state, j, VS=True) for j in range(len(np.unique(state['partition'])))]

	state['a'] = args.a; state['c'] = args.c
	#state['A'] = ref_ld_mat
	#state['B'] = [np.dot(ref_ld_mat[j], ref_ld_mat[j]) for j in range(len(ld_boundaries))]
	#state['A'] = [np.linalg.solve(ref_ld_mat[j]+state['a']*np.identity(ref_ld_mat[j].shape[0]), ref_ld_mat[j]) for j in range(len(ld_boundaries))]
	#state['B'] = [np.dot(ref_ld_mat[j], state['A'][j]) for j in range(len(ld_boundaries))]
	#state['L'] = [ref_ld_mat[j]**2 for j in range(len(ld_boundaries))]
	
	sz = np.array(sst_dict['N'])
	ss2 = pd.read_csv('/ysm-gpfs/pi/zhao-data/gz222/SDPR_revise/HDL/summ_stats/Mc_HDL.txt', '\t')
	snp2 = np.array(ss2['rsid']) 

	ref_ld_mat2 = []
	for i in range(len(ld_boundaries)):
	    start = ld_boundaries[i][0]
	    end = ld_boundaries[i][1]
	    N = sz[start:end]
	    tmp = np.minimum(np.array([N]*len(N)), np.array([N]*len(N)).T)*1.0 / (1.1*np.array([N]*len(N)) * np.array([N]*len(N)).T)
	    np.fill_diagonal(tmp, 1.0/N)
	    ref_ld_mat2.append(ref_ld_mat[i]*tmp)
	    SNP_sub = SNP[start:end]; sz_sub = sz[start:end]
	    idx1 = np.where(np.logical_and(np.isin(SNP_sub, snp2, invert=True), sz_sub<95e3))[0]
	    idx2 = np.where(np.logical_and(np.isin(SNP_sub, snp2), sz_sub<150e3))[0]
	    for m in range(len(idx1)):
		for n in range(len(idx2)):
		    ref_ld_mat2[i][idx1[m], idx2[n]] = 0
		    ref_ld_mat2[i][idx2[n], idx1[m]] = 0
	
	for i in range(len(ld_boundaries)):
	    eig = linalg.eigvalsh(ref_ld_mat2[i])
	    if (min(eig) < 0):
		ref_ld_mat2[i] -= 1.1*min(eig)*np.identity(ref_ld_mat2[i].shape[0])
	
	state['a'] = 0; state['c'] = 1
	state['A'] = [np.linalg.solve(ref_ld_mat2[j]+state['a']*np.diag(ref_ld_mat2[j]), ref_ld_mat[j]) for j in range(len(ld_boundaries))]
	state['B'] = [np.dot(ref_ld_mat[j], state['A'][j]) for j in range(len(ld_boundaries))]

	#l = np.concatenate([np.diag(np.dot(ref_ld_mat[j], ref_ld_mat[j])) for j in range(len(ld_boundaries))])
	#state['c'] = np.sqrt(stats.linregress(l*state['N_']/len(beta_margin), beta_margin**2*state['N_'])[1])

        for i in range(mcmc_samples):
            gibbs_stick_break(state=state, VS=True, ld_boundaries=ld_boundaries, ref_ld_mat=ref_ld_mat, threads=args.threads)
            trace['beta'][i,] = state['beta']*state['eta']
            trace['h2'].append(state['h2']*state['eta']**2)
	    #state['A'] = [np.linalg.solve(ref_ld_mat[j]+state['a']*np.identity(ref_ld_mat[j].shape[0]), ref_ld_mat[j]) for j in range(len(ld_boundaries))]
	    #state['B'] = [np.dot(ref_ld_mat[j], state['A'][j]) for j in range(len(ld_boundaries))]
	    #state['a'], state['c'], tmp = sample_a(state, ld_boundaries=ld_boundaries, ref_ld_mat=ref_ld_mat)
	    state['beta_margin_'] = beta_margin / state['c']
	    #accept.append(tmp)
    #         trace['pi'][i,:] = np.array(state['pi'].values())
        
        burnin = 200
        final_h2.append(np.mean(trace['h2'][burnin:]))
        final_postermean.append(np.mean(trace['beta'][burnin:], axis=0))
    #     final_pi.append(np.mean(trace['pi'][burnin:], axis=0))

    print "h2: " + str(np.sum(final_h2)) + "\n"
    print "max poster mean: " + str(np.max(final_postermean)) + "\n"
    print "a: " + str(state['a']) + "\n"
    print "c:" + str(state['c']) + "\n"
    print "Num of components: " + str(np.sum(np.array(state['suffstats'][0].values()) > 0)) + "\n"
    print "accept rate: " + str(np.mean(accept)) + "\n"

    res = pd.DataFrame({'SNP':SNP, 'A1':A1, 'beta':np.concatenate(final_postermean)})
    res.to_csv(args.out, index=False, 
           columns=['SNP','A1','beta'], sep="\t")

