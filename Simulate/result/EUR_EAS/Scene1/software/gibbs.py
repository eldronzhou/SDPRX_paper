#!/usr/bin/env python

import numpy as np
from scipy import stats, special, linalg 
from collections import Counter

# Initialize state 
def initial_state(data, num_clusters, N, alpha=1.0, a0k=.01, b0k=.01):
    cluster_ids = range(num_clusters)

    state = {
        'cluster_ids_': cluster_ids, 
        'beta_margin_': np.array(data),
        'b': None,
        'N_': N,
        'beta': np.zeros(len(data)),
        'num_clusters_': num_clusters,
        'num_points_': len(data),
        'alpha': alpha,
        'hyperparameters_': {
            "a0k": a0k,
            "b0k": b0k,
            "a0": 0.1,
            "b0": 0.1,
        },
        'suffstats': {cid: None for cid in cluster_ids},
        'assignment': np.array([np.random.choice(cluster_ids) for _ in data]),
        'pi': {cid: alpha / num_clusters for cid in cluster_ids},
        'V': None,
        'cluster_var': None,
	'h2': 0,
	'eta': 1
    }

    return state

# Calculate number of observations per cluster
def update_suffstats(state):
    suff_stats = dict(Counter(state['assignment']))
    for i in range(state['num_clusters_']):
        if i not in suff_stats.keys():
            suff_stats[i] = 0
    return suff_stats

# Sample cluster assignment
def sample_sigma2(state, VS=True):
    b = np.zeros(state['num_clusters_'])
    a = np.array(state['suffstats'].values() ) / 2.0 + state['hyperparameters_']['a0k'] 
    for cluster_id in state['cluster_ids_']:
        if (cluster_id == 0) and VS is True:
            continue
        b[cluster_id] = np.sum(state['beta'][np.where(state['assignment'] == cluster_id)]**2) / 2.0 + \
        state['hyperparameters_']['b0k']

    if VS is True:
        out = dict(zip(range(1, state['num_clusters_']), stats.invgamma(a=a[1:], scale=b[1:]).rvs()))
        out[0] = 0
    else: 
        out = dict(zip(range(0, state['num_clusters_']), stats.invgamma(a=a, scale=b).rvs()))
    return out

# Vetorized version of np.random.choice
def vectorized_random_choice(prob_matrix, items):
    # https://stackoverflow.com/a/34190035/8834998
    s = prob_matrix.cumsum(axis=0)
    r = np.random.rand(prob_matrix.shape[1])
    k = (s < r).sum(axis=0)
    k[np.where(k == len(items))] = len(items) - 1
    return items[k]

# Calculate b_j for j th block 
def calc_b(j, state, ld_boundaries, ref_ld_mat):
    start_i = ld_boundaries[j][0]
    end_i = ld_boundaries[j][1]
    ref_ld = ref_ld_mat[j]
    shrink_ld = ref_ld
    # shrink_ld = .9*ref_ld  + (1-.9)*np.identity(ref_ld.shape[0])
    b = state['eta'] * state['beta_margin_'][start_i:end_i] - \
        state['eta']**2 * (np.dot(shrink_ld, state['beta'][start_i:end_i]) - state['beta'][start_i:end_i])
    return b

# Sample z
def sample_assignment(j, state, ld_boundaries, VS=True):
    start_i = ld_boundaries[j][0]
    end_i = ld_boundaries[j][1]
    if VS is True:
        m = state['num_clusters_']; N = state['N_']
	b = state['b'][start_i:end_i]
        cluster_var = np.array( state['cluster_var'].values() )
	cluster_var = cluster_var*1.0
        C = -.5 * np.log(state['eta']**2 * N*cluster_var + 1) + np.log( np.array(state['pi'].values()) + 1e-40 )
        a = np.tile( (N*b )**2 , (m-1, 1)).T / (2* ( state['eta']**2 * N + 1.0/cluster_var[1:]) ) 
        log_prob_mat = np.insert(a, 0, 0, axis=1) + C
        logexpsum = special.logsumexp(log_prob_mat, axis=1)
        prob_mat = np.exp(log_prob_mat - np.tile(logexpsum, (m,1)).T)
    else:
        m = state['num_clusters_']
        sd = np.sqrt(state['cluster_var'].values())
        mat = np.tile(state['beta'], (m, 1)).T / sd
        prob_mat = state['pi'].values() * stats.norm.pdf(mat) / sd
        prob_mat /= prob_mat.sum(axis=1)[:,np.newaxis]
    return vectorized_random_choice(prob_mat.T, np.array(state['cluster_ids_']))

# Sample V
def sample_V(state):
    m = state['num_clusters_']; suffstats = np.array(state['suffstats'].values())
    a = 1 + suffstats[:-1]
    b = state['alpha'] + np.cumsum(suffstats[::-1])[:-1][::-1]
    sample_val = stats.beta(a=a, b=b).rvs()
    if 1 in sample_val:
        idx = np.argmax(sample_val == 1)
        sample_val[idx+1:] = 0
        sample_return = dict(zip(range(m-1), sample_val))
        sample_return[m-1] = 0
    else:
        sample_return = dict(zip(range(m-1), sample_val))
        sample_return[m-1] = 1
    return sample_return

# Compute pi
def update_p(state):
    m = state['num_clusters_']
    pi = {cid: None for cid in state['cluster_ids_']}
    pi[0] = state['V'][0]; V = state['V'].values()
    for i in range(1, m-1):
        pi[i] = np.prod( 1 - np.array(V[0:i]) ) * state['V'][i]
    pi[m-1] = 1 - np.sum(pi.values()[0:(m-1)])

    # last p may be less than 0 due to rounding error
    if pi[m-1] < 0: 
        pi[m-1] = 0
    return pi

# Sample alpha
def sample_alpha(state):
    m = np.size(np.where( np.array(state['V'].values()) != 0)); V = state['V'].values()
    a = state['hyperparameters_']['a0'] + m - 1
    b = state['hyperparameters_']['b0'] - np.sum( np.log( 1 - np.array(V[0:(m-1)]) ) )
    return stats.gamma(a=a, scale=1.0/b).rvs()

# Sample X ~ MVN(mu, cov) based on Cholesky factorization (scipy used SVD, robust but slower)
def sample_MVN(mu, cov):
    rv = stats.norm.rvs(size=mu.shape[0])
    C = linalg.cholesky(cov, lower=True)
    return np.dot(C, rv) + mu

# Sample beta
def sample_beta(j, state, ld_boundaries, ref_ld_mat, VS=True):
    start_i = ld_boundaries[j][0]
    end_i = ld_boundaries[j][1]
    ref_ld = ref_ld_mat[j]
    shrink_ld = ref_ld
    # shrink_ld = .9*ref_ld  + (1-.9)*np.identity(ref_ld.shape[0])
    cluster_var = np.array(state['cluster_var'].values()); inv_s2 = state['N_']
    beta_margin = state['beta_margin_']
    cluster_var = cluster_var*1.0

    if VS is True:
        idx = state['assignment'][start_i:end_i] != 0
        beta = np.zeros(end_i - start_i)

        if not any(idx):
            # all SNPs in this block are non-causal
            pass
        elif sum(idx) == 1:
            # only one SNP in this block is causal
            var_k = cluster_var[ state['assignment'][start_i:end_i][ np.where(idx == 1)] ]
            const = var_k / (var_k*state['eta']**2 + 1.0/inv_s2)
            bj = state['b'][start_i:end_i][ np.where(idx == 1) ]
            beta[idx == 1] = np.sqrt(const*1.0/inv_s2)*stats.norm.rvs() + const*bj
        else:
            # only sample SNPs that are causal
            shrink_ld = shrink_ld[idx,:][:,idx]
            mat = state['eta']**2*inv_s2*shrink_ld + np.diag(1.0 / cluster_var[state['assignment'][start_i:end_i][idx]])
            chol, low = linalg.cho_factor(mat, overwrite_a=False)
            cov_mat = linalg.cho_solve((chol, low), np.eye(chol.shape[0])) 
            mu = state['eta']*inv_s2*np.dot(cov_mat, beta_margin[start_i:end_i][idx])
            beta[idx == 1] = sample_MVN(mu, cov_mat)
        return beta

    else:    
        mat = inv_s2*shrink_ld + np.diag(1.0 / cluster_var[state['assignment'][start_i:end_i]])
        chol, low = linalg.cho_factor(mat, overwrite_a=False)
        cov_mat = linalg.cho_solve((chol, low), np.eye(chol.shape[0])) 
        mu = inv_s2*np.dot(cov_mat, beta_margin[start_i:end_i])
        return sample_MVN(mu, cov_mat)

def compute_varg(j, state, ld_boundaries, ref_ld_mat):
    start_i = ld_boundaries[j][0]
    end_i = ld_boundaries[j][1]
    ref_ld = ref_ld_mat[j]
    shrink_ld = ref_ld
    # shrink_ld = .9*ref_ld  + (1-.9)*np.identity(ref_ld.shape[0])
    idx = state['assignment'][start_i:end_i] != 0
    if not any(idx):
        return 0
    else:
        beta = state['beta'][start_i:end_i]
        shrink_ld = shrink_ld
        return np.dot(beta.T, np.dot(shrink_ld, beta))

def sample_eta(state):
    mu = np.sum(state['beta_margin_']*state['beta']) / (state['h2']+1e-6/state['N_'])
    var = 1.0/(state['N_']*state['h2']+1e-6)
    return np.sqrt(var)*stats.norm.rvs() + mu

# Wrap up gibbs sampler
def gibbs(state, ld_boundaries, ref_ld_mat, VS=True):
    #state['cluster_var'] = sample_sigma2(state, VS)
    state['b'] = np.concatenate([calc_b(j=j, state=state, ld_boundaries=ld_boundaries, ref_ld_mat=ref_ld_mat) 
                                for j in range(len(ld_boundaries))])
    state['assignment'] = np.concatenate([sample_assignment(j, state, ld_boundaries, VS) for j in range(len(ld_boundaries))])
    #state['suffstats'] = update_suffstats(state)
    #state['V'] = sample_V(state)
    #state['pi'] = update_p(state)
    #state['alpha'] = sample_alpha(state)
    state['beta'] = np.concatenate([sample_beta(j=j, state=state, ld_boundaries=ld_boundaries, ref_ld_mat=ref_ld_mat, VS=True)
                                   for j in range(len(ld_boundaries))])
    state['h2'] = np.sum([compute_varg(j=j, state=state, ld_boundaries=ld_boundaries, ref_ld_mat=ref_ld_mat) for j in range(len(ld_boundaries))])
    #state['eta'] = sample_eta(state)

# def log_likelihood(state):
#     m = state['num_clusters_']
#     sd = np.sqrt(state['cluster_var'].values())
#     a = stats.norm.pdf(np.tile(state['beta'], (m, 1)).T / sd)  / sd
#     return np.sum(np.log(a.dot(np.array(state['pi'].values()))))


