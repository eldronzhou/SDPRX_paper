import numpy as np
from scipy import stats, special, linalg
from collections import Counter
from joblib import Parallel, delayed

def initial_state(data1, data2, M, N1, N2, ld_boundaries, alpha=1.0, a0k=.5, b0k=.5):
    num_clusters = sum(M) 
    cluster_ids = range(num_clusters)

    state = {
        'cluster_ids_': cluster_ids, 
        'beta_margin1_': data1,
        'beta_margin2_': data2,
        'b1': np.zeros(len(data1)), 
        'b2': np.zeros(len(data2)), 
        'N1_': N1,
        'N2_': N2,
        'beta1': np.zeros(len(data1)),
        'beta2': np.zeros(len(data2)),
        'num_clusters_': num_clusters,
        'alpha': np.array([alpha]*4),
        'hyperparameters_': {
            "a0k": a0k,
            "b0k": b0k,
            "a0": 0.1,
            "b0": 0.1,
        },
        'suffstats': np.array([0]*(num_clusters)),
        'assignment': np.random.randint(num_clusters, size=len(data1)),
        'population': np.array([0]*4),
        'pi': np.array([alpha / num_clusters]*num_clusters),
        'pi_pop': np.array([.25, .25, .25, .25]),
        'pi_cluster': [np.array(1.0/M[i]) for i in range(len(M))],
        'V': [np.array([0]*M[i]) for i in range(len(M))],
        'cluster_var': np.array([0]*(num_clusters)),
        'varg1': np.array([0.0]*len(ld_boundaries)),
        'varg2': np.array([0.0]*len(ld_boundaries)),
        'h2_1': 0,
        'h2_2': 0,
        'eta': 1
    }
    
    # define indexes 
    state['population'][0] = 1 # null
    state['population'][1] = M[1] + 1 # pop 1 specific
    state['population'][2] = M[1] + M[2] + 1 # pop 2 specific
    state['population'][3] = num_clusters # shared with correlation
    return state


def calc_b(j, state, ld_boundaries, ref_ld_mat1, ref_ld_mat2):
    start_i = ld_boundaries[j][0]
    end_i = ld_boundaries[j][1]
    ref_ld1 = ref_ld_mat1[j]
    ref_ld2 = ref_ld_mat2[j]
    shrink_ld1 = ref_ld1; shrink_ld2 = ref_ld2
#     shrink_ld = .999*ref_ld  + (1-.999)*np.identity(ref_ld.shape[0])
#     b = state['eta'] * state['beta_margin_'][start_i:end_i] - \
#         state['eta']**2 * (np.dot(shrink_ld, state['beta'][start_i:end_i]) - state['beta'][start_i:end_i])
    b1 = state['eta']*np.dot(state['A1'][j], state['beta_margin1_'][start_i:end_i]) - state['eta']**2 * \
    (np.dot(state['B1'][j], state['beta1'][start_i:end_i]) - np.diag(state['B1'][j])*state['beta1'][start_i:end_i])
    b2 = state['eta']*np.dot(state['A2'][j], state['beta_margin2_'][start_i:end_i]) - state['eta']**2 * \
    (np.dot(state['B2'][j], state['beta2'][start_i:end_i]) - np.diag(state['B2'][j])*state['beta2'][start_i:end_i])
    state['b1'][start_i:end_i] = b1
    state['b2'][start_i:end_i] = b2


def vectorized_random_choice(prob_matrix, items):
    s = prob_matrix.cumsum(axis=0)
    r = np.random.rand(prob_matrix.shape[1])
    k = (s < r).sum(axis=0)
    k[np.where(k == len(items))] = len(items) - 1
    return items[k]

def sample_assignment(j, ld_boundaries, ref_ld_mat1, ref_ld_mat2, state, VS, rho):
    start_i = ld_boundaries[j][0]
    end_i = ld_boundaries[j][1]
    m = state['num_clusters_']; N1 = state['N1_']; N2 = state['N2_']
    b1 = state['b1'][start_i:end_i].reshape((end_i-start_i, 1))
    b2 = state['b2'][start_i:end_i].reshape((end_i-start_i, 1))
    B1 = state['B1'][j]; B2 = state['B2'][j]

    log_prob_mat = np.zeros((len(b1), m))
    
    # null or population 1 specific
#     idx = np.where(state['population'] <= 1)[0]
    idx = range(0, state['population'][1])
    cluster_var = state['cluster_var'][idx]
    pi = np.array(state['pi'])[idx]
    C = -.5 * np.log(state['eta']**2*N1*np.outer(np.diag(B1), cluster_var) + 1) + \
        np.log( pi + 1e-40 )
    a = (N1*b1)**2 / (2 * np.add.outer(state['eta']**2 * N1 * np.diag(B1),  1.0/cluster_var[1:]) )
    log_prob_mat[:,idx] = np.insert(a, 0, 0, axis=1) + C
    
    # population 2 specific
#     idx = np.where(state['population'] == 2)[0]
    idx = range(state['population'][1], state['population'][2])
    cluster_var = state['cluster_var'][idx]
    pi = np.array(state['pi'])[idx]
    C = -.5 * np.log(state['eta']**2*N2*np.outer(np.diag(B2), cluster_var) + 1) + \
        np.log( pi + 1e-40 )
    a = (N2*b2)**2 / (2 * np.add.outer(state['eta']**2 * N2 * np.diag(B2),  1.0/cluster_var) )
    log_prob_mat[:,idx] = a + C
    
    # shared with correlation
#     idx = np.where(state['population'] == 3)[0]
    idx = range(state['population'][2], state['population'][3])
    cluster_var = state['cluster_var'][idx]
    pi = np.array(state['pi'])[idx]
    
    ak1 = np.add.outer(.5*N1*state['eta']**2*np.diag(B1), .5/((1-rho**2)*cluster_var))
    ak2 = np.add.outer(.5*N2*state['eta']**2*np.diag(B2), .5/((1-rho**2)*cluster_var))
    ck = rho / ((1-rho**2)*cluster_var)
    mu1 = (2*b1 + (N2*1.0/N1)*b2*(ck/ak2)) / (4*ak1/N1 - ck**2/(ak2*N1))
    mu2 = (2*b2 + (N1*1.0/N2)*b1*(ck/ak1)) / (4*ak2/N2 - ck**2/(ak1*N2))
        
    C = -.5*np.log(4*ak1*ak2-ck**2) - .5*np.log(1-rho**2) - np.log(cluster_var) + np.log( pi + 1e-40 )
        
    a = ak1*mu1*mu1 + ak2*mu2*mu2 - ck*mu1*mu2
    log_prob_mat[:,idx] = a + C
    
    logexpsum = special.logsumexp(log_prob_mat, axis=1).reshape((len(b1), 1))
    prob_mat = np.exp(log_prob_mat - logexpsum)
    
    assignment = vectorized_random_choice(prob_mat.T, np.array(state['cluster_ids_']))
    
    return assignment


def sample_beta(j, state, ld_boundaries, ref_ld_mat1, ref_ld_mat2, rho, VS=True):
    start_i = ld_boundaries[j][0]
    end_i = ld_boundaries[j][1]
    N1 = state['N1_']; N2 = state['N2_']
    beta_margin1 = state['beta_margin1_'][start_i:end_i]; beta_margin2 = state['beta_margin2_'][start_i:end_i]
    A1 = state['A1'][j]; B1 = state['B1'][j]
    A2 = state['A2'][j]; B2 = state['B2'][j]
    
    assignment = state['assignment'][start_i:end_i]
    cluster_var = state['cluster_var'][state['assignment'][start_i:end_i]]
    
    beta1 = np.zeros(len(beta_margin1)); beta2 = np.zeros(len(beta_margin2))
    # null
    idx = state['assignment'][start_i:end_i] == 0
     
    # pop1 specific
#     idx1 = state['assignment'][start_i:end_i] == 1
    idx1 = (state['assignment'][start_i:end_i] >= 1) & (state['assignment'][start_i:end_i] < state['population'][1])
    
    # pop2 speicifc 
#     idx2 = state['assignment'][start_i:end_i] == 2
    idx2 = (state['assignment'][start_i:end_i] >= state['population'][1]) \
        & (state['assignment'][start_i:end_i] < state['population'][2])
    
    # shared with correlation
#     idx3 = state['assignment'][start_i:end_i] == 3
    idx3 = (state['assignment'][start_i:end_i] >= state['population'][2]) \
        & (state['assignment'][start_i:end_i] < state['population'][3])
    
    idx_pop1 = np.logical_or(idx1, idx3)
    idx_pop2 = np.logical_or(idx2, idx3)

    if all(idx):
        # all SNPs in this block are non-causal
        pass
    elif sum(idx1) == 1 and sum(idx3) == 0:
        # only one SNP in this block is causal
        var_k = cluster_var[idx1]
        const = var_k / (var_k*state['eta']**2*np.squeeze(B1[idx1,:][:,idx1]) + 1.0/N1)
        bj = state['b1'][start_i:end_i][idx1]
#             beta[idx == 1] = np.sqrt(var_k*1.0/inv_s2)*stats.norm.rvs() + const*bj
        beta1[idx1 == 1] = np.sqrt(const*1.0/N1)*stats.norm.rvs() + const*bj
    elif sum(idx2) == 1 and sum(idx3) == 0:
        var_k = cluster_var[idx2]
        const = var_k / (var_k*state['eta']**2*np.squeeze(B2[idx2,:][:,idx2]) + 1.0/N2)
        bj = state['b2'][start_i:end_i][idx2]
        beta2[idx2 == 1] = np.sqrt(const*1.0/N2)*stats.norm.rvs() + const*bj
    else:
        # two population LD matrix
        shrink_ld = np.block([[N1*B1[idx_pop1,:][:,idx_pop1], np.zeros((sum(idx_pop1), sum(idx_pop2)))],
             [np.zeros((sum(idx_pop2), sum(idx_pop1))), N2*B2[idx_pop2,:][:,idx_pop2]]])
        
        # variance covariance matrix for beta
        idx_cor1 = np.where(state['assignment'][start_i:end_i][idx_pop1] >= state['population'][2])[0]
#         idx_cor1 = np.where(state['assignment'][start_i:end_i][idx_pop1] > 2)[0]
        idx_cor2 = np.where(state['assignment'][start_i:end_i][idx_pop2] >= state['population'][2])[0]
#         idx_cor2 = np.where(state['assignment'][start_i:end_i][idx_pop2] > 2)[0]

        diag1 = np.diag(1.0/cluster_var[idx_pop1])
        cor1 = np.zeros((sum(idx_pop1), sum(idx_pop2)))
        diag2 = np.diag(1.0/cluster_var[idx_pop2])
        cor2 = np.zeros((sum(idx_pop2), sum(idx_pop1)))
        
#         rho = 0
        for i in range(len(idx_cor1)):
            cor1[idx_cor1[i],idx_cor2[i]] = -rho/(1-rho**2)*diag1[idx_cor1[i],idx_cor1[i]]
            cor2[idx_cor2[i],idx_cor1[i]] = -rho/(1-rho**2)*diag1[idx_cor1[i],idx_cor1[i]]
            diag1[idx_cor1[i],idx_cor1[i]] = 1.0/(1-rho**2)*diag1[idx_cor1[i],idx_cor1[i]]
            diag2[idx_cor2[i],idx_cor2[i]] = 1.0/(1-rho**2)*diag2[idx_cor2[i],idx_cor2[i]]
        
        var_mat = np.block([[diag1, cor1],
                    [cor2, diag2]])
        
        mat = state['eta']**2*shrink_ld + var_mat
        
        chol, low = linalg.cho_factor(mat, overwrite_a=False)
        cov_mat = linalg.cho_solve((chol, low), np.eye(chol.shape[0])) 
        
        # A matrix
        A_gamma = np.concatenate([N1*np.dot(A1[idx_pop1,:], state['beta_margin1_'][start_i:end_i]), 
               N2*np.dot(A2[idx_pop2,:], state['beta_margin2_'][start_i:end_i])])
        
        mu = state['eta']*np.dot(cov_mat, A_gamma)
        beta_tmp = sample_MVN(mu, cov_mat)
        beta1[idx_pop1] = beta_tmp[0:sum(idx_pop1)]
        beta2[idx_pop2] = beta_tmp[sum(idx_pop1):]
        
    state['beta1'][start_i:end_i] = beta1
    state['beta2'][start_i:end_i] = beta2


def sample_MVN(mu, cov):
    rv = stats.norm.rvs(size=mu.shape[0])
    C = linalg.cholesky(cov, lower=True)
    return np.dot(C, rv) + mu


def compute_varg(j, state, ld_boundaries, ref_ld_mat1, ref_ld_mat2):
    start_i = ld_boundaries[j][0]
    end_i = ld_boundaries[j][1]
    ref_ld1 = ref_ld_mat1[j]
    ref_ld2 = ref_ld_mat2[j]
    state['varg1'][j] = np.sum(state['beta1'][start_i:end_i] * np.dot(ref_ld1, state['beta1'][start_i:end_i]))
    state['varg2'][j] = np.sum(state['beta2'][start_i:end_i] * np.dot(ref_ld2, state['beta2'][start_i:end_i]))


def calc_num(j, state, ld_boundaries):
    start_i = ld_boundaries[j][0]
    end_i = ld_boundaries[j][1]
    A1 = state['A1'][j]
    A2 = state['A2'][j]
    return state['N1_']*np.dot(state['beta_margin1_'][start_i:end_i], np.dot(A1, (state['beta1'][start_i:end_i]))) + \
        state['N2_']*np.dot(state['beta_margin2_'][start_i:end_i], np.dot(A2, (state['beta2'][start_i:end_i])))

def calc_denum(j, state, ld_boundaries):
    start_i = ld_boundaries[j][0]
    end_i = ld_boundaries[j][1]
    B1 = state['B1'][j]
    B2 = state['B2'][j]
    beta1 = state['beta1'][start_i:end_i]
    beta2 = state['beta2'][start_i:end_i]
    return state['N1_']*np.dot(beta1, np.dot(B1, beta1)) + state['N2_']*np.dot(beta2, np.dot(B2, beta2))

def sample_eta(state, ld_boundaries):
    num = np.sum([calc_num(j, state, ld_boundaries) for j in range(len(ld_boundaries))])
    denum = np.sum([calc_denum(j, state, ld_boundaries) for j in range(len(ld_boundaries))])
    mu = num / (denum + 1e-6)
    var = 1.0 / (denum+1e-6)
    return np.sqrt(var)*stats.norm.rvs() + mu



def sample_sigma2(state, rho, VS=True):
    b = np.zeros(state['num_clusters_'])
    a = np.array(state['suffstats'].values() ) / 2.0 + state['hyperparameters_']['a0k'] 
#     for cluster_id in state['cluster_ids_']:
#         if (cluster_id == 0) and VS is True:
#             continue
#         b[cluster_id] = np.sum(state['beta'][np.where( (state['assignment'] == cluster_id) & 
#                         (state['partition'] == j) )]**2) / 2.0 + state['hyperparameters_']['b0k']
    table = [[] for i in range(state['num_clusters_'])]
    for i in range(len(state['assignment'])):
        table[state['assignment'][i]].append(i)
    
    # pop1 specifc
    for i in range(1, state['population'][1]):
        b[i] = np.sum(state['beta1'][table[i]]**2) / 2.0 + state['hyperparameters_']['b0k']
    
    # pop2 specifc
    for i in range(state['population'][1], state['population'][2]):
        b[i] = np.sum(state['beta2'][table[i]]**2) / 2.0 + state['hyperparameters_']['b0k']
    
    # shared with correlation
#     rho = 0
    for i in range(state['population'][2], state['population'][3]):
        a[i] += state['suffstats'][i] / 2.0
        beta1 = state['beta1'][table[i]]
        beta2 = state['beta2'][table[i]]
        b[i] = np.sum( (beta1**2 + beta2**2 - 2*rho*beta1*beta2) / 2*(1-rho**2) ) + state['hyperparameters_']['b0k']
    
#     beta_assn = [state['beta'][table[i]] for i in range(state['num_clusters_'])]
#     b = np.array([np.sum(beta_assn[i]**2) for i in range(state['num_clusters_'])]) / 2.0 + state['hyperparameters_']['b0k']

    out = np.array([0.0]*state['num_clusters_'])
    if VS is True:
        out[1:] = stats.invgamma(a=a[1:], scale=b[1:]).rvs()
        out[0] = 0
    else: 
        out = dict(zip(range(0, state['num_clusters_']), stats.invgamma(a=a, scale=b).rvs()))
    return out


def update_suffstats(state):
    suff_stats = dict(Counter( state['assignment'] ))
#     for i in range(state['num_clusters_']):
#         if i not in suff_stats.keys():
#             suff_stats[i] = 0
    suff_stats.update(dict.fromkeys(np.setdiff1d(range(state['num_clusters_']), suff_stats.keys()), 0))
    return suff_stats

def sample_V(state):
    for j in range(1,4):
        m = len(state['V'][j])
        suffstats = np.array(state['suffstats'].values()[state['population'][j-1]:state['population'][j]])
        a = 1 + suffstats[:-1]
        b = state['alpha'][j] + np.cumsum(suffstats[::-1])[:-1][::-1]
        sample_val = stats.beta(a=a, b=b).rvs()
        if 1 in sample_val:
            idx = np.argmax(sample_val == 1)
            sample_val[idx+1:] = 0
            sample_return = dict(zip(range(m-1), sample_val))
            sample_return[m-1] = 0
        else:
            sample_return = dict(zip(range(m-1), sample_val))
            sample_return[m-1] = 1
        state['V'][j] = sample_return.values()

def sample_pi_pop(state):
    m = np.array([0.0]*4)
    m[0] += state['suffstats'][0]
    m[1] += np.sum(state['suffstats'].values()[1:state['population'][1]])
    m[2] += np.sum(state['suffstats'].values()[state['population'][1]:state['population'][2]])
    m[3] += np.sum(state['suffstats'].values()[state['population'][2]:state['population'][3]])
    state['suff_pop'] = m
    state['pi_pop'] = dict(zip(range(0, 4), stats.dirichlet(m+1).rvs()[0]))
        
# Compute pi
def update_p(state):
    state['pi'][0] = state['pi_pop'][0]
    for j in range(1,4):
        m = len(state['V'][j])
        V = state['V'][j]
    #     for i in range(1, m-1):
    #         pi[i] = np.prod( 1 - np.array(V[0:i]) ) * state['V'][j][i]
        a = np.cumprod(1-np.array(V)[0:(m-2)])*V[1:(m-1)]
        pi = dict(zip(range(1, m), a))
        pi[0] = state['V'][j][0]
        pi[m-1] = 1 - np.sum(pi.values()[0:(m-1)])

        # last p may be less than 0 due to rounding error
        if pi[m-1] < 0: 
            pi[m-1] = 0
        state['pi_cluster'][j] = pi.values()
        idx = range(state['population'][j-1], state['population'][j])
        state['pi'][idx] = np.array(state['pi_cluster'][j])*state['pi_pop'][j]

# Sample alpha
def sample_alpha(state):
    for j in range(1,4):
        m = np.size(np.where( np.array(state['V'][j]) != 0)); V = state['V'][j]
        a = state['hyperparameters_']['a0'] + m - 1
        b = state['hyperparameters_']['b0'] - np.sum( np.log( 1 - np.array(V[0:(m-1)]) ) )
        state['alpha'][j] = stats.gamma(a=a, scale=1.0/b).rvs()


def gibbs_stick_break(state, rho, ld_boundaries, ref_ld_mat1, ref_ld_mat2, n_threads, VS=True):
    state['cluster_var'] = sample_sigma2(state, rho, VS)

    state['a'] = 0.1; state['c'] = 1
    state['A1'] = [np.linalg.solve(ref_ld_mat1[j]+state['a']*np.identity(ref_ld_mat1[j].shape[0]), ref_ld_mat1[j]) for j in range(len(ld_boundaries))]
    state['B1'] = [np.dot(ref_ld_mat1[j], state['A1'][j]) for j in range(len(ld_boundaries))]
    state['A2'] = [np.linalg.solve(ref_ld_mat2[j]+state['a']*np.identity(ref_ld_mat2[j].shape[0]), ref_ld_mat2[j]) for j in range(len(ld_boundaries))]
    state['B2'] = [np.dot(ref_ld_mat2[j], state['A2'][j]) for j in range(len(ld_boundaries))]
    
    for j in range(len(ld_boundaries)):
        calc_b(j, state, ld_boundaries, ref_ld_mat1, ref_ld_mat2)
        
    state['assignment'] = np.concatenate(Parallel(n_jobs=n_threads, require='sharedmem')(delayed(sample_assignment)(j=j, 
                                ld_boundaries= ld_boundaries, ref_ld_mat1=ref_ld_mat1, 
                                ref_ld_mat2=ref_ld_mat2, state=state, rho=rho, VS=True)
                                for j in range(len(ld_boundaries))))
    state['suffstats'] = update_suffstats(state) 
    sample_pi_pop(state)
    sample_V(state) 
    update_p(state) 
    sample_alpha(state) 

    for j in range(len(ld_boundaries)):
        sample_beta(j, state, ld_boundaries=ld_boundaries, ref_ld_mat1=ref_ld_mat1, 
            ref_ld_mat2=ref_ld_mat2, rho=rho, VS=True)
        compute_varg(j, state, ld_boundaries, ref_ld_mat1, ref_ld_mat2) 
        
    state['h2_1'] = np.sum(state['varg1'])
    state['h2_2'] = np.sum(state['varg2'])
    
    state['eta'] = sample_eta(state, ld_boundaries)














