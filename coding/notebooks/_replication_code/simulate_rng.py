## Script to run simulations for cluster paper
#
#
# Input requirements:
#
#   - "population" data frame with potential outcomes calculated and
#       columns as described in the one_simulation() function
#   - calculated values of mu_l and sigma_l from the same population as the one used
#       to generaet the dataframe.
#

import sys
import logging
import random
from datetime import datetime
from os.path import join
from time import time

import time
import random
import warnings
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import ShuffleSplit
from sklearn.utils import shuffle

# Setting random seed 
#random.seed(int(sys.argv[7]))
seeds = np.loadtxt("random_seeds.txt")
seed = int(seeds[int(sys.argv[7])])
rng = np.random.RandomState(seed)

# Logging and time check.
start_time = time.time()
logging.basicConfig(level=logging.INFO)

# Args
pk = float(sys.argv[1])
qk = float(sys.argv[2])
mu = float(sys.argv[3]) # these could be computed here but for ease of computation we pass them
std = float(sys.argv[4])
B = int(sys.argv[5])
nm = sys.argv[6]

# By default, creates a file called 'result.txt'.
# or appends to one already existing.
filename = "results_" +str(pk) + "_" + str(qk) + "_" + str(B) + "_" + nm + ".txt"

logging.info("Received parameters: " + str(sys.argv[1:]))

# Functions
## Functions for cluster replication

def FE_estimator(df, nmy, nmx, nm_cluster):
    # first we get tau fixed and state averages
    state_means = df.groupby(nm_cluster, as_index=False).agg({nmx: 'mean',
                                                  nmy: 'mean'})
    state_means['some_college_mean'] = state_means[nmx]
    state_means['learn_mean'] = state_means[nmy]
    
    lean_df2 = df.merge(state_means[[nm_cluster, 'some_college_mean', 'learn_mean']], 
                 how = 'left', on = nm_cluster)
    
    lean_df2['diff_some_college'] = lean_df2[nmx] - lean_df2['some_college_mean']
    lean_df2['diff_learn'] = lean_df2[nmy] - lean_df2['learn_mean']
    
    tau_fixed = np.sum(lean_df2[nmy]*lean_df2['diff_some_college']) / np.sum(lean_df2[nmx]*lean_df2['diff_some_college'])
    
    return tau_fixed

def compute_CCV(df, u, nmx, pk):
    """ @return: CCV variance for the case qk==1
        @param df dataframe with columns
            - 'u' indicator for random split (Z in the paper)
            - nmx name of the binary indicator W
            - 'state' name of the cluster indicator
            - 'Y' outcome variable
        @param u binary variable indicating which split to use for computing the model and which for evaluation
        @param nmx string with the name of the treatment indicator W
        @param pk double in (0,1] as described in the paper
    """
    
    u1 = df['u']==u
    u0 = df['u']==1-u
    w1 = df[nmx]==1
    w0 = df[nmx]==0
    
    # compute alpha, tau using first split
    alpha = df[u1 & w0]['Y'].mean()
    tau = df[u1 & w1]['Y'].mean() - df[u1 & w0]['Y'].mean()
    tau_full = df[w1]['Y'].mean() - df[w0]['Y'].mean()
        
    # compute for each m
    tau_ms = {}
    nm = 'tau_'
    pk_term = 0
    for m in df['state'].unique():
        ind_m = df['state']==m
        aux1 = df[u1 & ind_m & w1]['Y'].mean()
        aux0 =  df[u1 & ind_m & w0]['Y'].mean()
        
        aux1_full = df[ind_m & w1]['Y'].mean()
        aux0_full =  df[ind_m & w0]['Y'].mean()
        aux_tau = aux1 - aux0
        aux_tau_full = aux1_full - aux0_full
        
        if (np.isnan(aux1)) or (np.isnan(aux0)):
            aux_tau = tau
        
        aux_nm = nm + str(m)
        tau_ms[aux_nm] = aux_tau
        
        # compute the pk term in u0
        Nm = df[ind_m].shape[0] 
        aux_pk = Nm*((aux_tau_full-tau)**2)
        pk_term = pk_term + aux_pk
    
    # compute the residuals
    df['resU'] = df['Y'] - alpha - df[nmx]*tau
    
    # Wbar
    Wbar = df[u1][nmx].mean()
    #Wbar = df[nmx].mean() # to match Guido
    
    # pk term
    pk_term = pk_term*(1.-pk)/df.shape[0]
    
    # compute avg Z
    Zavg = np.float(np.sum(df['u']==1-u))/df.shape[0]

    # compute the normalized CCV using second split
    n = (df.shape[0]*(Wbar**2)*((1.-Wbar)**2))
    sum_CCV = 0
    for m in df['state'].unique():
        ind_m = df['state']==m
        df_m = df[u0 & ind_m]
        aux_nm = nm + str(m)

        # tau term
        tau_term = (tau_ms[aux_nm] - tau)*Wbar*(1.-Wbar)

        # Residual term
        res_term = (df_m[nmx] - Wbar)*df_m['resU']

        # square of sums
        sq_sum = np.sum(res_term - tau_term)**2
        
        # sum of squares
        sum_sq = np.sum((res_term - tau_term)**2)
        
        # compute CCV
        sum_CCV = sum_CCV + (1./(Zavg**2))*sq_sum - ((1.-Zavg)/(Zavg**2))*sum_sq + n*pk_term
    
    # normalize
    V_CCV = sum_CCV / n
    
    return V_CCV

# Sample according to 2 step procedure
def get_sample(df, state_df, nmx, random_seed):
    w1 = df[nmx]==1
    w0 = df[nmx]==0
    
    sample_df = pd.DataFrame()

    for m in df['state'].unique():
        ind_m = df['state']==m
        df_m = df[ind_m]

        # draw a Wbar
        Wbar = rng.choice(state_df[nmx])

        # draw treated units Nm Wbar
        num_T = int(np.rint(state_df.loc[state_df['state']==m,'Y']*Wbar))
        ind_T = df_m[w1].sample(num_T, replace=True, axis=0, random_state=random_seed) #.index.to_list()

        # draw untreated units Nm (1- Wbar)
        num_U = int(np.rint(state_df.loc[state_df['state']==m,'Y']*(1.- Wbar)))
        ind_U = df_m[w0].sample(num_U, replace=True, axis=0, random_state=random_seed) #.index.to_list()
        
        aux_sample = pd.concat([ind_T, ind_U])
        sample_df = pd.concat([sample_df, aux_sample], axis=0)
        
    return sample_df
    
# Bootstrap
def boostrap_taus(B, df, state_df, nmx, model):
    
    taus_OLS = []
    taus_FE = []
    
    for b in range(B):
        random_seed = b + seed
        sample = get_sample(df, state_df, nmx, random_seed)
        
        # OLS
        y = sample['Y']
        X = np.array(sample[nmx]).reshape(sample.shape[0],1)
        
        model.fit(X,y)
        
        # FE
        fe = FE_estimator(sample, 'Y', nmx, 'state')

        taus_OLS.append(model.coef_[0])
        taus_FE.append(fe)
        
        #if b%100 == 0:
            #pd.DataFrame(taus).to_csv('bootstrap_taus.csv')
            #print(b)
        
    return (taus_OLS, taus_FE)


# sample according to sampling process
# First sample clusters and units
# Then assign treatment

def sample_pq(pk, qk, df):
    """ returns indices of the generated sample for the original df
    """
    
    if (pk==1) & (qk==1):
        return df.index
    
    # sample clusters and units
    aux = pd.DataFrame()
    idx_sample = aux.index
    for m in df['state'].unique():
        sample_m = rng.binomial(1, qk)
        if sample_m == 1:
            ind_m = df['state']==m
            df_m = df[df['state']==m]
            units = rng.binomial(1, pk, size=df_m.shape[0])
            mask = np.multiply(list(df_m.index), units)
            mask = mask[mask!=0]
            idx_sample = idx_sample.union(pd.Index(mask))
            
    return idx_sample

def assign(mu, s, smpl_df):
    """ returns assignment vector for smpl_df
    """
    Wkis = []
    for m in smpl_df['state'].unique():
        odds = rng.normal(mu,s)
        Akm = np.exp(odds) / (1.+ np.exp(odds))

        ind_m = smpl_df['state']==m
        df_m = smpl_df[smpl_df['state']==m]
                
        Wki_m = rng.binomial(1, Akm, size=df_m.shape[0])
        Wkis = np.concatenate([Wkis, Wki_m])
        
    return Wkis


def one_sample(pk, qk, mu, s, df):
    """ returns one sample with W, Y and state from df
    """

    idx_sample = sample_pq(pk,qk,df)
    smpl_df = df.iloc[idx_sample]
    
    new_df = pd.DataFrame()
    
    new_df['state'] = smpl_df['state']
    new_df['W'] = assign(mu, s, smpl_df)
    new_df['Y'] = smpl_df['y1']*new_df['W'] + smpl_df['y0']*(1. - new_df['W'])
    new_df = new_df.reset_index()
    
    return new_df

def one_simulation(df, pk, qk, mu, s, B):
    # simulate sample
    smpl_df = one_sample(pk, qk, mu, s, df)

    # simulate FE
    fe = FE_estimator(smpl_df, 'Y', 'W', 'state')

    # simulate OLS
    y = smpl_df['Y']
    X = np.array(smpl_df['W']).reshape(smpl_df.shape[0],1)

    model = LinearRegression()
    model.fit(X,y)
    ols = model.coef_[0]

    Wbar = smpl_df['W'].mean()
    smpl_df['Uhat'] = smpl_df['Y'] - model.predict(smpl_df['W'].values.reshape(-1,1))
    
    # robust SE OLS
    wbar_factor = ((Wbar**2)*(1.-Wbar)**2)
    r_V = np.mean(smpl_df['Uhat']**2 * (smpl_df['W'] - Wbar)**2)/wbar_factor

    # cluster SE OLS
    cluster_SE = 0
    for m in smpl_df['state'].unique():
        dfm = smpl_df[smpl_df['state'] == m]
        err = np.sum(dfm['Uhat']*(dfm['W'] - Wbar))**2
        cluster_SE = cluster_SE + err

    cluster_SE = cluster_SE/(smpl_df.shape[0]*wbar_factor)

    # robust and cluster SE FE
    sum_tildeU = 0
    sum_tildeW = 0

    sum_tildeU_FE = 0

    num_lambdak = 0
    den_lambdak = 0

    for m in smpl_df['state'].unique():
        dfm = smpl_df[smpl_df['state'] == m]
        Ym = dfm['Y'].mean()
        Wmbar = dfm['W'].mean()
        dfm['Wtilde'] = dfm['W'] - Wmbar
        dfm['Utilde'] = (dfm['Y'] - Ym) - fe*dfm['Wtilde']

        sum_tildeU = sum_tildeU + np.sum((dfm['Wtilde']**2)*(dfm['Utilde']**2))
        sum_tildeU_FE = sum_tildeU_FE + np.sum((dfm['Wtilde']*dfm['Utilde']))**2
        sum_tildeW = sum_tildeW + np.sum((dfm['Wtilde']**2))

        num_lambdak = num_lambdak + Wmbar*(1. - Wmbar)
        den_lambdak = den_lambdak + (Wmbar**2)*((1. - Wmbar)**2)

    rFE_V = smpl_df.shape[0]*(sum_tildeU / (sum_tildeW **2))
    clusterFE_V = smpl_df.shape[0]*(sum_tildeU_FE / (sum_tildeW **2))

    Mk = smpl_df['state'].nunique()
    lambdak = 1. - qk*((num_lambdak/Mk)**2/(den_lambdak/Mk))

    # CCV FE
    CCV_FE_V = lambdak*clusterFE_V + (1. - lambdak)*rFE_V

    # CCV OLS
    # generate splits
    #sss = ShuffleSplit(n_splits=2, test_size=0.5)
    #sss.get_n_splits(smpl_df)

    #smpl_df['u'] = np.where(smpl_df.index.isin(next(sss.split(smpl_df))[0]),1,0)
    #V_CCV1 = compute_CCV(smpl_df, 1, 'W', pk)
    #V_CCV2 = compute_CCV(smpl_df, 0, 'W', pk)
    #se_avg = 0.5*(V_CCV1+ V_CCV2)

    # CCV new method
    sum_CCV = 0
    for i in range(4):
        random_seed = i + seed
        sss = ShuffleSplit(n_splits=2, test_size=0.5, random_state=random_seed)
        sss.get_n_splits(smpl_df)
        
        smpl_df['u'] = np.where(smpl_df.index.isin(next(sss.split(smpl_df))[0]),1,0)
        V_CCV = compute_CCV(smpl_df, 1, 'W', pk)
        sum_CCV = sum_CCV + V_CCV
    se_avg = 0.25*sum_CCV

    # bootstrap
    grouped_df = smpl_df.groupby('state', as_index=False).agg({'W': 'mean',
                                                               'Y': 'count'})
    tausOLS, tausFE = boostrap_taus(B, smpl_df, grouped_df,
                                    'W',  LinearRegression())
    
    return [ols, fe, smpl_df.shape[0], 
            cluster_SE, r_V, 
            rFE_V, clusterFE_V,
            CCV_FE_V, se_avg, 
            np.std(tausOLS), np.std(tausFE)]

def sim(i):
	return one_simulation(df, pk, qk, mu, std, B)

logging.info("Loading data.")

## Load data
df = pd.read_csv(nm)

# ignore warnings
warnings.filterwarnings('ignore')

# simulate
logging.info("Simulating...")
# ols, fe, Nk, cluster_V_OLS, r_V_OLS, r_V_FE, cluster_V_FE,ccv_V_FE, ccv_V_OLS, tscb_std_OLS, tscb_std_FE
results = one_simulation(df, pk, qk, mu, std, B)

# Append results to existing file, or create one if none exists.
logging.info("Done. Storing.")
with open(filename, "a") as f:
    for res in results:
        print >> f, res

# More logging as appropriate.
end_time = time.time()
logging.info("Finished simulation after seconds:")
logging.info(end_time - start_time)
