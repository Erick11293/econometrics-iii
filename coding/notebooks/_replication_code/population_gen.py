## generate population designs
#

import sys
import logging
import random
from datetime import datetime
from os.path import join
from time import time
from scipy.io import loadmat 

import warnings
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import ShuffleSplit
from sklearn.utils import shuffle

## Functions for asymptotic variances
#


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

def vk(df, mu, s, p, q):
    """ returns the asymptotic variance
        df must have cols u1, u0 and state
    """
    u1 = df['u1'].values
    u0 = df['u0'].values
    
    t1 = np.mean((1./mu)*(u1 **2) + (1./(1.-mu))*(u0 **2))
    t2 = -p*np.mean((u1 - u0)**2) - p*(s**2)*np.mean(((1./mu)*u1 + (1./(1.-mu))*u0)**2)
    
    df['diff'] = df['u1'] - df['u0']
    t3 = p*(1.-q)*np.sum(df.groupby('state')['diff'].sum() ** 2)/df.shape[0]
    
    df['sum'] = (1./mu)*df['u1'] + (1./(1.-mu))*df['u0']
    t4 = p*(s**2)*np.sum(df.groupby('state')['sum'].sum()**2)/df.shape[0]
    
    return np.sum([t1, t2, t3, t4])

def vk_cluster(df, mu, s, p, q):
    """ returns the asymptotic variance
        df must have cols u1, u0 and state
    """
    u1 = df['u1'].values
    u0 = df['u0'].values
    
    t1 = np.mean((1./mu)*(u1 **2) + (1./(1.-mu))*(u0 **2))
    t2 = -p*np.mean((u1 - u0)**2) - p*(s**2)*np.mean(((1./mu)*u1 + (1./(1.-mu))*u0)**2)
    
    df['diff'] = df['u1'] - df['u0']
    t3 = p*np.sum(df.groupby('state')['diff'].sum() ** 2)/df.shape[0]
    
    df['sum'] = (1./mu)*df['u1'] + (1./(1.-mu))*df['u0']
    t4 = p*(s**2)*np.sum(df.groupby('state')['sum'].sum()**2)/df.shape[0]
    
    return np.sum([t1, t2, t3, t4])

def tilde_vk(df, p, q, tauk, nmx, mu, s, Eaa, Ea2a, Eaa2, Ea2a2):
    """ df needs to have a nmx column, state and e1, e0
    """
    
    # get Akm population means
    grouped_df = lean_df2.groupby('state').agg({
            nmx: 'mean',
            'tau_m': 'mean',
            'e0': 'count' # use as count variable
        })
    #Eaa2 = np.mean(grouped_df[nmx]*(1. - grouped_df[nmx])**2)
    #Ea2a = np.mean((grouped_df[nmx]**2)*(1. - grouped_df[nmx]))
    #Ea2a2 = np.mean((grouped_df[nmx]**2)*(1. - grouped_df[nmx])**2)
    #Eaa = np.mean(grouped_df[nmx]*(1. - grouped_df[nmx]))
    
    t1 = Eaa2*np.mean(df['e1']**2) + Ea2a*np.mean(df['e0']**2)
    t2 = -p*Ea2a2*np.mean((df['e1'] - df['e0'])**2)
    
    aux = ((grouped_df['tau_m'] - tauk)**2)/df.shape[0]
    
    t3 = (Eaa - (5.+p)*Ea2a2 + 2.*q*(Eaa)**2)*np.sum((grouped_df['e0'])*aux)
    t4 = (p*Ea2a2 -p*q*(Eaa)**2)*np.sum((grouped_df['e0']**2)*aux)
    
    t5 = (mu*(1.-mu) - s**2)**2
    
    return np.sum([t1, t2, t3, t4]/t5)

def get_Akm_moments(mu_l, s_l, num):
    """ returns the mean and std of Akm over 
        num draws
    """
    akms = []
    for i in range(num):
        odds = np.random.normal(mu_l,s_l)
        Akm = np.exp(odds) / (1.+ np.exp(odds))
        akms.append(Akm)
    
    akms = np.array(akms)
    
    Eaa2 = np.mean(akms*(1. - akms)**2)
    Ea2a = np.mean((akms**2)*(1. - akms))
    Ea2a2 = np.mean((akms**2)*(1. - akms)**2)
    Eaa = np.mean(akms*(1. - akms))
    
    return (np.mean(akms), np.std(akms),
           Eaa, Ea2a, Eaa2, Ea2a2)

# set the seed
np.random.seed(2022)

# Args
pk = float(sys.argv[1])
qk = float(sys.argv[2])
design= int(sys.argv[3])
cx = float(sys.argv[4])
calpha = float(sys.argv[5])
ctau = float(sys.argv[6])

## Load data
mat = loadmat('census.mat') 
df = pd.DataFrame()
for m in mat.keys()[3:]:
    df[m] = mat[m].flatten()


# Prepare data to run functions
df['Y'] = df['learn']
df['W'] = df['educ']>12
lean_df = df[['Y', 'W', 'state']]

## Generate potential outcomes

# first we get tau fixed and state averages
state_means = df.groupby('state', as_index=False).agg({'W': 'mean',
                                              'Y': 'mean'})
state_means['W_mean'] = state_means['W']
state_means['Y_mean'] = state_means['Y']

# get the mu, s for the population
# log odds ratios
state_means['lm_W'] = np.log(state_means['W_mean']/(1.-state_means['W_mean']))
mu = state_means['lm_W'].mean()
stdev = np.std(state_means['lm_W'])*cx


lean_df2 = lean_df.merge(state_means[['state', 'W_mean', 'Y_mean']], 
                         how = 'left',
                         on = 'state')

lean_df2['diff_W'] = lean_df2['W'] - lean_df2['W_mean']
lean_df2['diff_Y'] = lean_df2['Y'] - lean_df2['Y_mean']

# get residuals OLD Version
tau_fixed = FE_estimator(df, 'Y', 'W', 'state')
lean_df2['resU'] = lean_df2['diff_Y'] - tau_fixed*lean_df2['diff_W']

# compute the residuals and
# get tau_m and alpha_m
# compute for each m
y0_ms = {}
y1_ms = {}
tau_kms = {}

w1 = (df['W']==1)
w0 = (df['W']==0)

for m in lean_df2['state'].unique():
    ind_m = (lean_df2['state']==m)
    aux1 = lean_df2[ind_m & w1]['Y'].mean()
    aux0 = lean_df2[ind_m & w0]['Y'].mean()
    
    y = lean_df2[ind_m]['Y']
    X = np.array(lean_df2[ind_m]['W']).reshape(lean_df2[ind_m].shape[0],1)
    
    # do a linear reg and get tau_kms
    model = LinearRegression()
    model.fit(X,y)
    tau_km = model.coef_[0]

    tau_kms[m] = tau_km
    y0_ms[m] = aux0
    y1_ms[m] = aux1

# compute the avg taus and alphas
tau_mean = np.array(list(tau_kms.values())).mean()
alpha_mean = np.array(list(y0_ms.values())).mean()

lean_df2['resU'] = lean_df2['diff_Y'] - lean_df2['state'].map(tau_kms)*lean_df2['diff_W']

# correct values according to factor
for m in lean_df2['state'].unique():
    tau_kms[m] = (tau_kms[m] - tau_mean)*ctau + tau_mean
    y0_ms[m] = (y0_ms[m] - alpha_mean)*calpha + alpha_mean

lean_df2['y0'] = lean_df2['state'].map(y0_ms) + lean_df2['resU']
lean_df2['y1'] = lean_df2['state'].map(y0_ms) + lean_df2['state'].map(tau_kms) + lean_df2['resU']

# save the prepped simulation dataframe
name = 'prepped_df_design' + str(design) + '.csv'
lean_df2.to_csv(name)

# compute ate for coverage
ate = np.mean(lean_df2['y1'] - lean_df2['y0'])

# Akm moments
Akm_mu, Akm_std, Eaa, Ea2a, Eaa2, Ea2a2 = get_Akm_moments(mu, stdev, 10000000)

# Asymptotic variances

# compute true residuals
alpha_k = lean_df2['y0'].mean()
tau_k = np.mean(lean_df2['y1'] - lean_df2['y0'])
lean_df2['u0'] = lean_df2['y0'] - alpha_k
lean_df2['u1'] = lean_df2['y1'] - (alpha_k + tau_k)

# compute asymptotic variances
vk_std = np.sqrt(vk(lean_df2, Akm_mu, Akm_std, pk, qk))

lean_df2['e0'] = lean_df2['y0'] - lean_df2['state'].map(y0_ms) 
lean_df2['e1'] = lean_df2['y1'] - lean_df2['state'].map(y0_ms) - lean_df2['state'].map(tau_kms)

# get tauk and tau_m, alpha_m
lean_df2['tau_m'] = lean_df2['state'].map(tau_kms)
lean_df2['alpha_m'] = lean_df2['state'].map(y0_ms) 

tilde_vk_std =np.sqrt(tilde_vk(lean_df2, pk, qk, tau_k, 'W', Akm_mu, Akm_std, Eaa, Ea2a, Eaa2, Ea2a2))

# store outputs
out = np.array([ate, mu, stdev, vk_std, tilde_vk_std])
name_out = 'output_design' + str(design) + '.csv'
np.savetxt(name_out, out, delimiter=',')



