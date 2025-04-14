## Generate Tables 2/3 from results
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

def check_coverage(est, var, n, z, true):
    """ For a given estimator and variance estimator, returns
        if true value covered.
        @param z is the critical value for two sided interval
    """
    lb = est - z*np.sqrt(var)/np.sqrt(n)
    ub = est + z*np.sqrt(var)/np.sqrt(n)
    
    if (true >= lb) & (true <= ub):
        return 1
    else:
        return 0

# Args
pk = float(sys.argv[1])
qk = float(sys.argv[2])
design= int(sys.argv[3])
B = int(sys.argv[4])
ate = float(sys.argv[5])
v_k = float(sys.argv[6])
tilde_v_k = float(sys.argv[7])

nm = "results_" + str(pk) + "_" + str(qk) + "_" + str(B) + "_prepped_df_design" + str(design) + ".csv.txt"
dat = np.loadtxt(nm)
dat = dat.reshape(dat.shape[0]/11, 11)

# column names
nms = ["ols", "fe", "Nk",
       "cluster_V_OLS", "r_V_OLS",
       "r_V_FE", "cluster_V_FE","ccv_V_FE",
       "ccv_V_OLS", "tscb_std_OLS", "tscb_std_FE"]
res = pd.DataFrame(dat, columns=nms)

## At the moment results only work for qk =1.

# Generate Table 2
table2 = pd.DataFrame()

table2['Design'] = [design,design]
Nk = res['Nk'].unique()[0]
table2['Nk_sd'] = [res['ols'].std()*np.sqrt(Nk),
                  res['fe'].std()*np.sqrt(Nk)]

table2['v_k'] = [v_k, tilde_v_k]

table2['robust'] = [np.sqrt(res['r_V_OLS']).mean(),
                    np.sqrt(res['r_V_FE']).mean()]

table2['cluster'] = [np.sqrt(res['cluster_V_OLS']).mean(),
                     np.sqrt(res['cluster_V_FE']).mean()]

table2['CCV'] = [np.sqrt(res['ccv_V_OLS']).mean(),
                 np.sqrt(res['ccv_V_FE']).mean()]

table2['TSCB'] = np.mean(res[['tscb_std_OLS', 'tscb_std_FE']]).values*np.sqrt(Nk)

tab_nm = "table2_design_" + str(design) + ".csv"
table2.to_csv(tab_nm)


## Generate Table 3
res[['tscb_var_OLS','tscb_var_FE']] = (res[['tscb_std_OLS', 'tscb_std_FE']]*np.sqrt(Nk))**2

# coverage
cover_ols_r = 0.
cover_fe_r = 0.
cover_ols_cluster = 0.
cover_fe_cluster = 0.
cover_ols_ccv = 0.
cover_fe_ccv = 0.
cover_ols_tscb = 0.
cover_fe_tscb = 0.
cover_vk = 0.
cover_tilde_vk = 0.

for i,row in res.iterrows():
    cover_ols_r = cover_ols_r + check_coverage(row['ols'],
                                               row['r_V_OLS'],
                                               Nk,
                                               1.96,
                                               ate)
    
    cover_fe_r = cover_fe_r + check_coverage(row['fe'],
                                               row['r_V_FE'],
                                               Nk,
                                               1.96,
                                               ate)
    
    cover_ols_cluster = cover_ols_cluster + check_coverage(row['ols'],
                                               row['cluster_V_OLS'],
                                               Nk,
                                               1.96,
                                               ate)
    
    cover_fe_cluster = cover_fe_cluster + check_coverage(row['fe'],
                                               row['cluster_V_FE'],
                                               Nk,
                                               1.96,
                                               ate)
    
    cover_ols_ccv = cover_ols_ccv + check_coverage(row['ols'],
                                               row['ccv_V_OLS'],
                                               Nk,
                                               1.96,
                                               ate)
    
    cover_fe_ccv = cover_fe_ccv + check_coverage(row['fe'],
                                               row['ccv_V_FE'],
                                               Nk,
                                               1.96,
                                               ate)  
    
    cover_vk = cover_vk + check_coverage(row['ols'],
                                         v_k**2,
                                         Nk,
                                         1.96,
                                         ate)  
    
    cover_tilde_vk = cover_tilde_vk + check_coverage(row['fe'],
                                                     tilde_v_k**2,
                                                     Nk,
                                                     1.96,
                                                     ate)  
# for all except Design 1
res1 = res
for i,row in res1.iterrows():
    cover_ols_tscb = cover_ols_tscb + check_coverage(row['ols'],
                                               row['tscb_var_OLS'],
                                               Nk,
                                               1.96,
                                               ate)
    
    cover_fe_tscb = cover_fe_tscb + check_coverage(row['fe'],
                                               row['tscb_var_FE'],
                                               Nk,
                                               1.96,
                                               ate)  

table3 = pd.DataFrame()

table3['Design'] = [design,design]
n = res.shape[0]

table3['v_k'] = [np.float(cover_vk)/n,
                 np.float(cover_tilde_vk)/n]

table3['robust'] = [np.float(cover_ols_r)/n,
                    np.float(cover_fe_r)/n]

table3['cluster'] = [np.float(cover_ols_cluster)/n,
                     np.float(cover_fe_cluster)/n]

table3['CCV'] = [np.float(cover_ols_ccv)/n,
                 np.float(cover_fe_ccv)/n]

table3['TSCB'] = [np.float(cover_ols_tscb)/res1.shape[0],
                  np.float(cover_fe_tscb)/res1.shape[0]]

nm3 = "table3_design_" + str(design) + ".csv"
table3.to_csv(nm3)
