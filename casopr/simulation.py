import os
import sys
import getopt

import parse_genet
import vi

import importlib

import simulate
import torch
import matplotlib.pyplot as plt
import numpy as np

import pyro.distributions as dist

import scipy.stats

chrom = 22

param_dict = {
    'ref_dir' : "/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD_PRScs/ldblk_ukbb_eur/", 
    'bim_prefix' : "test_data/test", 
    'sst_file' : "test_data/sumstats.txt", 
    'n_gwas' : 200000, 
    'out_dir' : "test_data",
    "seed" : 42, 
    "beta_std" : "False", 
    "n_iter" : 1000
}

## 'ref_dir' : "../ld/ldblk_1kg_eur", 

if '1kg' in os.path.basename(param_dict['ref_dir']):
    ref_df = parse_genet.parse_ref(param_dict['ref_dir'] + '/snpinfo_1kg_hm3')
elif 'ukbb' in os.path.basename(param_dict['ref_dir']):
    ref_df = parse_genet.parse_ref(param_dict['ref_dir'] + '/snpinfo_ukbb_hm3')
ref_df = ref_df[ref_df.CHR == chrom]

vld_df = parse_genet.parse_bim(param_dict['bim_prefix'] + ".bim")
vld_df = vld_df[vld_df.CHR == chrom]

sst_dict = parse_genet.parse_sumstats(ref_df, vld_df, param_dict['sst_file'], param_dict['n_gwas'])

ld_blk, ld_blk_sym, blk_size = parse_genet.parse_ldblk(param_dict['ref_dir'], sst_dict, chrom)

#mcmc_gtb.mcmc(param_dict['a'], param_dict['b'], param_dict['phi'], sst_dict, param_dict['n_gwas'], ld_blk, blk_size, param_dict['n_iter'], param_dict['n_burnin'], param_dict['thin'], int(chrom), param_dict['out_dir'], param_dict['beta_std'], param_dict['seed'])

beta_true, beta_mrg, annotations = simulate.simulate_sumstats(ld_blk, blk_size, param_dict['n_gwas'], p = len(sst_dict))

sst_dict["BETA"] = beta_mrg

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

importlib.reload(vi)
one = torch.tensor(1., device=device)
losses, beta, phi_est, stats = vi.vi(sst_dict, param_dict['n_gwas'], ld_blk, blk_size, device = device, annotations = annotations, max_iterations = param_dict['n_iter'], max_particles=4, desired_min_eig = 1e-3, min_iterations = 200, stall_window = 30, phi_as_prior = False, lr = 0.03, constrain_sigma = True)

# observations: 
# collapsed >> uncollapsed (deleted uncollapsed code)
# phi_as_prior = False seems to give better results, but phi_as_prior = True gives more robust training (less noisy loss), probably because constraint on psi is more useful. 
# Does constrain_sigma help? I think constrain_psi does. 
# Seem to need desired_min_eig=1e-3 at least to avoid numerical issues. Could go lower if using float64 it seemed. 

plt.plot(losses); plt.show()

plt.scatter(beta_true, beta)
plt.xlabel("True beta")
plt.ylabel("Infered beta")
plt.title("Pearson R=%1.3f" % scipy.stats.pearsonr(beta_true,beta)[0] )

stats["annotation_weights"] # ideally should show large positive value for second entry, and close to 0 for third (noise)

# convert standardized beta to per-allele beta
if param_dict["beta_std"] == 'False':
    beta /= np.sqrt(2.0*sst_dict['MAF']*(1.0-sst_dict['MAF']))
    
sst_dict["beta_shrunk"] = beta
sst_dict.to_csv("beta_shrunk.tsv", sep = "\t")

