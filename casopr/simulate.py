import numpy as np
import pyro.distributions as dist
import torch

def simulate_sumstats(ld_blk, blk_size, n_gwas, p, prop_nz = 0.2, beta_sd = 0.1, sigma_noise = 1. ): 
    sigma_over_sqrt_n = sigma_noise / torch.sqrt(torch.tensor(n_gwas))
    nz = torch.rand(p) < prop_nz
    beta_true = torch.where(nz, beta_sd * torch.randn(p), torch.zeros(p))

    annotations = torch.stack([torch.ones(p),nz,torch.randn(p)]).T # intercept, useful annotation, random annotation

    beta_mrg = torch.zeros(p)
    mm = 0
    for kk in range(len(ld_blk)):
        idx_blk = torch.arange(mm,mm+blk_size[kk])
        ld_torch = torch.tensor(ld_blk[kk], dtype = torch.float)
        L, V = torch.linalg.eigh(ld_torch)
        L[L < 0.] = 0.

        beta_mrg[idx_blk] = ld_torch @ beta_true[idx_blk] + sigma_over_sqrt_n * (V @ torch.diag(L.sqrt())) @ torch.randn(blk_size[kk])
        #ld_torch @ beta_true[idx_blk], 
        # covariance_matrix = ld_torch * sigma_over_sqrt_n**2).rsample()
        mm += blk_size[kk]

    return beta_true, beta_mrg, annotations


def stratified_LDSC(annotations, beta_mrg, ld_blk, blk_size):
    [p,k] = annotations.shape
    ldscore = torch.zeros([p,k])
    
    mm = 0
    for kk in range(len(ld_blk)):
        idx_blk = torch.arange(mm,mm+blk_size[kk])
        ld_torch = torch.tensor(ld_blk[kk], dtype = torch.float)
        ldscore[idx_blk,:] = ld_torch**2 @ annotations[idx_blk,:]
        mm += blk_size[kk]
    # min_iterations = 100, max_iterations = 1000, min_particles = 1, max_particles = 32, stall_window = 10, use_renyi = False, lr = 0.03

    chi2 = beta_mrg**2 
    tau = torch.linalg.solve(ldscore.T @ ldscore, ldscore.T @ chi2)

    return tau