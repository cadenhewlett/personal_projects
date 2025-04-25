import numpy as np
import pickle
import matplotlib.pyplot as plt
from sampler import *
from scipy.stats import laplace

np.random.seed(1924)
# 1. Simulate:
true_centers = [-5.0, 0.0, 5.0]
y_data = np.hstack([
  np.random.uniform(c-0.3, c+0.3, size=30)
  for c in true_centers
])

# 2. Base-sampler & log-prior:
draw_G0 = lambda: np.random.uniform(-10,10)
log_prior = lambda theta: 0 if -10 <= theta <= 10 else -999
# 3. Uniform likelihood:
def log_likelihood(y, theta):
    return ( -np.log(0.6)
             if (theta-0.3) <= y <= (theta+0.3)
             else -999 )

# 4. Run dpmm_aux_gibbs with, e.g.:
traces = dpmm_aux_gibbs(
  y=y_data,
  alpha=1.0,
  draw_G0=draw_G0,
  log_prior=log_prior,
  log_likelihood=log_likelihood,
  m=3,
  num_iters=2000,
  proposal_std=0.2,
  REPORT=True
)

c_final      = traces['c'][-1]
theta_final  = traces['theta'][-1]
unique = np.unique(c_final)
for lbl in unique:
    idx = (c_final == lbl)
    print(f"Cluster {lbl}: size={idx.sum()}, mean(y)={y_data[idx].mean():.3f}, theta_est={theta_final[lbl]:.3f}")

# cluster_counts = [ len(np.unique(c_arr)) for c_arr in traces['c'] ]
# Save output
output = {"traces": traces, "y_data": y_data}
with open('aux_gibbs_output.pkl', 'wb') as f:
    pickle.dump(output, f)