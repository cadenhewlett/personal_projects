import numpy as np

# Neal, R. M. (2000) JCGS, vol.9(2), pp.249–265
# Section 6: Gibbs Sampling with Auxiliary Parameters
# Algorithm 8 (p.261):
# Let the state of the Markov chain consist of c = (c1,...,cn) and θ = (θc : c ∈ {c1,...,cn})."

def dpmm_aux_gibbs(y, alpha, draw_G0, log_prior, log_likelihood,
                   m, num_iters, proposal_std,
                   burn_in=None, REPORT=True):
    """
    Dirichlet Process Mixture Model sampler using Neal's Algorithm 8.

    Parameters
    ----------
    y : array-like, shape (n,)
        Observations.
    alpha : float
        Concentration parameter of the Dirichlet Process.
    draw_G0 : callable
        Function returning an independent sample from the base distribution G0.
    log_prior : callable
        Function computing log density of G0 at theta. Signature: log_prior(theta).
    log_likelihood : callable
        Function computing log F(y_i | theta). Signature: log_likelihood(y_i, theta).
    m : int
        Number of auxiliary draws per update of each c_i.
    num_iters : int
        Number of full Gibbs sweeps to perform.
    proposal_std : float
        Standard deviation for Gaussian random-walk proposal in MH updates of theta_c.
    burn_in : int or None, optional
        Number of iterations to treat as burn-in. Defaults to one-third of num_iters.
    REPORT : bool, default True
        Whether to display a single-line updating progress.

    Returns
    -------
    traces : dict
        Contains:
          - 'c': list of arrays of cluster assignments at each saved iteration
          - 'theta': list of dicts mapping cluster label -> theta at each saved iteration
    """
    n = len(y)
    if burn_in is None:
        burn_in = num_iters // 3

    # 0. Initialization: each i in its own cluster
    c = np.arange(n, dtype=int)
    theta = {i: draw_G0() for i in range(n)}
    traces = {'c': [], 'theta': []}

    for it in range(num_iters):
        # Progress reporting on single console line
        if REPORT:
            pct = int(100 * (it + 1) / num_iters)
            if it + 1 == burn_in:
                print("\nBurn-in done!")
            print(f"\rRunning sampler: {pct}%", end='', flush=True)

        # === Step 1: update c[i] with auxiliary parameters ===
        for i in range(n):
            mask = np.arange(n) != i
            ex_labels, ex_counts = np.unique(c[mask], return_counts=True)
            count_dict = dict(zip(ex_labels, ex_counts))
            k_minus = len(ex_labels)
            current = c[i]
            is_singleton = current not in count_dict

            # Draw auxiliaries
            theta_aux = []
            if is_singleton:
                theta_aux.append(theta[current])
                for _ in range(m-1):
                    theta_aux.append(draw_G0())
            else:
                for _ in range(m):
                    theta_aux.append(draw_G0())

            # Compute weights following Neal, p.261
            w_ex = [count_dict[lbl] * np.exp(log_likelihood(y[i], theta[lbl]))
                    for lbl in ex_labels]
            w_aux = [(alpha/m) * np.exp(log_likelihood(y[i], th))
                     for th in theta_aux]
            weights = np.array(w_ex + w_aux, dtype=float)
            total = weights.sum()
            if total > 0:
                weights /= total
            else:
                weights = np.ones_like(weights) / len(weights)

            # Sample new label
            choice = np.random.choice(len(weights), p=weights)
            if choice < k_minus:
                c[i] = ex_labels[choice]
            else:
                aux_idx = choice - k_minus
                if is_singleton and aux_idx == 0:
                    c[i] = current
                else:
                    new_lbl = max(theta.keys(), default=-1) + 1
                    theta[new_lbl] = theta_aux[aux_idx]
                    c[i] = new_lbl

        # === Step 2: prune empty clusters & relabel ===
        active = np.unique(c)
        label_map = {old: new for new, old in enumerate(active)}
        c = np.array([label_map[ci] for ci in c], dtype=int)
        theta = {label_map[old]: theta[old] for old in active}

        # === Step 3: update theta_c via MH  ===
        for lbl in list(theta.keys()):
            y_clust = y[c == lbl]
            old = theta[lbl]
            prop = old + np.random.normal(scale=proposal_std)
            logp_old = log_prior(old) + np.sum([log_likelihood(yi, old) for yi in y_clust])
            logp_prop = log_prior(prop) + np.sum([log_likelihood(yi, prop) for yi in y_clust])
            if np.log(np.random.rand()) < (logp_prop - logp_old):
                theta[lbl] = prop

        # Save trace after burn-in
        if it >= burn_in:
            traces['c'].append(c.copy())
            traces['theta'].append(theta.copy())

    # Final newline after progress
    if REPORT:
        print()
    return traces
