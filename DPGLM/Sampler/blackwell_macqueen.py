import numpy as np

def blackwell_macqueen(N, alpha, g0 = lambda: np.random.normal()) :
    """
    Generate samples from a Dirichlet Process using the Blackwell-MacQueen Polya Urn scheme.

    This function sequentially samples from a Dirichlet Process (DP) with concentration parameter `alpha`
    and base distribution `G0`, implemented via the Blackwell-MacQueen Polya urn construction.
    The returned samples correspond to a marginal draw from the DP prior.

    Parameters
    ----------
    N : int
        The number of samples to generate.
    alpha : float
        The concentration parameter of the Dirichlet Process. Higher values lead to more unique clusters.
    g0 : callable, optional
        A function with no arguments that generates i.i.d. samples from the base distribution G0.
        Defaults to `lambda: np.random.normal()`, i.e., standard normal base distribution.

    Returns
    -------
    theta : np.ndarray
        A NumPy array of size `N` containing the sampled values from the DP prior.
    Examples
    --------
    >>> import numpy as np
    >>> np.random.seed(1924)
    >>> blackwell_macqueen(5, alpha=1.0)
    array([ 1.26815852 -1.08636907 -0.43013678  1.26815852 -1.08636907])

    >>> import numpy as np
    >>> np.random.seed(1928)
    >>> blackwell_macqueen(10, alpha=0.5, g0=lambda: np.random.exponential())
    array([1.57322682 1.94654845 0.09618378 1.57322682 1.57322682 1.57322682
           1.94654845 1.94654845 1.94654845 1.57322682])
    """
    theta = np.zeros(N)
    theta[0] = g0() # theta_1 ~ G0
    # then sample iteratively by Polya Urn scheme
    for n in range(1, N):
        # probability of drawing a new base measure sample
        p_new = alpha / (alpha + n - 1)
        if np.random.uniform() < p_new :
            # sample a new value from G0
            theta[n] = g0()
        else:
            # otherwise, resample from theta_{1:(n-1)}
            theta[n] = np.random.choice(theta[:n], size=1)[0]
    return(theta)
