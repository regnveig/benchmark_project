import sys
from typing import Optional, Iterable, List, Tuple
import hicreppy.utils.mat_process as cu
import numpy as np
import scipy.stats as ss
from scipy.sparse import SparseEfficiencyWarning
import warnings

def get_scc(
    mat1: 'scipy.sparse.csr_matrix', mat2: 'scipy.sparse.csr_matrix', max_bins: int
) -> float:
    """
    Compute the stratum-adjusted correlation coefficient (SCC) between two
    Hi-C matrices up to max_dist. A Pearson correlation coefficient is computed
    for each diagonal in the range of 0 to max_dist and a weighted sum of those
    coefficients is returned.
    Parameters
    ----------
    mat1 : scipy.sparse.csr_matrix
        First matrix to compare.
    mat2 : scipy.sparse.csr_matrix
        Second matrix to compare.
    max_bins : int
        Maximum distance at which to consider, in bins.
    Returns
    -------
    scc : float
        Stratum adjusted correlation coefficient.
    """
    corr_diag = np.zeros(len(range(max_bins)))
    weight_diag = corr_diag.copy()
    for d in range(max_bins):
        d1 = mat1.diagonal(d)
        d2 = mat2.diagonal(d)
        # Silence NaN warnings: this happens for empty diagonals and will
        # not be used in the end.
        with warnings.catch_warnings():
            # Warning does not exist in older scipy versions (<=1.2)
            try:
                warnings.filterwarnings(
                    "ignore", category=ss.PearsonRConstantInputWarning
                )
            except AttributeError:
                pass
            # Compute raw pearson coeff for this diag
            corr_diag[d] = ss.pearsonr(d1, d2)[0]
        # Compute weight for this diag
        r2k = cu.vstrans(d1, d2)
        weight_diag[d] = len(d1) * r2k
    # Normalize weights
    weight_diag /= sum(weight_diag)

    # Weighted sum of coefficients to get SCCs
    scc = np.sum(corr_diag * weight_diag)

    return scc