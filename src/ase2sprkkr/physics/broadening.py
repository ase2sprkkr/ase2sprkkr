from __future__ import annotations
import numpy as np
import scipy
from typing import Optional, Sequence, Tuple, Union
try:
    import scipy.sparse as sp
except Exception:
    sp = None  # sparse support optional

def create_gaussian_broadener(
        x_in: Sequence[float],
        sigma: Optional[Union[float, Sequence[float]]] = None,
        x_out: Optional[Sequence[float]] = None,
        fwhm: Optional[Union[float, Sequence[float]]] = None,
        n_sigma: float = 6.0,
        method: str = "auto",
        return_sparse: bool = True,
        include_dx: bool = True,
        discrete_normalize: bool = True,
        ) -> callable:

    def _ensure_1d_array(arr, name: str) -> np.ndarray:
        a = np.asarray(arr)
        if a.ndim != 1:
            raise ValueError(f"{name} must be a 1D array-like")
        return a.astype(float)


    def _fwhm_to_sigma(fwhm: float) -> float:
        return fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))


    def _compute_bin_widths(x: np.ndarray) -> np.ndarray:
        """
        Compute representative bin widths (Δx) for grid points x (uneven).
        Interior point: (x[i+1] - x[i-1]) / 2
        Endpoints: nearest neighbor spacing.
        """
        n = x.size
        if n == 1:
            return np.array([1.0], dtype=float)
        dx = np.empty(n, dtype=float)
        dx[0] = x[1] - x[0]
        dx[-1] = x[-1] - x[-2]
        if n > 2:
            dx[1:-1] = 0.5 * (x[2:] - x[:-2])
        if np.any(dx <= 0):
            raise ValueError("x must be strictly increasing.")
        return dx


    def gaussian_weight_matrix(
    ) -> Tuple[Union[np.ndarray, "sp.csr_matrix"], np.ndarray]:
        """
        Build weight matrix W such that y_out = W @ y_in (if include_dx=True).

        Parameters (high level)
        - x_in, x_out: input/output grids (1D). If x_out is None, x_out := x_in.
        - sigma/fwhm: scalar or array of length len(x_in).
        - include_dx: if True, W columns already include Δx_in factor.
        - discrete_normalize: if True, normalize each kernel column using Δx_out so
                             discrete sums approximate continuous integrals.
        - method: "auto", "dense", or "sparse".
        - return_sparse: if True and scipy.sparse is available, returns csr matrix in sparse mode.
        """
        x = _ensure_1d_array(x_in, "x_in")
        if x_out is None:
            xout = x.copy()
        else:
            xout = _ensure_1d_array(x_out, "x_out")

        # sort for stable building
        sort_in = np.argsort(x)
        x_sorted = x[sort_in]

        sort_out = np.argsort(xout)
        xout_sorted = xout[sort_out]

        dx_in = _compute_bin_widths(x_sorted)
        dx_out = _compute_bin_widths(xout_sorted)

        N = x_sorted.size
        M = xout_sorted.size

        # sigma array
        if sigma is None and fwhm is None:
            raise ValueError("Either sigma or fwhm must be provided.")
        if fwhm is not None:
            if np.ndim(fwhm) == 0:
                sigma_vals = np.full(N, _fwhm_to_sigma(float(fwhm)))
            else:
                fwhm_arr = np.asarray(fwhm, dtype=float)
                if fwhm_arr.shape != (N,):
                    raise ValueError("fwhm must be scalar or same length as x_in.")
                sigma_vals = _fwhm_to_sigma(fwhm_arr)
        else:
            if np.ndim(sigma) == 0:
                sigma_vals = np.full(N, float(sigma))
            else:
                sigma_arr = np.asarray(sigma, dtype=float)
                if sigma_arr.shape != (N,):
                    raise ValueError("sigma must be scalar or same length as x_in.")
                sigma_vals = sigma_arr

        if np.any(sigma_vals < 0):
            raise ValueError("sigma must be non-negative.")

        approx_ops = int(N) * int(M)
        if method == "auto":
            method_use = "dense" if approx_ops <= 10_000_000 else "sparse"
        elif method in ("dense", "sparse"):
            method_use = method
        else:
            raise ValueError("method must be 'auto', 'dense', or 'sparse'")

        zero_mask = sigma_vals == 0.0

        # Dense: build K then optionally renormalize columns w.r.t. dx_out
        if method_use == "dense":
            Xout = xout_sorted[:, None]    # (M,1)
            Xin = x_sorted[None, :]        # (1,N)
            Sigma = sigma_vals[None, :]    # (1,N)
            with np.errstate(divide="ignore", invalid="ignore"):
                diff = (Xout - Xin) / Sigma
                K = np.exp(-0.5 * diff * diff) / (np.sqrt(2.0 * np.pi) * Sigma)
                if np.any(zero_mask):
                    K[:, zero_mask] = 0.0

            # discrete normalize each column if requested:
            if discrete_normalize:
                # column_norm = sum_j K[j,i] * dx_out_j  (approx. continuous integral)
                col_norm = (K * dx_out[:, None]).sum(axis=0)  # shape (N,)
                # avoid division by zero: for sigma==0 keep col_norm=1 (we will assign delta)
                col_norm_safe = col_norm.copy()
                col_norm_safe[col_norm_safe == 0] = 1.0
                K = K / col_norm_safe[None, :]

            # Now include dx_in if requested. If include_dx=False, W contains only kernel values K.
            if include_dx:
                W_sorted = K * dx_in[None, :]
            else:
                W_sorted = K

            # handle sigma==0 deltas: assign mass to nearest output bin
            if np.any(zero_mask):
                cols = np.nonzero(zero_mask)[0]
                for col in cols:
                    xi = x_sorted[col]
                    j = np.searchsorted(xout_sorted, xi)
                    j = np.clip(j, 0, M - 1)
                    if j > 0 and abs(xout_sorted[j - 1] - xi) <= abs(xout_sorted[j] - xi):
                        j -= 1
                    # zero the column and place delta mass:
                    W_sorted[:, col] = 0.0
                    if include_dx:
                        W_sorted[j, col] = dx_in[col]
                    else:
                        # kernel-only matrix: place a 1.0 at that bin so user can multiply by dx_in later
                        W_sorted[j, col] = 1.0

            # restore output order
            inv_sort_out = np.argsort(sort_out)
            W = W_sorted[inv_sort_out, :]
            return W, xout.copy()

        # Sparse: similar logic but construct lists for COO format
        else:
            rows = []
            cols = []
            data = []
            for i in range(N):
                si = sigma_vals[i]
                xi = x_sorted[i]
                dxi = dx_in[i]
                if si == 0.0:
                    # delta -> nearest output row
                    j = np.searchsorted(xout_sorted, xi)
                    j = np.clip(j, 0, M - 1)
                    if j > 0 and abs(xout_sorted[j - 1] - xi) <= abs(xout_sorted[j] - xi):
                        j -= 1
                    rows.append(j); cols.append(i)
                    data.append(dxi if include_dx else 1.0)
                    continue

                cutoff = n_sigma * si
                xmin = xi - cutoff; xmax = xi + cutoff
                j0 = np.searchsorted(xout_sorted, xmin, side="left")
                j1 = np.searchsorted(xout_sorted, xmax, side="right")
                if j0 >= j1:
                    continue
                xseg = xout_sorted[j0:j1]
                diff = (xseg - xi) / si
                Kseg = np.exp(-0.5 * diff * diff) / (np.sqrt(2.0 * np.pi) * si)

                if discrete_normalize:
                    # normalization scalar for this column
                    norm = (Kseg * dx_out[j0:j1]).sum()
                    if norm == 0:
                        norm = 1.0
                    Kseg = Kseg / norm

                if include_dx:
                    wseg = Kseg * dxi
                else:
                    wseg = Kseg

                rows.extend((np.arange(j0, j1)).tolist())
                cols.extend([i] * (j1 - j0))
                data.extend(wseg.tolist())

            if sp is None:
                W_sorted = np.zeros((M, N), dtype=float)
                for r, c, v in zip(rows, cols, data):
                    W_sorted[r, c] += v
                inv_sort_out = np.argsort(sort_out)
                W = W_sorted[inv_sort_out, :]
                return W, xout.copy()
            else:
                W_sorted_csr = sp.csr_matrix((data, (rows, cols)), shape=(M, N))
                inv_sort_out = np.argsort(sort_out)
                W = W_sorted_csr[inv_sort_out, :]
                if return_sparse:
                    return W, xout.copy()
                else:
                    return W.toarray(), xout.copy()

    W = gaussian_weight_matrix()[0]

    def broaden(y_orig):
        """
        Apply Gaussian broadening to F0 using the precomputed matrix.

        Parameters
        ----------
        y_orig : array_like, shape (NV, Nx_orig)
            Functions to broaden along the last axis (energy)

        Returns
        -------
        y_broadened : ndarray, shape (NV, Nx_orig)
            Broadened functions
        """
        out=W @ y_orig.reshape(y_orig.shape[:1]+(-1,) )
        return out.reshape(y_orig.shape)

    return broaden






def create_lorentz_broadener(x_orig, gamma, wlortab=None, nelag=5):
    """
    Lorentz broadening matrix L for uneven grid x_orig with optional energy-dependent width.

    Parameters
    ----------
    x_orig : array_like, shape (N,)
        Monotonically increasing grid
    gamma : float
        Base Lorentz HWHM
    wlortab : array_like, shape (N,), optional
        Energy-dependent correction to gamma
    nelag : int
        Number of low-energy tail intervals

    Returns
    -------
    L : ndarray, shape (N, N)
        Broadening matrix such that F_broadened = L @ F0
    """
    x = np.asarray(x_orig)
    N = len(x)

    if wlortab is None:
        wlortab = np.zeros(N)
    wlortab = np.asarray(wlortab)

    L = np.zeros((N, N), dtype=float)
    de = max(np.diff(x).min(), 0.5)

    # --- Low-energy tail (vektorizovaně) ---
    x1_low = x[0] - np.arange(1, nelag+1)[:, None]*de   # shape (nelag,1)
    x2_low = x1_low + de                                 # shape (nelag,1)
    G_low = gamma + wlortab[0]
    W_low = (np.arctan((x2_low - x[None,:])/G_low) - np.arctan((x1_low - x[None,:])/G_low)) / np.pi
    L[:,0] += W_low.sum(axis=0)

    # --- Main contribution ---
    x_edges_left = x[:-1]
    x_edges_right = x[1:]
    G_main = gamma + 0.5*(wlortab[:-1] + wlortab[1:])

    atan_right = np.arctan((x_edges_right[:, None] - x[None,:])/G_main[:, None])
    atan_left  = np.arctan((x_edges_left[:, None]  - x[None,:])/G_main[:, None])
    W_main = (atan_right - atan_left)/np.pi

    # Trapezoid split: half to left, half to right
    L[:, :-1] += W_main.T / 2
    L[:, 1:] += W_main.T / 2

    # --- High-energy tail (vektorizovaně) ---
    x1_high = x[-1] + np.arange(0, 4)[:, None]*de
    x2_high = x1_high + de
    G_high = gamma + wlortab[-1]
    W_high = (np.arctan((x2_high - x[None,:])/G_high) - np.arctan((x1_high - x[None,:])/G_high)) / np.pi
    L[:,-1] += W_high.sum(axis=0)

    def broaden_y(y_orig):
        """
        Applies the precomputed Lorentzian broadening to y_orig using matrix multiplication.
        This nested function has access to the precomputed weight_matrix.
        """
        return np.dot(L, y_orig.reshape(y_orig.shape[:1]+(-1,) )).reshape(y_orig.shape)

    return broaden_y



