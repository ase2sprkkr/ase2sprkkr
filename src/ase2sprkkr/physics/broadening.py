from __future__ import annotations
import numpy as np
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
            raise ValueError("x must contain unique values.")
        return dx

    x = _ensure_1d_array(x_in, "x_in")
    if x_out is None:
        xout = x.copy()
    else:
        xout = _ensure_1d_array(x_out, "x_out")

    sort_in = np.argsort(x)
    x_sorted = x[sort_in]
    sort_out = np.argsort(xout)
    xout_sorted = xout[sort_out]

    N = x_sorted.size
    M = xout_sorted.size

    if sigma is None and fwhm is None:
        raise ValueError("Either sigma or fwhm must be provided.")
    if fwhm is not None:
        if np.ndim(fwhm) == 0:
            sigma_vals = np.full(N, _fwhm_to_sigma(float(fwhm)))
        else:
            fwhm_arr = np.asarray(fwhm, dtype=float)
            if fwhm_arr.shape != (N,):
                raise ValueError("fwhm must be scalar or same length as x_in.")
            sigma_vals = _fwhm_to_sigma(fwhm_arr[sort_in])
    else:
        if np.ndim(sigma) == 0:
            sigma_vals = np.full(N, float(sigma))
        else:
            sigma_arr = np.asarray(sigma, dtype=float)
            if sigma_arr.shape != (N,):
                raise ValueError("sigma must be scalar or same length as x_in.")
            sigma_vals = sigma_arr[sort_in]

    if np.any(sigma_vals < 0):
        raise ValueError("sigma must be non-negative.")

    dx_in = _compute_bin_widths(x_sorted)
    dx_out = _compute_bin_widths(xout_sorted)

    def gaussian_weight_matrix() -> Tuple[Union[np.ndarray, "sp.csr_matrix"], np.ndarray]:
        """
        Build weight matrix W such that y_out = W @ y_in_sorted (if include_dx=True).

        Parameters (high level)
        - x_in, x_out: input/output grids (1D). If x_out is None, x_out := x_in.
        - sigma/fwhm: scalar or array of length len(x_in), associated point-wise with x_in.
        - include_dx: if True, W columns already include Δx_in factor.
        - discrete_normalize: if True, normalize each kernel column using Δx_out so
                             discrete sums approximate continuous integrals.
        - method: "auto", "dense", or "sparse".
        - return_sparse: if True and scipy.sparse is available, returns csr matrix in sparse mode.
        """
        approx_ops = int(N) * int(M)
        if method == "auto":
            method_use = "dense" if approx_ops <= 10_000_000 else "sparse"
        elif method in ("dense", "sparse"):
            method_use = method
        else:
            raise ValueError("method must be 'auto', 'dense', or 'sparse'")

        zero_mask = sigma_vals == 0.0

        if method_use == "dense":
            Xout = xout_sorted[:, None]
            Xin = x_sorted[None, :]
            Sigma = sigma_vals[None, :]
            with np.errstate(divide="ignore", invalid="ignore"):
                diff = (Xout - Xin) / Sigma
                K = np.exp(-0.5 * diff * diff) / (np.sqrt(2.0 * np.pi) * Sigma)
                if np.any(zero_mask):
                    K[:, zero_mask] = 0.0

            if discrete_normalize:
                col_norm = (K * dx_out[:, None]).sum(axis=0)
                col_norm_safe = col_norm.copy()
                col_norm_safe[col_norm_safe == 0] = 1.0
                K = K / col_norm_safe[None, :]

            W_sorted = K * dx_in[None, :] if include_dx else K

            if np.any(zero_mask):
                cols = np.nonzero(zero_mask)[0]
                for col in cols:
                    xi = x_sorted[col]
                    j = np.searchsorted(xout_sorted, xi)
                    j = np.clip(j, 0, M - 1)
                    if j > 0 and abs(xout_sorted[j - 1] - xi) <= abs(xout_sorted[j] - xi):
                        j -= 1
                    W_sorted[:, col] = 0.0
                    W_sorted[j, col] = dx_in[col] if include_dx else 1.0

            inv_sort_out = np.argsort(sort_out)
            W = W_sorted[inv_sort_out, :]
            return W, xout.copy()

        rows = []
        cols = []
        data = []
        for i in range(N):
            si = sigma_vals[i]
            xi = x_sorted[i]
            dxi = dx_in[i]
            if si == 0.0:
                j = np.searchsorted(xout_sorted, xi)
                j = np.clip(j, 0, M - 1)
                if j > 0 and abs(xout_sorted[j - 1] - xi) <= abs(xout_sorted[j] - xi):
                    j -= 1
                rows.append(j)
                cols.append(i)
                data.append(dxi if include_dx else 1.0)
                continue

            cutoff = n_sigma * si
            xmin = xi - cutoff
            xmax = xi + cutoff
            j0 = np.searchsorted(xout_sorted, xmin, side="left")
            j1 = np.searchsorted(xout_sorted, xmax, side="right")
            if j0 >= j1:
                continue
            xseg = xout_sorted[j0:j1]
            diff = (xseg - xi) / si
            Kseg = np.exp(-0.5 * diff * diff) / (np.sqrt(2.0 * np.pi) * si)

            if discrete_normalize:
                norm = (Kseg * dx_out[j0:j1]).sum()
                if norm == 0:
                    norm = 1.0
                Kseg = Kseg / norm

            wseg = Kseg * dxi if include_dx else Kseg
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

        W_sorted_csr = sp.csr_matrix((data, (rows, cols)), shape=(M, N))
        inv_sort_out = np.argsort(sort_out)
        W = W_sorted_csr[inv_sort_out, :]
        if return_sparse:
            return W, xout.copy()
        return W.toarray(), xout.copy()

    W = gaussian_weight_matrix()[0]

    def broaden(y_orig):
        """
        Apply Gaussian broadening along the last (energy) axis.

        The values on that axis are associated point-wise with ``x_in`` and are
        therefore reordered together with an unsorted input grid before applying
        the precomputed matrix.

        Parameters
        ----------
        y_orig : array_like, shape (..., len(x_in))
            Function values to broaden along the last axis.

        Returns
        -------
        y_broadened : ndarray, shape (..., len(x_out))
            Broadened values in the original ``x_out`` order.
        """
        y = np.asarray(y_orig)
        if y.ndim == 0 or y.shape[-1] != N:
            raise ValueError("The last axis of y_orig must have the same length as x_in.")
        y_sorted = y[..., sort_in]
        flat = y_sorted.reshape((-1, N)).T
        out = W @ flat
        return np.asarray(out).T.reshape(y.shape[:-1] + (M,))

    return broaden


def create_lorentz_broadener(x_orig, gamma, wlortab=None, nelag=5, ylag=True):
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
        Number of grid points available for YLAG extrapolation (Fortran NELAG=5)
    ylag : bool, optional
        If True (default), extrapolate F below x[0] for the low-energy tail using
        linear Lagrange interpolation of the first two grid points, matching Fortran
        VECLORBRD/YLAG (N1=2).  If False, assume constant F = F[x[0]] below x[0].

    Returns
    -------
    callable
        Function ``broaden_y(y)`` that applies the precomputed broadening matrix.
    """
    x = np.asarray(x_orig)
    N = len(x)

    if wlortab is None:
        wlortab = np.zeros(N)
    wlortab = np.asarray(wlortab)

    L = np.zeros((N, N), dtype=float)
    de = max(x[1] - x[0], 0.5)  # Fortran: DE = MAX(E(2)-E(1), 0.5)

    # --- Low-energy tail: 4 intervals below x[0], matching Fortran DO J=1,4 ---
    j_arr = np.arange(1, 5)[:, None]  # shape (4,1)
    x1_low = x[0] - j_arr * de  # shape (4,1)
    x2_low = x1_low + de  # shape (4,1)
    G_low = gamma + wlortab[0]
    W_low = (np.arctan((x2_low - x[None, :]) / G_low) - np.arctan((x1_low - x[None, :]) / G_low)) / np.pi
    # W_low shape: (4, N)

    if ylag and N >= 2:
        # YLAG with N1=2: linear extrapolation using x[0] and x[1].
        # For interval j the Fortran evaluates F at X1=x[0]-j*de and X2=x[0]-(j-1)*de,
        # then accumulates W*(F(X1)+F(X2)).
        # Decomposed into the L matrix:
        #   F(X) = F[0]*(x[1]-X)/dx01 + F[1]*(X-x[0])/dx01
        # so F(X1)+F(X2) = F[0]*(2+(2j-1)*de/dx01) + F[1]*(-(2j-1)*de/dx01)
        dx01 = x[1] - x[0]
        total_coef0 = 2 + (2 * j_arr - 1) * de / dx01  # shape (4, 1)
        total_coef1 = -(2 * j_arr - 1) * de / dx01  # shape (4, 1)
        L[:, 0] += (W_low * total_coef0).sum(axis=0)
        L[:, 1] += (W_low * total_coef1).sum(axis=0)
    else:
        # Constant extrapolation: F(X) = F[x[0]] for all X < x[0]
        L[:, 0] += W_low.sum(axis=0)

    # --- Main contribution ---
    x_edges_left = x[:-1]
    x_edges_right = x[1:]
    G_main = gamma + 0.5 * (wlortab[:-1] + wlortab[1:])

    atan_right = np.arctan((x_edges_right[:, None] - x[None, :]) / G_main[:, None])
    atan_left = np.arctan((x_edges_left[:, None] - x[None, :]) / G_main[:, None])
    W_main = (atan_right - atan_left) / np.pi

    # Trapezoid split: half to left, half to right
    L[:, :-1] += W_main.T / 2
    L[:, 1:] += W_main.T / 2

    # --- High-energy tail (vektorizovaně) ---
    x1_high = x[-1] + np.arange(0, 4)[:, None] * de
    x2_high = x1_high + de
    G_high = gamma + wlortab[-1]
    W_high = (np.arctan((x2_high - x[None, :]) / G_high) - np.arctan((x1_high - x[None, :]) / G_high)) / np.pi
    L[:, -1] += W_high.sum(axis=0)

    def broaden_y(y_orig):
        """
        Applies the precomputed Lorentzian broadening to y_orig using matrix multiplication.
        This nested function has access to the precomputed weight_matrix.
        """
        return np.dot(L, y_orig.reshape(y_orig.shape[:1] + (-1,))).reshape(y_orig.shape)

    return broaden_y


def compute_wlortab(energies, n_valence):
    """
    Compute energy-dependent Lorentz broadening correction WLORTAB(E) as in
    Fortran INILORBRD (xband plotrxas).

    The total Lorentz half-width at energy E is: G(E) = G0 + WLORTAB(E).

    Parameters
    ----------
    energies : array_like
        Energy grid (eV), relative to E_Fermi.
    n_valence : int
        Number of valence electrons.  Special case: n_valence == 11 is
        treated as n_valence == 1 (matching Fortran logic).
        Returns zeros if n_valence is outside [1, 11].

    Returns
    -------
    wlortab : ndarray, shape (len(energies),)
    """
    CF1 = np.array(
        [
            0.001766920594,
            -0.184509337930,
            0.409983271398,
            0.067459351773,
            -0.610733613427,
            0.451717419499,
            -0.124878965835,
            0.014987644299,
            -0.000661473792,
        ]
    )
    CF11 = np.array(
        [
            0.048696144050,
            -0.537083543337,
            0.492816604919,
            0.398174247185,
            -0.444398434205,
            0.238062954735,
            -0.061552896506,
            0.007317769806,
            -0.000324445515,
        ]
    )

    nval = 1 if n_valence == 11 else n_valence
    if nval < 1 or nval > 10:
        return np.zeros(len(energies))

    w1 = (11 - nval) / 10.0
    w11 = (nval - 1) / 10.0
    CFZ = w1 * CF1 + w11 * CF11

    energies = np.asarray(energies, dtype=float)
    wlortab = np.zeros(len(energies))

    mask = energies > 2.5
    if mask.any():
        loge = np.log(energies[mask])
        pwrs = np.cumprod(np.outer(loge, np.ones(len(CFZ) - 1)), axis=1)
        pwrs = np.hstack([np.ones((pwrs.shape[0], 1)), pwrs])
        wlortab[mask] = pwrs @ CFZ

    idx25 = np.searchsorted(energies, 2.5, side="right")
    idx25 = np.clip(idx25, 1, len(energies) - 1)
    e_lo, e_hi = energies[idx25 - 1], energies[idx25]
    if e_hi > e_lo:
        t = (2.5 - e_lo) / (e_hi - e_lo)
        y25 = (1 - t) * wlortab[idx25 - 1] + t * wlortab[idx25]
    else:
        y25 = wlortab[idx25]

    mask_low = (energies > 0) & (energies <= 2.5)
    if mask_low.any():
        wlortab[mask_low] = y25 * energies[mask_low] / 2.5

    wlortab = np.maximum(0.0001, wlortab)
    wlortab[energies <= 0] = 0.0

    return wlortab
