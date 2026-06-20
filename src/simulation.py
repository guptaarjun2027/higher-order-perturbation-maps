"""
simulation.py

Solo reconstruction of the perturbed quadratic complex map simulation.

Original map (n=2):
    F(z) = e^(i*theta) * z  +  lambda * z^2  +  epsilon * z^(-2)

This file implements the three core components needed to reproduce
and extend the original research:

    1. Iteration engine     — apply F to a grid of starting points z0
    2. Escape detection     — |z| >= 20  (Growth-Valid Criterion)
    3. Pitch measurement    — Delta_n = d(phi)/d(log r) in the asymptotic tail

Parameters confirmed from original research:
    theta   = 0.7
    lambda  = 0.5
    epsilon = 0.2
    escape_radius = 20.0
"""

import numpy as np


# ---------------------------------------------------------------------------
# PARAMETERS
# ---------------------------------------------------------------------------

THETA   = 0.7
LAM     = 0.5
EPS     = 0.2
R_ESC   = 20.0      # escape radius — the Growth-Valid Criterion
MAX_ITER = 500      # hard cap per trajectory


# ---------------------------------------------------------------------------
# COMPONENT 1: ITERATION ENGINE
# ---------------------------------------------------------------------------

def iterate_grid(theta=THETA, lam=LAM, eps=EPS, n=2,
                  r_esc=R_ESC, max_iter=MAX_ITER,
                  resolution=1200, extent=3.0):
    """
    Apply F_n repeatedly to a dense grid of starting points z0.

    The general map is:
        F_n(z) = e^(i*theta) * z  +  lam * z^2  +  eps * z^(-n)

    n=2 is the original study. n=3,4,5 are the extension.

    Parameters
    ----------
    resolution : int
        Side length of the square grid. 1200 gives 1,440,000 trajectories,
        matching the original study exactly.
    extent : float
        Grid runs from -extent to +extent on both axes.

    Returns
    -------
    escape_iters : 2D int array, shape (resolution, resolution)
        Iteration count at which each point escaped.
        -1 means the point never escaped (trapped, Crash Set K).
    escaped : 2D bool array
        True where the point escaped (Escape Set E).
    z0_grid : 2D complex array
        The starting points, for reference.
    """
    # Build the starting-point grid
    axis = np.linspace(-extent, extent, resolution)
    Re, Im = np.meshgrid(axis, axis)
    z0_grid = Re + 1j * Im

    # Working copy — we iterate this in place
    z = z0_grid.astype(complex)

    # Guard: z=0 causes division blow-up in the z^(-n) term
    z[z == 0] = 1e-12

    escape_iters = np.full((resolution, resolution), -1, dtype=int)
    escaped      = np.zeros((resolution, resolution), dtype=bool)
    active       = np.ones( (resolution, resolution), dtype=bool)  # not yet decided

    rot = np.exp(1j * theta)   # precompute the rotation factor

    for i in range(1, max_iter + 1):
        z_a = z[active]
        z_a = np.where(z_a == 0, 1e-12, z_a)   # safety guard

        # Apply the map F_n
        z[active] = rot * z_a  +  lam * z_a**2  +  eps * z_a**(-n)

        # Check escape condition: |z| >= escape radius
        just_escaped        = active & (np.abs(z) >= r_esc)
        escape_iters[just_escaped] = i
        escaped             |= just_escaped
        active              &= ~just_escaped      # remove newly escaped points

        if not active.any():   # early exit if everything is decided
            break

    return escape_iters, escaped, z0_grid


def iterate_single(z0, theta=THETA, lam=LAM, eps=EPS, n=2,
                    r_esc=R_ESC, max_iter=MAX_ITER):
    """
    Iterate F_n from a single starting point z0, recording the full trajectory.
    Needed for pitch analysis (which requires the angular history of each orbit).

    Returns
    -------
    trajectory : 1D complex array  (z0, z1, z2, ..., up to escape or max_iter)
    escaped    : bool
    iter_count : int or None  (iteration at which escape occurred)
    """
    rot = np.exp(1j * theta)
    traj = [z0]
    z = z0 if z0 != 0 else 1e-12

    for i in range(1, max_iter + 1):
        if z == 0:
            z = 1e-12
        z = rot * z  +  lam * z**2  +  eps * z**(-n)
        traj.append(z)
        if abs(z) >= r_esc:
            return np.array(traj), True, i

    return np.array(traj), False, None


# ---------------------------------------------------------------------------
# COMPONENT 2: ESCAPE DETECTION
# (handled inline in the iteration functions above,
#  but exposed here explicitly for clarity and testing)
# ---------------------------------------------------------------------------

def has_escaped(z, r_esc=R_ESC):
    """
    Growth-Valid Criterion: a point z is considered to have escaped once
    its magnitude reaches or exceeds the escape radius.
    This isolates the asymptotic regime where the lambda*z^2 term dominates.
    """
    return abs(z) >= r_esc


# ---------------------------------------------------------------------------
# COMPONENT 3: PITCH MEASUREMENT
# ---------------------------------------------------------------------------

def measure_pitch(trajectory, tail_fraction=0.4):
    """
    Measure the spiral pitch kappa of an escaped trajectory using the
    stability metric:

        Delta_n = d(phi) / d(log r)

    where phi = unwrapped angular coordinate, r = |z|.

    If the orbit is converging to a logarithmic spiral phi = kappa * log(r),
    then Delta_n should stabilize to a constant kappa. The original study
    found this does NOT happen — Delta_n remains scattered and kappa_hat
    is statistically indistinguishable from zero, confirming decoupling.

    The analysis is performed only on the TAIL of the trajectory
    (the final `tail_fraction` of iterates after escape) to isolate true
    asymptotic behavior from the near-origin transient noise.

    Parameters
    ----------
    trajectory : complex array, the full orbit for one starting point
    tail_fraction : float, fraction of trajectory to treat as the tail

    Returns
    -------
    kappa_hat : float
        Best-fit spiral pitch (regression slope of phi vs log r over the tail).
        Returns np.nan if the tail is too short to measure reliably.
    delta_n : 1D float array
        Local pitch estimates at each step in the tail (finite differences).
    """
    r   = np.abs(trajectory)
    phi = np.unwrap(np.angle(trajectory))   # unwrap to make phi continuous
    log_r = np.log(r + 1e-12)              # log(|z|), small guard avoids log(0)

    # Use only the tail — the asymptotic window past the origin transient
    n_total    = len(trajectory)
    tail_start = max(int(n_total * (1 - tail_fraction)), 1)

    log_r_tail = log_r[tail_start:]
    phi_tail   = phi[tail_start:]

    if len(log_r_tail) < 3:    # need enough points to measure a trend
        return np.nan, np.array([])

    # Local pitch at each step: Delta_n = d(phi)/d(log r)
    delta_n = np.diff(phi_tail) / (np.diff(log_r_tail) + 1e-12)

    # Regression: fit phi = kappa * log(r) + const over the tail
    # kappa_hat is the slope — the best estimate of spiral pitch
    A = np.vstack([log_r_tail, np.ones(len(log_r_tail))]).T
    kappa_hat = np.linalg.lstsq(A, phi_tail, rcond=None)[0][0]

    return kappa_hat, delta_n


def pitch_statistics(kappa_list, n_windows=11):
    """
    Compute summary statistics for a collection of per-trajectory kappa_hat
    values, using the n=11 independent windows approach from the original study.

    The key test: does the 95% confidence interval cross zero?
    If yes => measured pitch is statistically indistinguishable from zero
           => angular decoupling confirmed.

    Returns a dict with mean, SE, 95% CI, and whether CI crosses zero.
    """
    from scipy import stats

    kappas = np.array([k for k in kappa_list if not np.isnan(k)])
    if len(kappas) < n_windows:
        return None

    # Shuffle and split into n_windows non-overlapping groups
    rng = np.random.default_rng(42)
    rng.shuffle(kappas)
    wsize = len(kappas) // n_windows
    window_means = [
        np.mean(kappas[w * wsize : (w+1) * wsize])
        for w in range(n_windows)
    ]
    window_means = np.array(window_means)

    mean_k = np.mean(window_means)
    se     = np.std(window_means, ddof=1) / np.sqrt(n_windows)
    t_crit = stats.t.ppf(0.975, df=n_windows - 1)
    ci_hw  = t_crit * se

    return {
        "n_trajectories" : len(kappas),
        "n_windows"      : n_windows,
        "mean_kappa"     : mean_k,
        "se"             : se,
        "ci_halfwidth"   : ci_hw,
        "ci_lower"       : mean_k - ci_hw,
        "ci_upper"       : mean_k + ci_hw,
        "ci_crosses_zero": (mean_k - ci_hw) < 0 < (mean_k + ci_hw),
    }