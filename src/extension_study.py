"""
extension_study.py

THE CORE EXTENSION — run this after validate_baseline.py passes.

Research question:
    Does the decoupling of radial escape and spiral pitch generalize
    across the perturbation family:

        F_n(z) = e^(i*theta) * z  +  lambda * z^2  +  eps * z^(-n)

    for n = 2, 3, 4, 5?

    And if so, does the degree of angular scarring scale predictably
    with the perturbation order n?

Three hypotheses (stated BEFORE seeing results — do not change after):

    A. MONOTONIC STRENGTHENING
       Decoupling intensifies with n. Higher-order inverse terms create
       more severe near-origin interference, leaving deeper angular scars.

    B. SATURATION
       Decoupling is present for all n >= 2 but does not monotonically
       increase — the outer escape regime is dominated by lambda*z^2
       regardless of perturbation order.

    C. WEAKENING / REVERSAL
       Higher-order perturbations decay faster (|eps*z^-n| -> 0 faster
       for large n), so orbits recover angular structure for large n,
       and decoupling weakens or disappears.

Usage:
    python extension_study.py              # quick run (300x300, ~seconds)
    python extension_study.py --full       # full run (1200x1200, ~minutes)
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from simulation import (
    iterate_grid, iterate_single, measure_pitch, pitch_statistics,
    THETA, LAM, EPS, R_ESC
)


# ---------------------------------------------------------------------------
# EPSILON SCALING
# ---------------------------------------------------------------------------
#
# Critical methodological decision:
#
# At a fixed epsilon, the perturbation term eps*z^(-n) grows MUCH stronger
# near the origin as n increases (since |z^(-n)| -> inf faster for larger n).
# This means a naive comparison at fixed epsilon confounds two things:
#   (a) the effect of perturbation ORDER (what we want to study)
#   (b) the effect of perturbation STRENGTH (a confounder)
#
# To isolate (a), we scale epsilon so that the perturbation term has the
# same magnitude at a reference radius r_ref:
#
#     |eps_n * r_ref^(-n)| = |eps_2 * r_ref^(-2)|  for all n
#     =>  eps_n = eps_2 * r_ref^(-(2-n)) = eps_2 * r_ref^(n-2)
#
# This means higher-order terms use a SMALLER epsilon, keeping the
# perturbation strength constant at r_ref.
#
# We use r_ref = 0.5 (inside the near-origin perturbation zone, where
# the angular scar is created, matching the physical motivation).
#
# This scaling is a scientific choice that MUST be stated clearly in the
# paper's methodology section. It makes the n comparison fair.

R_REF = 0.5   # reference radius for epsilon scaling

def scaled_eps(n, eps2=EPS, r_ref=R_REF):
    """Epsilon value for perturbation order n, scaled so that the
    perturbation strength at r_ref is the same as the original eps at n=2."""
    return eps2 * (r_ref ** (n - 2))


# ---------------------------------------------------------------------------
# SINGLE-N STUDY
# ---------------------------------------------------------------------------

def run_one_n(n, resolution, n_pitch_sample=3000, seed=42, verbose=True):
    """Run the full study for one value of perturbation order n."""
    eps_n = scaled_eps(n)
    if verbose:
        print(f"\n  n={n}  (eps scaled to {eps_n:.5f} at r_ref={R_REF})")

    # Grid iteration
    escape_iters, escaped, z0_grid = iterate_grid(
        theta=THETA, lam=LAM, eps=eps_n, n=n, r_esc=R_ESC, resolution=resolution
    )
    n_escaped = int(escaped.sum())
    n_trapped = int((~escaped).sum())

    if n_escaped == 0:
        if verbose:
            print(f"    WARNING: no escaped points — adjust epsilon scaling")
        return None

    # Sample escaped trajectories for pitch analysis.
    # Filter to slow escapers only (>= 15 iterations) — fast-escaping
    # trajectories have too few tail points for reliable pitch measurement.
    rng = np.random.default_rng(seed)
    slow_mask  = escaped & (escape_iters >= 15)
    escaped_z0 = z0_grid[slow_mask]
    if verbose:
        print(f"    Slow escapers (>= 15 iters): {len(escaped_z0):,}")
    sample = escaped_z0[rng.choice(
        len(escaped_z0), size=min(n_pitch_sample, len(escaped_z0)), replace=False
    )]

    kappa_hats = []
    for z0 in sample:
        traj, did_escape, _ = iterate_single(z0, theta=THETA, lam=LAM, eps=eps_n, n=n)
        if did_escape:
            kh, _ = measure_pitch(traj)
            if not np.isnan(kh):
                kappa_hats.append(kh)

    stats = pitch_statistics(kappa_hats)
    kappa_pred = THETA / np.log(max(LAM, 1 + 1e-9))
    pass_rate  = 100 * np.mean(np.abs(np.array(kappa_hats) - kappa_pred) <= 0.1)

    result = {
        "n"              : n,
        "eps_n"          : eps_n,
        "n_escaped"      : n_escaped,
        "n_trapped"      : n_trapped,
        "escape_iters"   : escape_iters,
        "n_kappa"        : len(kappa_hats),
        "pitch_pass_rate": pass_rate,
        "stats"          : stats,
    }

    if verbose and stats:
        print(f"    Escaped: {n_escaped:,}   Trapped: {n_trapped:,}")
        print(f"    mean kappa = {stats['mean_kappa']:.6f}   "
              f"CI = [{stats['ci_lower']:.6f}, {stats['ci_upper']:.6f}]")
        print(f"    CI crosses zero: {stats['ci_crosses_zero']}   "
              f"Pitch pass rate: {pass_rate:.2f}%")

    return result

def save_results_csv(results, filepath="results/comparison_table.csv"):
    """Export the comparison table to CSV for the paper and repository."""
    import csv
    import os

    os.makedirs(os.path.dirname(filepath), exist_ok=True)

    rows = []
    for n, r in results.items():
        s = r["stats"]
        rows.append({
            "n"              : n,
            "eps_n"          : r["eps_n"],
            "n_escaped"      : r["n_escaped"],
            "n_trapped"      : r["n_trapped"],
            "n_kappa_samples": r["n_kappa"],
            "mean_kappa"     : round(s["mean_kappa"], 6) if s else "N/A",
            "se"             : round(s["se"], 6) if s else "N/A",
            "ci_lower"       : round(s["ci_lower"], 6) if s else "N/A",
            "ci_upper"       : round(s["ci_upper"], 6) if s else "N/A",
            "ci_crosses_zero": s["ci_crosses_zero"] if s else "N/A",
            "pitch_pass_rate": round(r["pitch_pass_rate"], 2),
        })

    fieldnames = list(rows[0].keys())
    with open(filepath, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nSaved: {filepath}")
# ---------------------------------------------------------------------------
# FULL COMPARISON: n = 2, 3, 4, 5
# ---------------------------------------------------------------------------

def run_full_comparison(resolution=300, n_pitch_sample=3000):
    print("=" * 65)
    print("EXTENSION STUDY: Higher-Order Perturbation Family")
    print(f"F_n(z) = e^(i*{THETA})*z + {LAM}*z^2 + eps_n * z^(-n)")
    print(f"Resolution: {resolution}x{resolution}  ({resolution**2:,} trajectories)")
    print("=" * 65)

    results = {}
    for n in [2, 3, 4, 5]:
        r = run_one_n(n, resolution=resolution, n_pitch_sample=n_pitch_sample)
        if r:
            results[n] = r

    print_comparison_table(results)
    plot_all_fractals(results)
    plot_kappa_trend(results)
    save_results_csv(results)

    return results


def print_comparison_table(results):
    print("\n")
    print("=" * 90)
    print("COMPARISON TABLE")
    print("=" * 90)
    header = f"{'n':<5} {'eps_n':<8} {'escaped':>10} {'trapped':>10} {'mean_kappa':>12} {'CI_lower':>10} {'CI_upper':>10} {'CI_cross_0':>11} {'pass%':>7}"
    print(header)
    print("-" * 90)
    for n, r in results.items():
        s = r["stats"]
        if s:
            print(f"{n:<5} {r['eps_n']:<8.5f} {r['n_escaped']:>10,} {r['n_trapped']:>10,} "
                  f"{s['mean_kappa']:>12.6f} {s['ci_lower']:>10.6f} {s['ci_upper']:>10.6f} "
                  f"{str(s['ci_crosses_zero']):>11} {r['pitch_pass_rate']:>7.2f}%")
    print("=" * 90)
    print("\nConclusion guide:")
    print("  Pitch pass rate = 0% for all n  ->  decoupling generalizes")
    print("  mean_kappa trend across n        ->  tests Hypothesis A (monotonic), B (flat), C (decreasing)")
    print("  CI crosses zero                  ->  angular decoupling statistically confirmed for that n")


def plot_all_fractals(results):
    n_vals = list(results.keys())
    fig, axes = plt.subplots(1, len(n_vals), figsize=(5 * len(n_vals), 5))
    if len(n_vals) == 1:
        axes = [axes]

    for ax, n in zip(axes, n_vals):
        r = results[n]
        grid = r["escape_iters"].astype(float)
        grid[grid < 0] = np.nan
        ax.imshow(grid, cmap="inferno", origin="lower",
                   extent=[-3, 3, -3, 3],
                   norm=mcolors.LogNorm(vmin=1, vmax=500))
        s = r["stats"]
        kappa_str = f"κ={s['mean_kappa']:.4f}" if s else "κ=N/A"
        ax.set_title(f"n={n}  (eps={r['eps_n']:.4f})\n{kappa_str}  pass={r['pitch_pass_rate']:.0f}%",
                      fontsize=10)
        ax.set_xlabel("Re(z0)", fontsize=8)
        ax.set_ylabel("Im(z0)", fontsize=8)

    plt.suptitle("Escape-Time Fractals: F_n(z) = e^(i*theta)*z + lambda*z^2 + eps_n*z^(-n)",
                  fontsize=12, y=1.02)
    plt.tight_layout()
    plt.savefig("fractal_comparison.png", dpi=130, bbox_inches="tight")
    plt.close()
    print("\nSaved: fractal_comparison.png")


def plot_kappa_trend(results):
    """Plot mean kappa vs n — this is the key figure showing whether
    decoupling strengthens, saturates, or weakens with perturbation order."""
    ns = []
    means, lowers, uppers = [], [], []

    for n, r in results.items():
        s = r["stats"]
        if s:
            ns.append(n)
            means.append(s["mean_kappa"])
            lowers.append(s["ci_lower"])
            uppers.append(s["ci_upper"])

    if not ns:
        return

    ns = np.array(ns)
    means = np.array(means)
    lowers = np.array(lowers)
    uppers = np.array(uppers)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(ns, means, "o-", color="#4A90D9", linewidth=2, markersize=8,
             label="mean kappa (measured)")
    ax.fill_between(ns, lowers, uppers, alpha=0.25, color="#4A90D9",
                     label="95% CI")
    ax.axhline(0, color="gray", linestyle="--", linewidth=1, label="kappa = 0")

    ax.set_xlabel("Perturbation order n", fontsize=12)
    ax.set_ylabel("Mean spiral pitch (kappa)", fontsize=12)
    ax.set_title("Does angular decoupling scale with perturbation order?\n"
                  "F_n(z) = e^(i*theta)*z + lambda*z^2 + eps_n*z^(-n)", fontsize=11)
    ax.set_xticks(ns)
    ax.legend(fontsize=10)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("kappa_vs_n.png", dpi=130)
    plt.close()
    print("Saved: kappa_vs_n.png")


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--full", action="store_true",
                         help="Full 1200x1200 run (~minutes)")
    args = parser.parse_args()

    resolution = 1200 if args.full else 300
    run_full_comparison(resolution=resolution)