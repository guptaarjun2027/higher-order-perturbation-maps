"""
validate_baseline.py

Before running the extension (n=3,4,5), verify that the n=2 baseline
reproduces the original study's qualitative findings:

    - Pitch pass rate:  0%       (no trajectory stabilizes to the predicted spiral)
    - CI crosses zero:  True     (mean kappa statistically indistinguishable from zero)
    - 100% wedge detection       (outer escape region exists for all escaped points)

If these hold, the codebase is a valid foundation for the extension.
Run this FIRST before touching extension_study.py.

Usage:
    python validate_baseline.py
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from simulation import (
    iterate_grid, iterate_single, measure_pitch, pitch_statistics,
    THETA, LAM, EPS, R_ESC
)


def run_baseline_validation(resolution=300, n_pitch_sample=2000):
    print("=" * 60)
    print("BASELINE VALIDATION  (n=2, original map)")
    print(f"theta={THETA}, lambda={LAM}, epsilon={EPS}, R_esc={R_ESC}")
    print("=" * 60)

    # --- Step 1: Grid iteration ---
    print("\n[1] Running grid iteration...")
    escape_iters, escaped, z0_grid = iterate_grid(n=2, resolution=resolution)
    n_total   = resolution * resolution
    n_escaped = int(escaped.sum())
    n_trapped = n_total - n_escaped
    print(f"    Escaped (Set E): {n_escaped:,}  ({100*n_escaped/n_total:.1f}%)")
    print(f"    Trapped (Set K): {n_trapped:,}  ({100*n_trapped/n_total:.1f}%)")
    print(f"    Wedge detection: {'PASS — escaped region exists' if n_escaped > 0 else 'FAIL'}")

    # --- Step 2: Pitch analysis on a sample of escaped trajectories ---
    print(f"\n[2] Sampling {n_pitch_sample} escaped trajectories for pitch analysis...")
    rng = np.random.default_rng(42)
    escaped_z0 = z0_grid[escaped]
    sample     = escaped_z0[rng.choice(len(escaped_z0), size=min(n_pitch_sample, len(escaped_z0)), replace=False)]

    kappa_hats = []
    for z0 in sample:
        traj, did_escape, _ = iterate_single(z0, n=2)
        if did_escape:
            kh, _ = measure_pitch(traj)
            if not np.isnan(kh):
                kappa_hats.append(kh)

    print(f"    Valid kappa estimates: {len(kappa_hats):,}")

    # --- Step 3: Statistics ---
    print("\n[3] Computing statistics...")
    stats = pitch_statistics(kappa_hats)
    if stats:
        print(f"    mean kappa:        {stats['mean_kappa']:.6f}")
        print(f"    SE:                {stats['se']:.6f}")
        print(f"    95% CI:            [{stats['ci_lower']:.6f}, {stats['ci_upper']:.6f}]")
        print(f"    CI crosses zero?   {stats['ci_crosses_zero']}")
    
    # Pitch pass rate: fraction of trajectories whose kappa_hat is
    # "close" to the naive spiral prediction (theta/log(lambda)).
    # The original study found this to be 0%.
    kappa_pred  = THETA / np.log(max(LAM, 1 + 1e-9))
    pass_rate   = 100 * np.mean(np.abs(np.array(kappa_hats) - kappa_pred) <= 0.1)
    print(f"    Pitch pass rate:   {pass_rate:.2f}%  (target: 0%)")

    # --- Step 4: Quick fractal plot for visual confirmation ---
    print("\n[4] Generating fractal visualization for visual check...")
    _plot_fractal(escape_iters, resolution)

    # --- Summary ---
    print("\n" + "=" * 60)
    print("VALIDATION SUMMARY")
    print("=" * 60)
    checks = [
        ("Escaped region exists (wedge detection)", n_escaped > 0),
        ("Trapped basin exists",                    n_trapped > 0),
        ("Pitch pass rate = 0%",                    pass_rate == 0.0),
        ("CI crosses zero (decoupling confirmed)",  stats["ci_crosses_zero"] if stats else False),
    ]
    all_pass = True
    for name, result in checks:
        status = "PASS" if result else "FAIL"
        if not result:
            all_pass = False
        print(f"  [{status}]  {name}")
    print()
    if all_pass:
        print("All checks passed. Baseline is solid.")
        print("Proceed to extension_study.py")
    else:
        print("Some checks failed. Review before running the extension.")
    print("=" * 60)

    return all_pass


def _plot_fractal(escape_iters, resolution):
    grid = escape_iters.astype(float)
    grid[grid < 0] = np.nan   # trapped points -> white

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.imshow(grid, cmap="inferno", origin="lower",
               extent=[-3, 3, -3, 3],
               norm=mcolors.LogNorm(vmin=1, vmax=500))
    ax.set_title("Baseline (n=2): Escape-Time Fractal\n"
                  f"theta={THETA}, lambda={LAM}, eps={EPS}", fontsize=11)
    ax.set_xlabel("Re(z0)")
    ax.set_ylabel("Im(z0)")
    plt.tight_layout()
    plt.savefig("baseline_fractal.png", dpi=130)
    plt.close()
    print("    Saved: baseline_fractal.png  — compare visually to poster Fig. 6")


if __name__ == "__main__":
    run_baseline_validation()