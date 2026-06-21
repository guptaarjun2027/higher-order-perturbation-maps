# Research Hypotheses

**Project:** Higher-Order Perturbations in the Decoupled Escape Family  
**Researcher:** [Arjun Gupta]  
**Date:** [June 19th, 2026]  
**Status:** Written prior to running extension_study.py

---

## Background

The original research investigated the map F(z) = e^(iθ)z + λz² + εz⁻² (the n=2 case of the general family F_n(z) = e^(iθ)z + λz² + εz⁻ⁿ) and found that among trajectories that crossed the outer escape radius R₀, beyond which the quadratic term λz² dominates and the perturbation εz⁻² decay, certain phenomena occured. First, outwards growth was monotonic: after each iteration the point moved further from the diagram. Second, radial growth and angular change decopuled, this was due to angular scarring as the orbits escaped to infinity without converging to the predicted logarithmic spiral geometry φ ≈ κ log r. These findings led to the question whether similar phenomena would occur with higher order perturbations where epsilon is divided by z^(n) for n > 2.
[2-3 sentences summarizing the original finding — the decoupling result,
what angular scarring means, and what the natural open question is.]

## Map Under Study

F_n(z) = e^(i*theta) * z  +  lambda * z^2  +  eps_n * z^(-n)

Parameters: theta=0.7, lambda=0.5, eps_2=0.2, escape radius r >= 20
Perturbation orders studied: n = 2 (baseline), 3, 4, 5

## Research Question

Does the decoupling of radial escape and spiral pitch observed in the
n=2 case generalize across the perturbation family

    F_n(z) = e^(iθ)z + λz² + εz⁻ⁿ

for n = 3, 4, 5 — and if so, does the degree of angular scarring scale
systematically with the perturbation order n?

## Epsilon Scaling Rationale

At fixed ε, the perturbation term ε|z|⁻ⁿ grows by a factor of |z|⁻¹
for each unit increase in n near the origin, which would confound
perturbation order with perturbation strength across the comparison.
To isolate the effect of n alone, ε is scaled so that the perturbation
magnitude is equal across all n at a reference radius r_ref = 0.5,
which sits within the near-origin zone where the angular scar is formed:

    ε_n = ε_2 * r_ref^(n-2)

giving ε_2=0.200, ε_3=0.100, ε_4=0.050, ε_5=0.025.

## Hypotheses

### Hypothesis A — [Monotonic Strengthening (Scarring Dominates)]
The angular scar deepens with n because the near-origin perturbation grows more violent faster than the outer-regime decay provides healing. Mean kappa stays persistently nonzero and either increases or remains large across n. The trapped basin shrinks (higher-order poles destabilize more of the origin neighborhood) and the fractal boundary becomes more intricate. The perturbation magnitude near a reference radius r₀ < 1 scales as ε/r₀ⁿ, which grows without bound as n increases even after epsilon scaling, because the scaling only fixes the magnitude at exactly r_ref — everywhere else inside r_ref the perturbation is stronger.

### Hypothesis B — [Saturation (Rotation Dominates)]
The decoupling persists for all n but the degree of angular scarring doesn't systematically change, because the structural barrier is the rotational term e^(iθ)z, not the perturbation term. The perturbation creates the initial scar but the rotation prevents convergence regardless of how deep the scar is. Mean kappa is statistically similar across n=2,3,4,5. Pitch pass rate stays 0% throughout. The rotation term is invariant across the family. The poster's own conclusion identified the linear rotation term as the "structural barrier" — if that's truly the dominant mechanism, then changing n should have second-order effects at most.

### Hypothesis C — [Weakening (Recovery Dominates)]
Higher-order perturbations decay fast enough in the outer escape regime that orbits partially recover angular structure. Mean kappa trends toward zero as n increases, the CI becomes more likely to cross zero for larger n, and the trapped basin grows (higher-order poles create a larger region of origin-dominated dynamics that traps more points).  At r=20 (the escape radius), εz⁻ⁿ = ε/20ⁿ. For n=2 this is ε/400. For n=5 this is ε/3,200,000 — eight thousand times smaller. The angular scarring mechanism requires the perturbation to influence the trajectory; if it's negligible throughout the asymptotic window, it can't disrupt spiral formation.

## Pre-stated Prediction

I predict Hypothesis A will hold: as n increases, the more violent near-origin perturbation will accelerate angular scarring, causing decoupling to manifest earlier in each trajectory's escape path and producing a stronger suppression of spiral convergence. Concurrently, the intensifying pole at the origin is predicted to destabilize an increasing proportion of the near-origin neighborhood, causing the trapped basin (Crash Set K) to contract substantially as n increases toward n=5. This would suggest that perturbation order is the dominantfactor governing angular scar severity, rather than the fixed rotational term alone.