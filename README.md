# Higher-Order Perturbations in the Decoupled Escape Family

Independent extension of collaborative research on perturbed quadratic
complex maps. Original work: "Decoupling of Radial Escape and Spiral Pitch
in a Perturbed Quadratic Complex Mapping" (S-MATH-003, SEFMD 2025,
with Himanish Pasula & Rohan Nagaram).

## Research Question
Does the decoupling of radial escape and spiral pitch generalize across
the perturbation family F_n(z) = e^(iθ)z + λz² + εz⁻ⁿ for n = 2, 3, 4, 5?

## Map
F_n(z) = e^(i*0.7)*z + 0.5*z² + eps_n * z^(-n)

Parameters: θ=0.7, λ=0.5, ε₂=0.2, escape radius r≥20

## Quick Start
pip install -r requirements.txt
python src/validate_baseline.py      # confirm n=2 baseline first
python src/extension_study.py        # quick run (300x300)
python src/extension_study.py --full # full run (1200x1200)

## Status
- [x] Baseline validated (n=2)
- [x] Extension complete (n=3,4,5)
- [ ] Paper submitted to Journal of Emerging Investigators

## Building On
[https://github.com/guptaarjun2027/spiral_fractals]
