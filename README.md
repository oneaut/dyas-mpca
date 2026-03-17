# DyAS-MPCA: Stability-Guaranteed Dynamic Active Subspace Torque Allocation for Over-Actuated Hybrid Neuroprostheses

This repository contains the simulation code and parameters for:

> **M. Singh and N. Sharma**, "Stability-Guaranteed Dynamic Active Subspace Torque Allocation for Over-Actuated Hybrid Neuroprostheses," *IEEE Control Systems Letters*, 2025.

---

## Overview

Hybrid neuroprostheses (HNs) combine powered motors with functional electrical stimulation (FES) to assist or restore movement in patients with motor impairments. Because HNs have more actuators than joint degrees of freedom, an allocation problem must be solved in real time to distribute joint torques across the actuator set.

This work extends the **Dynamic Active Subspace (DyAS)** framework for model-predictive control allocation (MPCA) with three contributions:

1. **Analytical gradient** — replaces finite-difference rollout in online GCM estimation, reducing computation from *O(M_s n_w H d)* to *O(M_s n_w²)* with no plant evaluations.
2. **Approximation certificate** — a computable bound on the allocation cost loss from subspace truncation: δ_sub ≤ (1−η) C_J tr(**C**).
3. **Stability guarantees** — terminal-set invariance, Lyapunov decrease, ISS, and recursive feasibility for the reduced-order slack-QP allocator.

**System:** 3-DOF sit-to-stand (STS) hybrid neuroprosthesis, 8 actuators (3 motors + 5 FES channels)  
**Session:** 20 cycles × 7 s/cycle = 140 s total  
**Key result:** Active dimension locks at *l\** = 3 throughout; identical joint tracking to full-order allocation; up to 6.4 pp improvement in residual FES capacity.

---

## Repository Structure

```
dyas-mpca/
│
├── README.md
├── LICENSE
│
├── main/
│   ├── run_dyas_sim.m          % Top-level entry point — runs all three fatigue scenarios
│   ├── run_full_sim.m          % Full-order allocation baseline
│   └── compare_results.m       % Generate all manuscript figures and tables
│
├── core/
│   ├── dyas_update.m           % Algorithm 1: online DyAS identification (analytical gradient)
│   ├── slack_qp_allocate.m     % Inner-loop reduced slack-QP solver
│   ├── ctpd_outer_loop.m       % Outer computed-torque PD controller
│   ├── fatigue_update.m        % FES fatigue dynamics (Euler integration)
│   └── sts_dynamics.m          % 3-DOF STS plant dynamics (RK4)
│
├── params/
│   ├── sts_params.m            % All STS system parameters (see Table below)
│   ├── fatigue_scenarios.m     % Case A / B / C scenario definitions
│   └── dyas_params.m           % DyAS hyperparameters (M_s, η, μ, Δt_DyAS)
│
├── figures/
│   ├── plot_eigenvalue_evolution.m
│   ├── plot_fatigue_cycles.m
│   ├── plot_fatigue_timeseries.m
│   ├── plot_rmse_per_cycle.m
│   └── plot_endurance_bar.m
│
└── results/
    └── (auto-generated .mat files and figures saved here)
```

---

## Requirements

| Toolbox | Version | Purpose |
|---------|---------|---------|
| MATLAB | R2024a | Core simulation environment |
| Optimization Toolbox | any recent | `quadprog` (interior-point-convex) for allocation QP |

No other toolboxes are required. The simulation does not use CasADi, Simulink, or any third-party solvers.

---

## Quick Start

```matlab
% 1. Clone the repo and add all folders to your MATLAB path
addpath(genpath('dyas-mpca'))

% 2. Run the full 20-cycle, 3-scenario simulation
run_dyas_sim

% 3. Generate all figures from the manuscript
compare_results
```

Results are saved to `results/` as `.mat` files. Figures are saved as `.png` and `.pdf`.

To run a single scenario:

```matlab
% Available: 'A_mild', 'B_rapid', 'C_prefatigued'
run_dyas_sim('scenario', 'B_rapid')
```

To compare DyAS vs Full allocation for one scenario:

```matlab
run_full_sim('scenario', 'A_mild')
compare_results('scenario', 'A_mild')
```

---

## System Parameters

All parameters below correspond exactly to those used to produce the results in the paper.

### 3-DOF STS Neuroprosthesis

| Parameter | Symbol | Value |
|-----------|--------|-------|
| **Anthropometrics** (Ashby et al., 2004) | | |
| Segment masses (foot / shank / thigh / HAT) | m₁₋₄ | [1.2, 3.7, 8.6, 40.0] kg |
| Segment lengths | l₁₋₄ | [0.26, 0.43, 0.42, 0.55] m |
| **FES Channel Parameters** | | |
| Max isometric forces (TA / SOL / QUAD / HAMS / GLUT) | F_max | [800, 1200, 1500, 1000, 1400] N |
| Baseline fatigue rates w_f⁰ (TA / SOL / QUAD / HAMS / GLUT) | | [0.018, 0.015, 0.025, 0.020, 0.022] s⁻¹ |
| Baseline recovery rates w_r⁰ | | [0.009, 0.008, 0.013, 0.010, 0.011] s⁻¹ |
| Minimum active fraction | φ_min | 0.20 |
| **Motor and Actuator Bounds** | | |
| Peak motor torques — ankle / knee / hip | \|τ\|_max | [200, 250, 200] Nm |
| Motor activation bounds | [a_lb, a_ub] | [−1, 1] |
| FES activation bounds | [a_lb, a_ub] | [0, 1] |
| **Controller Gains** | | |
| Outer PD — proportional | K_p | diag(500, 450, 400) |
| Outer PD — derivative | K_d | diag(45, 42, 38) |
| Torque tracking weight | W_τ | diag(10, 10, 10) |
| Effort weight — motors / FES | R_u | diag(0.05, 0.05, 0.05, 2, 2, 2, 2, 2) |
| **DyAS Parameters** | | |
| Sample count | M_s | 20 |
| EMA weight | μ | 0.30 |
| Variance threshold | η | 0.95 |
| DyAS update period | Δt_DyAS | 0.10 s |
| Control timestep | Δt | 0.01 s (100 Hz) |
| Session length | | 140 s (20 cycles × 7 s) |

### Fatigue Scenarios

| Scenario | φ₀ | w_f / w_f⁰ | w_r / w_r⁰ | Clinical Description |
|----------|----|-----------|-----------|---------------------|
| A — Mild gradual | 1.00 | 0.7 | 1.3 | Maintenance therapy, rested patient |
| B — Rapid exhaustion | 1.00 | 2.5 | 0.4 | Intensive therapy, near-depletion |
| C — Pre-fatigued start | 0.45 | 1.5 | 0.8 | Patient arrives with residual fatigue |

---

## Algorithm: Online DyAS (Analytical Gradient)

The core computational improvement over prior work ([Singh et al., 2024]((https://ieeexplore.ieee.org/abstract/document/10354450))) is replacing finite-difference rollout with an exact closed-form gradient.

**Prior method (rollout):**
- Required 2 × M_s × n_w × H plant integrations per GCM update
- For n_w = 8, M_s = 40, H = 3, d = 4: **7680 RK4 evaluations per update** (~50 ms)

**This work (analytical):**
- H = **B_φᵀ W_τ B_φ + R_u** computed once per update
- r = **B_φᵀ W_τ τ_des + R_u a_nom** computed once per update  
- All M_s gradients: **G = 2H U_s − 2r 1ᵀ** (one matrix multiply)
- **0 plant evaluations, <0.1 ms per update**

```matlab
% Pseudocode for the analytical gradient step (see core/dyas_update.m)
B_phi = b * Phi_k;                          % Effective torque map
H     = B_phi' * W_tau * B_phi + R_u;       % Cost Hessian (once per update)
r     = B_phi' * W_tau * tau_des + R_u * a_nom;  % Gradient offset (once per update)
U_s   = sample_uniform(a_lb, a_ub, n_w, M_s);   % Random samples
G     = 2 * H * U_s - 2 * r * ones(1, M_s);     % All gradients, no plant calls
C_hat = (1/M_s) * G * G';                        % Batch GCM estimate
```

---

## Results Summary

| Metric | Full (n_w = 8) | DyAS (l* = 3) |
|--------|----------------|---------------|
| Ankle RMSE (deg) | 0.014 | **0.014** |
| Knee RMSE (deg) | 0.043 | **0.043** |
| Hip RMSE (deg) | 0.048 | **0.048** |
| Mean QP solve (ms) | 2.8 ± 0.3 | **0.4 ± 0.05** |
| QP speedup | — | **~7×** |
| DyAS update (ms) | — | **< 0.1** |
| Control budget at 100 Hz | 28% | **< 5%** |
| FES capacity gain — Case A | — | **+6.4 pp** |
| FES capacity gain — Case B | — | **+2.1 pp** |
| FES capacity gain — Case C | — | **+2.5 pp** |

---

## Citation

If you use this code or build on this work, please cite:

```bibtex
@article{singh2025dyas,
  author  = {Singh, Mayank and Sharma, Nitin},
  title   = {Stability-Guaranteed Dynamic Active Subspace Torque Allocation
             for Over-Actuated Hybrid Neuroprostheses},
  journal = {IEEE Control Systems Letters},
  year    = {2025},
  volume  = {X},
  pages   = {XX--XX},
  doi     = {10.1109/LCSYS.XXXX}
}
```

If you also use the prior DyAS-MPCA framework, please additionally cite:

```bibtex
@article{singh2024dynamic,
  author  = {Singh, Mayank and Lambeth, Krysten and Iyer, Ashwin and Sharma, Nitin},
  title   = {Dynamic Active Subspaces for Model Predictive Allocation
             in Over-Actuated Systems},
  journal = {IEEE Control Systems Letters},
  volume  = {8},
  pages   = {1234--1239},
  year    = {2024},
  doi     = {10.1109/LCSYS.2023.3342094}
}
```

---

## License

This project is released under the [MIT License](LICENSE).

---

## Contact

**Mayank Singh** — ms7513@columbia.edu
Department of Mechanical Engineering  
Columbia University, Ney York, NY

**Nitin Sharma** — nsharm23@ncsu.edu
