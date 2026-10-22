# Dynamical low-rank approximations for the optimal control of matrix ODEs

## Abstract

We consider the open-loop optimal control of a dynamical system whose state is described by a time-dependent matrix $X(t)\in\mathbb{R}^{m\times n}$. Standard gradient-based optimization methods may incur substantial computational and memory costs, particularly for large-scale systems and fine temporal discretizations. To alleviate this computational burden, we develop a mathematical framework to use Dynamical Low Rank Approximations (DLRA) in optimal control problems. We study two distinct paradigms: the first one, called *Optimize-then-Project* (OtP), first derives the optimality system and then approximates both the state and adjoint equations using DLRA. The second one, called *Project-then-Optimize* (PtO), first replaces the state equation with its DLR approximation, and then derives the optimality conditions for the resulting optimization problem. The two approaches differ fundamentally concerning the role of the adjoint variable. In OtP, the adjoint variable evolves on the manifold of fixed rank matrices, whereas in PtO, the adjoint variable belongs to the tangent space to the fixed-rank manifold at the evolving state variable. This latter formulation leads to a covariant Ordinary Differential Equation (ODE) describing the dynamics of the adjoint variable as the underlying tangent space evolves in time.
 
For both approaches we establish rigorous error estimates for the approximated gradients, identifying how their accuracy depends on the low-rank approximations of initial and terminal conditions, the components of the state and adjoint vector fields normal to the relevant tangent spaces, the stability properties of the dynamics, and the manifold curvature.
The PtO approach further requires tailored time discretization strategies for the resulting adjoint covariant ODE, which we show to be robust against small singular values of the state variable.
To illustrate the proposed framework, we apply these methodologies to the optimal control of a parabolic partial differential equation with random coefficients and initial conditions.

## Authors

- [Tommaso Vanzan](https://sites.google.com/view/tommaso-vanzan/) (Politecnico di Torino)
- [Benjamin Carrel](https://bhgmc.com) (Independent researcher)

## Reference

To cite the paper, you can use the following Bibtex template:
```
@article{DLRA_OCP,
 title={Dynamical low-rank approximations for the optimal control of matrix ODEs}, 
 author={Benjamin Carrel and Tommaso Vanzan},
 year={2026},
 eprint={2610.29554},
 archivePrefix={arXiv},
 primaryClass={math.NA},
 url={https://arxiv.org/abs/2610.29554}, 
}

```

## Installation

First of all, clone this repository with

```
git clonehttps://github.com/vanzantom/DLRA-OCP
```

The code is written in MATLAB and there is nothing to compile or install. The experiments of the paper were run with **MATLAB R2023b**. Any release from R2020a onward should work, since that is the first release with `exportgraphics`, which the post-processing scripts use to export the figures.

No additional toolbox is needed to reproduce the experiments of the paper. The meshes used are already provided in `data_input/` (files `geo_<1/h>.mat` for $h = 1/8, \dots, 1/128$). The [Partial Differential Equation Toolbox](https://www.mathworks.com/products/pde.html) is only needed to generate a mesh for a new mesh size $h$ (through `initmesh`).

The finite element routines in `fem/` are based on the [GDGMatlab](https://github.com/vanzantom/GDGMatlab) library and are already included, so you do not need to download it separately.

The repository is organized as follows:

| Folder / file | Content |
| --- | --- |
| `GenerateData.m` | Runs all experiments of the paper in sequence |
| `AnalysisCovariantODE.m`, `ConvergenceGradients.m`, `GeometricAnalysis.m`, `ConvergenceOptimization.m`, `TimeOptimization.m` | One driver per experiment (see below) |
| `Setupfun.m` | Problem definition (model parameters, mesh, FEM matrices, target state) |
| `solvers/` | Time integrators: full-rank, BUG integrator for DLRA, projected Euler for the covariant adjoint ODE |
| `opt/` | Objective functional, gradient and gradient-descent routines |
| `lib/` | Utilities (projections, weighted QR/SVD, norms, curvature estimation, ...) |
| `fem/` | $\mathbb{P}^1$ finite element assembly |
| `data_input/` | Precomputed meshes; the FEM matrices are also written here at run time |
| `data_output/` | Raw results (`.mat`) of the experiments and the post-processing scripts |
| `fig_output/` | Figures (`.eps`) as they appear in the paper |


## How to reproduce the experiments

The raw data of all experiments in the paper is already provided in `data_output/`. You can therefore regenerate the figures directly, without rerunning any experiment.

### Regenerating the figures from the provided data

The post-processing scripts use relative paths, so run them from inside `data_output/`:

```matlab
cd data_output
AnalysisCovariantODEPostProcessing   % Figure 1
GradientsPostProcessing              % Figure 2
GeometricAnalysisPostProcessing      % Figure 3
OptimizationPostProcessing           % Figure 4
TimeOptimizationPostProcessing       % Figure 5
```

The figures are written to `fig_output/`.

### Rerunning the experiments

The drivers add the subfolders to the MATLAB path themselves, but they use relative paths, so run them from the **root of the repository**. Each driver saves its results as a `.mat` file in `data_output/`, overwriting the provided one. The random seed is fixed (`rng(0)`), so the results are reproducible.

To rerun everything at once:

```matlab
GenerateData
```

Or run each experiment separately:

| Figure | Experiment | Command | Output in `data_output/` |
| --- | --- | --- | --- |
| 1 | Time integration of the covariant ODE (projected Euler vs. integration of the factors) | `AnalysisCovariantODE` | `AnalysisCovariantODE.mat` |
| 2 | Convergence of the OtP and PtO gradients and adjoints w.r.t. $\Delta t$ | `ConvergenceGradients(1/32, 500, 1)` | `ConvergenceGradients_Nh_961_N_500.mat` |
| 3 | Effective curvature, $1/\sigma_r$ and normal components of the vector fields | `GeometricAnalysis(1/32, 500, 1, 1e-3)` | `EffectiveCurvature_Nh_961_N_500_Deltat_1.00e-03.mat` |
| 4 | Convergence of gradient descent: full-rank vs. OtP and PtO | `ConvergenceOptimization(1/32, 500, 1, 1e-3)` | `ConvergenceOptimization_Nh_961_N_500.mat` |
| 5 | Computational times vs. problem size | `TimeOptimization(2^-4, 2^7, 1, 1e-3)`<br>`TimeOptimization(2^-5, 2^8, 1, 1e-3)`<br>`TimeOptimization(2^-6, 2^9, 1, 1e-3)` | `TimeOptimization_Nh_225_N_128.mat`<br>`TimeOptimization_Nh_961_N_256.mat`<br>`TimeOptimization_Nh_3969_N_512.mat` |

The arguments of the PDE experiments are `(h, N, T, dt)`: mesh size, number of Monte Carlo samples, final time and time step. With mesh size $h$ the number of interior FEM nodes is $m = (1/h - 1)^2$, and $n = N$. After rerunning an experiment, regenerate the figure with the corresponding post-processing script as described above.

> **Warning: memory and computational time.** The full-rank reference solutions are stored as dense $m \times n$ matrices at every time step. For $\Delta t = 10^{-3}$, a single full-rank trajectory requires about 4 GB for $m = 961$, $n = 500$, and about 16 GB for $m = 3969$, $n = 512$. The full-rank state and adjoint are both needed. `ConvergenceOptimization` and, above all, the largest case of `TimeOptimization` are therefore very memory- and time-consuming. In practice we ran them on a cluster. Figures 1–3 can be reproduced on a standard workstation.

## License

This repository is copyright by Tommaso Vanzan and Benjamin Carrel, and
is distributed under the terms of the GNU General Public
License (GPL) version 3 (or later). See accompanying file `LICENSE` or
<http://www.gnu.org/licenses/gpl.html>.

In short, this means that everyone is free to use the source code, to modify it and
to redistribute it on a free basis, as long as the derivatives are also published under GPL3 license. It is not in the public domain; it is copyrighted and there are restrictions on its distribution (see the license). For example, you cannot integrate this version (in full
or in parts) in any closed-source software you plan to distribute
(commercially or not). Please contact us for more information.

### Third-party code

The file `lib/lgwt.m` (Legendre–Gauss nodes and weights) was written by Greg von Winckel and is distributed on [MATLAB Central File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/4540-legendre-gauss-quadrature-weights-and-nodes) under its own license, which is not covered by the GPL notice above. See the header of that file.

