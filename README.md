# Equality-Constrained Convex Optimization in MATLAB

This repository contains two small MATLAB exercises on equality-constrained optimization.

The first solves a quadratic program through its Karush–Kuhn–Tucker (KKT) system and compares the result with MATLAB's `quadprog`. The second applies an equality-constrained Newton method to a logarithmic objective.

## Files

- `constraint Qp opt/KKT_Solve.m` — solves the KKT equations by block elimination
- `constraint Qp opt/KKT_System.m` — compares the custom KKT solution with `quadprog`
- `constraint Qp opt/EQ_NM.m` — equality-constrained Newton iterations
- `constraint Qp opt/NQE.m` — runs and plots the logarithmic-objective example

## Running the examples

Open MATLAB and run either `KKT_System.m` or `NQE.m` from the `constraint Qp opt` directory.

`KKT_System.m` requires the Optimization Toolbox for `quadprog`. `NQE.m` also uses YALMIP's `sdpvar` and `optimize` functions for comparison.

## Purpose

This is educational code written to make the KKT conditions and constrained Newton steps easier to inspect. The examples use fixed, low-dimensional problems and are not presented as a general-purpose solver.
