Supplementing material for the paper

## "A Splitting Scheme for Dispersive Shallow Moment Equations"

Ullika Scholz, Department of Mathematics, RWTH Aachen University, Aachen, Germany, scholz@acom.rwth-aachen.de
Robin Paar, Research Training Group EDDy, RWTH Aachen University, Aachen, Germany
Manuel Torrilhon, Department of Mathematics, RWTH Aachen University, Aachen, Germany


#### Abstract: 

The well-known Shallow Water Equations (SWE) have a wide range of applications, from weather forecasting to debris flow. They can be used for modeling incompressible free-surface flows whenever the shallowness allows for a vertical-averaging; i.e., vertical effects are negligible in comparison to horizontal ones. One of their biggest advantages is the low computational cost of their numerical solution. But vertical averaging comes with the price of losing information along the vertical axis, e.g., the flows velocity profile. 
To overcome this limitation, the Shallow Moment Models were introduced. They contain information on the vertical direction despite being dimensionally reduced.
Shallow moment models incorporating a non-hydrostatic pressure have been introduced before as Dispersive Shallow Moment Models (DSM). However, no method for solving the non-stationary equations has been presented to this date. This is mainly because it was unclear how to compute the pressure equation in form of the divergence-free constraint.
We show how the pressure equations of the DSM models can be rewritten in the form of a Poisson-like problem to enable their solution with a projection-type splitting scheme. For the linear equations, we present the calculations for the generalized model and point out what are the difficulties in the non-linear case. We state the first two linear models and the corresponding nonlinear counterparts. Finally, we introduce a hybrid Finite-Volume Finite-Difference method and discuss the non-stationary numerical results for an experiment with periodic boundary and uneven bottom topography.

#### How to use (requires [Wolfram Mathematica](https://www.wolfram.com/mathematica/) or the [Free Wolfram Engine for Developers](https://www.wolfram.com/engine/))

The five .wls files are Wolfram Mathematica scripts that run the experiments from the numerical section of the paper. Each file represents a test case and can be executed individually. The .wl files are automatically included during execution of the .wls files and only serve as a library.
The code was written for the Mathematica 13.0 Kernel. For any further questions feel free to contact the authors.
