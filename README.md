# VacancyAssistedChargeTransport.jl -- Numerical examples for analyzing the role of vacancy dynamics in innovative semiconductor devices

`VacancyAssistedChargeTransport.jl` includes all the required scripts and materials related to the PhD thesis
> *Modeling and simulation of vacancy-assisted charge transport in innovative semiconductor devices*

by Dilara Abdel, which was submitted to the Freie Universität Berlin.


The provided examples in `VacancyAssistedChargeTransport.jl` depend on the package [ChargeTransport.jl](https://github.com/PatricioFarrell/ChargeTransport.jl), which solves the drift-diffusion charge transport equations using the Voronoi finite volume method. This method is implemented through [VoronoiFVM.jl](https://github.com/j-fu/VoronoiFVM.jl).

The material is organized into four separate folders.

* **LargeTimeBehavior**: Here are example files that allow us to explore the large time behavior of the relative entropy with respect to the steady state and the quadratic $L^2$ errors for drift-diffusion models. These models describe charge transport in both perovskite solar cells (PSC) and TMDC-based memristive devices. The findings for the perovskite solar cells have already been published in [2].

* **VolumeExclusion**: Within this directory, we can discover files that comprehensively explore volume exclusion effects and their impact on perovskite charge transport modeling. Specifically, we compare a PSC charge transport model's internal states and current-voltage curves that rely on two different ionic current density descriptions. One of these descriptions assumes constant mobility, while the other considers a constant diffusion coefficient, with the remaining variable being density-dependent. Additional details can be found in [3].

* **TMDCDynamics**: In this directory, we find example files that allow us to validate the vacancy-assisted charge transport model against experimental hysteresis and pulse measurement data for lateral 2D $\text{MoS}_2$-based memristive devices. These files support the significance of vacancy dynamics in TMDC devices. For a comprehensive discussion of the simulation results, we refer to [4].

* **parameters**: This directory summarizes all the essential physical parameters required for performing the simulations.


----------------
The example files provided here are capable of reproducing the numerical results and findings presented in the following research papers:

[1] D. Abdel, P. Vágner, J. Fuhrmann and P. Farrell. [Modelling charge transport in perovskite solar cells: Potential-based and limiting ion depletion.](https://www.sciencedirect.com/science/article/abs/pii/S0013468621009865) Electrochimica Acta 390 (2021).

[2] D. Abdel, C. Chainais-Hillairet, P. Farrell and M. Herda. [Numerical analysis of a finite volume scheme for charge transport in perovskite solar cells.](https://doi.org/10.1093/imanum/drad034) IMA Journal of Numerical Analysis (2023).

[3] D. Abdel, N. E. Courtier and P. Farrell. [Volume exclusion effects in perovskite charge transport modeling.](https://doi.org/10.1007/s11082-023-05125-9) Optical and Quantum Electronics **55**, 884 (2023).

[4] B. Spetzler, D. Abdel, F. Schwierz, M. Ziegler and P. Farrell. [The Role of Vacancy Dynamics in Two-Dimensional Memristive Devices.](https://doi.org/10.48550/arXiv.2304.06527) Advanced Electronic Materials (accepted) (2023).

