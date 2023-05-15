# QuadSEM3D

QuadSEM3D: Construction of 3D full Hessian kernels on the fly based on spectral-element and adjoint method. 

Hessian kernels have two main applications: 
  1) in inversion, they can improve the convergence rate and reduce the problem of trade-offs among inverted multi-parameters; 
  2) after inversion, they can be used to analyze the resolution of the inverted model.

We have incorporated the QuadSEM3D implementation into the SPECFEM3D Cartesian, which many users already use. The QuadSEM3D implementation is used only when simultaneously computing the Frechet and full Hessian kernels on the fly is required, or when computing the full Hessian kernels is required. For other purposes that do not require the full Hessian kernels, we suggest using the SPECFEM3D Cartesian, although the QuadSEM3D can do the same things.

Although the QuadSEM3D and SPECFEM3D Cartesian share the same code files, their coding structures are different. We set two models assigned for each GLL point of the mesher for QuadSEM3D, where two or four sets of associated wavefields are computed for each GLL point, while the SPECFEM3D Cartesian is implemented to compute one set of wavefields (i.e., forward simulation) or two sets of wavefields (i.e., adjoint simulation) for each GLL point.

Simultaneously computing the Frechet and full Hessian kernels on the fly would take about a 2-fold computational cost when compared to the computation of the Frechet kernel on the fly alone. We suggest using the full Hessian kernels for inverting models with many scatterers or with small-scale heterogeneities within the model since the full Hessian kernels account for the scattering fields, which may be used to improve convergence rate of the inversion and mitigate inter-parameter trade-offs for multi-parameter inversion for such models. For the tested models, we observed that the use of full Hessian kernels improves the convergence rate, and this improvement effectively compensates for the higher computational cost per iteration compared to the L-BFGS in Seisflows. This results in smaller computational costs for the same level of misfit for the tested models. The full Hessian kernel also recovers structure at a lower frequency than that required for the L-BFGS at higher frequency stages for the tested models, which could further reduce the computational costs. In some cases, although the full Hessian kernel inversion and L-BFGS inversion show numerically the same misfit reduction, the full Hessian kernel inversion shows better inverting results when compared to the L-BFGS inversion, when both are compared to the true model. This conclusion may apply only to models with many scatterers or small-scale heterogeneities. This study also reflects the uncertainty of the inversions under numerically the same misfit.

Welcome to contact the author for the codes.

References:
1. Yujiang Xie, Catherine A. Rychert and Nicholas Harmon. Elastic and anelastic adjoint tomography with Fréchet and full Hessian kernels, Geophysical Journal International, 2023, https://doi.org/10.1093/gji/ggad114
2. Yujiang Xie, Catherine A. Rychert, Nicholas Harmon, Qinya Liu and Dirk Gajewski. On‐the‐Fly Full Hessian Kernel Calculations Based upon Seismic‐Wave Simulations, Seismological Research Letters 2021,92,3832-3844, doi: https://doi.org/10.1785/0220200410
