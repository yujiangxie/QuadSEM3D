# QuadSEM3D

QuadSEM3D is designed to compute 3D full Hessian kernels on the fly using spectral-element and adjoint methods.

Hessian kernels have two main applications:
1) In inversion, they can improve the convergence rate and mitigate trade-offs among inverted multi-parameters.
2) After inversion, they can be used to analyze the resolution of the inverted model.

QuadSEM3D is developed as an extension of the SPECFEM3D Cartesian code. We have integrated the QuadSEM3D implementation into the existing SPECFEM3D Cartesian framework, which is widely used by many users. 

While QuadSEM3D and SPECFEM3D Cartesian share the same code files, their coding structures differ. In QuadSEM3D, we assign two models to each GLL point of the mesher, resulting in the computation of two or four sets of associated wavefields for each GLL point. On the other hand, SPECFEM3D Cartesian is designed to compute one set of wavefields (forward simulation) or two sets of wavefields (adjoint simulation) for each GLL point.

Simultaneous computation of the Frechet and full Hessian kernels on the fly incurs approximately twice the computational cost compared to computing the Frechet kernel alone on the fly. We recommend using the full Hessian kernels for inverting models with many scatterers or small-scale heterogeneities, as they account for scattering fields and can improve the convergence rate of the inversion while mitigating trade-offs among multiple parameters. In our tests, we observed that the use of full Hessian kernels improves the convergence rate, effectively compensating for the higher computational cost per iteration compared to the L-BFGS method in Seisflows. This leads to lower computational costs for achieving the same level of misfit in the tested models. Additionally, the full Hessian kernel recovers structures at lower frequencies than what is required by the L-BFGS method at higher frequency stages, further reducing computational costs. In some cases, the full Hessian kernel inversion outperforms the L-BFGS inversion in terms of model recovery when compared to the true model, despite both showing similar misfit reduction numerically. It is important to note that this conclusion is specific to models with many scatterers or small-scale heterogeneities, and it highlights the uncertainty of inversions with the same misfit numerically.

Feel free to contact the author for access to the codes. The codes will be uploaded once the manual is ready.

References:
1. Yujiang Xie, Catherine A. Rychert and Nicholas Harmon. Elastic and anelastic adjoint tomography with Fréchet and full Hessian kernels, Geophysical Journal International, 2023, https://doi.org/10.1093/gji/ggad114
2. Yujiang Xie, Catherine A. Rychert, Nicholas Harmon, Qinya Liu and Dirk Gajewski. On‐the‐Fly Full Hessian Kernel Calculations Based upon Seismic‐Wave Simulations, Seismological Research Letters 2021,92,3832-3844, doi: https://doi.org/10.1785/0220200410
