QuadSEM3D is the 3D version of QuadSEM (https://github.com/yujiangxie/QuadSEM), which is based on the SPECFEM3D package. I constructed QuadSEM3D by modifying the src folder of SPECFEM3D. In other words, QuadSEM3D can solve four wave equations simultaneously, whereas SPECFEM3D can only solve one or two wave equations at a time.

One can use it in essentially the same way as SPECFEM3D; I have kept it as similar to SPECFEM3D as possible. The key difference is that QuadSEM3D can work with two models simultaneously to compute the full Hessian kernels on the fly (for elastic cases), which can substantially reduce storage and I/O costs and make the use of full Hessian kernels feasible for large-scale simulations and imaging.


References:

Yujiang Xie, Catherine A. Rychert and Nicholas Harmon. Elastic and anelastic adjoint tomography with and full Hessian kernels, Geophysical Journal International, 234, 1205-1235, 2023, https://doi.org/10.1093/gji/ggad114

Yujiang Xie, Catherine A. Rychert, Nicholas Harmon, Qinya Liu and Dirk Gajewski. On‐the‐Fly Full Hessian Kernel Calculations Based upon Seismic‐Wave Simulations, Seismological Research Letters 2021,92, 3832-3844, doi: https://doi.org/10.1785/0220200410



