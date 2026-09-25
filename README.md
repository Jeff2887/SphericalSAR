<!-- PROJECT LOGO -->
<br />
<div align="center">
  <a href="https://github.com/hanshang/CLR_vs_CDF_transformation](https://github.com/Jeff2887/SphericalSAR">    
    <img src="MQ.png" alt="Logo" width="250" height="250">
  </a>

<h3 align="center">Spherical spatial autoregressive model for spherically
embedded spatial data</h3>

</div>

<!-- ABOUT THE PROJECT -->
## Abstract
Spherically embedded spatial data are spatially indexed observations whose values naturally reside on or can be equivalently mapped to the unit sphere. Such data are increasingly ubiquitous in fields ranging from geochemistry to demography. However, analyzing such data presents unique difficulties due to the intrinsic non-Euclidean nature of the sphere, and rigorous methodologies for statistical modeling, inference, and uncertainty quantification remain limited. We introduce a unified framework to address these three limitations for spherically embedded spatial data. We propose a novel spherical spatial autoregressive model that leverages optimal transport geometry and then extend it to accommodate real-valued exogenous covariates. For either scenario with or without covariates, we establish the asymptotic properties of the estimators, derive a distribution-free Wald test for spatial dependence, complemented by a residual bootstrap procedure to enhance finite-sample performance, and develop approaches to forecast uncertainty quantification. The practical utility of these methodological advances is illustrated through extensive Monte-Carlo simulations, applications to Spanish geochemical compositions, and the Japanese age distribution of deaths.

### Code
1. simulation_basic_function.R --- basic functions for simulation
2. simulation.R --- run simulation studies
3. real_basic_function.R --- basic functions for real data
4. real_data_geochemical.R --- real data analysis for geochemical data in Spain
5. real_data_LTDC.R --- real data analysis for LTDC in Japan


## Contact
arXiv link: [https://arxiv.org/abs/2601.16385](https://arxiv.org/abs/2601.16385)

Jiazhen Xu - jiazhen.xu@mq.edu.au


