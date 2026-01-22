[![DOI](https://zenodo.org/badge/1140012919.svg)](https://doi.org/10.5281/zenodo.18343838) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/license/mit)
# Multi-Sensor INLA-SPDE-Bayesian-Flood-Mapping

This repository implements a Bayesian spatial modeling framework for flood detection
based on **latent Gaussian fields** constructed via the SPDE approach and estimated
using Integrated Nested Laplace Approximation (INLA).

The project focuses on modeling spatial dependence explicitly, providing an interpretable
and computationally efficient alternative to purely data-driven deep learning approaches,
particularly in data-scarce or heterogeneous environments.


## Methodology

- **Bayesian inference** via Integrated Nested Laplace Approximation (INLA)
- **SPDE formulation** to approximate continuous Gaussian fields
- **Triangular mesh construction** to balance spatial resolution and computational cost

The latent field captures structured spatial variability not explained by covariates**
and should be interpreted as an inferential construct, not a directly observable
physical process.

---

## Requirements ![R](https://img.shields.io/badge/r-%23276DC3.svg?style=for-the-badge&logo=r&logoColor=white)
  - R ≥ 4.2 
  - Required R packages:
    - INLA
    - sf
    - sp
    - terra/raster
      
  ```
  install.packages(
  "INLA",
  repos = c(getOption("repos"),
            INLA = "https://inla.r-inla-download.org/R/stable")
)
```


## References
	-	Lindgren, F., Rue, H., & Lindström, J. (2011). An explicit link between Gaussian fields and Gaussian Markov random fields. Journal of the Royal Statistical Society: Series B
	-	Rue, H., Martino, S., & Chopin, N. (2009). Approximate Bayesian inference for latent Gaussian models using INLA. Journal of the Royal Statistical Society: Series B

## License

This project is released under the MIT License
