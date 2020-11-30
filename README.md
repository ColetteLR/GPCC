## GPCC

Algorithms were developed in MATLAB for the implementation of the Student's t Gaussian process conditional copula (GPCC).
The code is based on an R software package written by Hernandez-Lobato et al. 
(R code available at https://github.com/lopezpaz/gaussian\_process\_conditional\_copulas/tree/master/code)

## Application

Run the 'Experiment.m' file with input a text (.txt) file without variable headers. 
The dataset should contain the transformed data in the probability space.
In this example code, time is used as the conditioning variable (Z=(1:n)'), but other conditioning variables can be used instead if required.

## References

Hernández-Lobato, J. M., Lloyd, J. R., & Hernández-Lobato, D. (2013). Gaussian process conditional copulas with applications to financial time series. Advances in Neural Information Processing Systems, 1736–1744.

Minka, T. P. (2013). Expectation Propagation for approximate Bayesian inference. http://arxiv.org/abs/1301.2294

Naish-Guzman, A., & Holden, S. (2009). The generalized FITC approximation. Advances in Neural Information Processing Systems 20 - Proceedings of the 2007 Conference, Ivm, 1–8.
