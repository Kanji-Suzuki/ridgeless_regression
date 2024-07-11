# ridgeless_regression

This repository analyzes the empirical performance of the ridgeless regression in the finance applicaiton based on 

Kelly, Bryan T. and Malamud, Semyon and Zhou, Kangying, The Virtue of Complexity in Return Prediction (December 13, 2021). Swiss Finance Institute Research Paper No. 21-90, Journal of Finance, forthcoming, Available at SSRN: https://ssrn.com/abstract=3984925 or http://dx.doi.org/10.2139/ssrn.3984925.

The out-of-sample performance of the timing strategy computed from the ridgeless regression is reported with different model complexity and ridge regularization.
The model complexity here is defined as the number of random features used to predict the return.

In addition to replicating their analysis, different weight distributions and activation functions are examined.
Moreover, I report the performance of the kernel regression with the Neural Network Gaussian Prosess corresponding to each case.
The main result is illustrated in `report/plot_stats.pdf`.

# Code Structure
- The analysis is conducted under `src`
- Functions are defined under `R`
- Reports are generated under `report`

