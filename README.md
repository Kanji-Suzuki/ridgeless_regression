# ridgeless_regression

This repository analyzes the empirical performance of the ridgeless regression in the finance applicaiton based on 

KELLY, B., MALAMUD, S. and ZHOU, K. (2024), The Virtue of Complexity in Return Prediction. J Finance, 79: 459-503. https://doi.org/10.1111/jofi.13298

The out-of-sample performance of the timing strategy computed from the ridgeless regression is reported with different model complexity and ridge regularization.
The model complexity here is defined as the number of random features used to predict the return.

In addition to replicating their analysis, different weight distributions and activation functions are examined.
Moreover, I report the performance of the kernel regression with the Neural Network Gaussian Prosess corresponding to each case.
The main result is illustrated in `report/plot_stats.pdf`.

## Code Structure
Please download data from the paper and save it under `data/raw`.

- The analysis is conducted under `src/`
- Functions are defined under `R/`
- Reports are generated under `report/`
- Data is saved under `data/`
  

