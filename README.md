# ridgeless_regression

This repository analyzes the empirical performance of the ridgeless regression in the finance applicaiton based on Kelly, Bryan T., Semyon Malamud, and Kangying Zhou. "The virtue of complexity everywhere." Available at SSRN 4166368 (2022).
The performance of the timing strategy based on the ridgeless regression is reported with different model complexity and ridge regularization.
The model complexity here is defined as the number of random features used to predict the return.

In addition to replicating their analysis, different weight distributions and activation functions are examined.
Moreover, I report the performance of the timing strategy based on the kernel regression with the Neural Network Gaussian Prosess corresponding to each case.
The main result is illustrated in `report/plot_stats.pdf`.

# Code Structure
- The analysis is conducted under `src`
- Functions are defined under `R`
- Reports are generated under `report`

