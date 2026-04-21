# Code Description

This code was developed as part of the Convex Optimization course (MATH 563) taught at McGill University. It implements 4 different convex optimization algorithms applied to the non-blind image denoising problem. The following algorithms are implemented:
- Primal Douglas-Rachford Splitting
- Primal Dual Douglas-Rachford Splitting
- Alternating Direction Method of Multipliers
- Chambolle–Pock

# Running the code

All algorithms can be run using the optsolve.m function. The file runme.m contains a "standard" implementation of our software package, which includes padding of the corrupted image to minimize aliasing errors. The function DefaultParams.m performs a default initialization of the hyperparameters for all the algorithms. These default hyperparameter values were found to yield good results during the validation of our software package.

Although not explicitly included in the project submission, all tests conducted using our code can be found in the following github repository [CONVEX-DENOISING](https://github.com/OlivierLefevre494/CONVEX-DENOISING).
