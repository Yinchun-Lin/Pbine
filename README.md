# Pbine
This repository provides a tools for p-values (dependent) combination.
There are two versions of Pbine. 

# 1. Combine p-values in two dimension 
This script allows you to combine two p-values by using 3 methods, Pbine, Fisher's, and decorrelation methods.
In addition, Pbine is allowed an advanced option to use. You can set different weight to each p-values. (More precisely, it changes the statistics from <img src="https://latex.codecogs.com/svg.image?p_{1}p_{2}" title="p_{1}p_{2}" /> to <img src="https://latex.codecogs.com/svg.image?p_{1}^{w}p_{2}^{(1-w)}" title="p_{1}^{w}p_{2}^{(1-w)}" />.)

# 2. Combine p-values in higher dimension
This script only allows you to combine n-dim p-values by using Pbine method.
(Note: when you apply on the case of dimension greater than four, it's will take more computation time.)
