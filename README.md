# Pbine

## Introduction
This repository provides a tool for a meta-analysis based on a combination of p-values from multiple studies. In our method, p-values are combined based on the test statistic: Fisher's p-value combination (Fisher, 1932) or weighted Fisher's p-value combination (Good, 1955). The method calculates statistical significance of the test statistic by using a numerical integration procedure (Lin, Liang, and Yang, 2022). To account for the correlation of p-values, a correlation matrix of p-values (sigma) should be provided or estimated based on the paired of observed p-value sequences. In addition to the proposed numerical integration method, the conventional Fisher’s method (Fisher, 1932) that it assumes all p-values are independent and the decorrelation method (Zaykin, 2002) are also provided.

## Input Argument
 - **p** \: A matrix (n by m) of p-value. This matrix is composed by p-value sequences with a size of n in m studies. In each of the n rows, m p-values (i.e., p-values from m studies) will be combined. 
 - **method** \: Three methods of p-value combination are provided. "Int": the proposed numerical integration method; "Fisher": Fisher's method that it assumes all p-values are independent (Fisher, 1932); "Decor": Decorrelation method (Zaykin, 2002).
 - **sigma** \: A user-specified correlation matrix (m by m) of the combined p-values. This correlation matrix is used in the "Int" and "Decor" methods. If n = 1, sigma must be provided. If n > 1 and sigma is not specified by a user, then sigma will be estimated based on the pairs of p-value sequence with a size of n.

## Output Value
 - **p** \: A vector of output statistical significance (n by 1) of p-value combination. Each cell in this vector indicates the output statistical significance for a combination of p-values from m studies.

## Example 
### 1. In the case of n > 1.
set.seed(123) <br />
p1 = runif(20) <br />
p2 = (p1 + runif(20))/2 \# p2 is correlated with p1. <br />
pm = cbind(p1, p2) <br />
pn = Pbine(pm) \# sigma is estimated by cor(p1, p2). <br />
print(cbind(p1, p2, pn)) <br />

### 2. In the case of n = 1, sigma is needed to be provided..
r12 = 0.2 <br />
r13 = 0.3 <br />
r23 = 0.6 <br />
COV = matrix(c(1, r12, r13, r12, 1, r23, r13, r23, 1), 3, 3) <br />
p1 = p2 = p3 = 0.05 <br />
pm = t(c(p1, p2, p3)) <br />
pn = Pbine(pm, sigma = COV, method = "Int") \# If n = 1, sigma is needed to be provided. <br />
print(c(p1, p2, pn))

### 3. The case of m > 2.
r12 = 0.2 <br />
r13 = 0.3 <br />
r14 = 0.4 <br />
r23 = 0.6 <br />
r24 = 0.2 <br />
r34 = 0.1 <br />
COV = matrix(c(1, r12, r13, r14, r12, 1, r23, r24, r13, r23, 1, r34, r14, r24, r34, 1), 4, 4) <br />
p1 = p2 = p3 = p4 = 0.05 <br />
pm = t(c(p1, p2, p3, p4)) <br />
Pbine(pm, sigma = COV) <br />
\# 0.01404982 <br />
Pbine(pm, method = "Fisher") <br />
\# 0.0023 <br />
\# NOTE 1: In this example, the Fisher's method overestimates the significance because p-values are correlated. <br />
\# NOTE 2: Please wait a while because it needs more time when m > 3.

## Reference
 - Fisher, R.A., Statistical methods for research workers. 4th ed. Edinburgh: Oliver and Boyd, 1932.
 - Good, I., On the weighted combination of significance tests. Journal of the Royal Statistical Society: Series B (Methodological), 1955. 17: p. 264–265.
 - Zaykin, D.V., et al., *Truncated product method for combining P-values.* Genet Epidemiol, 2002. 22(2): p. 170-85.



