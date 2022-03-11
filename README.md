# Pbine

## Introduction
This repository provides a tool for a meta-analysis based on a combination of p-values from multiple studies. In our method, p-values are combined based on the test statistic: Fisher's p-value combination (Fisher, 1970) or weighted Fisher's p-value combination (Hou, 2005). The method calculates statistical significance of the test statistic by using a numerical integration procedure (Lin, Liang, and Yang, 2022). To account for the correlation of p-values, a correlation matrix of p-values (sigma) should be provided or estimated based on the paired of observed p-value sequences. In addition to the proposed numerical integration method, the conventional Fisher’s method (Fisher, 1970) that it assumes all p-values are independent and the decorrelation method (Zaykin, 2002) are also provided.

## Input Argument
 - **p** \: A matrix (n by m) of p-value. This matrix is composed by p-value sequences with a size of n in m studies. In each of the n rows, m p-values (i.e., p-values from m studies) will be combined. 
 - **method** \: Three methods of p-value combination are provided. "Int": the proposed numerical integration method; "Fisher": Fisher's method that it assumes all p-values are independent (Fisher, 1970); "Decor": Decorrelation method (Zaykin, 2002).
 - **sigma** \: A user-specified correlation matrix (m by m) of the combined p-values. This correlation matrix is used in the "Int" and "Decor" methods. If n = 1, sigma must be provided. If n > 1 and sigma is not specified by a user, then sigma will be estimated based on the pairs of p-value sequence with a size of n.

## Output Value
 - **p** \: A vector of output statistical significance (n by 1) of p-value combination. Each cell in this vector indicates the output statistical significance for a combination of p-values from m studies.

## Example 
### 1. In the case of n > 1.
set.seed(123) <br />
p1 = runif(20) <br />
p2 = (p1 + runif(20))/2 \# p2 is correlated with p1. <br />
pm = cbind(p1, p2) <br />
Pbine(pm) \# sigma is estimated by cor(p1, p2). <br />
\# 0.42049297 0.78909252 0.47575662 0.92912077 0.88866757 0.12490452 0.54885020 0.83891085 0.49508593 0.37767268  <br />
\# 0.97024711 0.57237059 0.70451327 0.64749652 0.07448837 0.81221096 0.35668056 0.06729487 0.32888367 0.77738220


### 2. In the case of n = 1, sigma is needed to be provided..
r12 = 0.2 <br />
r13 = 0.3 <br />
r23 = 0.6 <br />
COV = matrix(c(1, r12, r13, r12, 1, r23, r13, r23, 1), 3, 3) <br />
p1 = p2 = p3 = 0.05 <br />
pm = t(c(p1, p2, p3)) <br />
Pbine(pm, sigma = COV, method = "Int") \# If n = 1, sigma is needed to be provided. <br />
\# 0.02151095 <br />

### 3. Comparison on higher dimension (m) with Fisher's method.
r12 = 0.2 <br />
r13 = 0.3 <br />
r14 = 0.4 <br />
r23 = 0.6 <br />
r24 = 0.2 <br />
r34 = 0.1 <br />
COV = matrix(c(1, r12, r13, r14, r12, 1, r23, r24, r13, r23, 1, r34, r14, r24, r34, 1), 4, 4) <br />
p1 = p2 = p3 = p4 = 0.05 <br />
pm = t(c(p1, p2, p3, p4)) <br />
\# Pbine(pm, sigma = COV) <br />
\# 0.01404982 <br />
\# Pbine(pm, method = "Fisher") <br />
\# 0.0023 <br />
\# NOTE: 0.01404982 < 0.0023 means Fisher's method overestimates the significance (i.e. there's an inflated type I error) when p-values are correlated.

## Reference
 - Fisher, R.A., Statistical methods for research workers. 4th ed. Edinburgh: Oliver and Boyd, 1932.
 - Hou, C.-D., *A simple approximation for the distribution of the weighted combination of non-independent or independent probabilities.* Statistics & Probability Letters, 2005. 73(2): p. 179-187.
 - Zaykin, D.V., et al., *Truncated product method for combining P-values.* Genet Epidemiol, 2002. 22(2): p. 170-85.



