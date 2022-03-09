# Pbine

## Introduction
This repository provides a tools for p-values (dependent) combination. We proposed a numerical integration method to combine correlated p-values (especially two p-values). Given a correlation matrix of correlated p-values <img src="http://www.sciweavers.org/tex2img.php?eq=%20%5CSigma%20&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt=" \Sigma " width="14" height="14" />, our method calculate p-value of Fisher's (or weighted Fisher's) statistics by using their exact joint probability function under the transformation <img src="http://www.sciweavers.org/tex2img.php?eq=T%3AR%3D1-%20%5CPhi%20%5C%7BC%20%5CPhi%20%5E%7B-1%7D%281-R%5E%7B%2A%7D%29%5C%7D&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt="T:R=1- \Phi \{C \Phi ^{-1}(1-R^{*})\}" width="239" height="21" /> of independent p-values, where <img src="http://www.sciweavers.org/tex2img.php?eq=C&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt="C" width="17" height="15" /> is the Cholesky factor of <img src="http://www.sciweavers.org/tex2img.php?eq=%20%5CSigma%20&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt=" \Sigma " width="14" height="14" /> and <img src="http://www.sciweavers.org/tex2img.php?eq=R%5E%7B%2A%7D&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt="R^{*}" width="25" height="17" /> and <img src="http://www.sciweavers.org/tex2img.php?eq=R&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt="R" width="17" height="15" /> reprents independent and dependent p-values respectively

## Input Argument
 - **p** \: A matrix (n by m) of p-values. n row reprent how many pairs of p-values to be combined and m column reprent how many p-values (studies) to be combined. 
 - **method** \: Three method can be selected. "Int": our method proposed above, "Fisher": Fisher's method, and "Decor": Decorrelation method (proposed by Zaykin, 2002)
 - **sigma** \: Correlation matrix of p-values used in "Int" and "Decor" methods.

## Output Value
 - **p** \: a vector of combined p-values with length equal to n. 

## Example 
r12 = 0.2<br />
r13 = 0.3<br />
r23 = 0.6<br />
COV = matrix(c(1, r12, r13, r12, 1, r23, r13, r23, 1), 3, 3)<br />
p1 = p2 = p3 = c(0.04, 0.05)<br />
pm = cbind(p1, p2, p3)<br />
Pbine(p1, p2, p3, sigma = COV)<br />
\# 0.01518269 0.02151095
