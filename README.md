# Pbine

## Introduction
This repository provides a tools for p-values (dependent) combination. We proposed a numerical integration method to combine correlated p-values (especially two p-values). Given a correlation matrix of correlated p-values 
![](http://www.sciweavers.org/upload/Tex2Img_1646821533/render.png)
, our method calculate p-value of Fisher's (or weighted Fisher's) statistics by using their exact joint probability density function under the transformation 
![](http://www.sciweavers.org/upload/Tex2Img_1646821672/render.png)
of independent p-values, where <img src="http://www.sciweavers.org/tex2img.php?eq=C&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt="C" width="17" height="15" /> is the Cholesky factor of
![1](http://www.sciweavers.org/upload/Tex2Img_1646818234/render.png);
![](http://www.sciweavers.org/upload/Tex2Img_1646821467/render.png)
and 
![1](http://www.sciweavers.org/upload/Tex2Img_1646818234/render.png)
represent independent and dependent p-values respectively.


## Input Argument
 - **p** \: A matrix (n by m) of p-values. n row represent how many pairs of p-values to be combined and m column represent how many p-values (studies) to be combined. 
 - **method** \: Three method can be selected. "Int": our method proposed above, "Fisher": Fisher's method, and "Decor": Decorrelation method (proposed by Zaykin, 2002).
 - **sigma** \: Correlation matrix of p-values used in "Int" and "Decor" methods.

## Output Value
 - **p** \: a vector of combined p-values with length equals to n. 

## Example 
r12 = 0.2<br />
r13 = 0.3<br />
r23 = 0.6<br />
COV = matrix(c(1, r12, r13, r12, 1, r23, r13, r23, 1), 3, 3)<br />
p1 = p2 = p3 = c(0.04, 0.05)<br />
pm = cbind(p1, p2, p3)<br />
Pbine(pm, sigma = COV, method = "Int")<br />
\# 0.01518269 0.02151095
