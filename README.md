# R package ```ciperm```
Non-parametric confidence intervals using permutation tests

----
This package implements the methodology described in https://arxiv.org/pdf/2111.14966.pdf

The general univariate case (ie. a model and a test that satisifes the conditions of section 2.1 in the paper) can be handled by calling ```ciperm```/```ciperm0``` directly.
Similarly, the general multivariate case in handled by calling ```ciperm```/```ciperm0```. 


```ciperm0``` and ```ciperm0.multi```: The workhorses that perform the actual permutation scheme

```ciperm```: Computes the confidence interval 
(ie. by using quantiles from ```ciperm0```)

```ciperm.multi```: Computes the confidence interval. 
Optionally calculates the joint confidence level and adjusted confidence intervals.

```ciperm.twosample```: 
User-friendly function for computing the two-sample confidence interval. 

```ciperm.linreg```: 
User-friendly function for computing confidence interval for the slope in linear regression. 

```alpha.multi```: Computes the joint confidence level (from output of ```ciperm0.multi```)

```adjusted_ci```: Simple bisection algorithm that computes adjusted confidence intervals (from output of ```ciperm0.multi```)

# Installation 
The recommended way is to use the devtools package, e.g. run `devtools::install_github("naolsen/ciperm")` from the R interface.
As the package is only based on R code, no special compilers are needed. 

# Remarks
Not that the package has not been optimised for speed, but uses "crude" tools like ```uniroot```. 
Furthermore, following the remarks in the article, there is an O(2^K) cost of calculating the adjusted confidence level where K is the number of dimensions/coordinates.

Please fell free to contribute to the package. 
