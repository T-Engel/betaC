
# betaC

The r package “betaC” provides code for the calculation of beta\_C, a
metric to quantify the non-random component in beta-diversity for a
given coverage.

## Installation

You can install the development version of betaC from
[GitHub](https://github.com/T-Engel/betaC) with:

``` r
devtools::install_github("T-Engel/betaC")
```

Please, also install the package “vegan” from CRAN. Furthermore,
“tidyverse” is recommended but the main functions will work without
it.

## Quick guide

If you are familiar with the concept of beta\_C and you just want to
have the code, here is all you need to know. If you’re unfamiliar with
the method you might want to read the paper and go through the basic
examples that follow below.

  - The main function of this package is `beta_C(x, C)` It’s first
    argument `x`takes a site-by-species abundance matrix as a matrix
    object or data frame (sites=rows, species= rows). The second
    argument `C` is target coverage used for standardization. The
    function returns beta\_C for the target coverage as a numeric value.
    Type `?beta_C` to see the documentation.

  - The other important function is `C_target(x)`. Its only argument `x`
    is a site-by-species abundance matrix. It returns the maximum
    possible coverage value that can be used to calculate beta\_C for
    the community matrix `x`. When comparing beta\_C across multiple
    communities it is recommended to choose the `C` argument in
    `beta_C(x, C)`such that it corresponds to the smallest `C_target()`
    output of all the communities. Type \`?C\_target to see the
    documentation.

## Example: Beta-diversity of the BCI dataset from vegan

Let’s look at the beta\_diversity of the `BCI` dataset that comes with
the package vegan. It’s a site by species abundance matrix with 50 plots
and 225 species from the Barro Colorado Island in Panama. If you want to
know more about the dataset type `?BCI` after loading the the vegan
package. We are intereted in the beta-diversity of the island because it
tells us something about the spatial structure species diversity. First,
let’s calculate Whittaker’s multiplicative beta-diverity as
\[\beta=\frac{\gamma}{\overline{\alpha}}\], where \(\gamma\) is the
gamma species richness (i.e. all plots combined) and
\(\overline{\alpha}\) is the alpha species richness (i.e. the average
plot richness).

``` r
library(betaC)
library(vegan)
#> Loading required package: permute
#> Loading required package: lattice
#> This is vegan 2.5-6
data(BCI)

# Whittakers multiplicative beta-diversity 
gamma=specnumber(colSums(BCI))
alpha= mean(specnumber(BCI))
beta_BCI=gamma/alpha
beta_BCI
#> [1] 2.478519
```

The samples have a beta-diversity of 2.48 The unit of this value is
“effective number of distinct sampling units”. This means that it
takes 2.48 hypothetical plots with complete species turnover to produce
the same beta-diversity that we observe in this data-set. Note that the
maximum possible value of beta is 50 here. In that case all 50 plots
would be completely unique in their species identities. The minimum
value of beta is 1 in which case all the plots share the same species.
However, this assumes that the sample size is big enough to randomly
sample all the species in the area - in other words that the samples are
complete. While this may be true on the gamma scale, where we have a
really large sample size, the alpha scale has a way smaller number of
individuals. Remember that the gamma scale has all the individuals of
the 50 subplots combined, while the average plot has only a 50th of this
sample size. Let’s have a look at this difference:

``` r
N_alpha= mean(rowSums(BCI))
N_gamma= sum(BCI)
barplot(c("gamma"=N_gamma,"alpha"= N_alpha), ylab = "Number of individuals", xlab = "Scale")
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="150%" />

To controll for this vast difference in sample size between the scales
we use individual based rarefaction (IBR). IBR lets you ask the question
“How many species can I expect to find if I sampled *n* Individuals
rather than the actual sample size of *N*?” This enables us to rescale
our gamma diversity estimate to the sample size observed at the alpha
scale. The relationship between the rarefied richness and the sample
size *n* used for rarefaction is called a individual based rarefaction
curve.

To illustrate this let’s have a look at the two-scale rarefaction curve.

``` r
library(tidyverse)

# alpha 
alphas_curves<-rarecurve(BCI)
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="150%" />

``` r
N_min= min(rowSums(BCI))
mean_alpha_curve= sapply(alphas_curves,function(x) return(as.numeric(x[1:N_min])) ) %>% rowMeans()

#gamma
gamma_curve= as.numeric(rarefy(colSums(BCI),sample = 1:N_gamma))

# plot them
plot(gamma_curve, xlab = "Number of individuals", ylab= "Rarefied richness", type = "l", lwd=2)
lines(mean_alpha_curve, col =3, lwd=2)
```

<img src="man/figures/README-unnamed-chunk-3-2.png" width="150%" />

## Functions

The main function `beta_C` uses the following helper functions to jump
between sample size and coverage. Both have been adapted from iNEXT.

  - `Chat`
  - `invChat`

Furthermore, there is a function that helps to determine the target
coverage value:

  - `C_target`

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/T-Engel/betaC.svg?branch=master)](https://travis-ci.org/T-Engel/betaC)
<!-- badges: end -->
