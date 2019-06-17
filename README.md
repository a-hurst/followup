followup: Follow-up Analyses for ANOVAs
---

**followup** is an R package designed for use with [afex](https://github.com/singmann/afex) that makes some common follow-up analyses for ANOVAs easy to do, including linear contrasts, polynomial trend analyses, and simple effects analysis. Contrasts and trends are powered internally by the `contrast()` function from the [emmeans](https://github.com/rvlenth/emmeans) package.

**Features include**:

- Computes effect sizes (partial or contrast eta-square, depending on the model) for trends and contrasts.
- Allows pooling of error for increased statistical power when splitting analyses on a between-subjects variable.
- Allows effortless splitting and pooling of error for contrasts and trends.
- Supports both ANCOVAs and regular ANOVAs.

**Disclaimer**: Although I have tested this package extensively with different types of datasets and models, it should still be considered beta software, and any results obtained should be treated with caution. In addition, since **followup** depends on a number of internal attributes from 'afex_aov' models in order to work correctly, any significant future internal changes to the afex package may break things.


## Installation

The easiest way to install **followup** is to use the 'devtools' package to install it directly from GitHub:

```r
require(devtools)
devtools::install_github('a-hurst/followup')
```
**followup** is not yet available on CRAN, so 'install.packages' unfortunately won't work.


## Functions

```r
# Performs linear contrasts 
getContrasts(mod, factor, contr, split = NULL, pooled = FALSE, adj = 'none')

# Performs polynomial trend analyses
getTrends(mod, factor, split = NULL, pooled = FALSE, adj = 'none')

# Performs simple effects analyses
simpleEffects(mod, split, pooled = FALSE, returns = 'summary', correction = 'GG')
```

You can read the documentation for all of these by accessing their documentation in R or RStudio (e.g. `?getContrasts`).


## Examples

First, we'll start off by creating an ANOVA model with afex:

```r
library(afex)
library(followup)

dat <- ToothGrowth
dat$dose <- as.factor(dat$dose)
dat$id <- as.factor(1:nrow(dat))

model <- aov_car(len ~ supp * dose + Error(id), data=dat)
nice(model)

# Anova Table (Type 3 tests)
# 
# Response: len
#      Effect    df   MSE         F ges p.value
# 1      supp 1, 54 13.19 15.57 *** .22   .0002
# 2      dose 2, 54 13.19 92.00 *** .77  <.0001
# 3 supp:dose 2, 54 13.19    4.11 * .13     .02
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘+’ 0.1 ‘ ’ 1
```

Now, let's do a follow-up simple effects analysis on our model, splitting on the 'supp' factor (tooth growth supplement, VC or OJ). To conserve statistical power, we'll also use pooled error and degrees of freedom from the omnibus model:

```r
simpleEffects(model, 'supp', pooled=TRUE)

# $VC
# Anova Table (Type 3 tests)
# 
# Response: len
#      num Df den Df    MSE      F     ges    Pr(>F)
# dose      2     54 13.187 62.541 0.83245 8.753e-15 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $OJ
# Anova Table (Type 3 tests)
# 
# Response: len
#      num Df den Df    MSE      F     ges    Pr(>F)
# dose      2     54 13.187 33.565 0.69961 3.363e-10 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

Next, let's do a linear contrast comparing tooth growth in the high dose condition to the two lower dose conditions:

```r
clist <- list(two.vs.others = c(-1/2, -1/2, 1))
getContrasts(model, 'dose', clist)

# NOTE: Results may be misleading due to involvement in interactions
#  contrast      estimate df      MSE        F     etaSq p.value
#  two.vs.others    10.93 54 13.18715 120.7892 0.6564634  <.0001
#
# Results are averaged over the levels of: supp
```

Finally, we'll do a polynomial trend analyses for levels of dosage, splitting by type of supplement:

```r
getTrends(model, 'dose', split='supp')

# supp = OJ:
#  contrast  estimate df      MSE            F        etaSq p.value
#  linear       12.83 27 12.29633 134.09916235 0.9996600967  <.0001
#  quadratic    -6.11 27 12.29633   0.04559625 0.0003399033  0.8325
# 
# supp = VC:
#  contrast  estimate df      MSE            F        etaSq p.value
#  linear       18.16 27 14.07796  58.46332329 0.9297157460  <.0001
#  quadratic     0.58 27 14.07796   4.41968535 0.0702842540  0.0450

```

I'll hopefully add a vignette or two for this package in the future, once I figure out how those work.
