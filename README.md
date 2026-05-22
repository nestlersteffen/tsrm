# tsrm

Maximum Likelihood and Bayesian Estimation of the Parameters of the Triadic and the Standard Social Relations Model. Allows univariate and bivariate round-robin data and allows to examine the influence of person-level, dyad-level, and triad-level predictors (in case of the TSRM). Finally, the classic ANOVA estimator for the SRM and the TSRM is also implemented.

## Installation

``` r
# install.packages("devtools")
devtools::install_github("nestlersteffen/tsrm")
```

To update, simply rerun the installation command.

## Examples

A triadic social relations model is fitted with REML:

``` r
data(CurryEmerson)
colnames( CurryEmerson ) <- c("Group","Judge","Actor","Partner","y")
fit <- tsrm(y ~ 1, 
            p_var = c("Actor","Partner","Judge"), 
            g_var = "Group",
	    control=list(large=TRUE), # use faster sparse matrix operations
            data = Kenzer)
summary(fit)
```

A social relations model with an predictor of the perceiver effects is fitted with REML:

``` r
data(Kenzer)
fit <- srm(y ~ 1, 
            p_var = c("actor.id","partner.id"), 
            g_var = "rrgroup.id",
	    control=list(random_group=TRUE), # random effects for groups
            data = Kenzer)
summary(fit)
```

## Contributing

Issues and pull requests are not actively monitored. For questions, 
suggestions, or bug reports, please contact me directly via mail.

## Status

Work in progress.