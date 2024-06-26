# synthdid: Synthetic Difference in Differences Estimation

This package implements the synthetic difference in difference estimator (SDID) for the average treatment effect in panel data,
as proposed in Arkhangelsky et al (2019). We observe matrices of outcomes Y and binary treatment indicators W
that we think of as satisfying Y<sub>ij</sub> = L<sub>ij</sub> + &tau;<sub>ij</sub> W<sub>ij</sub> + &epsilon;<sub>ij</sub>.
Here &tau;<sub>ij</sub> is the effect of treatment on the unit i at time j, and we estimate the average effect of
treatment when and where it happened: the average of &tau;<sub>ij</sub> over the observations with W<sub>ij</sub>=1.
All treated units must begin treatment simultaneously, so W is a block matrix: W<sub>ij</sub> = 1 for i > N<sub>0</sub> and j > T<sub>0</sub>
and zero otherwise, with N<sub>0</sub> denoting the number of control units and T<sub>0</sub> the number of observation times
before onset of treatment. This applies, in particular, to the case of a single treated unit or treated period.

This package is currently in beta and the functionality and interface is subject to change.

Some helpful links for getting started:

- The [R package documentation](https://synth-inference.github.io/synthdid/) contains usage examples and method reference.
- The [online vignettes](https://synth-inference.github.io/synthdid/articles/more-plotting.html) contains a gallery of plot examples.
- For community questions and answers around usage, see [Github issues page](https://github.com/synth-inference/synthdid/issues).

### Installation

```R
devtools::install_github("susanathey/MCPanel")
devtools::install_github("apoorvalal/ebal")
devtools::install_github("apoorvalal/synthdid")
```

### Example

#### Omnibus function

The `panel_estimate` function accepts a data frame, unit_id, time_id, treatment, outcome name, and optionally an inference method (jackknife/bootstrap/placebo) and internally performs the reshaping to $N \times T$ matrices and calls several estimators to produce a large table with multiple estimates and corresponding standard errors.

![image](https://user-images.githubusercontent.com/12086926/232068341-f0bd39fe-26f4-4ebd-a948-4ad0572bf265.png)

#### SDID call

```R
library(synthdid)

# Estimate the effect of California Proposition 99 on cigarette consumption
data('california_prop99')
setup = panel.matrices(california_prop99)
tau.hat = synthdid_estimate(setup$Y, setup$N0, setup$T0)
se = sqrt(vcov(tau.hat, method='placebo'))
sprintf('point estimate: %1.2f', tau.hat)
sprintf('95%% CI (%1.2f, %1.2f)', tau.hat - 1.96 * se, tau.hat + 1.96 * se)
plot(tau.hat)
```



#### References
Dmitry Arkhangelsky, Susan Athey, David A. Hirshberg, Guido W. Imbens, and Stefan Wager.
<b>Synthetic Difference in Differences</b>, 2019.
[<a href="https://arxiv.org/abs/1812.09970">arxiv</a>]
