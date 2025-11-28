Client asks for a pre/post analysis with differential followup -- the time between initial and final scores was not the same across patients.

A simple paired analysis will then distort the reported mean.  If time is measured in months, $m$, then the paired estimate is

$$ E[\hat \delta] = \sum_{m \in \mathcal{M}} E[\hat\delta \mid M=m] \times \Pr(M=m) $$
and those patients who go a  long time between measurements could potentially up weight their contribution to the paired difference estimate.  Additionally, it may be the case that those subjects who are particularly bad at baseline may be treated longer, thereby going even longer between measurements.

What we want to do is report an estimate of $\hat \delta$ conditional on treatment time but also adjusting for baseline.  If $y$ is the final measurement, $x$ is the initial measurement, $\delta=y-x$ is the change score, and $m$ is the time between measurements, then we want to fit a model

$$ \delta_i = \gamma_0 + \gamma_1 x_i + \gamma_2 m_i + u_i $$

This will allow us to present an estimate of the difference between final and initial measurements at any point in time, hence making the estimate a bit more "apples to apples".  Here is an example in R

``` r
library(tidyverse)
library(marginaleffects)


cleaned_data <- tibble::tribble(
  ~initial, ~final,          ~months,
  1.9,   1.32,               28,
  1.9,    1.6, 21.1612903225806,
  2.4,      2, 20.0333333333333,
  2.08,   1.46, 37.9032258064516,
  1.6,    1.6, 52.5161290322581,
  1.8,   1.33, 17.9642857142857,
  1.5,   1.51, 16.0357142857143,
  2.1,    1.9, 29.3548387096774,
  1.8,    1.6, 71.9354838709677,
  2.21,    2.2, 2.96774193548387,
  2.35,   1.54,             18.4,
  1.95,    1.6, 15.1290322580645,
  1.72,   1.51,               26,
  1.4,   1.24, 11.4666666666667,
  1.9,    1.4, 26.5806451612903,
  2.4,   1.83, 56.3666666666667,
  1.9,    1.4, 18.2258064516129,
  2,    1.4,  49.741935483871,
  1.6,    1.5, 19.1666666666667,
  2.3,    1.8, 22.0666666666667,
  1.4,    1.4, 21.4193548387097,
  2.7,   1.38, 62.9677419354839,
  2.8,    1.5, 72.8709677419355,
  1.55,   1.26, 22.5161290322581,
  1.85,   1.67, 42.0666666666667,
  1.8,    1.7,  31.741935483871,
  1.8,    1.6,  37.258064516129,
  1.8,    1.3, 31.9677419354839,
  2.1,   1.89, 35.8709677419355,
  1.9,    1.4, 44.6451612903226,
  1.58,    1.1,               10,
  1.61,    1.4, 12.5333333333333,
  1.6,   1.49, 16.0019305019305,
  1.31,   1.08, 37.3333333333333,
  1.8,    1.5,             20.8,
  1.8,    1.6, 32.6666666666667
) %>% 
  mutate(
    difference = final - initial
  )



fit <- lm(difference ~ initial + months, data = cleaned_data)

# Now, just average over predictions

avg_predictions(
  fit, 
  variables = list(months=c(12, 24, mean(cleaned_data$months), 36))
)
#> 
#>  months Estimate Std. Error      z Pr(>|z|)    S  2.5 % 97.5 %
#>    12.0   -0.305     0.0537  -5.69   <0.001 26.2 -0.410 -0.200
#>    24.0   -0.345     0.0373  -9.27   <0.001 65.5 -0.418 -0.272
#>    30.4   -0.367     0.0344 -10.66   <0.001 85.7 -0.434 -0.299
#>    36.0   -0.386     0.0366 -10.52   <0.001 83.6 -0.457 -0.314
#> 
#> Type: response
```

Note that the averaged prediction at the empirical mean for `months` matches the empirical mean for `differences`.

There is yet another way to do this with the regular ANCOVA model.  Here it is

```r
# Or, with ANCOVA
fit <- lm(final ~ initial + months, data = cleaned_data)

# 12 momnth estimate
J <- c(1, mean(cleaned_data$initial), 12)
V <- vcov(fit)

(point_estimate <- as.numeric(J %*% coef(fit) ) - mean(cleaned_data$initial))
#> [1] -0.3050546


se <- sqrt(as.numeric(J %*% V %*% J))

#Matches up to rounding
(confint <- point_estimate + c(-1, 1)*qnorm(1-0.025) *se)
#> [1] -0.4102181 -0.1998911
```