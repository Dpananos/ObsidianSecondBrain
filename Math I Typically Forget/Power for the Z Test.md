
## One Tailed Test

Suppose we collect a sample $X_1, \cdots, X_n$ and are willing to assume $E[X] = \mu_a$ and $\operatorname{Var}(X) \leq \infty$.  We want to test the null hypothesis

$$ \begin{align} H_0:& \mu \leq \mu_0 \\ H_a:& \mu>\mu_0 \end{align} $$
If the sample is sufficiently large we can leverage the CLT to justify a Z-test.  The probability we reject the null hypothesis (assuming the alternative is the case) is called the "Statistical Power", and is denoted

$$ 1 - \beta  = \Pr \left( \dfrac{\bar X - \mu_0}{\sigma/\sqrt{n}} \gt z_{1-\alpha} \Bigg \vert \mu = \mu_a\right) $$

We can solve for $n$, the required sample size to achieve the desired power.  First, note the test statistic is not a standard normal.  Let's do some algebra to turn this into a standard normal random variable.  Denote

$$ Z = \dfrac{\bar X - \mu_a}{\sigma / \sqrt{n}} \>. $$
If $\Delta = \mu_a - \mu_0$ then

$$ 1 - \beta = \Pr \left( Z > z_{1-\alpha} - \dfrac{\Delta}{\sigma / \sqrt{n}} \right) \>. $$
Now that we have expressed the probability in terms of a standard normal random variable, the above expression can be written as

$$ 1-\beta = 1 - \Phi\left( z_{1-\alpha} - \dfrac{\Delta}{\sigma / \sqrt{n}}\right) = \Phi \left(\dfrac{\Delta}{\sigma / \sqrt{n}} - z_{1-\alpha}\right) \>. $$

Here, $\Phi(x)$ is the CDF of a standard normal distribution, and we have used the identity $\Phi(x) = 1-\Phi(x)$ to obtain the last equality.  Finally, this implies

$$ z_{1-\beta} =  \dfrac{\Delta}{\sigma / \sqrt{n}} - z_{1-\alpha} $$
from which we can solve for $n$, yielding

$$ n = \dfrac{(z_{1-\beta} + z_{1-\alpha})^2 \sigma^2}{\Delta^2} $$
## Two Tailed Test

In the event we want the power for a two tailed test, the power is 

$$1-\beta = 1 - \Phi\left(z_{1-\alpha/2} - \frac{\Delta\sqrt{n}}{\sigma}\right) + \Phi\left(-z_{1-\alpha/2} - \frac{\Delta\sqrt{n}}{\sigma}\right)$$

Typically, one of these terms is going to be _very small_ and the sample size for the two tailed test is the same as the one tailed test substituting $1-\alpha/2$ where appropriate.

## Simulation

While these formulae are nice, it can be hard to obtain closed form expressions for most problems of interest, and so resulting to simulation is probably the best way to estimate power.  Shown below is an example in R


``` r
# Sampling dist of true mean
estimate_power <- function(n, mua=0.2, mu0=0, sigma=1, alpha=0.05, nsims=1e6){
  
  # True sampling dist
  sample_mean <- rnorm(as.integer(nsims), mua, sigma/sqrt(n))
  # Alternatively, sample a standard normal with this mean
  test_stat <- (sample_mean - mu0) / (sigma / sqrt(n))
  
  mean(abs(test_stat) > qnorm(1-alpha/2))
  
}


ns <- seq(25, 350, 25)

# From Simulation
(power_curve <- sapply(ns, estimate_power) |> round(digits=3))
#>  [1] 0.170 0.293 0.410 0.516 0.608 0.687 0.753 0.807 0.851 0.886 0.913 0.934
#> [13] 0.950 0.962

# Actual power
(actual_power <- sapply(ns, \(x) powertools::ztest.1samp(N=x, delta=0.2)) |> round(digits=3))
#>  [1] 0.170 0.293 0.410 0.516 0.609 0.688 0.754 0.807 0.851 0.885 0.913 0.934
#> [13] 0.950 0.963
```

