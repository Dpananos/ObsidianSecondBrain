## The Problem with Peeking

Suppose you're running an A/B test and you check the results every day. Each day, if the p-value drops below 0.05, you call it significant and stop. This feels reasonable — why wait if the answer is already clear? — but it inflates your false positive rate badly. Each additional look is another chance to get lucky, and over many looks you'll exceed your nominal $\alpha$ by a wide margin. This is the "peeking problem."

The naive fix is to commit to a single look at a pre-specified sample size. But this is operationally painful: you have to ignore all data until the experiment ends, and if you mis-specified the sample size you either under-power yourself or run the experiment longer than needed.

Sequential testing offers a principled middle ground. You can look as often as you like, at any time, and the false positive guarantee still holds — as long as you use the right threshold.

## The Always-Valid Confidence Interval

The key idea is to replace the fixed critical value $z_{1-\alpha/2}$ (typically 1.96) with a *data-dependent* bound $B(n)$ that grows slowly as sample size $n$ increases. The bound is chosen so that, no matter when you stop looking, the probability that zero is ever excluded from the interval when there is no true effect is at most $\alpha$.

The specific bound used here comes from a normal mixture martingale argument. For a given sample size $n$ and significance level $\alpha$, the bound is:

$$B(n) = \sqrt{\frac{\tau}{\,n\,}\left(\log\frac{\tau}{\rho} - 2\log\alpha\right)}$$

where $\rho$ is a tuning parameter (related to a lower bound on the sample sizes you expect to observe) and $\tau = n + \rho$. This bound is always larger than 1.96 — that's the price of the sequential guarantee. You need a wider interval at any given moment to compensate for the fact that you could have stopped at any earlier moment and gotten lucky.

Once you have $B(n)$, the always-valid confidence interval at time $n$ is just:

$$\hat{\mu} \pm B(n) \cdot \text{SE}(n)$$

where $\hat\mu$ is your current point estimate and $\text{SE}(n) = \sigma / \sqrt{n}$. Reject the null when zero falls outside this interval — and you're allowed to check this continuously.

## The Sequential P-Value

In a classical fixed-sample test, the p-value has a clean closed form. Your test statistic is $T = \hat\mu / \text{SE}$, and under the null $T \sim N(0,1)$, so:

$$p_{\text{fixed}} = 2\,\Phi(-|T|)$$

That's just the probability of seeing something at least as extreme as $T$ by chance. Simple.

In the sequential setting, the p-value is still well-defined, but it can't be read off a normal table. Instead, it's defined by *inverting the always-valid CI*: the p-value is the smallest $\alpha$ at which the CI would exclude zero. In other words, it's the $\alpha_0$ that solves:

$$B(n,\, \alpha_0) = |T|$$

Intuitively, $\alpha_0$ is the "just barely significant" level — tighten the CI any further and zero falls back inside. There's no closed form for this because $B$ depends on $\alpha$ in a nonlinear way (through $\rho$), so you find it numerically with a root-finder.

The crucial point is that this p-value is *always valid*: if you compute it at every look and reject the first time it drops below $\alpha$, your overall type I error rate is still bounded by $\alpha$. You don't need to correct for multiple looks. A classical p-value computed at interim peeks doesn't have this property.

One way to build intuition: a fixed-sample p-value is the tail probability under the null. A sequential p-value is more like asking, "at what significance level would my current CI boundary happen to land exactly on zero?" Both answer "how surprising is this result?", but the sequential version accounts for the fact that you could have looked earlier and gotten a different answer.

## Power and the Minimum Detectable Effect

Now suppose you want to know: at a fixed sample size $n$, how powerful is this sequential test? That is, if the true lift is $\delta$, how often will the always-valid interval exclude zero?

The estimate $\hat\mu$ is approximately $N(\delta, \text{SE}^2)$, so the test statistic $\hat\mu / \text{SE}$ is $N(\delta/\text{SE},\, 1)$. We reject when this exceeds $B$ in absolute value, so:

$$\text{Power}(\delta) = \Phi\!\left(\frac{\delta}{\text{SE}} - B\right) + \Phi\!\left(-\frac{\delta}{\text{SE}} - B\right)$$

The second term is negligible for any effect worth caring about, so power is essentially determined by how far $\delta/\text{SE}$ is above $B$.

Setting $\text{Power} = 1 - \beta$ (say 80%) and solving for $\delta$ gives the minimum detectable effect (MDE):

$$\text{MDE} = (B + z_\beta) \cdot \text{SE}$$

where $z_\beta = \Phi^{-1}(1-\beta)$, e.g. $z_{0.80} \approx 0.842$.

This is exactly analogous to the fixed-sample MDE:

$$\text{MDE}_{\text{fixed}} = (z_{1-\alpha/2} + z_\beta) \cdot \text{SE}$$

The only change is that $z_{1-\alpha/2}$ (1.96) is replaced by the sequential bound $B(n)$, which is larger. So the sequential test needs a bigger effect to achieve the same power at the same $n$, which makes sense — you're buying flexibility to peek at the cost of some sensitivity.

## A Common Mistake

There is a tempting but wrong shortcut for computing the sequential MDE: just take the fixed-sample MDE and scale it by $B/z$:

$$\text{MDE}_{\text{wrong}} = \underbrace{(z + z_\beta) \cdot \text{SE}}_{\text{fixed MDE}} \times \frac{B}{z} = \left(B + z_\beta \cdot \frac{B}{z}\right) \cdot \text{SE}$$

Comparing to the correct formula:

| Term | Correct | Wrong |
|---|---|---|
| Threshold part | $B \cdot \text{SE}$ | $B \cdot \text{SE}$ |
| Power part | $z_\beta \cdot \text{SE}$ | $z_\beta \cdot (B/z) \cdot \text{SE}$ |

The threshold term is right, but the power term is scaled up by $B/z > 1$. The result is an MDE that is too large — the test is actually more sensitive than advertised. If you run a simulation at the stated MDE, you'll observe power closer to 95% than 80%, because the MDE you computed was conservative. The correct multiplier to go from a fixed-sample MDE to a sequential MDE is $(B + z_\beta)\,/\,(z + z_\beta)$, not $B/z$.

## Verifying it in R

The script below implements the sequential bound, the always-valid p-value (via numerical inversion), simulates the power curve directly, and overlays the analytical formula. The simulated points and the analytical line should agree closely — if they don't, something is wrong with one or the other.

### Algorithm

Before reading the code, here is exactly how the interval $\hat\mu \pm B(n)\cdot\text{SE}$ is computed, step by step. The one trick worth flagging up front: the code never carries $\alpha$ around directly. It carries the quantity

$$\lambda := -2\log\alpha$$

which the source calls `negative_2_log_alpha`. Everything — the tuning parameter, the bound, the p-value root-find — is expressed in terms of $\lambda$. This is what "log space" in the function name refers to, and it is why the same bound function serves both the CI (where $\lambda$ comes from the confidence level) and the p-value (where $\lambda$ is the unknown being solved for).

**Inputs.** Point estimate $\hat\mu$, standard error $\text{SE} = \sigma/\sqrt{n}$, current sample size $n$, confidence level $1-\alpha$, and a planned sample-size scale $n_0$ (`sample_size_lower_bound`, default $10^4$).

1. **Convert the confidence level to log space.** Set $\alpha = 1 - \texttt{confidence\_level}$ and

   $$\lambda = -2\log\alpha.$$

   In `sequential_bound_multiplier` this is the literal expression `-2 * log(1 - confidence_level)`, passed straight into the bound function.

2. **Pick the tuning parameter $\rho$.** This is the only non-obvious line. $\rho$ is chosen by a closed-form approximation that makes the confidence sequence *tightest near the planned scale* $n_0$:

   $$\rho = \frac{n_0}{\lambda + \log(1+\lambda)} - 1.$$

   In code: `rho <- sample_size_lower_bound / (negative_2_log_alpha + log1p(negative_2_log_alpha)) - 1`, where `log1p(`$\lambda$`)` $= \log(1+\lambda)$. Intuitively, $\rho$ sets the "center of mass" of the normal mixture; tuning it to $n_0$ means you pay the smallest width penalty around the sample size you actually expect to be looking at.

3. **Form $\tau$.** 

   $$\tau = n + \rho \qquad(\texttt{tau <- observations + rho}).$$

4. **Evaluate the bound.** 

   $$B(n) = \sqrt{\frac{\tau}{n}\Big(\log\tfrac{\tau}{\rho} + \lambda\Big)} = \sqrt{\frac{\tau}{n}\Big(\log\tfrac{\tau}{\rho} - 2\log\alpha\Big)}.$$

   The second form is the $B(n)$ written earlier in this note — the code *adds* `negative_2_log_alpha` rather than *subtracting* $2\log\alpha$, which is the same number. This is `sequential_normal_mixture_bound_log_space`.

5. **Scale and center.** The half-width is $B(n)\cdot\text{SE}$, so

   $$\text{CI} = \hat\mu \pm B(n)\cdot\text{SE}.$$

   In `sequential_confidence_interval`: `bound <- multiplier * standard_error`, then `c(lower = point_estimate - bound, upper = point_estimate + bound)`. Reject the null when this interval excludes zero — i.e. `lower > 0 || upper < 0`.

**The p-value reuses steps 1–4.** Because $B$ is written as a function of $\lambda$, inverting $B(n,\alpha_0) = |T|$ for the test statistic $T = \hat\mu/\text{SE}$ is just a one-dimensional root find over $\alpha_0$. As $\alpha_0$ shrinks, $\lambda = -2\log\alpha_0$ grows and $B$ grows monotonically, so the root is unique. `sequential_p_value` calls `uniroot` on $\alpha_0 \in [10^{-12},\, 1-10^{-10}]$ and returns the crossing; if no crossing exists in range (the result is nowhere near significant) it falls back to $1.0$.

```r
DEFAULT_SAMPLE_SIZE_LOWER_BOUND <- 1e4
.MIN_P_VALUE <- 1e-12

sequential_normal_mixture_bound_log_space <- function(
    observations,
    negative_2_log_alpha,
    sample_size_lower_bound = DEFAULT_SAMPLE_SIZE_LOWER_BOUND
) {
  rho <- sample_size_lower_bound /
    (negative_2_log_alpha + log1p(negative_2_log_alpha)) - 1
  tau <- observations + rho
  sqrt((tau / observations) * (log(tau / rho) + negative_2_log_alpha))
}

sequential_bound_multiplier <- function(
    observations,
    confidence_level,
    sample_size_lower_bound = DEFAULT_SAMPLE_SIZE_LOWER_BOUND
) {
  sequential_normal_mixture_bound_log_space(
    observations,
    -2 * log(1 - confidence_level),
    sample_size_lower_bound
  )
}

sequential_confidence_interval <- function(
    point_estimate,
    standard_error,
    sample_size,
    confidence_level = 0.95,
    sample_size_lower_bound = DEFAULT_SAMPLE_SIZE_LOWER_BOUND
) {
  bound <- sequential_bound_multiplier(
    sample_size,
    confidence_level,
    sample_size_lower_bound
  ) * standard_error

  c(lower = point_estimate - bound, upper = point_estimate + bound)
}

# Always-valid p-value: smallest alpha at which the sequential CI excludes zero.
# Found by numerically inverting B(n, alpha) = |T|.
sequential_p_value <- function(
    point_estimate,
    standard_error,
    sample_size,
    sample_size_lower_bound = DEFAULT_SAMPLE_SIZE_LOWER_BOUND
) {
  t_stat <- abs(point_estimate / standard_error)

  f <- function(alpha) {
    sequential_normal_mixture_bound_log_space(
      sample_size,
      -2 * log(alpha),
      sample_size_lower_bound
    ) - t_stat
  }

  tryCatch(
    uniroot(f, interval = c(.MIN_P_VALUE, 1 - 1e-10))$root,
    error = function(e) 1.0
  )
}

simulate_sequential_detection_rate <- function(
    true_lift,
    n = 10000,
    sigma = 1,
    nsim = 1000,
    confidence_level = 0.95,
    sample_size_lower_bound = DEFAULT_SAMPLE_SIZE_LOWER_BOUND,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  se <- sigma / sqrt(n)

  detected <- replicate(nsim, {
    estimate <- rnorm(1, mean = true_lift, sd = se)

    ci <- sequential_confidence_interval(
      point_estimate = estimate,
      standard_error = se,
      sample_size = n,
      confidence_level = confidence_level,
      sample_size_lower_bound = sample_size_lower_bound
    )

    ci["lower"] > 0 || ci["upper"] < 0
  })

  mean(detected)
}

analytical_sequential_power <- function(
    true_lift,
    n = 10000,
    sigma = 1,
    confidence_level = 0.95,
    sample_size_lower_bound = DEFAULT_SAMPLE_SIZE_LOWER_BOUND
) {
  se <- sigma / sqrt(n)
  B <- sequential_bound_multiplier(n, confidence_level, sample_size_lower_bound)
  pnorm(true_lift / se - B) + pnorm(-true_lift / se - B)
}

true_lift <- seq(0.0, 0.1, 0.005)
pwr <- sapply(true_lift,
       \(x)
       simulate_sequential_detection_rate(
         true_lift = x,
         n = 10000,
         sigma = 1,
         nsim = 1000,
         confidence_level = 0.95,
         seed = 123
       ))

pwr_analytical <- sapply(true_lift, analytical_sequential_power)

plot(true_lift, pwr, ylim = c(0, 1), ylab = "Power", xlab = "True Lift")
lines(true_lift, pwr_analytical, col = "red")
legend("topleft", legend = c("Simulated", "Analytical"), col = c("black", "red"), pch = c(1, NA), lty = c(NA, 1))
```
