
## I.I.D. Case

The $p^{th}$ percentile is the point $q$ such that $q = F^{-1}(p)$.  Here, $F(\cdot)$ is the CDF of some distribution and therefore $F^{-1}(\cdot)$ is the inverse CDF (or percentile function).

A percentile metric then seeks an estimate of the $p^{th}$ percentile, $\hat q$ , and the associated standard error (or confidence interval).  Given i.i.d. data $X_1 \leq X_2 \leq \cdots \leq X_n$ (sorted in ascending order so their index is also their rank), the $p^{th}$ percentile is estimated via $\hat q = X_{\lfloor np \rfloor}$.  The standard error for the percentile is provided by a Central Limit Theorem for percentiles, stating

$$\sqrt{n}\left\{X_{\lfloor n p\rfloor}-q\right\} \rightarrow N\left[0, \frac{p(1-p)}{f\left\{q\right\}^2}\right]$$
Note that the standard error involves an estimate of the probability density function evaluated at the true -- and unknown -- quantile.  Well...shit.  How can we get a standard error for the estimate $\hat q$ in light of this?

The trick is to compute $Y_i = I(X_i \leq X_{\lfloor np \rfloor})$ and note that $S = \sum_i Y_i$ is binomial 

$$ S \sim \operatorname{Binomial}(p;n) \>. $$
We can estimate $\hat p = S/n$ and construct a Wald confidence interval for $\hat p$

$$(L, U) =  \hat p \pm z_{1-\alpha/2} \dfrac{\sqrt{\hat p (1-\hat p)}}{n}$$

The points $(L, U)$ can be interpreted as ranks, which can be used to obtain a confidence interval for $\hat q$

$$ \Pr(q \in (X_L, X_U)) = 1-\alpha $$
We can verify this with simulation.

``` r
n <- 1000
p <- 0.75
a <- 0.05
z <- qnorm(1-a/2)
true_q <- qnorm(p)

r <- replicate(50000, {
  vals <- rnorm(n)
  X <- sort(vals)
  qhat <- X[floor(n*p)]
  Y <- as.numeric(X < qhat)
  phat <- mean(Y)
  s <- sqrt(phat*(1-phat)/n)
  l <- min(phat - z*s, 1)
  u <- max(0, phat + z*s)
  
  qs <- X[floor(n*c(l, u))]
  dplyr::between(true_q, min(qs), max(qs))
})

coverage <- scales::percent(mean(r))

glue::glue("Estimated Coverage for Percentile: {coverage}")
#> Estimated Coverage for Percentile: 95%
```

This assumes that the $X_i$ are all i.i.d.  This is not always the case (e.g. for example, in core web vitals, we may get several latency measures from the same user).  

## Clustered Case

Now, suppose the observations are clustered.  Let $c = 1, \cdots, K$ be the cluster index, and let $i = 1, \cdots, N_c$ index the number of observations within that cluster.  The observations are then $X_{c, i}$.  We can take a similar approach and let $Y_{c, i} = I(X_{c, i} \leq \hat q)$.  Now, let $S_c = \sum_i^{n_c} Y_{c, i}$ and note that 

$$\hat p = \dfrac{\sum_{c=1}^{K} \sum_{i=1}^{N_c} Y_{c, i}}{\sum_{c=1}^{K} N_c} = \dfrac{\sum_{c=1}^K S_c}{\sum_{c=1}^{K} N_c} = \dfrac{\bar S}{\bar N} \>.$$Hence in the clustered case, $\hat p$ is the ratio of two means.  We can use the delta method to obtain an estimate of the variance for this quantity.  The delta method will provide the sampling variance for the ratio of means (for more on how to do this, see [[Ratio Metrics]]).  Then, we can follow our procedure above and fetch the ranks which provide us a confidence interval for $\hat q$.  Here is the new procedure in R

``` r
n <- 1000
K <- 50            # Number of clusters
p <- 0.75
a <- 0.05
z <- qnorm(1-a/2)
mu_global <- 10
sd_cluster <- 2    
sd_error <- 1      
cluster_id <- sample(1:K, size=n, replace = T)

r <- replicate(50000, {
  
  u_c <- rnorm(K, mean = 0, sd = sd_cluster)
  obs_means <- mu_global + u_c[cluster_id]
  vals <- obs_means + rnorm(n, mean = 0, sd = sd_error)
  
  X <- sort(vals)
  q_hat <- X[floor(n * p)]
  
  Y <- as.numeric(vals<q_hat)
  S <- tapply(Y, cluster_id, sum)
  N <- tapply(Y, cluster_id, length)
  
  Sbar <- mean(S)
  Nbar <- mean(N)
  Svar <- var(S)/K
  Nvar <- var(N)/K
  covar <- cov(S, N)/K
  
  grad <- c(1/Nbar, -Sbar / (Nbar*Nbar))
  covmat <- matrix(c(Svar, covar, covar, Nvar), nrow = 2)
  delta_method <- grad %*% covmat %*% grad
  s <- sqrt(delta_method)
  
  l <- max(0, p -z * s)
  u <- min(1, p + z * s)
  
  qs <- X[floor(n*c(l, u))]
  
  # Check coverage against the TRUE Marginal Percentile
  # The marginal distribution is N(mu, sqrt(sd_cluster^2 + sd_error^2))
  true_q <- qnorm(p, mean = mu_global, sd = sqrt(sd_cluster^2 + sd_error^2))
  
  dplyr::between(true_q, min(qs), max(qs))
})

coverage <- scales::percent(mean(r), accuracy=0.001)

glue::glue("Estimated Coverage for Percentile: {coverage}")
#> Estimated Coverage for Percentile: 94.438%
```

## More Efficient Approach

Note how in either approach above we needed to make 3 passes over the data: one pass to get $X_{\lfloor np \rfloor}$ , one to get $X_L$, and a final pass to get $X_U$.  We can do this in fewer passes over the data, as Deng explains in [@dengApplyingDeltaMethod2018].  This requires adjusting the naive quantiles by some adjustment factor

``` r
n <- 10000
K <- 200            # Number of clusters
p <- 0.75
a <- 0.05
z <- qnorm(1-a/2)
mu_global <- 10
sd_cluster <- 2    
sd_error <- 1      
cluster_id <- sample(1:K, size=n, replace = T)

r <- replicate(50000, {
  
  u_c <- rnorm(K, mean = 0, sd = sd_cluster)
  obs_means <- mu_global + u_c[cluster_id]
  vals <- obs_means + rnorm(n, mean = 0, sd = sd_error)
  
  X <- sort(vals)
  s <- sqrt(p*(1-p)/n)
  l <- max(0, p - z*s)
  u <- min(1, p + z*s)
  
  # Get the quantiles in a single pass.
  qs <- X[floor(n * c(l, p, u))]
  q_hat <- qs[2]
  XL <- qs[1]
  XU <- qs[3]
  
  Y <- as.numeric(vals<=q_hat)
  S <- tapply(Y, cluster_id, sum)
  N <- tapply(Y, cluster_id, length)
  
  Sbar <- mean(S)
  Nbar <- mean(N)
  Svar <- var(S)/K
  Nvar <- var(N)/K
  covar <- cov(S, N)/K
  
  grad <- c(1/Nbar, -Sbar / (Nbar*Nbar))
  covmat <- matrix(c(Svar, covar, covar, Nvar), nrow = 2)
  delta_method <- grad %*% covmat %*% grad
  se <- sqrt(delta_method)
  
  correction_factor <- se/s
  
  corrected_qs <- c(
    q_hat - correction_factor*(q_hat - XL),
    q_hat + correction_factor*(XU - q_hat)
  )
  
  
  # Check coverage against the TRUE Marginal Percentile
  # The marginal distribution is N(mu, sqrt(sd_cluster^2 + sd_error^2))
  true_q <- qnorm(p, mean = mu_global, sd = sqrt(sd_cluster^2 + sd_error^2))
  
  dplyr::between(true_q, min(corrected_qs), max(corrected_qs))
})

coverage <- scales::percent(mean(r), accuracy=0.001)

glue::glue("Estimated Coverage for Percentile: {coverage}")
#> Estimated Coverage for Percentile: 94.122%
```

Here is the code Giorgio wrote, which returns the standard error instead of the CI directly.

```python
# /// script
# requires-python = ">=3.14"
# dependencies = [
#     "numpy",
#     "pandas",
#     "scipy",
#     "tqdm",
# ]
# ///
import numpy as np
import pandas as pd
from scipy.stats import norm
import tqdm

def percentile_delta_method(df, metric_col_name, cluster_col_name, p, alpha=0.05):
    """
    Estimate quantile and its confidence interval using the Delta method with clustering.

    This implements "outer confidence interval with post-adjustment" algorithm
    at the end of section 4.2 in Deng, Knoblich, Lu (2018), or equivalently
    Algorithm 1 in Yao, Li, and Lu (2014).
    
    Parameters:
    df: pd.DataFrame - Fact-level DataFrame with columns for cluster ID and the metric of interest
    metric_col_name: str - Name of the column containing the metric of interest (e.g., 'page_load_time')
    cluster_col_name: str - Name of the column containing the cluster ID (e.g., 'user_id')
    p: float - Desired percentile (between 0 and 1)
    alpha: float - Significance level for confidence interval (default 0.05)
    
    Returns:
    tuple - (percentile estimate, standard error)
    """
    n = len(df)
    z_critical = norm.ppf(1 - alpha/2)
    
    # Step 1: Compute lower and upper percentiles for confidence interval
    L_p = p - z_critical * np.sqrt(p*(1-p)/n)
    U_p = p + z_critical * np.sqrt(p*(1-p)/n)

    L = int(n * L_p) # convert to ranks
    U = int(n * U_p) + 1 # take the ceiling

    # Step 2: Fetch 3 percentiles. Here, exact rank because the data is small. In production,
    #        could be replaced with approximate percentiles (but keeping high precision).
    sorted_values = df[metric_col_name].sort_values()
    X_np = sorted_values.iloc[int(n*p)]
    X_L = sorted_values.iloc[L]
    X_U = sorted_values.iloc[U]
    
    # Step 3: Compute indicator variable (0/1 "is this observation below the target percentile?")
    df['I'] = (df[metric_col_name] <= X_np).astype(int)
    
    # Step 4: Compute variance of I using Delta method
    K = df[cluster_col_name].nunique()
    N = df.groupby(cluster_col_name).size()
    S = df.groupby(cluster_col_name)['I'].sum()
    
    S_bar = S.mean()
    N_bar = N.mean()
    
    var_S = S.var()
    var_N = N.var()
    cov_SN = np.cov(S, N)[0, 1]
    
    var_I = (1 / (K * N_bar**2)) * (var_S - 2*(S_bar/N_bar)*cov_SN + (S_bar/N_bar)**2 * var_N)
    
    # Step 5: Compute correction factor
    c = np.sqrt(var_I / (p * (1 - p) / n))
    
    # Step 6: Compute standard error
    se = c * (X_U - X_L) / (2 * z_critical)
    
    return X_np, se

def between(x:float, l:float, u:float):
    return (x<u)&(l<x)


if __name__ == "__main__":

    nsims = 50_000
    covered = 0
    
    # Parameters from your R Code
    n = 1000
    K = 50
    p = 0.75
    mu_global = 10
    sd_cluster = 2
    sd_error = 1
    
    # The true marginal distribution is N(mu, sqrt(sd_cluster^2 + sd_error^2))
    # True quantile calculation
    total_sd = np.sqrt(sd_cluster**2 + sd_error**2)
    true_ptile = norm.ppf(p, loc=mu_global, scale=total_sd)

    for i in tqdm.tqdm(range(nsims)):
        # Generate Cluster IDs (0 to K-1)
        cluster_ids = np.random.choice(np.arange(K), size=n, replace=True)
        
        # Generate Cluster Effects (random intercept)
        u_c = np.random.normal(0, sd_cluster, size=K)
        
        # Map cluster effects to observations
        # R equivalent: obs_means <- mu_global + u_c[cluster_id]
        obs_means = mu_global + u_c[cluster_ids]
        
        # Add observation error
        # R equivalent: vals <- obs_means + rnorm(n, mean = 0, sd = sd_error)
        vals = obs_means + np.random.normal(0, sd_error, size=n)
        
        df = pd.DataFrame({
            'user_id': cluster_ids,
            'outcome': vals
        })

        X_np, se = percentile_delta_method(df, metric_col_name='outcome', cluster_col_name='user_id', p=p, alpha=0.05)

        covered += between(true_ptile, X_np-2*se, X_np+2*se)

    print(f"Coverage: {covered/nsims:.2%}")
```

## The Density Estimate Revisited

Originally, I claimed that the asymptotic sampling distribution of the quantile satisfied

$$\sqrt{n}\left\{X_{\lfloor n p\rfloor}-q\right\} \rightarrow N\left[0, \frac{p(1-p)}{f\left\{q\right\}^2}\right]$$

Which should imply (in the clustered case)

$$ X_{\lfloor np \rfloor} \sim N \left[ q, \dfrac{\operatorname{Var}(I)}{nf\{q\}^2}\right] $$
The algorithm described in the steps above is an alternate way to estimate the sampling variability of the quantile, with sampling variability

$$ \sigma^2_{OCI} = \dfrac{(X_U-X_L)^2}{4z^2_{1-\alpha/2}} \dfrac{n \operatorname{Var}(I)}{p(1-p)} $$
in which $\operatorname{Var}(I) = \nabla \hat p^\top \Omega \nabla \hat p$ (i.e. the result from the delta method).  These two results should provide similar estimates, and so it makes sense to equate them

$$\frac{\operatorname{Var}(I)}{f(q)^2} = \frac{(X_U-X_L)^2}{4z^2} \cdot \frac{n \operatorname{Var}(I)}{p(1-p)}$$

solving for the density yields

$$f(q) = \frac{2z_{1-\alpha/2}}{X_U - X_L} \sqrt{\frac{p(1-p)}{n}}$$

which can be interpreted as a finite difference formula for the slope of the CDF evaluated at $q$.  Hence, we can use this term in place of the quantiles completely.

```python

# /// script
# requires-python = ">=3.14"
# dependencies = [
#     "numpy",
#     "pandas",
#     "scipy",
#     "tqdm",
# ]
# ///
import numpy as np
import pandas as pd
from scipy.stats import norm
import tqdm

def percentile_w_density(df, metric_col_name, cluster_col_name, p, alpha=0.05):
    """
    Estimate quantile and its confidence interval using the Delta method with clustering.

    This implements "outer confidence interval with post-adjustment" algorithm
    at the end of section 4.2 in Deng, Knoblich, Lu (2018), or equivalently
    Algorithm 1 in Yao, Li, and Lu (2014).
    
    Parameters:
    df: pd.DataFrame - Fact-level DataFrame with columns for cluster ID and the metric of interest
    metric_col_name: str - Name of the column containing the metric of interest (e.g., 'page_load_time')
    cluster_col_name: str - Name of the column containing the cluster ID (e.g., 'user_id')
    p: float - Desired percentile (between 0 and 1)
    alpha: float - Significance level for confidence interval (default 0.05)
    
    Returns:
    tuple - (percentile estimate, standard error)
    """
    n = len(df)
    z_critical = norm.ppf(1 - alpha/2)
    
    # Step 1: Compute lower and upper percentiles for confidence interval
    L_p = p - z_critical * np.sqrt(p*(1-p)/n)
    U_p = p + z_critical * np.sqrt(p*(1-p)/n)

    L = int(n * L_p) # convert to ranks
    U = int(n * U_p) + 1 # take the ceiling

    # Step 2: Fetch 3 percentiles. Here, exact rank because the data is small. In production,
    #         could be replaced with approximate percentiles (but keeping high precision).
    sorted_values = df[metric_col_name].sort_values()
    X_np = sorted_values.iloc[int(n*p)]
    X_L = sorted_values.iloc[L]
    X_U = sorted_values.iloc[U]
    
    # Step 3: Compute indicator variable (0/1 "is this observation below the target percentile?")
    df['I'] = (df[metric_col_name] <= X_np).astype(int)
    
    # Step 4: Compute variance of I using Delta method
    K = df[cluster_col_name].nunique()
    N = df.groupby(cluster_col_name).size()
    S = df.groupby(cluster_col_name)['I'].sum()
    
    S_bar = S.mean()
    N_bar = N.mean()
    
    var_S = S.var()
    var_N = N.var()
    cov_SN = np.cov(S, N)[0, 1]
    
    var_I = (1 / (K * N_bar**2)) * (var_S - 2*(S_bar/N_bar)*cov_SN + (S_bar/N_bar)**2 * var_N)

    # Skip the correction factor, just estimate the density directly using X_U and X_L
    # Step 6: Compute standard error
    density = 2*z_critical / (X_U - X_L) * np.sqrt(p*(1-p)/n)
    se = np.sqrt(var_I)/ (density)
    
    return X_np, se

def between(x:float, l, u):
    return (x<u)&(l<x)


if __name__ == "__main__":

    nsims = 50_000
    covered = 0
    true_ptile = norm().ppf(0.975)

    for i in tqdm.tqdm(range(nsims)):
        df = pd.DataFrame({
            'user_id' : np.random.choice(np.arange(100), size=1000),
            'outcome': np.random.normal(size=1000)
        })

        X_np, se = percentile_delta_method(df, metric_col_name='outcome', cluster_col_name='user_id', p=0.975, alpha=0.05)

        covered += between(true_ptile, X_np-2*se, X_np+2*se)

    print(f"Coverage: {covered/nsims:.2f}%")
```