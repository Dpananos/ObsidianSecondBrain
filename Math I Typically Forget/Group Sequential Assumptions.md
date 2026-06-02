The canonical group sequential assumptions state that the vector of test statistics across looks is multivariate normal. For any two looks $k$ and $j$ (where $k < j$), the joint distribution of $(Z_k, Z_j)^T$ is bivariate normal:

$$E[Z_k] = \theta \sqrt{\mathcal{I}_k}$$
$$\operatorname{Var}(Z_k) = 1 \quad \text{and} \quad \operatorname{Cov}(Z_k, Z_j) = \sqrt{\frac{\mathcal{I}_k}{\mathcal{I}_j}} = \sqrt{r}$$

Where $r = \frac{\mathcal{I}_k}{\mathcal{I}_j}$ is the information fraction between the two looks.

To simulate on the estimated effect scale ($\hat{\theta}$), we use the identity $\hat{\theta}_k = \frac{Z_k}{\sqrt{\mathcal{I}_k}}$. If we fix the initial sampling variance at the first look, $SE_1^2 = \frac{1}{\mathcal{I}_1}$, the subsequent variances scale dynamically with the information fraction:

$$\mathcal{I}_j = \frac{\mathcal{I}_k}{r} \implies SE_j^2 = r \cdot SE_k^2$$

Taking the square root reveals the true relationship:

$$SE_j = \sqrt{r} \cdot SE_k$$

Thus, the standard error at a later look shrinks precisely by a factor of $\sqrt{r}$, which is exactly equal to the covariance ($\operatorname{Cov}(Z_k, Z_j)$) between the standardized test statistics.

# Simulating in R

Here is a quick and dirty way to simulate the the $Z$ and the $\hat \theta$.

```r
library(MASS)

nsims <- 1000
theta <- rnorm(nsims, 0.0, 0.05)
se <- 0.05 # Initial standard error (at r = 0.2)

r <- c(0.2, 0.5, 1.0)

# 1. Covariance matrix (Correct)
Sigma <- outer(r, r, FUN = function(a, b) sqrt(pmin(a, b) / pmax(a, b)))

# 2. Maximum information scale based on initial look
I_max <- (1 / se^2) / min(r)   # (1 / 0.0025) / 0.2 = 2000
I <- I_max * r                # c(400, 1000, 2000)
sqrt_I <- sqrt(I)             # c(20, 31.62, 44.72)

# 3. Simulate Z-statistics with the correct sqrt(I) mean shift
z_null <- mvrnorm(nsims, rep(0, length(r)), Sigma)
mean_shift <- outer(theta, sqrt_I)

z <- z_null + mean_shift

# 4. Broadcast back to the estimator scale (theta_hat)
theta_hat <- sweep(z, MARGIN = 2, STATS = sqrt_I, FUN = "/")
```