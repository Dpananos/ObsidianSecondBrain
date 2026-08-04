
Consider a random variable, $\theta$, which has Student's T distirbution with location parameter 0 and scale parameter $\sigma$, namely

$$ \theta \sim t_{\nu}(0, \sigma) \>. $$


We can leverage the fact that the Student's T distribution is a mixture of a normal and a scaled  gamma distribution to approximate a the Student's T with a series of normal distributions.

Write

$$ \lambda \sim \operatorname{Gamma} \left( \dfrac{\nu}{2}, \dfrac{\nu}{2}  \right) \qquad \theta \mid \lambda \sim \operatorname{Normal} \left( 0, \dfrac{\sigma^2}{\lambda}\right)$$



We can discretize the mixture distribution in the following way

```r
df <- 3
scale <- 0.05
a <- 1e-4
rng <- scale * qt(c(a, 1-a), df=df)
x <- seq(min(rng), max(rng), 0.01)


# Numer of components I want to use to approximate

K <- 1000
k <- 1:K
uk <- (k-0.5)/K
lambda_k <- qgamma(uk, shape=df/2, rate=df/2)
tau_k <- scale / sqrt(lambda_k)

mu_k <- 0


dens <- sapply(tau_k, \(z) dnorm(x, mu_k, z)) 
wt <- rep(1/K, K)


plot(x, log(dens %*% wt))
lines(x, log(dt(x/scale, df=df)/scale), col='red')
```


