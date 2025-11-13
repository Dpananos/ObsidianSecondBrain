Suppose we randomly sample companies (clusters), and then sample employees within those companies.  Let $y$ be whatever outcome we are measuring.  If there are $g=1, \cdots, G$ clusters, in which there are $i=1, \cdots, m_g$ employees then we can model the outcome for individual $i$ in cluster $g$ as 

$$ y_{i, g} = \mu + a_g + \varepsilon_{i, g} $$

in which we assume $a_g$ are i.i.d with mean 0 and variance $\sigma^2_a$, $\varepsilon_{i, g}$ are i.id with mean 0 and variance $\sigma^2_e$ , and $a_g \perp \varepsilon_{i, g}$ (the two are statistically independent of one another).  By applying properties of the variance operator, it can then by shown that

$$ \operatorname{Var}(y_{i, g}) = \operatorname{Var}(a_g) + \operatorname{Var}(\varepsilon_{i, g}) = \sigma^2_a + \sigma^2_{\varepsilon} $$
owing to the fact that the two are independent  (hence their covariance is 0).  The covariance between to subjects from the same cluster is then

$$ \operatorname{Cov}(y_{i, g}, y_{j, g}) = \sigma^2_a $$
and the covariance of subjects across clusters is 0 (again, owing to independence between the $a$ and the $\varepsilon$)

$$ \operatorname{Cov}(y_{i, g}, y_{j, h}) = 0 $$
## Why Does The Sample Mean Have Larger Variance Under Clustered Sampling?

Consider the sample mean of a single cluster, $g$, which is

$$ \bar{y}_g = \dfrac{1}{m_g} \sum_{i=1}^{m_g} y_{i, g}$$
It can be shown that the variance of the sample mean is then

$$ \operatorname{Var}(\bar{y}_g) = \dfrac{1}{m^2_g} \left[ \sum_i \operatorname{Var}(y_{i, g}) + 2\sum_{l \lt k} \operatorname{Cov}(y_{l, g}, y_{k, g}) \right] $$

The second sum looks hairy, but remember the index is over subjects and we assume that the covariance between any two subjects is the same. Hence, this is the sum of the covariance a total of $\binom{m_g}{2}$ times (you can verify this by drawing it out.  There are $m_g-1$ factors, then $m_g-2$, then $m_g-3$ and so on).  After substitution, the variance is

$$  \operatorname{Var}(\bar{y}_g) = \dfrac{1}{m^2_g} \left[ m_g(\sigma^2_a + \sigma^2_\varepsilon) + 2 \binom{m_g}{2}\sigma^2_a \right] $$
The grand mean, over all clusters, is then

$$ \bar{y} = \dfrac{\sum_g m_g \bar{y}_g}{n} $$
where $n = \sum_g m_g$.  Hence, the sampling variance for the grand mean is

$$ \operatorname{Var}(\bar y) = \dfrac{1}{n^2} \sum_g\left[  m_g(\sigma^2_a + \sigma^2_\varepsilon) + 2 m_g \binom{m_g}{2}\sigma^2_a \right] $$
owing to the fact that sample means across companies are independent.  Further simplification is possible by noting that $\binom{m_g}{2} = m_g(m_g-1)/2$  and so

$$  \operatorname{Var}(\bar y) = \dfrac{1}{n^2} \sum_g\left[  m_g \sigma^2_\varepsilon + m_g^2\sigma^2_a \right] = \dfrac{\sigma^2_\varepsilon}{n} + \dfrac{\sigma^2_a}{n} \left[\sum_g\dfrac{m^2_g}{n}\right] $$

Hence, the clustered standard error is larger than the naive standard error whenever at least one cluster has multiple observations.  We need to account for this, else our confidence intervals will be wrong.  Using the naive approach which ignores the clustering would instead target the estimand

$$ V_{naive} =  \dfrac{\sigma^2_\varepsilon}{n} + \dfrac{\sigma^2_a}{n} $$

which is smaller than $\operatorname{Var}(\bar y)$ whenever $m_g>1$ for at least one $g \in G$.
## The Intraclass Correlation (ICC)

The relationship between the variance components is often summarized by the Intraclass Correlation Coefficient (ICC), denoted $\rho$. This measures the proportion of the total variance that is attributable to the between-cluster variance component $\sigma^2_a$:

$$\rho = \dfrac{\sigma^2_a}{\sigma^2_a + \sigma^2_\varepsilon}$$

The ICC $\rho$ is also the theoretical correlation between any two individuals $i$ and $j$ from within the same cluster $g$.

$$\operatorname{Corr}(y_{i,g}, y_{j,g}) = \dfrac{\operatorname{Cov}(y_{i,g}, y_{j,g})}{\sqrt{\operatorname{Var}(y_{i,g}) \operatorname{Var}(y_{j,g})}} = \dfrac{\sigma^2_a}{\sigma^2_a + \sigma^2_\varepsilon} = \rho$$

If $\rho=0$, then $\sigma^2_a = 0$, all clusters are identical, and observations are independent. In this case, the naive variance $V_{naive}$ is correct.

## The Design Effect

The "Design Effect" ($Deff$) quantifies the inflation of variance due to clustering. It is defined as the ratio of the true (clustered) variance to the naive variance:

$$Deff = \dfrac{\operatorname{Var}(\bar y)}{V_{naive}} = \dfrac{\dfrac{\sigma^2_\varepsilon}{n} + \sigma^2_a\sum_g\dfrac{m^2_g}{n^2}}{\dfrac{\sigma^2_\varepsilon + \sigma^2_a}{n}}$$

This ratio shows precisely how much the naive approach understates the true sampling variance.

### Balanced Case

In the simpler **balanced case**, where $m_g = m$ for all $G$ clusters, we have $n = G \cdot m$. The term in the true variance simplifies:

$$\sum_g\dfrac{m^2_g}{n^2} = \sum_{g=1}^G \dfrac{m^2}{(Gm)^2} = G \cdot \dfrac{m^2}{G^2 m^2} = \dfrac{1}{G}$$

The true variance for the balanced case is $\operatorname{Var}(\bar y) = \frac{\sigma^2_\varepsilon}{Gm} + \frac{\sigma^2_a}{G}$. By substituting the definitions of $\rho$, $\sigma^2_a$, and $\sigma^2_\varepsilon$, the Design Effect simplifies to the well-known formula:

$$Deff = 1 + \rho(m-1)$$

This shows that the variance inflation is a function of the within-cluster correlation $\rho$ and the cluster size $m$.

