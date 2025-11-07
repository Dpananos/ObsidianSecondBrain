Let $X$ and $Y$ be two random variables so that $E[X] = \mu_X$, $E[Y] = \mu_Y$, $\operatorname{Var}[X] = \sigma_X^2 \leq \infty$, $\operatorname{Var}[Y] = \sigma_Y^2 \leq \infty$ and $\operatorname{Cov}(X, Y) = \sigma_{XY}$, and consider the ratio of these two random variables

$$ R(X,Y) = \dfrac{X}{Y} $$
Then, the variance of $R$ can be approximated using the Delta Method

$$ 
\begin{align}
\operatorname{Var}(R(X, Y)) &\approx \dfrac{\mu_X^2}{\mu_Y^2} \left( \dfrac{\sigma^2_{X}}{\mu_X^2} + \dfrac{\sigma^2_{Y}}{\mu_Y^2} - 2 \dfrac{\sigma_{XY}}{\mu_X \mu_Y} \right) \\
&\approx \nabla R(\mu_X, \mu_Y) ^\top \> \Omega \> \nabla R(\mu_X, \mu_Y)\\
\end{align}$$
where $\nabla R(\mu_X, \mu_Y) = \begin{bmatrix} \dfrac{1}{\mu_Y}, \dfrac{-\mu_X}{\mu_Y^2}\end{bmatrix}^\top$ and 
$$\Omega = \begin{bmatrix} \sigma^2_X & \sigma_{XY} \\ \sigma_{XY} & \sigma^2_Y \end{bmatrix}$$
is the covariance matrix for the $X$ and $Y$. 

Ratio metrics are _wildly_ useful.
1.  For clustered experiments, the standard error obtained by the delta method is equivalent to a [[Cluster Robust Standard Errors]] in the case where one performs a regression on a treatment indicator and clusters the standard error [@dengEquivalenceDeltaMethod2021]
2. Percentiles can be estimated via a ratio metric, even in the case where clustering must be adjusted for [@dengApplyingDeltaMethod2018].
