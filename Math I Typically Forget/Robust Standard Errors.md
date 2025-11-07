Consider the OLS problem

$$ 
\hat{\beta} = \underset{\beta \in \mathbb{R}^p}{\arg \min} \left \lVert y - X\beta \right \rVert \>.
$$

The solution is 

$$ \hat{\beta} = (X^\top X)^{-1} X^\top y $$
and the variance/covariance matrix for the coefficients is estimated by

$$ \widehat{\Sigma} = \dfrac{\varepsilon \varepsilon ^ 
\top}{n-p} (X^\top X)^{-1} \>. $$

In computing $\widehat{\Sigma}$, we make the assumption of homogeneity of variance. To see this, assume that the true model is

$$ y = X \beta  + \varepsilon  \>.$$
This implies

$$ \hat{\beta} =  (X^\top X)^{-1} X^\top (X\beta + \varepsilon) = \beta +  (X^\top X)^{-1} X^\top \varepsilon $$

We can now apply the variance operator to this quantity and leverage [[Statistics on Matrices and Vectors]] to obtain

$$ \operatorname{Var}(\hat{\beta} \mid X) = (X^\top X)^{-1} X^\top \operatorname{Var}(\varepsilon \mid X) \, X\, (X^\top X)^{-1} $$
Under homogeneity of variance $\operatorname{Var}(\varepsilon \mid X) = \sigma^2 I$ and hence we obtain the usual covariance formula.  To relax this assumption, we assume $\operatorname{Var}(\varepsilon \mid X) =  \Omega$ and hence obtain

$$ \operatorname{Var}(\hat{\beta} \mid X) = (X^\top X)^{-1} X^\top \Omega \, X\, (X^\top X)^{-1} \>.$$

To obtain EHW standard errors, use

$$ \widehat{\Omega} = \operatorname{diag}(\hat{\varepsilon}^2_1, \cdots, \hat{\varepsilon}^2_n) $$
where $\varepsilon_i = y_i - X \hat{\beta}$.

## A Useful Example

Consider a regression on a binary treatment variable $d_i$ so that $\sum_i d_i = n_t$ and $n - n_t = n_c$.

Then

$$ X^\top X = 
\begin{bmatrix} 
\sum_i 1 & \sum_i d_i \\
\sum_i d_i & \sum_i d_i
\end{bmatrix}=
\begin{bmatrix} 
n & n_t \\
n_t & n_t
\end{bmatrix} $$

$$ (X^\top X)^{-1} = \dfrac{1}{n_tn_c}
\begin{bmatrix} 
n_t & -n_t \\
-n_t & n
\end{bmatrix} $$

Note that $n \cdot n_t - n_t^2 = n_t(n-n_t)$.  The matrix $X^\top \widehat{\Omega} X$ is a weighted version of $X^\top X$.  This, along side the fact that $\sum \varepsilon^2_i = \displaystyle\sum_{i:d_i=1} \varepsilon^2_i + \displaystyle\sum_{i:d_i=0} \varepsilon^2_i$  yields

$$ X^T \widehat{\Omega}X = 
\begin{bmatrix} 
\sum_i \varepsilon^2_i & \sum_i \varepsilon^2_i d_i \\
\sum_i \varepsilon^2_i d_i & \sum_i \varepsilon^2_i d_i
\end{bmatrix}
=
\begin{bmatrix}
S_t^2 + S_c^2 & S_t^2 \\
S_t^2 & S_t^2
\end{bmatrix}
$$
The resulting covariance matrix has some nice cancelation happen so as to obtain

$$(X^\top X)^{-1} X^\top \widehat{\Omega} \, X\, (X^\top X)^{-1} =
\begin{bmatrix}
\dfrac{S_c^2}{n^2_c} & -\dfrac{S_c^2}{n^2_c} \\
-\dfrac{S_c^2}{n^2_c} & \dfrac{S_c^2}{n^2_c} + \dfrac{S_t^2}{n_t^2}
\end{bmatrix}$$


Note:
 1. The intercept and treatment effect are always negatively correlated.  Makes sense, larger the intercept the smaller the treatment effect is going to be.
 2. The standard error for the treatment effect is the same estimator we would get had we taken the textbook formula approach.  If we had not done this, then our covariance matrix would have been

$$ (X^\top X)^{-1} = \dfrac{s^2}{n_tn_c}
\begin{bmatrix} 
n_t & -n_t \\
-n_t & n
\end{bmatrix} $$

and our treatment effect standard error would be $\dfrac{n}{n_tn_c} \dfrac{S_t^2 + S_c^2}{n-2} \approx (S_t^2 + S_c^2) \left( \dfrac{1}{n_t} + \dfrac{1}{n_c} \right)$ which pools the variance. As a consequence, it is justifiable to always apply EHW to these types of regressions.