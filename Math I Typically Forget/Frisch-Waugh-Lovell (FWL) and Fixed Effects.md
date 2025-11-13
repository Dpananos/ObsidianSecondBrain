## FWL Theorem

FWL states that the estimate of $\gamma$ from the model $Y=X\beta + D\gamma + u$ will be the same as the estimate obtained from $MY = MD\gamma + Mu$ , where

$$ M = (I-X(X^\top X)^{-1}X^\top) = I-H $$
Note that $H$ is the "hat" matrix because it "puts the hat on $Y$" in regression.  $H$ is also a projection matrix -- meaning that $\hat Y = HY$ is the projection of $Y$ onto the column space of $X$, $\mathcal{C} (X)$ .  Makes sense, since the prediction from the linear model is, literally,  a linear combination of the columns of $X$.

This also means $I-H$ is a projection matrix, but left multiplication by this matrix projects the vector onto $\mathcal{C}(X)^\perp$, the orthogonal compliment of the column space of $X$.  Note that $(I-H)X = \mathbf{0}$ for this reason.  We'll come back to this after a brief interlude of algebra.

## Why is $\gamma$ Left Unchanged?

It would first behoove us to understand what the estimate of $\gamma$ is under each model.  Under the second model

$$ \hat \gamma = ((MD)^\top MD)^{-1} (MD)^\top MY \>.$$
Since $M$ is a projection matrix it is idempotent and symmetric, and so this expression simplifies to

$$ \hat \gamma = (D^\top M D)^{-1} D^\top M Y \>.$$
Under the first model, the estimate of $\gamma$ is more involved.  First, concatenate $X$ and $D$ column-wise so that $Z = \begin{bmatrix} X & D \end{bmatrix}$.  We can then express the estimates of $\beta$ and $\gamma$ in terms of inverses of block matrices

$$ \begin{bmatrix} \beta \\ \gamma \end{bmatrix} = (Z^\top Z)^{-1} Z^\top Y $$
Note that $Z^\top Z$ is a block matrix with the following form, so

$$ Z^\top Z = \begin{bmatrix} X^\top X & X^\top D \\ D^\top X  & D^\top D\end{bmatrix} $$
and therefore
$$ \begin{bmatrix} \beta \\ \gamma \end{bmatrix} = \begin{bmatrix} X^\top X & X^\top D \\ D^\top X  & D^\top D\end{bmatrix}^{-1} \begin{bmatrix} X^\top Y \\ D^\top Y \end{bmatrix} \>.$$

From this form, we can see we need the bottom row.  The inverse of a 2x2 block matrix is a fucking hairy thing, here it is

$$
\mathbf{P}^{-1} = \begin{bmatrix}
\mathbf{A} & \mathbf{B} \\
\mathbf{C} & \mathbf{F}
\end{bmatrix}^{-1} = \begin{bmatrix}
(\mathbf{A} - \mathbf{B} \mathbf{F}^{-1} \mathbf{C})^{-1} & - (\mathbf{A} - \mathbf{B} \mathbf{F}^{-1} \mathbf{C})^{-1} \mathbf{B} \mathbf{F}^{-1} \\
- \mathbf{F}^{-1} \mathbf{C} (\mathbf{A} - \mathbf{B} \mathbf{F}^{-1} \mathbf{C})^{-1} & \mathbf{F}^{-1} + \mathbf{F}^{-1} \mathbf{C} (\mathbf{A} - \mathbf{B} \mathbf{F}^{-1} \mathbf{C})^{-1} \mathbf{B} \mathbf{F}^{-1}
\end{bmatrix}
$$

and so the inverse of our block matrix is 

$$
\begin{bmatrix} X^\top X & X^\top D \\ D^\top X & D^\top D \end{bmatrix}^{-1} = 
\begin{bmatrix}
(X^\top X)^{-1} + (X^\top X)^{-1} X^\top D \Delta^{-1} D^\top X (X^\top X)^{-1} & -(X^\top X)^{-1} X^\top D \Delta^{-1} \\
-\Delta^{-1} D^\top X (X^\top X)^{-1} & \Delta^{-1}
\end{bmatrix}
$$

$$ \text{where } \Delta = D^\top D - D^\top X (X^\top X)^{-1} X^\top D $$
The estimate for $\gamma$ is then

$$ \hat \gamma  = -\Delta^{-1}D^\top X(X^\top X)^{-1}X^\top Y + \Delta^{-1}D^\top Y $$
Note that the hat matrix appears in $\Delta$ and $\hat \gamma$. After some elbow grease, it can be shown

$$ 
\begin{align}
\hat \gamma &= \Delta^{-1}\bigl(D^\top Y - D^\top X (X^\top X)^{-1} X^\top Y\bigr)\\
&= \Delta^{-1} D^\top (I - X (X^\top X)^{-1} X^\top) Y\\
&= (D^\top M D)^{-1} D^\top M Y
\end{align}quart
$$

This proves that the application of the FWL theorem provides the same coefficients as a regression on the full matrix $Z$.

## Interpretation

Focusing on the modified regression, $MY = MD\gamma + Mu$ will be most helpful.  Recall that $M$ is a projection matrix which projects the right multiplied vector onto $\mathcal{C}(X)^\perp$.  So that means that $MY$ is the part of $Y$ which is uncorrelated with the columns of $X$, and likewise $MD$ is the part of $D$ which us uncorrelated with $X$. Then, when we perform the regression of $MY$ onto $MD$, we get the effect of $D$ on $Y$ independent of $X$.  That sounds a hell of a lot like $\gamma$ to me!
