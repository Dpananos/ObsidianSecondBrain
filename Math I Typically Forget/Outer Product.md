
## Motivation

Recall [[Projections]] and note that the projection is a linear operation on a vector.  Hence, it has matrix representation.  Suppose we have a vector, $\mathbf{u}$, which we want to project onto a vector, $\mathbf{v}$.  The projection is

$$ \operatorname{Proj}_{\mathbf{v}}(\mathbf{u}) =  \dfrac{\left< \mathbf{u}, \mathbf{v} \right>}{\left< \mathbf{v}, \mathbf{v}\right>} \mathbf{v} = \dfrac{\mathbf{u}^\top \mathbf{v}}{\mathbf{v}^\top \mathbf{v}} \mathbf{v} $$
How can we turn this into a matrix which acts on $\mathbf{u}$?  We can re-write this as 

$$ \dfrac{\mathbf{u}^\top \mathbf{v}}{\mathbf{v}^\top \mathbf{v}} \mathbf{v} = \dfrac{\mathbf{v}^\top \mathbf{u}}{\mathbf{v}^\top \mathbf{v}} \mathbf{v}$$

Note that $\mathbf{v}^T \mathbf{u}$ is a scalar, so this allows further re-writing of the projection as 

$$ \dfrac{\mathbf{v}^\top \mathbf{u}}{\mathbf{v}^\top \mathbf{v}} \mathbf{v} = \dfrac{ \mathbf{v} (\mathbf{v}^\top \mathbf{u})}{\mathbf{v}^\top \mathbf{v}}  $$
Which then places $u$ on the far right, appearing as if $\mathbf{v} \mathbf{v}^\top$ is a matrix.  If it were, then we could obtain the matrix representation of the projection.  Since $\mathbf{v}$ is a $n\times 1$  , then $\mathbf{v}^\top$ is a $1 \times n$ vector, and hence the product is an $n \times n$ matrix.

Naively following the matrix multiplication implied, this would mean that

$$ \mathbf{v} \mathbf{v}^\top = \begin{bmatrix} v_1 \mathbf{v}^\top  \\ v_2 \mathbf{v}^\top \\ \vdots \\ v_n \mathbf{v}^\top  \end{bmatrix}

=

\begin{bmatrix}

v_1v_1 & v_1v_2 & \cdots & v_1v_n \\

v_2v_1 & v_2v_2 & \cdots & v_2v_n \\

\vdots & \vdots & \ddots & \vdots \\

v_nv_1 & v_nv_2 & \cdots & v_nv_n

\end{bmatrix} $$
We have just derived the motivation and interpretation of the outer product.  If $\mathbf{v}$ is a unit vector, then $\mathbf{v} \mathbf{v}^\top$ is the matrix representation of the projection.

## Connections to OLS

We have seen in places like , [[Frisch-Waugh-Lovell (FWL) and Fixed Effects]], [[Projections]], and [[Robust Standard Errors]] that the the $X(X^\top X)^{-1}X^\top y$ is the prediction from a linear model.  The prediction also happens to be the projection of $y$ onto the column space of $X$.  What would happen of $X$ were a single column?  Then, it would be a vector and we would recover something very similar to what we see above!  

Hence, we have further intuition for what the outer product should look like and how it should be interpreted.  