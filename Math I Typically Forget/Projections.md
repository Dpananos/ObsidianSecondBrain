## Projection of $\mathbf{v}$ onto $\mathbf{u}$

Suppose we have vectors $\mathbf{u}, \mathbf{v} \in \mathbb{R}^n$. The projection of $\mathbf{v}$ onto $\mathbf{u}$ is written as $\mathbf{w} = \operatorname{proj}_{\mathbf{u}}(\mathbf{v})$ and is the vector in the direction of $\mathbf{u}$ such that $\langle \mathbf{v}-\mathbf{w}, \mathbf{u} \rangle = 0$. We know that $\mathbf{w}$ is collinear to $\mathbf{u}$, so we can write it as $a \cdot \mathbf{u}$, so we just need to determine $a$.  To find the appropriate scalar, recall that geometrically the projection satisfies

$$ \left| \mathbf{v} \right| \cos(\theta) = \left| \operatorname{proj}_{\mathbf{u}}(\mathbf{v}) \right|  \>. $$

We also know that

$$ \mathbf{u} \cdot \mathbf{v} = \left| \mathbf{u} \right| \left| \mathbf{v} \right| \cos(\theta). $$

Simple algebra yields

$$ \left| \operatorname{proj}_{\mathbf{u}}(\mathbf{v}) \right| = \frac{\mathbf{u} \cdot \mathbf{v}}{\left| \mathbf{u} \right|}. $$

This is the length of the projection. To create a vector with that length in the direction of $\mathbf{u}$, simply make $\mathbf{u}$ a unit vector, yielding

$$  
\mathbf{w} = \operatorname{proj}_{\mathbf{u}}(\mathbf{v})  
= \frac{\mathbf{u} \cdot \mathbf{v}}{|\mathbf{u}|^2} \mathbf{u}  
= \frac{\langle \mathbf{u}, \mathbf{v} \rangle}{\langle \mathbf{u}, \mathbf{u} \rangle} \mathbf{u}.  
$$

## Orthogonal Projection onto a Subspace

A single vector is itself a subspace of the vector space in which it is an element.  As such, it is reasonable to ask if the concept of vector projection can be extended to subspaces of dimension larger than one.

The "orthogonal projection" of a vector $\mathbf{v}$ onto the subspace spanned by the basis $\left\{ \mathbf{u}_1, \cdots, \mathbf{u}_n \right\}$ is the vector $\mathbf{w} = \sum_k a_k \mathbf{u}_k$ such that $\left< \mathbf{v} - \mathbf{w}, \mathbf{u}_k\right>=0 \quad \forall k = 1, \cdots, n$.   The coordinates with respect to the basis are typically of interest since, given the coordinates, we can construct the orthogonal projection via linear combination.

We can express the orthogonality with respect to the basis more succinctly by combining the basis vectors column wise into the matrix $U$.  Then, we want to solve the following equation for $a$

$$ U^\top(\mathbf{v} - U \mathbf{a}) = 0 \>. $$
Here, we have written $\mathbf{w} = U \mathbf{a}$, and left multiplication by $U^\top$ expresses the inner product among all basis vectors.  Simple matrix algebra yields the solution

$$ U^\top U \mathbf{a} = U^\top \mathbf{v{}} \>. $$
$$  \mathbf{a} = (U^\top U)^{-1} U^\top \mathbf{v} \>. $$
Here, we have used the fact that $\left\{ \mathbf{u}_1, \cdots, \mathbf{u}_n \right\}$ are a basis, hence $U$ is full rank, hence the Graham matrix $U^\top U$ is invertible. As a corollary, the matrix

$$ P = U (U^\top U)^{-1} U^\top  $$
is known as a "Projection matrix" since $P\mathbf{v}$ maps $\mathbf{v}$ to its projection in the subspace spanned by the columns of $U$.  This matrix has also been called the "hat matrix" in courses on Linear Regression since "it puts the hat on y".  This serves to highlight that the prediction from a model fit by OLS is the projection of the vector $\mathbf{y}$ onto the column space of the design matrix.

When the columns of $U$ are an orthogonal basis, a great simplification occurs

$$ U^\top U = \operatorname{diag}(||\mathbf{u}_1||^2, \cdots, ||\mathbf{u}_n||^2) $$
which is then easily invertible

$$ (U^\top U)^{-1} = \operatorname{diag}(||\mathbf{u}_1||^{-2}, \cdots, ||\mathbf{u}_n||^{-2}) $$which then means the coordinates for the orthogonal projection are

$$ a_j = \dfrac{\left<\mathbf{u_j}, \mathbf{v}\right>}{\left<\mathbf{u_j}, \mathbf{u_j}\right>} $$
which is _further_ simplified when the basis is an orthonormal basis.  This underscores the elegance of orthogonality.  Interestingly, it would seem that the orthogonal projection is the linear combination of projections onto each basis vector.

## Applications

Projections matrices have appeared in [[Frisch-Waugh-Lovell (FWL) and Fixed Effects]]