Consider a similar scenario to [[Robust Standard Errors]], in which we have an OLS problem

$$  
\hat{\beta} = \underset{\beta \in \mathbb{R}^p}{\arg \min} \left \lVert y - X\beta \right \rVert \>. $$
with solution 

$$ \hat{\beta} = (X^\top X)^{-1} X^\top y  \>.$$

We have seen that assuming $\varepsilon \varepsilon^\top \neq \sigma^2 I$ leads to the following estimator for the variance/covariance

$$  \operatorname{Var}(\hat{\beta} \mid X) = (X^\top X)^{-1} X^\top \Omega \, X\, (X^\top X)^{-1} \>.$$

In the event we sample clusters, $c \in C$, as opposed to individual units (e.g. sample companies and analyze at the employee level) we need to use cluster robust standard errors.  The estimator is similar, but now the "meat" of the sandwich is expressed as 

$$ X^\top \Omega X = \sum_{c \in C} X^\top_c \Omega_c X_c \>.$$

Although it is possible to understand $\Omega$ as a block diagonal matrix, it is much easier to consider $\Omega_c$.  Within each cluster, we estimate the variance/covariance matrix and then add these together.  The operation is even more easily understood in code as an increment of a matrix of 0s.

```python
import numpy as np

def clustered_standard_errors(X: np.ndarray, 
                             residuals: np.ndarray, 
                             clusters: np.ndarray) -> np.ndarray:
    """
    Computes the clustered (robust) covariance matrix (sandwich estimator).

    Args:
        X: Regressor matrix, shape (n_obs, n_features)
        residuals: Residuals from the OLS regression, shape (n_obs,)
        clusters: Array of cluster IDs for each observation, shape (n_obs,)

    Returns:
        The (n_features, n_features) clustered covariance matrix.
    """
    
    # 1. Get dimensions
    n_obs, p = X.shape
    
    # 2. Calculate the "Bread"
    # (X'X)^-1
    try:
        # Use pinv for numerical stability, as in your original code
        bread = np.linalg.pinv(X.T @ X)
    except np.linalg.LinAlgError:
        print("Error: X'X matrix is not invertible.")
        return np.full((p, p), np.nan)

    # 3. Calculate the "Meat"
    # Sum over clusters g: [ (X_g' u_g) (X_g' u_g)' ]
    meat = np.zeros((p, p))
    unique_clusters = np.unique(clusters)
    G = len(unique_clusters) # Number of clusters

    for c in unique_clusters:
        # Find all observations in this cluster
        ix = (clusters == c)
        
        # Select the X and residuals for this cluster
        Xc = X[ix, :]
        uc = residuals[ix] # This will be a 1D array
        
        # Calculate the cluster's contribution
        # Xc.T is (p, n_c), uc is (n_c,) -> v is (p,)
        v = Xc.T @ uc
        
        # Add to the "meat"
        # np.outer((p,), (p,)) creates the (p, p) matrix
        meat += np.outer(v, v)

    # 4. Apply Degrees of Freedom Correction
    # This matches the correction used by default in Stata
    # (G / (G-1)) * ((n-1) / (n-p))
    # Note: n_obs = n
    correction_factor = (G / (G - 1)) * ((n_obs - 1) / (n_obs - p))

    # 5. Assemble the "Sandwich"
    # (X'X)^-1 [ (Meat) ] (X'X)^-1
    cov = bread @ meat @ bread
    
    # Return the corrected covariance matrix
    return cov * correction_factor
```

## Bias Correction