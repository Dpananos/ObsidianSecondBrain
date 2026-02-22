
Let $\phi(x)$ denote the density function for a standard normal random variable, and let $\Phi(x) = \int_{-\infty}^x \phi(t) \, dt$ be the associated cumulative distribution function.  There are several properties of the density/cdf that are important to remember.

1) $\phi$ is an even function.  This means $\phi(x) = \phi(-x)$ which can be readily seen from the density. 
				$$ \phi(x) = \dfrac{1}{\sqrt{2\pi}} \exp(-x^2) $$
	Note that only the standard normal is even, and more generally, the pdf for a normal random variable is symmetric about it's mean, so $f(x+\mu) = f(-x+\mu)$, where

			$$ f(x) = \dfrac{1}{\sqrt{2\pi} \sigma} \exp\left(\dfrac{-(x-\mu)^2}{2\sigma^2}\right) $$
2) If $Z \sim \mathcal{N}(\mu, \sigma)$ then $(Z - \mu) / \sigma \sim \mathcal{N}(0, 1)$.  This procedure is said to be standardizing $Z$ -- that is transforming the random variable such that the density of the transformed variable is now standard normal.