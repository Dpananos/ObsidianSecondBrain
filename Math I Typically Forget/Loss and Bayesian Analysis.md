
## The Posterior
Consider the following Bayesian model

$$ \hat \theta \vert \theta \sim \mathrm{Normal}(\theta, s^2) $$
$$  \theta \sim \mathrm{Normal}(\mu, \tau^2) $$
If we assume $s$ is known, we can do inference on $\theta$ to obtain the posterior mean and variance as

$$ m = v^2 \left(\dfrac{\hat \theta}{s^2} + \dfrac{\mu}{\tau^2}\right) $$
$$ v^2 = \left(\dfrac{1}{s^2} + \dfrac{1}{\tau^2}\right)^{-1} $$

## Loss and Expected Loss

Let $a = \left\lbrace 0, 1\right\rbrace$ denote the action of shipping control ($a=0$) or treatment ($a=1$).  The loss associated with each is

$$ \mathcal{L}(\theta, 0) = \max \left\lbrace 0, \theta \right\rbrace $$
$$ \mathcal{L}(\theta, 1) = \max \left\lbrace 0, -\theta \right\rbrace $$
or

$$ \mathcal{L}(\theta, a) = \max \left\lbrace 0, (1-2a)\theta \right\rbrace   $$

The expected loss for each decision is

$$ E_\theta[\mathcal{L}(\theta, 0)] = E[\theta \mid \theta>0, \hat \theta]\Pr(\theta>0 \mid \hat \theta) $$

Note that 

$$ \Pr(\theta>0 \mid \hat \theta) = 1 -\Phi\left(\dfrac{-m}{v}\right) = \Phi\left(\dfrac{m}{v}\right) $$
and that the expectation above can be obtained from a truncated normal

$$ E[\theta \mid \theta>0, \hat \theta] = m + v \dfrac{\varphi(m/v)}{\Phi(m/v)} $$

Therefore, the expected loss is

$$ E_\theta[\mathcal{L}(\theta, 0)] = \Phi(m/v) m + v \varphi(m/v) $$
A similar calculation can be performed to calculate the expected loss for shipping treatment

$$ \Pr(\theta \lt 0 \mid \hat \theta) = \Phi\left(\dfrac{-m}{v}\right) $$

$$ E[\theta \mid \theta\lt 0, \hat \theta] = m - v \dfrac{\varphi(m/v)}{\Phi(-m/v)} $$

$$ E_\theta[\mathcal{L}(\theta, 1)] = -\Phi(-m/v) m + v \varphi(m/v) $$

Here, we have made judicious use of the fact that $\varphi$ is an even function.

## The Expected Value of Perfect Information

Let $a^\star(\theta)$ be the action we would take if we knew $\theta$.

$$a^\star (\theta) = \begin{cases} 1 & \theta>0 \\ 0 & \theta \lt 0 \end{cases} $$

The loss would then be

$$ \mathcal{L}(\theta, a^\star(\theta)) = \underset{a}{\min} \left\lbrace \mathcal{L}(\theta, a) \right\rbrace = 0 $$

We don't know $\theta$, so we choose the action which minimizes the expected loss loss

$$ \mathcal{L}^\star =  \min \left\lbrace E_\theta[\mathcal{L}(\theta, 1)], E_\theta[\mathcal{L}(\theta, 0)] \right\rbrace$$

The EVPI is 
$$ EVPI = \mathcal{L}^\star - E_\theta[\mathcal{L}(\theta, a^\star)] = \mathcal{L}^\star - 0 $$

so in this case, the EVPI is just the smallest expected loss.

## The Expected Value of Sample Information

Suppose we have observed $\hat\theta$, so that

$$ \theta\mid\hat\theta\sim\mathrm{Normal}(m,v^2), $$

and are considering collecting an additional, independent estimate $\hat\theta_{\text{new}}$ satisfying

$$ \hat\theta_{\text{new}}\mid\theta
\sim\mathrm{Normal}(\theta,s_{\text{new}}^2). $$

After observing $\hat\theta_{\text{new}}$, the new posterior variance and mean will be

$$
v_{\text{new}}^2
=\left(\dfrac{1}{v^2}+\dfrac{1}{s_{\text{new}}^2}\right)^{-1}
$$

and

$$
M
=v_{\text{new}}^2
\left(
\dfrac{m}{v^2}
+\dfrac{\hat\theta_{\text{new}}}{s_{\text{new}}^2}
\right).
$$

The capital letter $M$ emphasizes that the future posterior mean is random
before the new estimate is observed. The posterior predictive distribution of
the new estimate is

$$
\hat\theta_{\text{new}}\mid\hat\theta
\sim\mathrm{Normal}(m,v^2+s_{\text{new}}^2).
$$

To derive the distribution of $M$, first define the weight placed on the new
estimate as

$$
K
=\dfrac{v_{\text{new}}^2}{s_{\text{new}}^2}
=\dfrac{v^2}{v^2+s_{\text{new}}^2}.
$$

The future posterior mean can then be written as

$$
\begin{aligned}
M
&=v_{\text{new}}^2
\left(
\dfrac{m}{v^2}
+\dfrac{\hat\theta_{\text{new}}}{s_{\text{new}}^2}
\right)\\
&=(1-K)m+K\hat\theta_{\text{new}}\\
&=m+K(\hat\theta_{\text{new}}-m).
\end{aligned}
$$

Thus, $M$ is an affine transformation of the normally distributed posterior
predictive estimate and is itself normally distributed. Its preposterior mean
is

$$
E[M\mid\hat\theta]
=m+K\{E[\hat\theta_{\text{new}}\mid\hat\theta]-m\}
=m,
$$

and its preposterior variance is

$$
\begin{aligned}
\mathrm{Var}(M\mid\hat\theta)
&=K^2\mathrm{Var}(\hat\theta_{\text{new}}\mid\hat\theta)\\
&=K^2(v^2+s_{\text{new}}^2)\\
&=\dfrac{v^4}{v^2+s_{\text{new}}^2}.
\end{aligned}
$$

Meanwhile,

$$
\begin{aligned}
v^2-v_{\text{new}}^2
&=v^2-\dfrac{v^2s_{\text{new}}^2}{v^2+s_{\text{new}}^2}\\
&=\dfrac{v^4}{v^2+s_{\text{new}}^2}.
\end{aligned}
$$

Therefore,

$$
M\mid\hat\theta
\sim\mathrm{Normal}(m,w^2),
\qquad
w^2=v^2-v_{\text{new}}^2.
$$

The variance identity can also be obtained directly from the law of total
variance:

$$
\begin{aligned}
\mathrm{Var}(\theta\mid\hat\theta)
&=E\left[
\mathrm{Var}(\theta\mid\hat\theta,\hat\theta_{\text{new}})
\mid\hat\theta
\right]\\
&\quad+\mathrm{Var}\left(
E[\theta\mid\hat\theta,\hat\theta_{\text{new}}]
\mid\hat\theta
\right)\\
&=v_{\text{new}}^2+\mathrm{Var}(M\mid\hat\theta).
\end{aligned}
$$

Hence, the variance of the future posterior mean is exactly the reduction in
posterior variance obtained from the new sample.

After seeing the new sample, we choose the action having the smaller posterior
expected loss. Thus, the expected value of sample information is the reduction
in optimal expected loss:

$$
\mathrm{EVSI}
=\mathcal{L}^\star
-E_{\hat\theta_{\text{new}}\mid\hat\theta}
\left[
\min_{a\in\{0,1\}}
E_{\theta\mid\hat\theta,\hat\theta_{\text{new}}}
\{\mathcal{L}(\theta,a)\}
\right].
$$

To obtain a closed form, define

$$
g(x,\sigma)
=E[|X|],
\qquad
X\sim\mathrm{Normal}(x,\sigma^2),
$$

where

$$
g(x,\sigma)
=2\sigma\varphi(x/\sigma)
+x\left\lbrace2\Phi(x/\sigma)-1\right\rbrace.
$$

Using $x_+=(|x|+x)/2$ and $(-x)_+=(|x|-x)/2$, the optimal expected loss under
a $\mathrm{Normal}(x,\sigma^2)$ posterior is

$$
\min_a E[\mathcal{L}(\theta,a)]
=\dfrac{g(x,\sigma)-|x|}{2}.
$$

Consequently, the current optimal expected loss is

$$
\mathcal{L}^\star
=\dfrac{g(m,v)-|m|}{2}.
$$

Given the future posterior mean $M$, the future optimal expected loss is

$$
\dfrac{g(M,v_{\text{new}})-|M|}{2}.
$$

By the law of iterated expectations,

$$
E_{M\mid\hat\theta}[g(M,v_{\text{new}})]
=E_{\theta\mid\hat\theta}[|\theta|]
=g(m,v).
$$

Because $M\mid\hat\theta\sim\mathrm{Normal}(m,w^2)$,

$$ E_{M\mid\hat\theta}[|M|]=g(m,w). $$

Therefore,

$$
\boxed{
\mathrm{EVSI}
=\dfrac{g(m,w)-|m|}{2}
},
\qquad
w=\sqrt{v^2-v_{\text{new}}^2}.
$$

This expression has the expected limiting behavior. If the new sample contains
no information, then $v_{\text{new}}^2\to v^2$, $w\to0$, and
$\mathrm{EVSI}\to0$. If the new sample reveals $\theta$ perfectly, then
$v_{\text{new}}^2\to0$, $w\to v$, and
$\mathrm{EVSI}\to\mathrm{EVPI}$.
