## The Posterior

Consider the Bayesian model

$$
\hat\theta \mid \theta \sim \mathrm{Normal}(\theta,s^2)
$$

and

$$
\theta \sim \mathrm{Normal}(\mu,\tau^2).
$$

If $s$ is known, then the posterior distribution of $\theta$ is

$$
\theta\mid\hat\theta\sim\mathrm{Normal}(m,v^2),
$$

where

$$
m=v^2\left(
\frac{\hat\theta}{s^2}
+\frac{\mu}{\tau^2}
\right)
$$

and

$$
v^2=\left(
\frac{1}{s^2}
+\frac{1}{\tau^2}
\right)^{-1}.
$$

## Utility Loss and Posterior Expected Loss

Let $a\in\{0,1\}$ denote the action of shipping control ($a=0$) or treatment
($a=1$). Define the utility loss as

$$
\mathcal{L}_{\mathrm{utility}}(\theta,a)=-a\theta.
$$

Thus,

$$
\mathcal{L}_{\mathrm{utility}}(\theta,0)=0
$$

and

$$
\mathcal{L}_{\mathrm{utility}}(\theta,1)=-\theta.
$$

This formulation takes control as the baseline. If treatment is beneficial,
shipping treatment produces negative loss, representing a gain relative to
control. Equivalently, one could define the utility as
$U(\theta,a)=a\theta$ and minimize its negative.

The posterior expected loss, or posterior risk, of action $a$ is

$$
\begin{aligned}
\mathcal{R}(a\mid\hat\theta)
&=E_\theta[
\mathcal{L}_{\mathrm{utility}}(\theta,a)
\mid\hat\theta
] \\
&=-aE[\theta\mid\hat\theta] \\
&=-am.
\end{aligned}
$$

Consequently,

$$
\mathcal{R}(0\mid\hat\theta)=0
$$

and

$$
\mathcal{R}(1\mid\hat\theta)=-m.
$$

The Bayes action is therefore

$$
a^\star(\hat\theta)
=
\begin{cases}
1, & m>0, \\
0, & m<0.
\end{cases}
$$

Either action is optimal when $m=0$. The optimal current posterior expected
loss is

$$
\begin{aligned}
\mathcal{R}^\star
&=\min_{a\in\{0,1\}}\mathcal{R}(a\mid\hat\theta) \\
&=\min(0,-m) \\
&=-m_+,
\end{aligned}
$$

where $x_+=\max(0,x)$.

## The Expected Value of Perfect Information

Suppose we knew the true value of $\theta$ before choosing an action. The
optimal action would be

$$
a^\star(\theta)
=
\begin{cases}
1, & \theta>0, \\
0, & \theta<0.
\end{cases}
$$

The corresponding loss is

$$
\begin{aligned}
\min_{a\in\{0,1\}}
\mathcal{L}_{\mathrm{utility}}(\theta,a)
&=\min(0,-\theta) \\
&=-\theta_+.
\end{aligned}
$$

Therefore, the expected loss under perfect information is

$$
E_\theta\left[
\min_{a\in\{0,1\}}
\mathcal{L}_{\mathrm{utility}}(\theta,a)
\mid\hat\theta
\right]
=
-E[\theta_+\mid\hat\theta].
$$

For $\theta\mid\hat\theta\sim\mathrm{Normal}(m,v^2)$,

$$
\Pr(\theta>0\mid\hat\theta)
=
\Phi\left(\frac{m}{v}\right)
$$

and

$$
E[\theta\mid\theta>0,\hat\theta]
=
m+v
\frac{\varphi(m/v)}{\Phi(m/v)}.
$$

It follows that

$$
\begin{aligned}
E[\theta_+\mid\hat\theta]
&=
E[\theta\mid\theta>0,\hat\theta]
\Pr(\theta>0\mid\hat\theta) \\
&=
m\Phi(m/v)+v\varphi(m/v).
\end{aligned}
$$

The expected value of perfect information is the reduction in expected loss
obtained by learning $\theta$ before acting:

$$
\begin{aligned}
\mathrm{EVPI}
&=
\mathcal{R}^\star
-
E_\theta\left[
\min_{a\in\{0,1\}}
\mathcal{L}_{\mathrm{utility}}(\theta,a)
\mid\hat\theta
\right] \\
&=
-m_+ + E[\theta_+\mid\hat\theta].
\end{aligned}
$$

Therefore,

$$
\boxed{
\mathrm{EVPI}
=
m\Phi(m/v)
+v\varphi(m/v)
-m_+
}.
$$

An equivalent expression can be obtained by defining

$$
g(x,\sigma)
=
E[|X|],
\qquad
X\sim\mathrm{Normal}(x,\sigma^2),
$$

where

$$
g(x,\sigma)
=
2\sigma\varphi(x/\sigma)
+x\left\{2\Phi(x/\sigma)-1\right\}.
$$

Because

$$
x_+=\frac{|x|+x}{2},
$$

we have

$$
E[\theta_+\mid\hat\theta]
=
\frac{g(m,v)+m}{2}
$$

and

$$
m_+
=
\frac{|m|+m}{2}.
$$

Hence,

$$
\boxed{
\mathrm{EVPI}
=
\frac{g(m,v)-|m|}{2}
}.
$$

## The Expected Value of Sample Information

Suppose we have observed $\hat\theta$, so that

$$
\theta\mid\hat\theta
\sim\mathrm{Normal}(m,v^2),
$$

and are considering collecting an additional, independent estimate
$\hat\theta_{\mathrm{new}}$ satisfying

$$
\hat\theta_{\mathrm{new}}\mid\theta
\sim\mathrm{Normal}(\theta,s_{\mathrm{new}}^2).
$$

After observing $\hat\theta_{\mathrm{new}}$, the new posterior variance and
mean will be

$$
v_{\mathrm{new}}^2
=
\left(
\frac{1}{v^2}
+\frac{1}{s_{\mathrm{new}}^2}
\right)^{-1}
$$

and

$$
M
=
v_{\mathrm{new}}^2
\left(
\frac{m}{v^2}
+\frac{\hat\theta_{\mathrm{new}}}{s_{\mathrm{new}}^2}
\right).
$$

The capital letter $M$ emphasizes that the future posterior mean is random
before the new estimate is observed. The posterior predictive distribution of
the new estimate is

$$
\hat\theta_{\mathrm{new}}\mid\hat\theta
\sim
\mathrm{Normal}(m,v^2+s_{\mathrm{new}}^2).
$$

Define the weight placed on the new estimate as

$$
K
=
\frac{v_{\mathrm{new}}^2}{s_{\mathrm{new}}^2}
=
\frac{v^2}{v^2+s_{\mathrm{new}}^2}.
$$

The future posterior mean can then be written as

$$
\begin{aligned}
M
&=
v_{\mathrm{new}}^2
\left(
\frac{m}{v^2}
+\frac{\hat\theta_{\mathrm{new}}}{s_{\mathrm{new}}^2}
\right) \\
&=(1-K)m+K\hat\theta_{\mathrm{new}} \\
&=m+K(\hat\theta_{\mathrm{new}}-m).
\end{aligned}
$$

Thus, $M$ is an affine transformation of the normally distributed posterior
predictive estimate and is itself normally distributed. Its preposterior mean
is

$$
\begin{aligned}
E[M\mid\hat\theta]
&=
m+K\left\{
E[\hat\theta_{\mathrm{new}}\mid\hat\theta]-m
\right\} \\
&=m,
\end{aligned}
$$

and its preposterior variance is

$$
\begin{aligned}
\mathrm{Var}(M\mid\hat\theta)
&=
K^2
\mathrm{Var}(\hat\theta_{\mathrm{new}}\mid\hat\theta) \\
&=
K^2(v^2+s_{\mathrm{new}}^2) \\
&=
\frac{v^4}{v^2+s_{\mathrm{new}}^2}.
\end{aligned}
$$

Meanwhile,

$$
\begin{aligned}
v^2-v_{\mathrm{new}}^2
&=
v^2-
\frac{v^2s_{\mathrm{new}}^2}
{v^2+s_{\mathrm{new}}^2} \\
&=
\frac{v^4}{v^2+s_{\mathrm{new}}^2}.
\end{aligned}
$$

Therefore,

$$
M\mid\hat\theta
\sim\mathrm{Normal}(m,w^2),
\qquad
w^2=v^2-v_{\mathrm{new}}^2.
$$

The variance identity can also be obtained from the law of total variance:

$$
\begin{aligned}
\mathrm{Var}(\theta\mid\hat\theta)
&=
E\left[
\mathrm{Var}
(\theta\mid\hat\theta,\hat\theta_{\mathrm{new}})
\mid\hat\theta
\right] \\
&\quad+
\mathrm{Var}\left(
E[\theta\mid\hat\theta,\hat\theta_{\mathrm{new}}]
\mid\hat\theta
\right) \\
&=
v_{\mathrm{new}}^2
+\mathrm{Var}(M\mid\hat\theta).
\end{aligned}
$$

Hence, the variance of the future posterior mean is exactly the reduction in
posterior variance obtained from the new sample.

After observing the new sample, the posterior expected losses are

$$
\mathcal{R}_{\mathrm{new}}(0)=0
$$

and

$$
\mathcal{R}_{\mathrm{new}}(1)=-M.
$$

The future optimal posterior expected loss is therefore

$$
\min_{a\in\{0,1\}}
\mathcal{R}_{\mathrm{new}}(a)
=
\min(0,-M)
=
-M_+.
$$

The expected value of sample information is the reduction in optimal expected
loss:

$$
\begin{aligned}
\mathrm{EVSI}
&=
\mathcal{R}^\star
-
E_{M\mid\hat\theta}
\left[
\min_{a\in\{0,1\}}
\mathcal{R}_{\mathrm{new}}(a)
\right] \\
&=
-m_+
-
E_{M\mid\hat\theta}[-M_+] \\
&=
E_{M\mid\hat\theta}[M_+]-m_+.
\end{aligned}
$$

Since $M\mid\hat\theta\sim\mathrm{Normal}(m,w^2)$,

$$
E[M_+\mid\hat\theta]
=
m\Phi(m/w)+w\varphi(m/w).
$$

Thus,

$$
\boxed{
\mathrm{EVSI}
=
m\Phi(m/w)
+w\varphi(m/w)
-m_+
},
\qquad
w=\sqrt{v^2-v_{\mathrm{new}}^2}.
$$

Equivalently,

$$
E[M_+\mid\hat\theta]
=
\frac{g(m,w)+m}{2},
$$

so

$$
\boxed{
\mathrm{EVSI}
=
\frac{g(m,w)-|m|}{2}
},
\qquad
w=\sqrt{v^2-v_{\mathrm{new}}^2}.
$$

Expressions involving $m/w$ are interpreted by continuity when $w=0$.

This expression has the expected limiting behavior. If the new sample contains
no information, then $v_{\mathrm{new}}^2\to v^2$, $w\to0$, and

$$
\mathrm{EVSI}\to0.
$$

If the new sample reveals $\theta$ perfectly, then
$v_{\mathrm{new}}^2\to0$, $w\to v$, and

$$
\mathrm{EVSI}\to\mathrm{EVPI}.
$$
