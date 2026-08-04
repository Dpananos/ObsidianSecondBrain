
# Notation

Let $A$ be a binary treatment indicator, and let $M$ be mediation (i.e. $M$ is caused by $A$ and is a cause of the outcome), and let $Y$ be the outcome of interest.

* $M(a)$ is the potential value of the mediator when $A=a$
* $Y(a, m)$ is the potential outcome if the treatment were set to $A=a$ and the mediator were set to $M = m$.

There are then four primary observable and counterfactual states for a binary treatment

* $Y(1, M(1))$ : The potential outcome when treated, with the mediator taking its natural value under the treatment.
* $Y(0, M(0))$ : The potential outcome when untreated, with the mediator taking its natural value under the treatment.
* $Y(1, M(0))$ : A cross-world counterfactual.  The outcome if treated but the mediator is forced to take the value it would have naturally taken under control
* $Y(0, M(1))$ : A cross-world counterfactual.  The outcome if untreated but the mediator is forced to take the value it would have naturally taken under treatment.

# Estimands

## Total Effect

$$ TE = E[Y(1, M(1))] - E[Y(0, M(0))] $$
## Controlled Direct Effect

$$ CDE(m) = E[Y(1, m)] - E[Y(0, m)] $$
## Natural Direct Effect

This is the effect were we to disable the effect through the mediator.  We call this "natural" because the mediator is held at the level it would have naturally taken under control condition $A=0$.

$$ NDE = E[Y(1, M(0))] - E[Y(0, M(0))] $$

## Natural Indirect Effect

This captures the effect that happens solely through the mediator.

$$NIE = E[Y(1, M(1)] - E[Y(1, M(0)]$$

## Decomposition

Note, that we can decompose the total effect into the sum of the natural direct effect and natural indirect effect.


$$\begin{align}
NIE + NDE  &=  E[Y(1, M(1)] - E[Y(1, M(0)] + E[Y(1, M(0))] - E[Y(0, M(0))] \\
&= E[Y(1, M(1)] - \cancel{E[Y(1, M(0)]} + \cancel{E[Y(1, M(0)]} - E[Y(0, M(0))] \\ 
& = E[Y(1, M(1))] - E[Y(0, M(0))] \\
& = TE
\end{align}$$

