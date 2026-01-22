# An Energy-aware Skill Score for Earth science climate models

## Introduction

We need a score that allows researchers to compare different codes running on different machines.

A score based on kWh is useful as it accounts for different clock frequencies, memory bandwidths and any other hardware features that contribute to performance. Such a score is therefore hardware agnostic; it is also independent of the machine location as far as the financial cost of running a model is concerned.

Allowing for differences between code behaviour requires more thought. We must consider not only differences in functionality but also differences in input data. A code’s computational profile will depend on its starting parameters. The fact that we are considering Earth science (ES) codes only (e.g. ocean, atmosphere and climate models) at least makes the problem of finding a suitable skill score tractable.

First, we can assume that any ES code will output a certain amount of data every $n$ time steps, where the time step can represent years or hours depending on the model. The amount of data generated is a reasonable proxy for the number of properties being simulated. We refer to this as the Raw Scientific Output (RSO), measured in bytes. The associated performance score is simply the RSO per unit of energy consumed. We depart here from the more usual simulation time per unit of runtime as such a score does not give us an indication of how much work is being done per timestep.


## Weighting the score by difficulty of prediction

Next, we explore how to assess the level of achievement implied by the model prediction. We start with a ground truth that can be compared with computer-generated predictions. For a given property, $p$,  such as temperature or salinity, we can calculate the differences between a time series of $N$ observed values ($x_{t}$) and the corresponding forecasted values ($f_{t}$), summarising these as a root mean square relative error, $\sqrt{\varepsilon_{p}/N}$, where

```math
\varepsilon_{p} = \sum_{t=1}^{N} \left( \frac{x_{t} - f_{t}}{f_{t}} \right)^{2} \;.
```

However, given the distribution of $x$, how difficult is it to make a successful prediction? Mayer & Yang (2024) propose making use of the lag $h$ autocorrelation of the observations, $\gamma(h)$. The idea here is that $x_{t}$ (measured at time $t$) will exhibit some degree of correlation with $x_{t-h}$ (the same property measured at time $t-h$). The autocorrelation is (minus) one if the observations are fully (anti) correlated, and so, the difficulty of making a successful prediction is lowest when $\gamma = \pm1$. It follows that the difficulty is highest when the autocorrelation is zero. In fact, the difficulty is impossible for a randomly varying time series characterised by $\gamma = 0$. We now introduce a skill based score,

```math
\chi_{p,h} = \frac{1 - \sqrt{\varepsilon_{p,h} / n}}{\gamma^2(h)}
```
where
```math
\varepsilon_{p,h} = \sum_{t=1}^{N} \left( \frac{x_{t} - f_{t}}{f_{t}} \right)^2 \;,
```
```math
2\leq n \leq N \;,
```
```math
t_{i+1} - t_{i} = h \,\,\forall\,\, i \in \{1..n\} \;.
```

This is specific to a particular property and autocorrelation lag ($h$). In the numerator of the expression for $\chi_{p,h}\,$, we subtract the error term from one so that a low error corresponds to a high skill score. Further, the error term is also tied to a single $h$ value, which is now the time interval between each successive forecast ($f_{t}$) and observation ($x_{t}$), i.e. the summation for $\chi_{p,h}$ is over some subset of the original time series of size $N$.

The (discrete) set of relevant lag values will depend on the scope of the simulation, whether it is a climate model or a weather forecast. In general, the autocorrelation will tend to zero for smaller lag values: fluctuations are random once time scales are short enough, e.g. minute-to-minute variations in wind speed. It is necessary therefore to bound $h$ such that undefined skill scores are avoided. (The maximum $h$ value is limited by the simulation time.) Although we propose to discard low lag values when calculating a skill score, the impact of short term variability is known to influence longer term variations as a consequence of Hasselmann’s stochastic theory (Hasselmann 1976). We expect therefore that the skill scores associated with longer lag values will reflect the model’s success in allowing long-term phenomena to be influenced by continuous short-term random excitations. In other words, the fidelity of a model on timescales corresponding to small $h$ values is still being captured by the skill score. 

There are three more points worth mentioning. In some cases, $\sqrt{\varepsilon_{p,h}/n}>1$, which would cause the skill score to become more negative the closer the autocorrelation is to zero. This can be remedied by changing the  $\gamma^2(h)$ denominator in the equation above to $1-\gamma^2(h)$ whenever the error exceeds one. Thus, the skill score now becomes less negative for observations that show lower autocorrelation and are therefore harder to predict.

Secondly, as with the RSO, we can divide the skill score, $\chi$, by the unit of energy usage to obtain an energy-aware skill-based score. Thirdly, autocorrelations can be applied to spatially varying properties too.

## Visualisation

How best to visualise skill scores? We could draw [circular bar plots](https://python-graph-gallery.com/web-circular-barplot-with-matplotlib/), where the bars take the appearance of flower petals. The radial position of the bars corresponds to particular autocorrelation lags. The length and colour of the bars give the skill scores. Each such plot would be tied to a specific property.

## Conclusion

The RSO per kWh metric allows comparison between different ES codes running on different machines, but has only a loose association with skill. The amount of data generated by a simulation is related to the model resolution and/or the number of properties modelled. We assume that the volume of data being produced is a proxy for the complexity of a code’s computations and so is also related to the skill required to run the model. This assumption would be invalidated however if the data being computed are diagnostic values that are based on a relatively small number of prognostic fields.

The skill score derived from the autocorrelation of recorded observations is a more direct attempt to measure the level of the achievement represented by a model’s successful predictions. Making the skill score energy aware is simply a matter of dividing by the energy consumed: a high skill score is penalised if a model requires many kWhs to generate a prediction.


## References

Hasselmann, K. [Stochastic climate models Part I. Theory](https://onlinelibrary.wiley.com/doi/10.1111/j.2153-3490.1976.tb00696.x). Tellus 28, 473–485 (1976).

Mayer, M. J., Yang, D. [Potential root mean square error skill score](https://pubs.aip.org/aip/jrse/article-abstract/16/1/016501/3267300/Potential-root-mean-square-error-skill-score?redirectedFrom=fulltext). J. Renewable Sustainable Energy 16, 016501 (2024).
