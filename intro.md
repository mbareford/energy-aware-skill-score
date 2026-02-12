# An Energy-aware Skill Score for Earth science climate models

## Introduction

We need a score that allows researchers to compare different codes running on different machines.

A score based on kWh is useful as it accounts for different clock frequencies, memory bandwidths and any other hardware features that contribute to performance. Such a score is therefore hardware agnostic; it is also independent of the machine location as far as the financial cost of running a model is concerned.

Allowing for differences between code behaviour requires more thought. We must consider not only differences in functionality but also differences in input data. A code’s computational profile will depend on its starting parameters. The fact that we are considering Earth science (ES) codes only (e.g. ocean, atmosphere and climate models) at least makes the problem of finding a suitable skill score tractable.

First, we can assume that any ES code will output a certain amount of data every $n$ time steps, where the time step can represent years or hours depending on the model. The amount of data generated is a reasonable proxy for the number of properties being simulated. We refer to this as the Raw Scientific Output (RSO), measured in bytes. The associated performance score is simply the RSO per unit of energy consumed. We depart here from the more usual simulation time per unit of runtime as such a score does not give us an indication of how much work is being done per timestep.


## Weighting the score by difficulty of prediction

Next, we explore how to assess the level of achievement implied by the model prediction. We start with a ground truth that can be compared with computer-generated predictions. For a given property, $p$,  such as temperature or salinity, we can calculate the differences between a time series of $N$ observed values ($x_{p,i}$) and the corresponding forecasted values ($f_{p,i}$).

An obvious first approach is to summarise these differences as a root mean square relative error (RMSE), $\sqrt{\frac{\varepsilon_{p}}{N}}$, where

```math
\varepsilon_{p} = \sum_{i=1}^{N} \left( \frac{x_{p,i} - f_{p,i}}{x_{p,i}} \right)^{2} \;.
```

Basing a skill score on such an expression gives rise to two problems however. Firstly, it is possible for $x_{p,i}$ to take zero values (fairly common for measurements of temperature and precipitation) causing $\varepsilon_{p}$ to be undefined. Secondly, converting the RMSE into a skill score by, say, subtracting it from one, can result in a negative score whenever the error is sufficiently large.

A better technique is to calculate the mean absolute scaled error (MASE) as recommended by Hyndman & Koehler (2006). This method is based on the mean absolute point-to-point difference calculated over the full time series of observations.

```math
\omega_{p} = \frac{1}{N-1}\sum_{i=2}^{N} | x_{p,i} - x_{p,i-1} | \;.
```

The $\omega_{p}$ term provides the scaling for the MASE, which we denote by $\xi_{p}$,

```math
\xi_{p} = \frac{1}{N}\sum_{i=1}^{N} \frac{|x_{p,i} - f_{p,i}|}{\omega_{p}} \;.
```

This scaled error ($\xi_{p}$) is less than one for predictions that beat the average one-step naïve forecast computed in-sample. Conversely, $\xi_{p} > 1$ if the prediction is worse than the average one-step forecast. Further, an undefined $\xi_{p}$ can now only occur if $\omega_{p}=0$, meaning all meaurements in the time series are equal. Such a situation is easy to detect and its likelihood depends on the number of elements in the time series. An alternative error scaling would need to be considered for the case when $x_{p,i} = x_{p,i-1} \forall i \in \{2..N\}$. 

Knowing that the range of $\xi_{p}$ starts from zero and is unbounded, we need a way to map $\xi_{p}$ to the range $0-1$, allowing us to derive a skill score.

```math
\chi_{p} = 1 - \frac{\xi_{p}}{\xi_{p}+1} \;.
```

The expression above avoids negative skill scores, since $\xi_{p}=\infty$ yields a zero score and $\xi_{p}=0$ gives a score of one.

However, given the distribution of $x_{p,i}$, how difficult was it to make the successful prediction in the first place? Mayer & Yang (2024) propose making use of the lag $h$ autocorrelation of the observations, $\gamma(h)$. The idea here is that $x_{p,i}$ (measured at time $i$) will exhibit some degree of correlation with $x_{p,i-h}$ (the same property measured at time $i-h$). The autocorrelation is (minus) one if the observations are fully (anti) correlated, and so, the difficulty of making a successful prediction is lowest when $\gamma = \pm1$. It follows that the difficulty is highest when the autocorrelation is zero. In fact, the difficulty is impossible for a randomly varying time series characterised by $\gamma = 0$. We now modify the expression for $\chi_{p}$ by introducing the autocorrelation,

```math
\chi_{p,h} = \Bigg[1 - \frac{\xi_{p,h}}{\xi_{p,h}+1}\Bigg]\big(1-|\gamma(h)|\big)
```
where
```math
\xi_{p,h} = \frac{1}{n}\sum_{i=1}^{n} \frac{|x_{p,i} - f_{p,i}|}{\omega_{p}} 
```
and
```math
\omega_{p,h} = \frac{1}{n-1}\sum_{i=2}^{n} | x_{p,i} - x_{p,i-1} | 
```
for
```math
2\leq n \leq N \;,
```
```math
i_{t+1} - i_{t} = h \,\,\forall\,\, t \in \{1..n\} \;.
```

And so, $\chi_{p,h}$ is now specific to a particular property and autocorrelation lag ($h$). The $\xi_{p,h}$ and $\omega_{p,h}$ terms are tied to a single $h$ value, which is now the time interval between each successive forecast ($f_{p,i}$) and observation ($x_{p,i}$), i.e. the summations used to calculate $\xi_{p,h}$ and $\omega_{p,h}$ are over some subset of the original time series of size $N$. All values of $h$ compatible with the time series can be used with the formulation above. The unweighted skill score (i.e. the expression in square brackets that partly determines $\chi_{p,h}$) is preserved when $\gamma(h)=0$. At the opposite extreme, $\chi_{p,h}=0$ when $\gamma(h)=\pm1$. A zero skill score can occur if and only if the observed values are perfectly correlated or anti-correlated. Note, the unweighted skill score itself can only become zero if $\xi_{p,h}=\infty$, which could only happen if the difference between a forecast and observation is also infinite, i.e. the simulation making the forecast has obviously failed and so is not fit to be scored.

The (discrete) set of relevant lag values will depend on the scope of the simulation, whether it is a climate model or a weather forecast. The maximum $h$ value will of course be limited by the simulation time. For the smallest lag values, the autocorrelation may approach zero: fluctuations can appear random over short time scales, e.g. minute-to-minute strong variations in wind speed. Nevertheless, the impact of short term variability is known to influence longer term variations as a consequence of Hasselmann’s stochastic theory (Hasselmann 1976). We expect therefore that the skill scores associated with longer lag values will reflect the model’s success in allowing long-term phenomena to be influenced by continuous short-term random excitations. In other words, the fidelity of a model on timescales corresponding to small $h$ values is still being captured by the skill scores based on long lag times.

As with the RSO, we can divide a skill score ($\chi_{p,h}$) by the unit of energy usage to obtain an energy-aware skill-based score.

Lastly, all of the above can be applied to spatially varying properties too.

## Visualisation

How best to visualise skill scores? We could draw [circular bar plots](https://python-graph-gallery.com/web-circular-barplot-with-matplotlib/), where the bars take the appearance of flower petals. The radial position of the bars corresponds to particular autocorrelation lags. The length and colour of the bars give the skill scores. Each such plot would be tied to a specific property.

More conventially, we could simply plot a property skill score against autocorrelation lag. At smaller scales, randomness is more present, suppressing $\chi_{p,h}$ but also minimising the reduction in the score caused by the $1 - |\gamma(h)|$ component. We anticipate less randomness over longer lag values, giving rise to a higher unweighted skill score that is dampened by rising $|\gamma(h)|$. In other words, a plot of $\chi_{p,h}$ against $h$ may be close in form to a horizontal line that can be summarised as a mean $\chi_{p}$ over all lag values.

## Conclusion

The RSO per kWh metric allows comparison between different ES codes running on different machines, but has only a loose association with skill. The amount of data generated by a simulation is related to the model resolution and/or the number of properties modelled. We assume that the volume of data being produced is a proxy for the complexity of a code’s computations and so is also related to the skill required to run the model. This assumption would be invalidated however if the data being computed are diagnostic values that are based on a relatively small number of prognostic fields.

The skill score based on the autocorrelation of recorded observations is a more direct attempt to measure the level of the achievement represented by a model’s successful predictions. Making the skill score energy aware is simply a matter of dividing by the energy consumed: a high skill score is penalised if a model requires many kWhs to generate a prediction.


## References

Hasselmann, K. [Stochastic climate models Part I. Theory](https://onlinelibrary.wiley.com/doi/10.1111/j.2153-3490.1976.tb00696.x). Tellus 28, 473–485 (1976).

Hyndman, R. J., Koehler, A, B. [Another look at measures of forecast accuracy](https://www.sciencedirect.com/science/article/pii/S0169207006000239). International Journal of Forecasting, 22, 4, 679-6988, (2006). 

Mayer, M. J., Yang, D. [Potential root mean square error skill score](https://pubs.aip.org/aip/jrse/article-abstract/16/1/016501/3267300/Potential-root-mean-square-error-skill-score?redirectedFrom=fulltext). J. Renewable Sustainable Energy 16, 016501 (2024).
