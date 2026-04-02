# An Energy-aware Skill Score for Climate Models

Assume a set of time-averaged observations for some property distributed over a grid, $\bar{x}_{i}$, where $i$ is the grid cell index covering the range $1-N$. The corresponding forecasted values are denoted by $\bar{f}_{i}$.

As regards the error part of the skill score, all values are averaged over some time period no shorter than one month. The two sets of monthly averages representing the model  and ground truths can then be compared to derive an error term,

```math
\xi_{i} =  \frac{\sum_{i=1}^{N} (\bar{x}_{i} - \bar{f}_{i})^2} {\sum_{i=1}^{N} (\bar{x}_{i} - \bar{x}_{\mu})^2} \;,
```
where $\bar{x}_{\mu}$ is the spatial average of the temporally-averaged ground truth observations.

The expression for $\xi_{i}$ is based on the coefficient of determination. If the averages derived from the model truth are far from the corresponding ground truth, $\xi_{i}$ will be large, unless the ground truth data deviates substantially from the spatial mean. Conversely, $\xi_{i}$ will be small if the model and ground truths match closely.

Knowing that $\xi_{i}$ starts from zero and is unbounded, we need a way to map $\xi_{i}$ to the range $0-1$, allowing us to derive a skill score.

```math
\chi_{i} = 1 - \frac{\xi_{i}}{\xi_{i}+1} \;.
```

The expression above avoids negative skill scores, since $\xi_{i}=\infty$ yields a zero score and $\xi_{i}=0$ gives a score of one.

However, given the distribution of $x_{i}$ over time, how difficult was it to make the successful prediction in the first place? Mayer & Yang (2024) propose making use of the autocorrelation of the observations, $\gamma_{i}(h)$, where $h$ is the lag value. The idea here is that $x_{i,j}$ (measured at time $j$) will exhibit some degree of correlation with $x_{i,j-h}$ (the same property measured at time $j-h$). The autocorrelation is (minus) one if the observations are fully (anti) correlated, and so, the difficulty of making a successful prediction is lowest when $\gamma = \pm1$. It follows that the difficulty is highest when the autocorrelation is zero. In fact, the difficulty is impossible for a randomly varying time series characterised by $\gamma = 0$.

The autocorrelation lag ($h$) should be expressed in units of the finest timescale captured by the ground truth data, e.g. if using the ERA5 dataset (Soci et al. 2024), $h$ will be some multiple of one hour. In this way, the skill score will reflect the true volatility. For any climate model, there is a set of relevant lag values, where the maximum $h$ value is determined by the simulation time. For the smallest lag values, the autocorrelation may approach zero: fluctuations can appear random over short time scales, e.g. minute-to-minute strong variations in wind speed. Nevertheless, the impact of short term variability is known to influence longer term variations as a consequence of Hasselmann’s stochastic theory (Hasselmann 1976). We expect therefore that the skill scores associated with longer lag values will reflect the model’s success in allowing long-term phenomena to be influenced by continuous short-term random excitations. In other words, the fidelity of a model on timescales corresponding to small $h$ values is still being captured by the skill scores based on long lag times.

We now modify the expression for $\chi_{i}$ by introducing the autocorrelation.

```math
\chi_{i,h} = \Bigg[1 - \frac{\xi_{i}}{\xi_{i}+1}\Bigg]\big(1-|\gamma_{i}(h)|\big)
```

All values of $h$ compatible with the ground truth time series can be used with the formulation above. The unweighted skill score (i.e. the expression in square brackets that partly determines $\chi_{i,h}$) is preserved when $\gamma_{i}(h)=0$. At the opposite extreme, $\chi_{i,h}=0$ when $\gamma_{i}(h)=\pm1$. A zero skill score can occur if and only if the observed values are perfectly correlated or anti-correlated. Note, the unweighted skill score itself can only become zero if $\xi_{i}=\infty$, which could only happen if the difference between a forecast and observation is also infinite, i.e. the simulation making the forecast has obviously failed and so is not fit to be scored.

The autocorrelation and skill score are calculated for each spatial position $i$ defined by the grid. Once this is done the score can be averaged spatially to give $\chi_{h}$.


## Visualisation

How best to visualise skill scores? We could draw [grouped bar plots](https://matplotlib.org/stable/gallery/lines_bars_and_markers/barchart.html), where each group corresponds to a particular autocorrelation lag and the bars within a group refer to different properties. The higher the bar the higher the skill score.


## References

Hasselmann, K. [Stochastic climate models Part I. Theory](https://onlinelibrary.wiley.com/doi/10.1111/j.2153-3490.1976.tb00696.x). Tellus 28, 473–485 (1976).

Hyndman, R. J., Koehler, A, B. [Another look at measures of forecast accuracy](https://www.sciencedirect.com/science/article/pii/S0169207006000239). International Journal of Forecasting, 22, 4, 679-6988, (2006). 

Mayer, M. J., Yang, D. [Potential root mean square error skill score](https://pubs.aip.org/aip/jrse/article-abstract/16/1/016501/3267300/Potential-root-mean-square-error-skill-score?redirectedFrom=fulltext). J. Renewable Sustainable Energy 16, 016501 (2024).

Soci, C., Hersbach, H., Simmons, A., Poli, P., Bell, B., Berrisford, P., et al. (2024) The ERA5 global reanalysis from 1940 to 2022. Quarterly Journal of the Royal Meteorological Society, 150(764), 4014–4048. Available from: https://doi.org/10.1002/qj.4803 .