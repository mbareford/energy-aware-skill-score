# An Energy-aware Skill Score for Weather Forecasts and Climate Models

## Introduction

We need a score that allows researchers to compare different codes running on different machines. A common performance metric is simulation time per
unit of runtime (e.g. years per day), but such metrics say nothing about the accuracy of the simulation or even how much compute is being done per timestep.
And, in any case, a *skill* score must capture the difficulty of the forecast. This last condition requires the existence of a ground truth that can be
compared with computer-generated predictions.

For weather forecast codes, we can use the autocorrelation of a ground truth time series as a measure of the predictive difficulty, which can then
be combined with the error (the difference between the computed results and the observations) to give a skill score. This approach is made easier by
the fact that the ground truth dataset (e.g. ERA5) and the output from the forecast code will both follow the Gregorian calendar. 

Climate models on the other hand are driven by the observed radiative forcing and for that reason, hour-by-hour comparisons with Earth-based measurements
would not be appropriate. What can be compared is data that has been averaged over some suitably long timeframe, one that is no shorter than one month.
That qualification also deals with the common practice of running climate models over 30-day months: comparing monthly averages works round the fact
that ground truth datasets will of course feature months that have different numbers of days (28-31).

Hence, different skill score calculations are needed for [weather forecasts](weather-forecasts.md) and [climate models](climate-models.md).


## Energy Awareness

Making the skill score energy aware is simply a matter of dividing by the energy consumed by the model run: a high skill score is penalised if a model requires many kWhs to generate a prediction.

A score based on kWh is useful as it accounts for different CPU/GPU clock frequencies, memory bandwidths and any other hardware features that contribute to performance.
Such a score is therefore hardware agnostic and it is also independent of the machine location as far as the financial cost of running a model is concerned.