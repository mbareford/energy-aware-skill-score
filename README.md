# An Energy-aware Skill Score for Earth science climate models

## Introduction

We need a score that allows researchers to compare different codes running on different machines.

A score based on kWh is useful as it accounts for different clock frequencies, memory bandwidths and any other hardware features that contribute to performance.
Such a score is therefore hardware agnostic; it is also independent of the machine location as far as the financial cost of running a model is concerned.

Allowing for differences between code behaviour requires more thought. We must consider not only differences in functionality but also differences in input data.
A code’s computational profile will depend on its starting parameters. The fact that we are considering Earth science (ES) codes only (e.g. ocean, atmosphere and climate models)
at least makes the problem of finding a suitable skill score tractable.

First, we can assume that any ES code will output a certain amount of data every $n$ time steps, where the time step can represent years or hours depending on the model.
The amount of data generated is a reasonable proxy for the number of properties being simulated. We refer to this as the Raw Scientific Output (RSO), measured in bytes.
The associated performance score is simply the RSO per unit of energy consumed. We depart here from the more usual simulation time per unit of runtime as such a score
does not give us an indication of how much work is being done per timestep.


## Weighting the score by difficulty of prediction

Next, we explore how to assess the level of achievement implied by the model prediction. We start with a ground truth that can be compared with computer-generated predictions.
For a given property, $p$,  such as temperature or salinity, we can calculate the differences between a time series of $N$ observed values ($x_{t}$) and the corresponding forecasted values ($f_{t}$),
summarising these as a root mean square relative error, $p / N$, where

1. $$\varepsilon_{p} = \sum_{t=1}^{N} \left( \frac{x_{t} - f_{t}}{f_{t}} \right)^{2}$$ .