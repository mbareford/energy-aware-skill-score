# An Energy-aware Skill Score for Weather Forecasts and Climate Models

The motivation and mathematics behind the energy-aware skill score are covered in the [Introduction](intro.md).

Also stored here are two Python scripts, ([calc-skill-score-wf.py](calc-skill-score-wf.py)) and ([calc-skill-score-cm.py](calc-skill-score-cm.py)).
Those scripts are currently under development and are intended to calculate skill scores for weather forecast codes and climate models.

Essentially, the skill score calculation compares a "model" truth with a "ground" truth.
The latter is obtained from the ERA5 data archives held on the [JASMIN server](https://help.jasmin.ac.uk/docs/long-term-archive-storage/ceda-archive/).
Likewise, example model truths can be taken from CMIP model outputs also stored on JASMIN.

Lastly, there is [setup-pyenv.sh](setup-pyenv.sh), a script for installing a Python environment that will allow you to run the Python codes mentioned above.