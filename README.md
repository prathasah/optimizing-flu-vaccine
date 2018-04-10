# optimizing-flu-vaccine

## Age structured influenza transmission model with optimization procedure

**[Pratha Sah](https://prathasah.github.io/)
[\<pratha.sah@yale.edu\>](mailto:pratha.sah@yale.edu),
 Jan Medlock, Meagan C. Fitzpatrick, Burton H. Singer, Alison P. Galvani.
[The Center for Infectious Disease Modeling and Analysis](http://cidma.yale.edu/).**

**Copyright 2017-2019, Pratha Sah et al. All rights reserved.
Released under the [GNU AGPL 3](LICENSE).**

This repository contains the code and data used to simulate seasonal influenza in the United States and
optimize uptake of low efficacy vaccine across age-groups. Please cite the following paper, if you use our code in any
form or create a derivative work
> Pratha Sah, Jan Medlock, Meagan C. Fitzpatrick, Burton H. Singer, Alison P. Galvani.
> Optimizing the impact of low-efficacy influenza vaccines.
> PNAS. 2018.


Dependencies
===============================================
* [Python 2.7](http://python.org/)
* [NumPy & SciPy](http://www.scipy.org/)
* [matplotlib](http://matplotlib.org/)
* [pandas](http://pandas.pydata.org/)
* [ggplot2](http://ggplot2.org/)

Parameter data
===============================================

* [epidemiology_raw](epidemiology_raw) contains the data used to
  parametrize the model and references to the sources of these data.

* [sampled_parameter_set.csv](sampled_parameter_set.csv) contains 1000 random samples of the parameters that were used in our paper

