# Rankcluster

[![Travis build status](https://travis-ci.com/modal-inria/Rankcluster.svg?branch=master)](https://travis-ci.com/modal-inria/Rankcluster) [![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/modal-inria/Rankcluster?branch=master&svg=true)](https://ci.appveyor.com/project/modal-inria/Rankcluster)

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/Rankcluster)](https://cran.r-project.org/package=Rankcluster) [![Total Downloads](http://cranlogs.r-pkg.org/badges/grand-total/Rankcluster?color=blue)](http://cranlogs.r-pkg.org/badges/grand-total/Rankcluster) [![Downloads](https://cranlogs.r-pkg.org/badges/Rankcluster)](https://cran.rstudio.com/web/packages/Rankcluster/index.html)

The code was originally on an [R-forge repository](https://r-forge.r-project.org/projects/rankclust/).


This package proposes a model-based clustering algorithm for ranking data. 
Multivariate rankings as well as partial rankings are taken into account.
This algorithm is based on an extension of the Insertion Sorting Rank (ISR) model for ranking data, which is a meaningful
and effective model parametrized by a position parameter (the modal ranking, quoted by mu) and a dispersion parameter (quoted by pi).
The heterogeneity of the rank population is modelled by a mixture of ISR, whereas conditional independence assumption is considered for multivariate rankings.


## Installation

From github:
``` r
library(devtools)
install_github("modal-inria/Rankcluster", build_vignettes = TRUE)
```

From CRAN:
``` r
install.packages("Rankcluster", repos = "https://cran.rstudio.com")
```

## Vignettes

Once the package is installed, a vignette showing an example and one describing the data format are available using the R commands:

``` r
RShowDoc("Rankcluster", package = "Rankcluster")
RShowDoc("dataFormat", package = "Rankcluster")
```

## Credits

**Rankcluster** is developed by Quentin Grimonprez, Julien Jacques and Christophe Biernacki.

Copyright Inria - Université de Lille

## Licence

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
[GNU Affero General Public License](https://www.gnu.org/licenses/agpl-3.0.en.html) for more details.


## References

* J. Jacques and C. Biernacki (2012), Model-based clustering for multivariate partial ranking data, Inria Research Report n 8113. [link](https://hal.inria.fr/hal-00743384/document)
* C. Biernacki and J. Jacques (2013), A generative model for rank data based on sorting algorithm, Computational Statistics and Data Analysis, 58, 162-176. [link](https://www.sciencedirect.com/science/article/pii/S0167947312003118)
* J. Jacques, Q. Grimonprez and C. Biernacki (2014), Rankcluster: An R Package for Clustering Multivariate Partial Rankings, The R Journal 6:1, pages 101-110. [link](https://journal.r-project.org/archive/2014/RJ-2014-010/index.html)


