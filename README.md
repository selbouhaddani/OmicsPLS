[![GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
[![Build Status](https://travis-ci.org/selbouhaddani/OmicsPLS.svg?branch=master)](https://travis-ci.org/selbouhaddani/OmicsPLS)
[![codecov](https://codecov.io/gh/selbouhaddani/OmicsPLS/branch/master/graph/badge.svg)](https://codecov.io/gh/selbouhaddani/OmicsPLS)
[![Build Status](https://cranlogs.r-pkg.org/badges/last-month/OmicsPLS)](https://cran.r-project.org/package=OmicsPLS)

# OmicsPLS

Welcome to the OmicsPLS package page on Github. OmicsPLS contains an implementation of O2PLS.

O2PLS is designed for integrating two high dimensional datasets and the O2PLS package is an implementation of this method.

## Installing

Install the package with:

    install.packages("OmicsPLS")

Note that copying the quotation marks can sometimes go wrong, type them yourself to be sure.
If that doesn't work, ensure that the required packages are installed (`ggplot2` and `parallel`).
Otherwise try the development version on this GitHub repo:
    
    require(devtools)
    devtools::install_github("selbouhaddani/OmicsPLS")
    
Or download the source .tar.gz or the Windows binaries .zip at my ZippedPackages repo. 

## Loading

After installing the package, load the functions by typing

    library(OmicsPLS)
    
## Bugs, questions, complaints, etc

For questions, complaints, bugs, etc, you can e-mail me (<s.el_bouhaddani at lumc.nl>) or file an issue here at GitHub.

## Citing the package
When using the OmicsPLS R-package in your research, please cite the corresponding article by running command 

    citation("OmicsPLS")

or copy-paste:
Bouhaddani, S., Houwing-duistermaat, J., Jongbloed, G., Salo, P., Perola, M., & Uh, H.-W. (2016). Evaluation of O2PLS in Omics data integration. BMC Bioinformatics BMTL Supplement. DOI: 10.1186/s12859-015-0854-z

## Atlas of science publication
Also please see http://atlasofscience.org/simultaneous-integrated-analysis-of-biological-datasets-an-evaluation-of-o2pls/ for a gentle explanation and illustration of the just mentioned article.

# Updates

- We are working on a software tutorial paper and vignette to be published soon! Please see the **OmicsPLS_vignette.pdf** for a walkthrough of the OmicsPLS package with real data. In the meantime just email me for more information on the package.
