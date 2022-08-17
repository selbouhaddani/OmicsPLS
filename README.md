[![GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
[![Build Status](https://travis-ci.org/selbouhaddani/OmicsPLS.svg?branch=master)](https://travis-ci.org/selbouhaddani/OmicsPLS)
[![codecov](https://codecov.io/gh/selbouhaddani/OmicsPLS/branch/master/graph/badge.svg)](https://codecov.io/gh/selbouhaddani/OmicsPLS)
[![Build Status](https://cranlogs.r-pkg.org/badges/last-month/OmicsPLS)](https://cran.r-project.org/package=OmicsPLS)
[![R-CMD-check](https://github.com/selbouhaddani/OmicsPLS/workflows/R-CMD-check/badge.svg)](https://github.com/selbouhaddani/OmicsPLS/actions)

# **Data integration with OmicsPLS**

Welcome to the OmicsPLS package page on Github. OmicsPLS contains an implementation of O2PLS.

OmicsPLS is designed for integrating two high dimensional datasets, visualise relationships between the features and inspect groupings across subjects. The OmicsPLS package can be applied in many fields, in particular life sciences. 

### Installing the package

Install the package with:

    install.packages("OmicsPLS")

Note that copying the quotation marks can sometimes go wrong, type them yourself to be sure.
If that doesn't work, ensure that the required packages are installed (`ggplot2` and `parallel`).
Otherwise try the development version on this GitHub repo:
    
    require(devtools)
    devtools::install_github("selbouhaddani/OmicsPLS")
    
Also, it is likely that in the GitHub version, bug fixes have been incorporated that are found in the CRAN version.

### Loading OmicsPLS

After installing the package, load the functions by typing

    library(OmicsPLS)
    
### Bugs, questions, complaints, etc

For questions, complaints, bugs, etc, you can e-mail me (<s.elbouhaddani at umcutrecht.nl>) or file an issue here at GitHub.



# **OmicsPLS software article**

[Our software article is published in BMC Bioinformatics!](https://doi.org/10.1186/s12859-018-2371-3) Please see the article for detailed explanation of the methods, functions and real data applications, including a vignette on omics data integration. 

### Citing the package

When using the OmicsPLS R package in your research, please cite the corresponding software article by running command 

    citation("OmicsPLS")

or copy-paste:
el Bouhaddani, S., Uh, H. W., Jongbloed, G., Hayward, C., Klarić, L., Kiełbasa, S. M., & Houwing-Duistermaat, J. (2018). Integrating omics datasets with the OmicsPLS package. BMC Bioinformatics, 19(1). https://doi.org/10.1186/s12859-018-2371-3

### Atlas of science publication

Also please see http://atlasofscience.org/simultaneous-integrated-analysis-of-biological-datasets-an-evaluation-of-o2pls/ for a gentle explanation and illustration of the just mentioned article.



# **Updates**

- Group Sparse O2PLS is implemented! Check out the newest version 2.0.2 of OmicsPLS for this penalized approach to O2PLS (both on CRAN and GitHub). See https://selbouhaddani.eu/2021-04-27-OmicsPLS-1.3.0-new/ for the updates. This version of OmicsPLS is on the main branch. 
- A *Probabilistic* version of O2PLS is available in the [PO2PLS GitHub](https://github.com/selbouhaddani/PO2PLS) repo. This approach is based on a rigorous statistical framework to integrate two datasets and obtain estimates of uncertainty. It has better performance than O2PLS - [see publication](https://doi.org/10.1111/rssc.12583) - in terms of prediction and interpretation if the number of components are correctly specified, and comes with a global test for the relationship between two datasets.
